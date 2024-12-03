"""
    Core Dispatcher for SurEVirus
"""
#!/usr/bin/env python3

import argparse
import os
import sys
import pysam
import pyfaidx
import max_is_calc

def eprint(*args, **kwargs):
    "Prints messages to stderr"
    print(*args, file=sys.stderr, **kwargs)

cmd_parser = argparse.ArgumentParser(description='SurVirus, a virus integration caller.')
cmd_parser.add_argument('input_files', help='Input files, separated by a comma.')
cmd_parser.add_argument('workdir', help='Working directory for SurVirus to use.')
cmd_parser.add_argument('host_reference', help='Reference of the host organism in FASTA format.')
cmd_parser.add_argument('virus_reference', help='References of a list of viruses in FASTA format.')
cmd_parser.add_argument('host_and_virus_reference', help='Joint references of host and viruses.')
cmd_parser.add_argument('--threads', type=int, default=1, help='Number of threads to be used.')
cmd_parser.add_argument('--bwa', default='bwa', help='BWA path.')
cmd_parser.add_argument('--bedtools', help='Bedtools path.', default='bedtools')
cmd_parser.add_argument('--samtools', help='Samtools path.', default='samtools')
cmd_parser.add_argument('--dust', help='Dust path.', default='dust')
cmd_parser.add_argument('--excluded_regions_bed', \
                        help=   'Can use to exclude repetetive regions or uninteresting'
                                'portions of the genome')
cmd_parser.add_argument('--minClipSize', type=int, default=20, \
                        help='Min size for a clip to be considered.')
cmd_parser.add_argument('--maxClipAlt', type=int, default=6,
                        help='Max number of alternative alignments for a \
                                clip to be considered well mapped')
cmd_parser.add_argument('--maxSCDist', type=int, default=10, help='Max SC distance.')
cmd_parser.add_argument('--fq', action='store_true', help='Input is in fastq format.')
cmd_args = cmd_parser.parse_args()

SURVIRUS_PATH = os.path.dirname(os.path.realpath(__file__))

def execute(cmd):
    "Logs and executes a command"
    eprint("Executing:", cmd)
    if not os.system(cmd) == 0:
        eprint("Execution Failed - Exiting")
        sys.exit(1)

input_names = cmd_args.input_files.split(',')
if cmd_args.fq:
    if len(input_names) != 2:
        eprint("Two colon-separated fastq files are required.")
        sys.exit(1)

# Create config file in workdir
config_file_name = cmd_args.workdir + "/config.txt"
with open(config_file_name, "w", encoding ="utf-8") as config_file:
    config_file.write(f"threads {cmd_args.threads}\n")
    config_file.write(f"min_sc_size {cmd_args.minClipSize}\n")
    config_file.write(f"max_sc_dist {cmd_args.maxSCDist}\n")
    config_file.write(f"max_clip_alt {cmd_args.maxClipAlt}\n")

# Generate general distribution of insert sizes
with open(f"{cmd_args.workdir}/contig_map", "w", encoding= "utf-8") as contig_map:
    reference_fa = pyfaidx.Fasta(cmd_args.host_and_virus_reference)
    i = 1
    for k in list(reference_fa.keys()):
        contig_map.write("{k} {i}\n")
        i += 1

# count viruses
reference_fa = pyfaidx.Fasta(cmd_args.virus_reference)
n_viruses = len(list(reference_fa.keys()))

#Store the virusIDs
with open(f"{cmd_args.workdir}/virusNames.list","w",encoding="utf-8") as virus_names:
    for key in reference_fa.keys():
        virus_names.write(f"{key}\n")
#Store the lengths of all contigs
with pyfaidx.Fasta(cmd_args.host_and_virus_reference) as joint_ref, \
     open(f"{cmd_args.workdir}/contig-lengths.tab","w",encoding="utf-8") as contig_lengths:
    for key,seq in joint_ref.items():
        contig_lengths.write(f"{key}\t{seq.len()}\n")

#Construct Bam workspace for fastq input
bam_workspaces = []
if cmd_args.fq:
    bam_workspace = f"{cmd_args.workdir}/bam_0/"
    if not os.path.exists(bam_workspace):
        os.makedirs(bam_workspace)
    bam_workspaces.append(bam_workspace)

    max_read_len, max_is = \
        max_is_calc.get_max_is_from_fq(cmd_args.workdir, input_names[0], input_names[1], \
        cmd_args.host_and_virus_reference, cmd_args.bwa, cmd_args.threads)
    with open(f"{bam_workspace}/stats.txt", "w", encoding="utf-8") as stat_file:
        stat_file.write(f"max_is {max_is}\n")
    with open(config_file_name, "a", encoding="utf-8") as config_file:
        config_file.write(f"read_len {max_read_len}\n")

    isolate_cmd = f"{SURVIRUS_PATH}/isolate_relevant_pairs_fq \
                    {input_names[0]} {input_names[1]} {cmd_args.host_reference} \
                    {cmd_args.virus_reference} {cmd_args.workdir} {bam_workspace}"
    execute(isolate_cmd)

#Define function to map clips
def map_clips(prefix, reference):
    "Maps clips with a given prefix to the given reference"
    bwa_aln_cmd = f"{cmd_args.bwa} aln -t {cmd_args.threads} {reference} {prefix}.fa \
                    -f {prefix}.sai"
    bwa_samse_cmd = f"{cmd_args.bwa} samse -n {cmd_args.maxClipAlt} {reference} \
                        {prefix}.sai {prefix}.fa | \
                        {cmd_args.samtools} view -b -F 2304 > {prefix}.full.bam"
    execute(bwa_aln_cmd)
    execute(bwa_samse_cmd)

    filter_unmapped_cmd = f"{cmd_args.samtools} view -b -F 4 {prefix}.full.bam > \
                                {prefix}.aln.bam"
    execute(filter_unmapped_cmd)

    dump_unmapped_fa = f"{cmd_args.samtools} fasta -f 4 {prefix}.full.bam > \
                            {prefix}.unmapped.fa"
    execute(dump_unmapped_fa)

    bwa_mem_cmd = f"{cmd_args.bwa} mem -t {cmd_args.threads} {reference} \
                    {prefix}.unmapped.fa | \
                    {cmd_args.samtools} view -b -F 2308 > {prefix}.mem.bam"
    execute(bwa_mem_cmd)

    cat_cmd = f"{cmd_args.samtools} cat {prefix}.aln.bam {prefix}.mem.bam  >| {prefix}.bam"
    execute(cat_cmd)

    pysam.sort("-@", str(cmd_args.threads), "-o", f"{prefix}.cs.bam", f"{prefix}.bam")

#Iterate over workspaces
for bam_workspace in bam_workspaces:
    bwa_cmd = f"{cmd_args.bwa} mem -t {cmd_args.threads} \
                {cmd_args.host_and_virus_reference} {bam_workspace}/retained-pairs_1.fq \
                {bam_workspace}/retained-pairs_2.fq | \
                {cmd_args.samtools} view -b -F 2304 | samtools sort -n - | \
                samtools fixmate -m - - | samtools sort - | \
                samtools markdup -Sr - {bam_workspace}/retained-pairs.remapped.bam"
    execute(bwa_cmd)

    samtools_sort_cmd = f"{cmd_args.samtools} sort -@ {cmd_args.threads} \
                          -o {bam_workspace}/retained-pairs.remapped.cs.bam \
                          {bam_workspace}/retained-pairs.remapped.bam"
    execute(samtools_sort_cmd)

    #Name sort the pairs for later processing
    samtools_sort_cmd = f"{cmd_args.samtools} sort -n -@ {cmd_args.threads} \
                        -o {bam_workspace}/retained-pairs.namesorted.bam \
                        {bam_workspace}/retained-pairs.remapped.bam"
    execute(samtools_sort_cmd)

    extract_clips_cmd = f"{SURVIRUS_PATH}/extract_clips {cmd_args.virus_reference} \
                            {cmd_args.workdir} {bam_workspace}"
    execute(extract_clips_cmd)

    # map virus clips
    map_clips(f"{bam_workspace}/virus-clips", cmd_args.host_reference)

    # map host clips
    map_clips("{bam_workspace}/host-clips", cmd_args.virus_reference)

#Extract candidate host and viral positions
extract_regions_cmd = f"{SURVIRUS_PATH}/extract_regions \
        {cmd_args.virus_reference} {cmd_args.workdir} {bam_workspace}"
execute(extract_regions_cmd)

## Filter out any junctions in annotated repeat regions
if cmd_args.exclude_regions_bed is not None:
    bedtools_cmd = f"{cmd_args.bedtools} subtract -A \
                    -a {cmd_args.workdir}/junction-candidates.bed \
                    -b {cmd_args.exclude_regions_bed} >| \
                    {cmd_args.workdir}/junction-candidates.cs.bed"
    os.replace(f"{cmd_args.workdir}/junction-candidates.cs.bed", \
                "{cmd_args.workdir}/junction-candidates.bed")

# Collapse Junctions together into regions
cluster_junctions_cmd = f"{SURVIRUS_PATH}/cluster_junctions \
        {cmd_args.virus_reference} {cmd_args.workdir} {bam_workspace}"
execute(cluster_junctions_cmd)

## Split Regions which appear to overlap with repeat regions
if cmd_args.exclude_regions_bed is not None:
    bedtools_cmd = f"{cmd_args.bedtools} subtract \
                    -a {cmd_args.workdir}/regions-candidates.bed \
                    -b {cmd_args.exclude_regions_bed} >| \
                    {cmd_args.workdir}/regions-candidates.cs.bed"
    os.replace(f"{cmd_args.workdir}/regions-candidates.cs.bed", \
                "{cmd_args.workdir}/regions-candidates.bed")

## Pair up host-viral regions and assign reads to each edge, filter edges
##  with too few reads
enum_edges_cmd = f"{SURVIRUS_PATH}/enumerate_edges.awk \
                    {cmd_args.workdir}/virusNames.list \
                    {cmd_args.workdir}/region-candidates.bed | \
                    sort -k3,3nr |> {cmd_args.workdir}/edges.tab"
execute(enum_edges_cmd)

#Get the sequences of the regions of interest
extract_regions_cmd = f"{SURVIRUS_PATH}/edge2edgeRegionBed.awk \
                            {cmd_args.workdir}/bam_0/stats.txt \
                            {cmd_args.workdir}/contig-lengths.tab \
                            {cmd_args.workdir}/edges.tab | \
                        {cmd_args.bedtools} getfasta -s -name \
                            -fi {cmd_args.host_and_virus_reference} -bed - >| \
                        {cmd_args.workdir}/regions.fna"
execute(extract_regions_cmd)

##Extract the fasta sequences of the reads on the edges
extract_reads_cmd = f"{cmd_args.samtools} view \
                            {cmd_args.workdir}/bam_0/retained-pairs.namesorted.bam | \
                        {SURVIRUS_PATH}/extractEdgeReads.awk \
                            {cmd_args.workdir}/virusNames.list \
                            {cmd_args.workdir}/edges.tab /dev/stdin >| \
                        edge_reads.fna"
execute(extract_reads_cmd)

readsx = cmd_args.workdir + "/readsx"
if not os.path.exists(readsx):
    os.makedirs(readsx)

##Edgemapper Command
edge_mapper_cmd = f"{SURVIRUS_PATH}/edge_mapper \
        {cmd_args.virus_reference} {cmd_args.workdir} {bam_workspace}"
execute(edge_mapper_cmd)

##Dust Filter the breakpoint sequences
dust_cmd = f"{cmd_args.dust} {cmd_args.workdir}/host_bp_seqs.fa > \
        {cmd_args.workdir}/host_bp_seqs.masked.bed"
execute(dust_cmd)

dust_cmd = f"{cmd_args.dust} {cmd_args.workdir}/virus_bp_seqs.fa > \
        {cmd_args.workdir}/virus_bp_seqs.masked.bed"
execute(dust_cmd)

####Retreive the breakpoint region coverage
#Pull out a bed formated version of the called breakpoints, separately for host and virus
split_cmd = f"{SURVIRUS_PATH}/splitCallBed.awk {cmd_args.workdir}/results.txt"
execute(split_cmd)
calc_coverage_cmd = f"{{ {cmd_args.bedtools} slop -s -l {max_is} -r 0 \
                                -i {cmd_args.workdir}/h.bed \
                                -g {cmd_args.workdir}/contig-lengths.tab; \
                            {cmd_args.bedtools} slop -s -l 0 -r {max_is} \
                                -i {cmd_args.workdir}/v.bed \
                                -g {cmd_args.workdir}/contig-lengths.tab; \
                        }} | \
                        {cmd_args.samtools} bedcov -d 1 - \
                            {cmd_args.workdir}/bam_0/retained-pairs.remapped.cs.bam | \
                        {SURVIRUS_PATH}/adjustCoverage.awk \
                            {cmd_args.workdir}/virusNames.list /dev/stdin \
                            {cmd_args.workdir}/results.txt >| \
                        results.remapped.txt"
execute(calc_coverage_cmd)

filter_cmd = f"{SURVIRUS_PATH}/filter {cmd_args.workdir} > \
                {cmd_args.workdir}/results.t1.txt"
execute(filter_cmd)

filter_cmd = f"{SURVIRUS_PATH}/filter {cmd_args.workdir} --remapped > \
                {cmd_args.workdir}/results.remapped.t1.txt"
execute(filter_cmd)

filter_cmd = f"{SURVIRUS_PATH}/filter {cmd_args.workdir} --print-rejected > \
                {cmd_args.workdir}/results.discarded.txt"
execute(filter_cmd)

eprint("Done")
