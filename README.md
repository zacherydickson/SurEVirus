# SurEVirus

## Description

SurEVirus (Survey of Edges supporting Virus integrations) is a modified reimplementation of SurVirus which is described in:

    "SurVirus: a repeat-aware virus integration caller" (2021). R. Rajaby, Y. Zhou, Y. Meng, X. Zeng, G. Li, P. Wu, and WK. Sung. Nucleic Acids Research (6). doi: 10.1093/nar/gkaa1237

The major underlying concepts which made SurVirus repeat aware are still present. This re-implementation aims to be more conservative and more efficient.

As compared to SurVirus an additional filter has been put into place which a read to be chimeric or split in the primary alignment, and all identified alternative alignments. The rationale being that viral integrations are rare events, if there is an alternate explanation for a read which does not involve an integration, it should not be used to support the existance of an integration.

Greater efficiency is acheived by considering the problem of calling integrations in a repeat aware setting to be a problem of identifying adequately supported edges in a bipartite graph. After regions of interest are identified in the host and viral genomes, chimeric and split reads lend support to an edge between host regions and viral regions. Processing these edges rather than regions can be done more efficiently and in parallel. This allows for a dramatic speedup relative to SurVirus in cases where there are many candidate regions in the host genome. This is often the case when considering viruses which have regions with high similarity to host regions; a common feature when integration is part of the the viral strategy.

The same requirements for supporting a junction are still present:
    1. All reads supporting a junction support the same orientation of host and virus sequences
    2. All reads supporting a junction are consistent (up to sequencing errors) with the consensus sequence of the junction breakpoints
    3. Any given template supports one and only junction

For the most part, SurEVirus can be used as a drop-in replacement for SurVirus, as input formats and output formats are matched. However, CRAM support has been removed on the input side, and the meaning of SPLIT\_READS in the output is subtly different.

## Compiling

SurEVirus has been compiled and tested with gcc 11.4.0, so we recommend this version or higher.

First of all, the required external libraries (downloaded with the source code) must be compiled with
```
utils/build_libs.sh
```

Then, run
```
cmake -DCMAKE_BUILD_TYPE=Release . && make
```

## Preparing the references

SurEVirus needs three references:
1) the host genome
2) the virus(es) reference, one fasta sequence per virus
3) a concatenation of the host genome and the viruses genomes

Each of the references should be indexed with bwa and samtools. For example, suppose the host genome is contained in a file host.fa, and the virus genomes are in virus.fa. You should run
```
bwa index host.fa
samtools faidx host.fa

bwa index virus.fa
samtools faidx virus.fa

cat host.fa virus.fa > host+virus.fa
bwa index host+virus.fa
samtools faidx host+virus.fa
```

## Required software

Python 3 and libraries NumPy (http://www.numpy.org/), PyFaidx (https://github.com/mdshw5/pyfaidx) and PySam (https://github.com/pysam-developers/pysam) are required. 

Recent versions of samtools, bwa and sdust are required. The latest versions at the moment of writing are 1.13 for samtools, 0.7.18 for bwa.
For sdust, we recommend the implementation at https://github.com/lh3/sdust

## Running

The bare minimum command for running SurEVirus is 
```
python surveyor input_files /path/to/empty/workdir /path/to/host/reference /path/to/virus/reference /path/to/host+virus/reference 
```

input_files can be a list of comma-separated bam_files, for example
```
python surveyor	input1.bam,input2.bam /path/to/empty/workdir /path/to/host/reference /path/to/virus/reference /path/to/host+virus/reference
```
Note that if multiple BAM files are present, they must all be relative to the same sample.

input_files can also be a pair of comma-separated fastq files, containing read pairs
```
python surveyor reads_1.fq,reads_2.fq /path/to/empty/workdir /path/to/host/reference /path/to/virus/reference /path/to/host+virus/reference --fq
```
Note that in this second case, the flag --fq must be provided.

If samtools, bwa or sdust are not in your PATH, or if you wish to provide an alternative location for either of them, you can do so with the --samtools, --bwa and --dust flags
```
python surveyor input_files /path/to/empty/workdir /path/to/host/reference /path/to/virus/reference /path/to/host+virus/reference --samtools /path/to/samtools --bwa /path/to/bwa --dust /path/to/sdust
```

Finally, an useful flag is --threads, which you can use to specify the number of threads you wish to assign to SurEVirus. The default is 1.

## Output

The final output will be placed in the workdir, under the name results.remapped.t1.txt.

The following line is an example of a predicted integration:
```
ID=0 chr13:-73788865 type16:-3383 SUPPORTING_PAIRS=700 SPLIT_READS=725 HOST_PBS=0.927897 COVERAGE=0.684993
```

The ID is simply a unique number. The second and the third fields are the host and the virus coordinates for the integration. In this example, the sample contains a junction between chr13:73788865, negative strand, and type16:3383, negative strand.
Supporting pairs is the number of read pairs that supports the integration.
While split reads is the subset of split reads that overlap the junction (i.e. they map partly to the host breakpoint and partly to the virus breakpoint).

HOST_PBS and COVERAGE are quality metrics. Since they are used in the filtering of the results, the user can safely ignore them most of the time. 
HOST_PBS is interpretable as the fraction of base matches between the supporting reads that are mapped to the host breakpoint and the reference sequence. In the example, a value of 0.927897 means that when performing a local alignment between the supporting reads and the host reference sequence near the breakpoint, we produce ~92.8% base matches (as opposed to mismatches and gaps). It is actually calculated based on alignment scores.
SurEVirus internally analyses the distribution of insert sizes in the input sample, and determines the maximum insert size (maxIS) that is considered not to be an artifact (i.e. what is usually referred to as discordant by most SV callers). When remapping reads, SurEVirus considers a maxIS bp-long window next to the suspected breakpoint, for both virus and host. The rationale behind this is that if a read is more than maxIS bp away from a breakpoint, its pair would not be able to be chimeric, as the mate would not be able to cross the breakpoint.
COVERAGE is the fraction of such maxIS window that is covered by reads, averaged between host and virus.
