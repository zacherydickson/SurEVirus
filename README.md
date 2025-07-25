# SurEVirus

## Description

SurEVirus (Survey of Edges supporting Virus integrations) is a modified reimplementation of SurVirus which is described in:

    "SurVirus: a repeat-aware virus integration caller" (2021). R. Rajaby, Y. Zhou, Y. Meng, X. Zeng, G. Li, P. Wu, and WK. Sung. Nucleic Acids Research (6). doi: 10.1093/nar/gkaa1237

The major underlying concepts which made SurVirus repeat aware are still present. This re-implementation aims to be more conservative and more efficient.

As compared to SurVirus an additional filter has been put into place: for chimeric or split read to be considered, it cannot be non-chimeric or non-split in any identified alternative alignments. The rationale being that viral integrations are rare events, if there is an alternate explanation for a read which does not involve an integration, it should not be used to support the existence of an integration.

Greater efficiency is achieved by considering the problem of calling integrations in a repeat aware setting to be a problem of identifying adequately supported edges in a bipartite graph. After regions of interest are identified in the host and viral genomes, chimeric and split reads lend support to an edge between host regions and viral regions. Processing these edges rather than regions can be done more efficiently and in parallel. This allows for a dramatic speedup relative to SurVirus in cases where there are many candidate regions in the host genome. This is often the case when considering viruses which have regions with high similarity to host regions; a common feature when integration is part of the viral strategy. Additional speedup is achieved by using a branched queue data structure when selecting the final set of edges.

The same requirements for supporting a junction are still present:
    1. All reads supporting a junction support the same orientation of host and virus sequences
    2. All reads supporting a junction are consistent (up to sequencing errors) with the consensus sequence of the junction breakpoints
    3. Any given template may be used to support one and only junction

If one was using fastq input, SurEVirus can be used almost as a drop-in replacement for SurVirus, as input formats and output formats are matched. However, BAM and CRAM support has been removed on the input side, and the meaning of SPLIT\_READS in the output is subtly different.

Another difference from SurVirus is that SurEVirus can be fully reproducible. SurVirus would allow for variability in the order of processing when performing operations in parallel which led to minor differences in bwa output. SurEVirus has reconfigured processing and added options (including specifying insert size parameters to prevent BWA from estimating them) which ensure that given identical input separate runs of SurEVirus will have identical output.

## Compiling

SurEVirus has been compiled and tested with gcc 11.4.0, so we recommend this version or higher.

First of all, the required external libraries (downloaded with the source code) must be compiled with
```
util/build_libs.sh
```

Then, run
```
cmake -DCMAKE_BUILD_TYPE=Release . && make
```

## Required software

Python 3 and libraries [NumPy](http://www.numpy.org), [PyFaidx](https://github.com/mdshw5/pyfaidx) and [PySam](https://github.com/pysam-developers/pysam>) are required. 

Recent versions of samtools, bwa and dust are required. The latest versions at the moment of writing are 1.13 for samtools, 0.7.18 for bwa.

For dust, we recommend [this](https://github.com/lh3/sdust) sdust implementation

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

## Preprocessing the input reads

SurEVirus does not perform any pre-processing on the input fastq files.
Inputs for for SurEVirus should first pre processed to remove adapters, low-quality bases and reads, and polynucleotide artifacts such as poly-A or Poly-Gs.

We suggest [fastp](https://github.com/OpenGene/fastp).

## Running

The bare minimum command for running SurEVirus is 
```
python surveyor reads_1.fq[.gz] reads_2.fq[.gz] /path/to/empty/workdir /path/to/host/reference /path/to/virus/reference /path/to/host+virus/reference 
```

reads 1 and 2 are fastq formatted fwd and reverse reads

## Options

If samtools, bwa or sdust are not in your PATH, or if you wish to provide an alternative location for either of them, you can do so with the `--samtools`, `--bwa` and `--dust` flags
```
python surveyor input_files /path/to/empty/workdir /path/to/host/reference /path/to/virus/reference /path/to/host+virus/reference --samtools /path/to/samtools --bwa /path/to/bwa --dust /path/to/sdust
```

A useful flag is --threads, which you can use to specify the number of threads you wish to assign to SurEVirus. The default is 1.

If there is a region of the viral or host region one does not wish to consider, these regions can be placed into a BED formatted file and provided to SurEVirus with the `--excluded_regions_bed` option.

If insert size and read length values are already known, SurEVirus does not need to estimate them, they can be provided with the `---isParams` option.
The parameters are provided as a comma separated list of 4 values:
- The Read Length
- Mean insert size
- Standard deviation from mean insert size
- The upper standard deviation from mean insert size ( std.dev calculated only for insert sizes larger than the
  mean )

Other options are provided for fine tuning of the insert sizes, clip sizes, allowed alternative alignments.
Use `surveyor --help` for more information.

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
