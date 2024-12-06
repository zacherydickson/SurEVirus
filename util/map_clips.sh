#!/bin/bash

if [ "$#" -lt 3 ]; then
	>&2 echo "Usage: $(basename $0) prefix reference max_alt_clip";
	exit 1;
fi

prefix=$1;shift;
ref=$1;shift;
max_alt=$1;shift

#Map with ALN for short reads
bwa aln -t 30 $ref $prefix.fa -f $prefix.sai;
bwa samse -n $max_alt $ref $prefix.sai $prefix.fa | samtools view -b -F 2304 >| $prefix.full.bam;
#Extract the mapped and unmapped reads
samtools view -b -F 4 $prefix.full.bam >| $prefix.aln.bam;
samtools fasta -f 4 $prefix.full.bam >| $prefix.unmapped.fa;
#Map with mem to recover more reads
bwa mem -t 30 $ref $prefix.unmapped.fa | 
	awk -v ma=$max_alt '
		{n=0} /XA:Z/ {match($0,"XA:Z:(.+;)",arr); n=split(arr[1],a,";");} (n<=ma+1)
	' | samtools view -b -F 2308 >| $prefix.mem.bam;
#Combine the results
samtools cat $prefix.aln.bam $prefix.mem.bam >| $prefix.bam;
#Sort the output
samtools sort -@ 30 -o $prefix.cs.bam $prefix.bam;


