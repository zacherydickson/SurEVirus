#!/bin/bash

execDir=$(dirname $(readlink -f $0));
adjCov="$execDir/adjustCoverage.awk"

if [ "$#" -lt 11 ]; then
    >&2 echo "Usage: $(basename $0) bedtools samtools threads tmpBed maxIS ctgLen hBed vBed bamFile vNames rawRes > adjRes";
    exit 1;
fi

bedtools=$1;shift
samtools=$1;shift
threads=$1;shift
tmpBed=$1;shift
maxIS=$1;shift
ctgLen=$1;shift
hBed=$1;shift
vBed=$1;shift
bamFile=$1;shift
vNames=$1;shift
rawRes=$1;shift


"$bedtools" slop -s -l "$maxIS" -r 0 -i "$hBed" -g "$ctgLen" >| "$tmpBed";
"$bedtools" slop -s -l 0 -r "$maxIS" -i "$vBed" -g "$ctgLen" >> "$tmpBed";
"$samtools" depth -@ "$threads" -b "$tmpBed" "$bamFile" |
    awk '($NF){print $1"\t"$2-1"\t"$2}' |
    "$bedtools" intersect -a "$tmpBed" -b /dev/stdin -wa -wb |
    awk '   function out(){if(!lastiv){return} print lastiv"\t"s}
            {iv=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}
            (iv != lastiv){out();s=0;lastiv=iv} {s++} END{out()}' |
    awk '   BEGIN {OFS=" "}
            (ARGIND == 1) {isVirus[$1]=1;next}
            (ARGIND == 2){v=0;if(isVirus[$1]){v=1}Cov[$4,v]=$NF/($3-$2);next}
            (Cov[$1,0] > $(NF-1)) { $(NF-1) = Cov[$1,0] }
            (Cov[$1,1] > $NF) { $NF=Cov[$1,1] } 
            1
        ' "$vNames" /dev/stdin "$rawRes"
