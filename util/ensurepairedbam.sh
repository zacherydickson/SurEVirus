#!/bin/bash
set -o pipefail

if [ "$#" -lt 1 ]; then
	>&2 echo -e "Usage: $(basename $0) inbam > outbam\n\tRemoves any read ids that appear only once";
	exit 1;
fi

function main {
	bamFile=$1;
	samtools view -h "$bamFile" | 
		awk '
			(ARGIND == 1){Exclude[$1]=1;next}
			/^@/{print;next}
			(!Exclude[$1])
		' <(getRemoveList "$bamFile") /dev/stdin |
	samtools view -O BAM
}

function getRemoveList {
	bamFile=$1;
	samtools view "$bamFile" |
		awk -F '\\t' '
			{ Count[$1]++ }
			END {
				for (k in Count){
					if(Count[k] == 1){
						print k
					}
				}
		} '
}

main "$@"
