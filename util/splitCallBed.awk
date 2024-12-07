#!/usr/bin/awk -f

BEGIN{
    OFS="\t"
    if(!outdir){
	outdir="."
    }
}

function printbed(s,id,of, a){
    split(s,a,":");
    print a[1],a[3]-1,a[4],id,".",a[2] > of
}
{
    printbed($2,$1,outdir"/h.bed");
    printbed($3,$1,outdir"/v.bed");
}

