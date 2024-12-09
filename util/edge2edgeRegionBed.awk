#!/usr/bin/awk -f

BEGIN {
    FS = "[ \t:]"
    OFS="\t"
}

(ARGIND == 1){
    maxIS=$2;
    next
}
(ARGIND == 2){
    cLen[$1]=$2;
    next
}

{
    for(i=1;i<=2;i++){
	split($i,a,",");
	s=a[2]-maxIS;
	if(s < 0){
	    s=0
	}
	e=a[3]+maxIS;
	if(e > cLen[a[1]]){
	    e = cLen[a[1]]
	}
	entry = a[1]"\t"s"\t"e"\t"$i"\t"".""\t"a[4]
	if(!Count[entry]++){
	    printf("%s\t%d\t%d\t%s\t%c\t%c\n",a[1],s,e,$i,".",a[4])
	}
    }
}
