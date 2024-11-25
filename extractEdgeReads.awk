#!/usr/bin/awk -f

BEGIN {
    FS="\t";
    IS_READ1=0x40;
    IS_READ2=0x80;
}
(ARGIND == 1){
    isViralContig[$1]=1;
    next
}
(ARGIND == 2) {
    n=split($2,readList,",");
    for(i=1;i<=n;i++){
	name=readList[i];
	if(substr(name,length(name)-1,1) == "_"){
	    char=substr(name,length(name));
	    name=substr(name,1,length(name)-2);
	    if(char=="1"){
		useR1[name]=1
	    }
	    if(char=="2"){
		useR2[name]=1
	    }
	} else {
	    isPaired[name]=1
	}
	InSet[name]=1;
    }
    next
}
(   InSet[$1] && \
    (	isPaired[$1] || \
	(useR1[$1] && and($2,IS_READ1)) || \
	(useR2[$1] && and($2,IS_READ2)) \
    ) \
) {
    suffix="1";
    if(isPaired[$1]){
	suffix="H";
	if(isViralContig[$3]){
	    suffix="V"
	}
    } else if(and($2,IS_READ2)){
	suffix="2"
    }
    print ">"$1"_"suffix;
    print $10;
}
