#!/usr/bin/awk -f

BEGIN {
    OFS=" "
}

(ARGIND == 1) {
    isVirus[$1]=1;
    next
}

(ARGIND == 2){
    v=0;
    if(isVirus[$1]) {
	v=1
    }
    Cov[$4,v]=$NF/($3-$2);
    next
}

(Cov[$1,0] > $(NF-1)) {
    $(NF-1) = Cov[$1,0]
}

(Cov[$1,1] > $NF) {
    $NF=Cov[$1,1]
}

1
