#!/usr/bin/awk -f

BEGIN{
    FS="\t";
    OFS="\t"
}

(ARGIND == 1){
    isVirus[$1]=1;
    next;
}

{
    reg=$1","$2","$3","$6;
    n=split($4,a,",");
    for(i=1;i<=n;i++){
	if(isVirus[$1]){
	    viralRegions[a[i],++nVirReg[a[i]]]=reg;
	}else{
	    hostRegions[a[i],++nHostReg[a[i]]]=reg;
	}
    }
}

END{
    for(read in nVirReg){
        for(i=1;i<=nHostReg[read];i++){
	   for(j=1;j<=nVirReg[read];j++){
    	       edge = hostRegions[read,i]":"viralRegions[read,j];
    	       edgeReads[edge,++nEdgeRead[edge]]=read
    	   }
        }
    }
    for(edge in nEdgeRead){
	if(nEdgeRead[edge] < 3){continue}
        s=edgeReads[edge,1];
        for(i=2;i<=nEdgeRead[edge];i++){
	   s=s","edgeReads[edge,i]
        }
        print edge,s,nEdgeRead[edge]
    }
}
