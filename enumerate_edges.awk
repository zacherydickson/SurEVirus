#!/usr/bin/awk -f

BEGIN{
    FS="\t";
    virus="NC_007605.1"
    OFS="\t"
}

{
    reg=$1","$2","$3","$6;
    n=split($4,a,",");
    isViral=0;
    if($1==virus){isViral=1}
    for(i=1;i<=n;i++){
	if(isViral){
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
