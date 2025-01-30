#include <iostream>
#include <stdexcept>
#include "sam_utils.h"
#include "str_utils.h"
#include <htslib/thread_pool.h>
#include <htslib/sam.h>
#include <string>
#include <thread>
#include <chrono>

//=== TYPE DECLARATIONS

//=== FUNCTION DECLARATIONS

std::string Bam2UID(const bam1_t* read, bam_hdr_t* hdr);
bool OutputThreadResult(hts_tpool_process* q, samFile * out, bam_hdr_t* hdr);
void* ProcessRead(void * arg);
void PrintUsage();


//=== MAIN

int main(int argc, const char* argv[]) {
    if(argc < 3){ //Check for valid Input
	PrintUsage();
	return 1;
    }
    std::string inFile(argv[1]);
    std::string outFile(argv[2]);
    int nThread = 2;
    if(argc >= 4){
	try {
	    nThread = std::stoi(argv[3]);
	} catch (std::invalid_argument & e) {
	    std::cerr << "[ERROR] nThreads must be a number\n";
	    return 1;
	}
    }
    //Turns out HTSlib has native support for "-" as stdin/stdout, just
    //pass through
    //Open the Input File
    open_samFile_t* in;
    try {
	in = open_samFile(inFile.c_str(),false,false);
    } catch(std::string & e) {
	std::cerr << "[ERROR] " << e << "\n";
	return 1;
    }
    //Open the output File
    samFile* out = sam_open(outFile.c_str(),"wb");
    if(sam_hdr_write(out,in->header) != 0) {
	std::cerr << "[ERROR] Could not write to " + outFile + "\n";
	return 1;
    }
    //Init the thread pool
    hts_tpool *p = hts_tpool_init(nThread);
    hts_tpool_process *q = hts_tpool_process_init(p,nThread*2,0);
    //Read all reads from the file, 
    bam1_t* curRead = bam_init1();
    while (sam_read1(in->file, in->header, curRead) >= 0) {
	int blk;
	bam1_t* inRead = bam_dup1(curRead);
	do{
	    blk = hts_tpool_dispatch2(p,q,ProcessRead,inRead,1);
	    if(!OutputThreadResult(q,out,in->header)){
	        ; // No result to output - wait maybe
	    }
	    if(blk == -1){ // The queue is full - wait a bit
		std::this_thread::sleep_for(std::chrono::milliseconds(1)); 
	    }
	} while(blk==-1);
    }

    //Wait for all jobs
    hts_tpool_process_flush(q);
    while(OutputThreadResult(q,out,in->header)) {
	; //Just keep processing results
    }
    //Clean up
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);
    bam_destroy1(curRead);
    sam_close(out);
    close_samFile(in);
}

//=== FUNCTION DEFINITIONS

//Given a bam entry construct a string which should uniquely identify it
//  Segment 0 if neither read1 or 2 is set, segment 3 if both are
//Input - a bam1_t* representing the read to make an id for
//	- a bam_hdr_t* to convert reference ids to strings
//Output - a string in the format: QName/Segment@Chr:Pos
std::string Bam2UID(const bam1_t* read, bam_hdr_t* hdr){
    if(!read) return "_NULL_";
    std::string qname = bam_get_qname(read);
    int segment = 0;
    if(read->core.flag & BAM_FREAD1) segment += 1;
    if(read->core.flag & BAM_FREAD2) segment += 2;
    std::string sID (sam_hdr_tid2name(hdr,read->core.tid));
    return  qname + '\\' + std::to_string(segment) + '@' + sID + ':' +
	    std::to_string(read->core.pos + 1);
}

//Checks if a process has returned a result and outputs the read if it has
//  Cleans up the read object
//Input a hts thread pool process
//Output - true if a result was processed
//	 - false if there was no result ready
bool OutputThreadResult(hts_tpool_process* q, samFile* out, bam_hdr_t* hdr){
    hts_tpool_result* r;
    if(!(r = hts_tpool_next_result(q))){
	return false;
    }
    bam1_t* resRead = (bam1_t*)hts_tpool_result_data(r);
    if(sam_write1(out,hdr,resRead) < 0){
        std::cerr << "[WARNING] Failed to write " + Bam2UID(resRead,hdr) << "\n"; 
    }
    bam_destroy1(resRead);
    hts_tpool_delete_result(r,0);
    return true;
}

//Given a bam entry determine the value of X2 and add/update the X2 tag
//This then calls the output function
//Input - a bam1_t read to process
//	- a samFile to which to write
//	- a header for the tid conversion
//Output- None, writes to samFile 
void* ProcessRead(void* arg){
    bam1_t* read = (bam1_t*) arg;
    if(!read) return read;
    //Remove existing X2
    uint8_t* pX2 = bam_aux_get(read,"X2");
    if(pX2) bam_aux_del(read,pX2);
    //Calc X2
    int x2 = 1;
    uint8_t* pXA = bam_aux_get(read,"XA");
    if(pXA) {
	std::string xaStr = bam_aux2Z(pXA);
	x2 += strsplit(xaStr,';').size();
    }
    //Update Read
    bam_aux_update_int(read,"X2",x2); 
    //Done
    return read; 
}

//Prints the Usage Message
//Inputs - None
//Outputs - None, writes to stderr
void PrintUsage() {
    std::cerr	<< "===Description\n"
	<< "\tGiven a [SB]AM file, add/update the X2 tag\n"
	<< "\tThe X2 tag is the integer number of alignments, determined as the length\n"
	<< "\t of the XA tag (if present) plus 1\n"
	<< "\tThe X2 tag will always be at the end of the output entries\n"
	<< "===Usage\n"
	<< "\taddX2tag in out [nThread=1]\n"
	<< "===ARGUMENTS\n"
	<< "\tin PATH\tPath to a [SB]AM file to process\n"
	<< "\tout PATH\tPath to an output file which will be BAM formatted\n"
	<< "\tNOTE: Use '-' to indicate that input/output are stdin/stdout\n"
	<< "\tnThread [1,∞)εZ = 2\tThe number of threads to use\n"
	<< "\t writing is the real bottleneck, so more than two is unlikely to be useful\n"
	<< "===Output\n"
	<< "\tOutput is BAM formatted\n"
	;
}
