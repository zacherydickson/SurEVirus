#include <iostream>
#include "sam_utils.h"
#include <htslib/sam.h>
#include <memory>
#include <unordered_set>
#include <vector>

//=== TYPE DECLARATIONS

typedef std::vector<bam1_t*> bamVec_t;

class CReadBlock {
    //TODO: Reimplement so that adding sequences does the processing
    //as it goes
    //	essentially a two segment approach
    //Members
    public:
	const std::string name;
    protected:
	std::vector<bam1_t*,2>
	bam1_t* segment1;
	bam1_t* segment2;
    //Con/Destruction
    public:
	CReadBlock(const std::string & name) :
	    name(name), segment1(nullptr),segment2(nullptr) {}
	~CReadBlock() {this->destroy();}
    //Accsessors
	bool 
	bamVec_t::const_iterator cbegin() const {return this->m_Buf.cbegin();}
	bamVec_t::const_iterator cend() const {return this->m_Buf.cend();}
	size_t size() const {return this->m_Buf.size();}
    //Public Methods
    public:
	void addRead(const bam1_t* read) {this->m_Buf.push_back(bam_dup1(read));}
	void destroy() {
	    if(this->m_Buf.size()){
		for(bam1_t* & read : this->m_Buf){
		    bam_destroy1(read);
		}
		this->m_Buf.clear();
	    }
	}
	void process();
};

void CReadBlock::process() {

    //TODO: Implement
}

//=== FUNCTION DECLARATIONS

std::string Bam2UID(const bam1_t* read, bam_hdr_t* hdr);
void OutputReadBlock(const CReadBlock & block, samFile * out, bam_hdr_t* hdr);
void PrintUsage();

//=== GLOBAL CONSTANTS

const uint16_t BAM_FNONPRIMARY = BAM_FSECONDARY & BAM_FSUPPLEMENTARY;

//=== MAIN

int main(int argc, const char* argv[]) {
    if(argc < 3){ //Check for valid Input
	PrintUsage();
	return 1;
    }
    std::string inFile(argv[1]);
    std::string outFile(argv[2]);
    //TODO: Handle '-' for stdin or stdout
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
    //Read all reads from the file, process all blocks of reads with the
    //same read id
    bam1_t* curRead = bam_init1();
    std::unique_ptr<CReadBlock> pCurBlock(nullptr);
    std::unordered_set<std::string> observedQNames;
    while (sam_read1(in->file, in->header, curRead)) {
	std::string qname = bam_get_qname(curRead);
	//Handle the First Read
	if(!pCurBlock){ 
	    pCurBlock = std::make_unique<CReadBlock>(qname);
	}
	//If properly sorted, then each qname will appear in only 
	//  one block
	if(observedQNames.count(qname)){
	    std::cerr << "[ERROR] Input file is not sorted by name\n";
	    return 1;
	}
	//Skip Secondary/Supplementary reads
	if(curRead->core.flag & BAM_FNONPRIMARY) continue;
	//Skip Reads which are marked neither/both Segment 1 or 2
	if( (curRead->core.flag & BAM_FREAD1) ==
	    (curRead->core.flag & BAM_FREAD2)) continue;
	//Check if we have entered a new read block
	if(qname != pCurBlock->name){
	    //Process the Read Block
	    pCurBlock->process();
	    //Output the remaining Reads in the Block
	    OutputReadBlock(*pCurBlock,out,in->header); 
	    //Track the Blocks we've seen as a check for proper sorting
	    observedQNames.insert(pCurBlock->name);
	    //Replace the current block, destroying the old one
	    pCurBlock = std::make_unique<CReadBlock>(qname);
	}
	//Store a copy of the current read in the block
	pCurBlock->addRead(bam_dup1(curRead));
    }
    //Clean up
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
    std::string qname = bam_get_qname(read);
    int segment = 0;
    if(read->core.flag & BAM_FREAD1) segment += 1;
    if(read->core.flag & BAM_FREAD2) segment += 2;
    std::string sID (sam_hdr_tid2name(hdr,read->core.tid));
    return  qname + '\\' + std::to_string(segment) + '@' + sID + ':' +
	    std::to_string(read->core.pos);
}

//Iterate over reads in a read block and output to the bam output
//Input - a CReadBlock reference to output
//	- a samFile* to which to write
//	- a bam_hdr_t* for Chr conversion
//Output - None, writes to samFile
void OutputReadBlock(const CReadBlock & block, samFile * out, bam_hdr_t* hdr){
    for(auto it = block.cbegin(); it != block.cend(); it++){
	if(sam_write1(out,hdr,*it) != 0){
	    std::cerr << "[WARNING] Failed to write " + Bam2UID(*it,hdr) << "\n"; 
    	}
    }
}

//Prints the Usage Message
//Inputs - None
//Outputs - None, writes to stderr
void PrintUsage() {
    std::cerr	<< "===Description\n"
	<< "\tGiven a name-sorted [SB]AM file, ensure each read name\n"
	<< "\t and has exactly one primary alignment for each segment.\n"
	<< "\tBWA mem occasionally will have more than one entry per segment,\n"
	<< "\t neither of which is marked as secondary or supplementary.\n"
	<< "\tThis program collapses such reads and adds other mappings to the\n"
	<< "\t supplementary alignment tag (XA).\n"
	<< "===Usage\n"
	<< "\tcleanBAM in out\n"
	<< "===ARGUMENTS\n"
	<< "\tin PATH=-\tPath to a [SB]AM file to process\n"
	<< "\tout PATH=-\tPath to an output file which will be BAM formatted\n"
	<< "\tNOTE: Use '-' to indicate that input/output are stdin/stdout\n"
	<< "===Output\n"
	<< "\tOutput is BAM formatted, with only primary alignments.\n"
	<< "\t As a result, Alignments matching flag 0xF00 are excluded.\n"
	;
}
