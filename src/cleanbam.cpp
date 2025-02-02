#include <algorithm>
#include <array>
#include <iostream>
#include "sam_utils.h"
#include <htslib/sam.h>
#include <memory>
#include <set>
#include <unordered_set>
#include <vector>
#include "utils.h"

//FUNCTION DECL
std::string Bam2UID(const bam1_t* read, bam_hdr_t* hdr);

//=== TYPE DECLARATIONS

typedef std::unordered_set<std::string> VirusNameSet_t;
typedef std::vector<bam1_t*> bamVec_t;
typedef std::unique_ptr<std::set<std::string>> XAStrSet_pt;
typedef std::array<XAStrSet_pt,2> XAStrSetArr_t;
typedef std::array<std::string,2> PrimaryXAStrArr_t;

class CReadBlock {
    //StaticMembers
    public:
	static const uint16_t BAM_FNONPRIMARY = BAM_FSECONDARY | BAM_FSUPPLEMENTARY;
	static const uint16_t BAM_FANYUNMAP = BAM_FUNMAP | BAM_FMUNMAP;
    //Members
    public:
	const std::string name;
    protected:
	bamVec_t m_Buf;
    //Con/Destruction
    public:
	CReadBlock(const std::string & name) :
	    name(name), m_Buf() {}
	~CReadBlock() {this->clear();}
    //Accsessors
	bamVec_t::const_iterator cbegin() const {return this->m_Buf.cbegin();}
	bamVec_t::const_iterator cend() const {return this->m_Buf.cend();}
	size_t size() const {return this->m_Buf.size();}
    //Public Methods
    public:
	void addRead(const bam1_t* read);
	void clear();
	bool process(bam_hdr_t* hdr,const VirusNameSet_t & virusNameSet);
    protected:
	bool selectBestMates();
    //Static Methods
    public:
	static void addBamAltstoXASet(bam1_t* read, std::set<std::string> & XASet);
	static int compareBams(bam1_t* & a, bam1_t* & b);
        static bool filterOnXASet(  XAStrSetArr_t & xaStrSetArr,
                                    const PrimaryXAStrArr_t & primeXAStrArr,
                                    const VirusNameSet_t & virusNameSet);
	static std::string primaryBam2xaStr(bam1_t* read, bam_hdr_t* hdr);
};

//Process the XA and SA tags of a read to extract the alignment strings
//Then adds these to a provided set of strings
//Input - a read object to process
//	- a reference set of strings
//Output- None, modifes the string set
void CReadBlock::addBamAltstoXASet(bam1_t* read, std::set<std::string> & xaSet){
    for (char tag1 : {'S','X'}){
        char tag[2] = {tag1,'A'};
        uint8_t* pTag = bam_aux_get(read,tag);
        if(!pTag) continue;
        for(auto & str : strsplit(bam_aux2Z(pTag),';')){
	    if(tag1 == 'S'){ //SA strings are formatted differently
		std::vector<std::string> s = strsplit(str,',');
		str =	s.at(0) + ',' + s.at(2) + s.at(1) + ',' +
			s.at(3) + ',' + s.at(5);
	    }
            xaSet.insert(str);
        }
    }
}

void CReadBlock::addRead(const bam1_t* read){
    //Refuse to add reads which are marked neither/both Segment 1 or 2
    if( (read->core.flag & BAM_FREAD1) ==
        (read->core.flag & BAM_FREAD2)) return;
    this->m_Buf.push_back(bam_dup1(read));
}

void CReadBlock::clear() {
    if(this->m_Buf.size()){
        for(bam1_t* & read : this->m_Buf){
            bam_destroy1(read);
        }
	this->m_Buf.clear();
    }
}

//Compares Two Bams in a way such that they can be sorted by 'goodness'
// for processing.
// -1: Better, 0: The Same, 1: Worse
// Primary alignments are better
// Alignments with both segments mapped are better
// Higher Scoring alignments are better
// Segment 1 is better than segment 2
// Earlier References (by header order) are better
// Lower mapping positions are better
// Otherwise equal
int CReadBlock::compareBams(bam1_t* & a, bam1_t* & b){
    uint16_t nonPrimeResA = a->core.flag & CReadBlock::BAM_FNONPRIMARY;
    uint16_t nonPrimeResB = b->core.flag & CReadBlock::BAM_FNONPRIMARY;
    //Compare primaryness
    if(nonPrimeResA != nonPrimeResB){
        //A is less if it has fewer NonPrime Flags
        return (nonPrimeResA < nonPrimeResB) ? -1 : 1;
    }
    //Compare mappedness
    uint16_t unmapResA = a->core.flag & CReadBlock::BAM_FANYUNMAP;
    uint16_t unmapResB = b->core.flag & CReadBlock::BAM_FANYUNMAP;
    if(unmapResA != unmapResB){
        //A is less if it has less Unmapped Flags
        return (unmapResA < unmapResB) ? -1 : 1;
    }
    //Compare Alignment Score
    uint8_t* asA = bam_aux_get(a,"AS");
    uint8_t* asB = bam_aux_get(b,"AS");
    if(asA || asB){ //At least one has an alignment score
	if(asA && asB){ //Both have alignment scores
	    int64_t scoreA = bam_aux2i(asA);
	    int64_t scoreB = bam_aux2i(asB);
	    if(scoreA != scoreB){
		//A is less if it has a higher score
		return (scoreA > scoreB) ? -1 : 1;
	    }
	} else { //Only one has an alignment score
	    //A is less if it is the only one with an alignment score
	    return (asA) ? -1 : 1;
	}
    }
    //Compare Segment order
    if((a->core.flag & BAM_FREAD1) != (b->core.flag & BAM_FREAD1)){
        //A is less if it is from read 1
        return (a->core.flag & BAM_FREAD1) ? -1 : 1;
    }
    //Compare Reference sequence
    if(a->core.tid != b->core.tid){
	//A is less if it is mapped to an earlier (by header) sequence
	return (a->core.tid < b->core.tid) ? -1 : 1;
    }
    //Compare Position in genome
    if(a->core.pos != b->core.tid){
	//A is less if it maps to an earlier position
	return (a->core.pos < b->core.pos) ? -1 : 1;
    }
    //Otherwise they are the same
    return 0;
}

//This is an effort to clean up the data a bit and deal witht eh potential
//issue of template switching
//For Cleaning:
//Filter non-primary XA strings that are double clipped (clipped on both
//sides)
//  if a primary XA string is double clipped, then the read block fails
//For Template Switching:
//  Any read which is clipped on one side which has an alternate
//  alignment to the same organism which is clipped on the other side
//  is considered a potential instance of template switching and should be
//  removed
//Input - A two element array of uniq pointers to a set of strings
//          representing the complete XA set
//      - A two element array of strings representing the primary
//          alignment string
//      - A set of virus names to differentiate host from virus
//Output - True if the readblock passes the filters
//       - False otherwise
//SideEffect: non-primary XA strings may be reomved from the provided set
//      if they fail 
bool CReadBlock::filterOnXASet( XAStrSetArr_t & xaStrSetArr,
                                const PrimaryXAStrArr_t & primeXAStrArr,
                                const VirusNameSet_t & virusNameSet)
{
    //Iterate over segments
    for(int segIdx = 0; segIdx < 2; segIdx++){
        const std::string & primeXAStr =primeXAStrArr[segIdx];
        CXA primeXAObj(primeXAStr);
        uint8_t primeClipFlag = primeXAObj.clipSide();
        bool primeIsHost = virusNameSet.count(primeXAObj.chr) == 0;
        fprintf(stderr,"PRIME: %s -> %u %d\n",primeXAStr.c_str(),primeClipFlag,primeIsHost);
        //If the primary alignment is double clipped it fails
        if(primeClipFlag == 0b11)
            return false;
        XAStrSet_pt & xaStrSet_p = xaStrSetArr[segIdx];
        for(auto it = xaStrSet_p->begin(); it != xaStrSet_p->end(); ){
            const std::string & xaStr = *it;
            if(xaStr == primeXAStr){ //Primary already processed
                it++;
                continue;
            }
            CXA xaObj(xaStr);
            uint8_t clipFlag = xaObj.clipSide();
            //Check for Double Clip
            if(clipFlag == 0b11){
                //Remove this xastr from the set
                it = xaStrSet_p->erase(it);
                continue;
            }
            //Need to check for Possible Template Switch
            if(primeClipFlag){ //At this point either left or right clipped
                bool isHost = virusNameSet.count(xaObj.chr) == 0;
                fprintf(stderr,"\t%s -> %u %d\n",xaStr.c_str(),clipFlag,isHost);
                //Internal Clip:
                // Prime clip and this clip are on opposite sides
                // Maps to same organism (isHost status is the same)
                if(primeClipFlag == (~clipFlag&3) && primeIsHost == isHost){
                    return false;
                }
            }
            it++;
        }
    }
    return true;
}

std::string CReadBlock::primaryBam2xaStr(bam1_t* read, bam_hdr_t* hdr){
    if(read->core.flag & BAM_FUNMAP)
	return "";
    std::string rname = sam_hdr_tid2name(hdr,read->core.tid);
    char strand = (read->core.flag & BAM_FREVERSE) ? '-' : '+';
    std::string cigar = get_cigar_code(read);
    uint8_t* nm = bam_aux_get(read,"NM");
    int editDist = (!nm) ? 0 : bam_aux2i(nm);
    std::string str = rname + ',' + strand + std::to_string(read->core.pos+1) + ',' +
	    cigar + ',' + std::to_string(editDist);
    return str;
}

//Processes the contents of the read block and collapsed down to a vector
//of only one entry each for the forward and reverse valid segments
//  Adds X2 tag: number of optimal and suboptimal alignments
//Input - A Header for tid to rname conversion purposes
//Output - true if a valid segment pair is found
//	 - false if not
bool CReadBlock::process(bam_hdr_t* hdr, const VirusNameSet_t & virusNameSet) {

    //Sort all reads in the block so the 'best' read will appear first
    std::sort(	this->m_Buf.begin(),this->m_Buf.end(),
		[] (bam1_t* & a, bam1_t* & b) {
		    return (CReadBlock::compareBams(a,b) == -1);
		});
    bam1_t* & read = this->m_Buf.front();
    //Clear the block and return if the best read is non primary
    if(read->core.flag & CReadBlock::BAM_FNONPRIMARY){
	//this->clear();
	return false;
    }
    if(read->core.flag & BAM_FUNMAP){
	//return if the best read is unmapped
	return false;
    }
    if(!selectBestMates()){
	//this->clear();
	return false;
    }
    //Load the pre-existing XA tag for the segments
    XAStrSetArr_t XAStrSet;
    PrimaryXAStrArr_t primaryXAStr;
    for(int i = 0; i <=1; i++){
	primaryXAStr[i] = CReadBlock::primaryBam2xaStr(this->m_Buf[i],hdr);
	XAStrSet[i].reset(new std::set<std::string>()); 
	//Insert this alignment into the XA set to prevent duplicates
	XAStrSet[i]->insert(primaryXAStr[i]);
	//Store the existing XA and SA tags
	CReadBlock::addBamAltstoXASet(this->m_Buf[i],*(XAStrSet[i]));
	//Delete the existing XA and SA tags
	for (char tag1 : {'S','X'}){
	    char tag[2] = {tag1,'A'};
	    uint8_t* pTag = bam_aux_get(this->m_Buf[i],tag);
	    if(!pTag) continue;
	    //Delete the existing XA string
	    bam_aux_del(this->m_Buf[i],pTag);
	}
    }
    //Collect all XA strings
    for(int i = this->m_Buf.size() - 1; i > 1; i--){
	bam1_t* & read = this->m_Buf[i];
	std::string xaStr = CReadBlock::primaryBam2xaStr(read,hdr);
	int segment = (read->core.flag & BAM_FREAD1) ? 0 : 1;
	XAStrSet[segment]->insert(xaStr);
	//Process the existing XA and SA tags
	CReadBlock::addBamAltstoXASet(read,*(XAStrSet[segment]));
	this->m_Buf.pop_back();
    }
    fprintf(stderr,"%s\n",this->name.c_str());
    //Filter The block on the basis of XA properties
    if(!CReadBlock::filterOnXASet(XAStrSet,primaryXAStr,virusNameSet)){
        return false;
    }
    for(int i = 0; i <=1; i++){
	int j = std::abs(i-1);
	bam1_t* & read = this->m_Buf[i];
	//Add an X2 tag (number of optimal and subobtimal alignments)
	uint8_t* pX2 = bam_aux_get(read,"X2");
	if(pX2) bam_aux_del(read,pX2);
	uint32_t x2 = XAStrSet[i]->size();
	bam_aux_update_int(read,"X2",x2);
	//Add Y2 tag (number of opt/subotp alignmtes in mate)
	uint8_t* pY2 = bam_aux_get(read,"Y2");
	if(pY2) bam_aux_del(read,pY2);
	uint32_t y2 = XAStrSet[j]->size();
	bam_aux_update_int(read,"Y2",y2);
	//Don't add the XA tag if only the decoy XAstr is present
	if(XAStrSet[i]->size() < 2) continue;
	std::string xaStr = "";
	for(const auto & str : *(XAStrSet[i])){
	    //Skip the decoy XAStr
	    if(str == primaryXAStr[i]) continue;
	    xaStr += str + ';';
	}
	bam_aux_update_str(read,"XA",xaStr.size(),xaStr.c_str());
    }
    return true;
}


bool CReadBlock::selectBestMates() {
    bam1_t* & read = this->m_Buf.front();
    //Identify the matching reads on the segments
    int seg1Idx = -1;
    int seg2Idx = -1;
    int* pSegIdx = nullptr;
    uint16_t oppSeg = 0;
    if(read->core.flag & BAM_FREAD1){ // Best is from segment 1
	seg1Idx = 0;
	pSegIdx = &seg2Idx;
	oppSeg = BAM_FREAD1;
    } else { // Best is from segment 2
	seg2Idx = 0;
	pSegIdx = &seg1Idx;
	oppSeg = BAM_FREAD2;
    }
    int32_t mtid = (read->core.tid & BAM_FMUNMAP) ? -1 : read->core.mtid;
    hts_pos_t mpos = read->core.mpos; 
    //Find matching segment 2 if it exists
    int firstIdx = -1;
    for(int i = 1; i < this->m_Buf.size() && *pSegIdx == -1; i++){
        bam1_t* & oppR = this->m_Buf[i];
        //Skip over other segment one entries
        if(oppR->core.flag & oppSeg) continue;
        //Record the first opposite segment
        if(firstIdx == -1) firstIdx = i;
        //Take the firts oppR if no target to find, otherwise go until
        // the target is found
        if(mtid == -1 || (mtid == oppR->core.tid && mpos == oppR->core.pos)) {
	   *pSegIdx=i;
	   continue;
        }
    }
    //If no exact match was found, but any mapping for the other segment
    //was found, use that
    if(*pSegIdx == -1) *pSegIdx = firstIdx;
    //No mapping found for the other segment
    if(seg1Idx == -1 || seg2Idx == -1) return false;
    //Reorder the buffer so that the first two entries are best segment1
    //and best segment2
    if(seg1Idx != 0){ // seg2Idx is at 0th potition
	if(seg1Idx != 1) { // move seg1 to 1th position if necessary
	    std::swap(this->m_Buf[seg1Idx],this->m_Buf[1]);
	}
	//swap seg1 and seg2
	std::swap(this->m_Buf[0],this->m_Buf[1]);
    } else if(seg2Idx != 1){
	// seg1Idx is at the 0th postion as desired
	// just move the seg2 to the 1th position if it isn't already
	// there
	std::swap(this->m_Buf[1],this->m_Buf[seg2Idx]);
    }
    //Set the mate tid and position so the segments match
    for(int i = 0, j = 1; i <= 1; i++, j--){
	bam1_t* & curRead = this->m_Buf[i];
	bam1_t* & curMate = this->m_Buf[j];
	//Test if the selected mate matches the nominal mate
	bool bDiff = false;
	bDiff = (curRead->core.mtid != curMate->core.tid ||
		 curRead->core.mpos != curMate->core.pos ||
		 (curRead->core.flag & BAM_FMREVERSE) != 
		    (curMate->core.flag & BAM_FREVERSE));
	if(!bDiff) continue;
	curRead->core.mtid = curMate->core.tid;
	curRead->core.mpos = curMate->core.pos;
	//Mask at desired bit
	uint16_t mateRevFlagMask = ~0 & BAM_FMREVERSE;
	//Unset the mate reverse flag
	curRead->core.flag &= ~mateRevFlagMask;
	//Set the bit if necessary
	if(curMate->core.flag & BAM_FREVERSE)
	    curRead->core.flag |= mateRevFlagMask;
	curRead->core.isize = 0;
    }
    //Unset any nonprimary flags
    for(int i = 0; i<=1; i++){
	bam1_t* & curRead = this->m_Buf[i];
	if(curRead->core.flag & CReadBlock::BAM_FNONPRIMARY){
	    uint16_t mask = ~0 & BAM_FNONPRIMARY;
	    curRead->core.flag &= ~mask;
	}
    }
    return true;
}

//=== FUNCTION DECLARATIONS

void OutputReadBlock(const CReadBlock & block, samFile * out, bam_hdr_t* hdr);
void PrintUsage();

//=== GLOBAL CONSTANTS - custom Types

///////const uint16_t BAM_FNONPRIMARY = BAM_FSECONDARY & BAM_FSUPPLEMENTARY;

//=== MAIN

int main(int argc, const char* argv[]) {
    if(argc < 4){ //Check for valid Input
	PrintUsage();
	return 1;
    }
    std::string viralFile(argv[1]);
    std::string inFile(argv[2]);
    std::string outFile(argv[3]);


    VirusNameSet_t vNameSet;
    LoadVirusNames(viralFile,vNameSet);

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
    //Read all reads from the file, process all blocks of reads with the
    //same read id
    bam1_t* curRead = bam_init1();
    std::unique_ptr<CReadBlock> pCurBlock(nullptr);
    std::unordered_set<std::string> observedQNames;
    while (sam_read1(in->file, in->header, curRead) >= 0) {
	std::string qname(bam_get_qname(curRead));
	//Handle the First Read
	if(!pCurBlock){ 
	    pCurBlock.reset(new CReadBlock(qname));
	}
	//If properly sorted, then each qname will appear in only 
	//  one block
	if(observedQNames.count(qname)){
	    std::cerr << "[ERROR] Input file is not sorted by name\n";
	    return 1;
	}
	///////Skip Secondary/Supplementary reads
	/////if(curRead->core.flag & BAM_FNONPRIMARY) continue;
	//Check if we have entered a new read block
	if(qname != pCurBlock->name){
	    //Process the Read Block
	    if(pCurBlock->process(in->header,vNameSet)) {
		//Output the remaining Reads in the Block
		OutputReadBlock(*pCurBlock,out,in->header); 
	    }
	    //Track the Blocks we've seen as a check for proper sorting
	    observedQNames.insert(pCurBlock->name);
	    //Replace the current block, destroying the old one
	    pCurBlock.reset(new CReadBlock(qname));
	}
	//Store a copy of the current read in the block
	pCurBlock->addRead(bam_dup1(curRead));
    }
    //Process the last read block
    if(pCurBlock && pCurBlock->process(in->header,vNameSet)) {
        //Output the remaining Reads in the Block
        OutputReadBlock(*pCurBlock,out,in->header); 
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
	    std::to_string(read->core.pos + 1);
}

//Iterate over reads in a read block and output to the bam output
//Input - a CReadBlock reference to output
//	- a samFile* to which to write
//	- a bam_hdr_t* for Chr conversion
//Output - None, writes to samFile
void OutputReadBlock(const CReadBlock & block, samFile * out, bam_hdr_t* hdr){
    for(auto it = block.cbegin(); it != block.cend(); it++){
	if(sam_write1(out,hdr,*it) < 0){
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
	<< "\tIf there are reads where the listed mate doesn't exist\n"
	<< "\t that mate alignment will be lost (Very rare edge case)\n"
        << "\tThis Does additional filters:\n"
        << "\t\tRemoves Reads where the best alignment is clipped on both sides\n"
        << "\t\tRemoves alternative alignments which are clipped on both sides\n"
        << "\t\tRemoves Reads which may be the result of template switching\n"
        << "\t\t\tThis is marked by having a clipped primary alignment and\n"
        << "\t\t\t an alt alignment mapping to the same organim which isclipped\n"
        << "\t\t\t on the opposite side: Ex. HostChrA:LClip and HostChrB:RClip\n"
	<< "\tAdds X2 and Y2 tag to entries indicating the total number of\n"
	<< "\t alignments (both primary and secondary) in the entry and its mate\n"
	<< "===Usage\n"
	<< "\tcleanBAM ViralGenome.fa in out\n"
	<< "===ARGUMENTS\n"
        << "\tViralGenome PATH\tA fasta formatted file which will be used to\n"
        << "\t differentiate host and viral reads\n"
	<< "\tin PATH\tPath to a [SB]AM file to process\n"
	<< "\t\tIt is expected that the input has MC tags from fixmate\n"
	<< "\tout PATH\tPath to an output file which will be BAM formatted\n"
	<< "\tNOTE: Use '-' to indicate that input/output are stdin/stdout\n"
	<< "===Output\n"
	<< "\tOutput is BAM formatted, with only primary alignments.\n"
	<< "\t Supplementary and Secondary alignments are added to XA tag\n"
	<< "\t of the primary alignment.\n"
	;
}
