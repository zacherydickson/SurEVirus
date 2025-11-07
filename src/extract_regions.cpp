#include <iostream>
#include <vector>
#include <array>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <stdexcept>

#include "sam_utils.h"
#include "config.h"
#include <cptl_stl.h>
#include "utils.h"

enum JunctionSide_t {JS_HOST, JS_VIRUS};
enum GoodClipType_t {
    GCT_R1L = 0, //good left clip from R1
    GCT_R1R = 1, //good right clip from R1
    GCT_R2L = 2, //good left clip from R2
    GCT_R2R = 3, //good right clip from R2
};

////Note Junctions fall into one of 8 situations:
////  0		1	2		3
////  (A) H+V+, (B) H+V-, (C) H-V+, and (D) H-V-
////  (a)	V-H-, (b) V+H-, (c) V-H+, and (d) V+H+ 
////  Note that the latter 4 are equivalent to the former 4
////  The reported junction strand will be from the former 4
////00 -> ++, 01 -> +-, 10 -> -+, 11 -> --
//uint8_t JunctionOrientation[32] = {
//    //	    Seq		Cli	Anc	aVir	Left
//    0b00, //   1	+	+	V	L
//    0b11, //   1 	+	+	V	R
//    0b11, //   1 	+	+	H	L
//    0b00, //   1 	+	+	H	R
//    0b11, //   1 	+	-	V	L
//    0b00, //   1 	+	-	V	R
//    0b00, //   1 	+	-	H	L
//    0b11, //   1 	+	-	H	R
//    0b10, //   1 	-	+	V	L
//    0b01, //   1 	-	+	V	R
//    0b10, //   1 	-	+	H	L
//    0b01, //   1 	-	+	H	R
//    0b01, //   1 	-	-	V	L
//    0b10, //   1 	-	-	V	R
//    0b10, //   1 	-	-	H	L
//    0b01, //   1 	-	-	H	R
//    0b11, //   2 	+	+	V	L
//    0b00, //   2 	+	+	V	R
//    0b00, //   2 	+	+	H	L
//    0b11, //   2 	+	+	H	R
//    0b00, //   2 	+	-	V	L
//    0b11, //   2 	+	-	V	R
//    0b11, //   2 	+	-	H	L
//    0b00, //   2 	+	-	H	R
//    0b01, //   2 	-	+	V	L
//    0b10, //   2 	-	+	V	R
//    0b10, //   2 	-	+	H	L
//    0b01, //   2 	-	+	H	R
//    0b10, //   2 	-	-	V	L
//    0b01, //   2 	-	-	V	R
//    0b01, //   2 	-	-	H	L
//    0b10  //   2 	-	-	H	R
//};

//Same 4(8) cases as previously but the information is more constrained
//only 3 tests needed
uint8_t PairedJunctionOrientation[8] = {
    //R1 Host
    //	R1 Fwd
    0b01, //R2 Fwd
    0b00, //R2 Rev
    // R1 Rev
    0b11, //R2 Fwd
    0b10, //R2 Rev
    //R1 Virus
    //	R1 Fwd
    0b01, //R2 Fwd
    0b11, //R2 Rev
    //	R1 Rev
    0b00, //R2 Fwd
    0b10, //R2 Rev
};

typedef std::array<bam1_t*,4> ClipArray_t;
typedef std::unordered_map<std::string, ClipArray_t> GoodClipMap_t;
typedef std::unordered_set<std::string> QNameSet_t;

std::unordered_set<std::string> VirusNameSet;
GoodClipMap_t GoodClipMap;
QNameSet_t GoodClipSet;
//std::mutex mtx;

//===== Function Declarations

std::string ConstructCandidateString();
void DestroyGoodClips();
//char DetermineClipJunctionStrand (uint8_t flag);
std::array<char,2> DetermineJunctionOrientation (   bool bViralAnchor,
		    bool isLeftClip, bool bAnchorRev, bool bClipRev, bool isR1);
std::array<char,2> DeterminePairedJunctionOrientation(bool r1Virus, bool r1Rev,
		    bool r2Rev);
bool IsLeftOfJunction(uint8_t flag);
void LoadAnchorOrientation(std::string fname);
void LoadGoodClips(std::string fname);
//void OutputBEDEntries(	std::ofstream & outbed, const bam1_t* read,
//			std::string cname, uint8_t clipflag = 0x0,
//			std::string qname = std::string());
void ParseReadXA(bam1_t *read, std::string primaryContig,std::vector<CXA> & out);
//void ProcessAnchor(bam1_t *read, std::string cname, std::ofstream & outbed);
//void ProcessClip(bam1_t *read, std::string cname, std::ofstream & outbed);
void ProcessPair(   bam1_t *r1, bam1_t *r2, std::string cname1,
		    std::string cname2, std::ofstream &outbed);
void ProcessPairs(std::string fname, std::ofstream & outbed);
void ProcessSplitRead(	bam1_t *anchor, bam1_t clip, int jSide,
			std::string primaryContig, std::string clipCName,
			std::ofstream & outbed);
void ProcessSplitReads(	std::string anchor_fname, std::string clip_fname, 
			int jSide, std::ofstream & outbed);

// ===== MAIN

//Process BAM files from the working environment to extract all
//candidate junctions based on BWA alignment
//Inputs - A path to the viral reference in fasta format
//	 - A path to the working directory
//	 - A path to the bam workspace
//Outputs - A bed6 formatted file detailing potential junction positions
//Each entry describes a side of a junction, whether host or virus can
//be determined by checing the contig id in column 1
//The indexes describe the interval in which the junction could be
//  For split reads this is a single position (Defined by the split)
//  For chimeric reads this ranges up to the maxIS from the junction
//  side of the mapped portion of the read (this is limited by 
//The name column is id of the fragment supporting the candidate
//  a single fragment can support multiple candidates
//The score column is unused 
//The strand column
int main(int argc, char* argv[]) {
    //##PARSE INPUTS
    std::string virus_names_file = argv[1];
    std::string workdir = argv[2];
    std::string workspace = argv[3];

    //##Files to be used from the workspace
    std::string bam_fname = workspace + "/retained-pairs.namesorted.bam";
    std::string clip_bam_fnames[2];
    //Host clips = virus-side, virus clips = host side
    clip_bam_fnames[JS_HOST] = workspace + "/host-clips.cs.bam";
    clip_bam_fnames[JS_VIRUS] = workspace + "/virus-clips.cs.bam";
    //Host anchors = host-side, virus anchors = virus side
    std::string anchor_bam_fnames[2];
    anchor_bam_fnames[JS_HOST] = workspace + "/host-anchors.bam";
    anchor_bam_fnames[JS_VIRUS] = workspace + "/virus-anchors.bam";

    //##Output file
    std::string bed_fname = workdir + "/junction-candidates.bed";
    std::ofstream outbed(bed_fname);

    /*//Set up the thread pool
    int nThread = parse_config_threads(workdir + "/config.txt");
    ctpl::thread_pool thread_pool(nThread); */

    //##LOAD DATA INTO GLOBAL VARIABLES
    //Load names of viral contigs
    LoadVirusNames(virus_names_file,VirusNameSet);
    //Load the ids and directions of clips which map properly
    for (int side = JS_HOST; side <= JS_VIRUS; side++){
        LoadGoodClips(clip_bam_fnames[side]);
        ProcessSplitReads(  anchor_bam_fnames[side],clip_bam_fnames[side],
			    side,outbed);
        DestroyGoodClips();
    }

    //Pass over the paired reads to find valid chimeras
    ProcessPairs(bam_fname,outbed);
}

//===== Function Defintions

//Given information to be printed, constructs a string describing the candidate breakpoint
//Inputs -
//Output - a string
std::string ConstructCandidateString(   std::string chr, size_t pos,
                                        std::string qname, char strand)
{
    std::string str =  chr + '\t' + std::to_string(pos) + '\t' +
	                std::to_string(pos+ 1) + '\t' +
                        qname + "\t.\t" + strand;
    return str;
}

//Looks up the case in a precalculated table based on 5 boolean values
//This table gives information on if the human and viral sides are in the
//+ or - orientation
//The result can then be matched to the anchor and clip as needed
//Inputs - 5 boolean values defining the orientation
//Output - a string 
std::array<char,2> DetermineJunctionOrientation (  bool bViralAnchor, bool isLeftClip,
					    bool bAnchorRev, bool bClipRev,
					    bool isR1) {
    char anchorStrand, clipStrand;
    ////Result is two bits where the one's bit is one if the virus is rev
    //// and the twos bit is one if the host is rev
    uint8_t result = 0;
    bool bRightClip = !isLeftClip;
    bool bIs5Prime = (bViralAnchor != bRightClip);
    if(!bIs5Prime) result |= 0b10; //Host is Inverted at the 3' junction
    //virus matches host in fwd insertions, and doesn't in reverse insertions
    result |= ((result >> 1) ^ bClipRev); // XOR can act as a negate
    //If the anchor and clip strand do not match ( 0b01=1 or 0b10=2 )
    //And the virus is the anchor, the strands need to be flipped
    if(bClipRev && bViralAnchor){
	result = (~result) & 0b11;
    }
    if(bViralAnchor){ // Anchor takes viral result
        anchorStrand = (result & 0b01) ? '-' : '+';
        clipStrand = (result & 0b10) ? '-' : '+';
    } else {
        anchorStrand = (result & 0b10) ? '-' : '+';
        clipStrand = (result & 0b01) ? '-' : '+';
    }
    return {anchorStrand,clipStrand};
}


std::array<char,2> DeterminePairedJunctionOrientation(bool r1Virus, bool r1Rev,
		    bool r2Rev) {
    char hostStrand, virStrand;
    uint8_t flag = 0;
    if(r1Virus) flag |= 0x4;
    if(r1Rev) flag |= 0x2;
    if(r2Rev) flag |= 0x1;
    uint8_t result = PairedJunctionOrientation[flag];
    hostStrand = (result & 0x2) ? '-' : '+';
    virStrand = (result & 0x1) ? '-' : '+';
    return {hostStrand,virStrand};
}

//Opens a bam file containing mapped clips.
//It is assumed that the bam file only contains primary mappped clips
//(no secondary/supplementary/unmapped)
//All reads in these files are assumed to define good clips
//Each alignment object is stored for later use (the entire clips file is
//loaded into memory)
//Input  - a string reperesenting a file name
//Output - none, modifies the global GoodClipMap
void LoadGoodClips(std::string fname){
    open_samFile_t* clips_file = open_samFile(fname.c_str(), false, false);
    bam1_t* read = bam_init1();

    while (sam_read1(clips_file->file, clips_file->header, read) >= 0) {
        std::string clip_name = bam_get_qname(read);
        std::string qname = clip_name.substr(0, clip_name.length()-4);
	GoodClipSet.insert(qname);
	if(!GoodClipMap.count(qname)){
	    GoodClipMap[qname] = {nullptr,nullptr,nullptr,nullptr};
	}
	bool isLeftClip(clip_name[clip_name.length()-3] == 'L');
	bool isR1(clip_name[clip_name.length()-1] == '1');
	int idx = 0;
	if(!isLeftClip) idx += 1;
	if(!isR1) idx += 2;
	GoodClipMap[qname][idx]	= bam_dup1(read);
    }
    close_samFile(clips_file);
    bam_destroy1(read);
}

void DestroyGoodClips(){
    for( auto & pair : GoodClipMap ){
	for( int i = 0; i < 4; i++){
	    if(pair.second[i]){
		bam_destroy1(pair.second[i]);
		pair.second[i] = nullptr;
	    }
	}
    }
    GoodClipMap.clear();
}

void ParseReadXA (  bam1_t *read, std::string primaryContig,
		    std::vector<CXA> & out){
    uint8_t * nm = bam_aux_get(read,"NM"); 
    int nmVal = (nm) ? bam_aux2i(nm) : 0;

    std::string xaStr = primaryContig + "," + 
			((read->core.flag & BAM_FREVERSE) ? "-" : "+")  +
			std::to_string(read->core.pos) + "," +
			"1M" + "," + std::to_string(nmVal);
    out.emplace_back(xaStr);
    out.front().nCigar = read->core.n_cigar;
    out.front().cigar = (uint32_t*) std::realloc(out.front().cigar,
					    sizeof(uint32_t) * (out.front().nCigar));
    memcpy( out.front().cigar,bam_get_cigar(read),
	    sizeof(uint32_t) * read->core.n_cigar);


    uint8_t * xa = bam_aux_get(read,"XA");
    if(!xa) return;
    std::string xaListStr = bam_aux2Z(xa);
    size_t pos = 0, prev = 0;
    while((pos = xaListStr.find(';',prev+1)) != std::string::npos){
	xaStr = xaListStr.substr(prev,pos-prev);
	prev=pos+1;
	out.emplace_back(xaStr);
    }
}

void ProcessPair(   bam1_t *r1, bam1_t *r2, std::string cname1,
		    std::string cname2, std::ofstream & outbed){
    std::string qname = bam_get_qname(r1);
    //Clipped Reads get priority
    if(GoodClipSet.count(qname)) return;
    //Small optimization to quickly exclude non-chimeric reads
    if(r1->core.tid == r2->core.tid) return;
    bool r1IsVirus = (VirusNameSet.count(cname1));
    bool r2IsVirus = (VirusNameSet.count(cname2));
    //Again can quickly skip non-chimerics (Double host or double virus)
    if(r1IsVirus == r2IsVirus) return;
    //Exclude Homopolymers
    if(is_poly_ACGT(r1) || is_poly_ACGT(r2)) return;



    std::unordered_set<std::string> potentialEntries;

    std::vector<CXA> r1Mappings;
    std::vector<CXA> r2Mappings;
    ParseReadXA(r1,cname1,r1Mappings);
    ParseReadXA(r2,cname2,r2Mappings);

    for(int i = 0; i < r1Mappings.size(); i++){
	const CXA & r1Map = r1Mappings[i];
	for(int j = 0; j < r2Mappings.size(); j++){
	    const CXA & r2Map = r2Mappings[j];
	    //All alignments for a segment mapping to host must be to host
	    //All alignments for a segment mapping to virus must be to virus
	    if(r1IsVirus != bool(VirusNameSet.count(r1Map.chr))) return;
	    if(r2IsVirus != bool(VirusNameSet.count(r2Map.chr))) return;
	    //Construct Entries
	    std::string hostChr = r1Map.chr;
	    std::string virChr = r2Map.chr;
	    hts_pos_t hostPos = r1Map.endpos();
	    hts_pos_t virPos = r2Map.pos;
	    if(r1IsVirus){ // r2 is Host Side
		std::swap(hostChr,virChr);
		hostPos = r2Map.endpos();
		virPos = r1Map.pos;
	    }
	    std::array<char,2> strands = DeterminePairedJunctionOrientation(
		    r1IsVirus,r1Map.bRev,r2Map.bRev);
	    //Store the unique entries
            potentialEntries.insert(
                    ConstructCandidateString(   hostChr,hostPos,qname,
                                                strands.front()));
            potentialEntries.insert(
                    ConstructCandidateString(   virChr,virPos,qname,
                                                strands.back()));
	}
    }

    //The pair passed all filters, can output the entries now
    for(auto entry : potentialEntries){
	outbed << entry << "\n";
    }

}

//Opens a given bam file and outputs all of the host/virus side junctions
//each properly mapped pair of reads suppports, only read pairs which
//only map in chimeric configurations are accepted
//Inputs - a path to a mapped clip file
//	 - an ofstream object to which to write
//	 - Also uses the Global Good Clips Set for filtering
//Output - None, writes to outbed
void ProcessPairs(std::string fname, std::ofstream & outbed){
    open_samFile_t* bam_file = open_samFile(fname.c_str(), false, false);
    bam1_t* read1 = bam_init1();
    bam1_t* read2 = bam_init1();


    while (sam_read1(bam_file->file, bam_file->header, read1) >= 0 &&
	   sam_read1(bam_file->file, bam_file->header, read2) >= 0 ) {
	if(read1->core.flag & (BAM_FUNMAP | BAM_FMUNMAP)) continue;
	std::string qname = bam_get_qname(read1);
	if(qname != bam_get_qname(read2)){
	    throw std::invalid_argument("Mates not adjacent in paired reads bam");
	}
	std::string cname1 = sam_hdr_tid2name(bam_file->header,read1->core.tid);
	std::string cname2 = sam_hdr_tid2name(bam_file->header,read2->core.tid);
	//Provide the pair in segment 1 segment 2 order
	if(read1->core.flag & BAM_FREAD1){
	    ProcessPair(read1,read2,cname1,cname2,outbed);
	} else {
	    ProcessPair(read2,read1,cname2,cname1,outbed);
	}
    }


    close_samFile(bam_file);
    bam_destroy1(read1);
    bam_destroy1(read2);
}

void ProcessSplitRead(	bam1_t *anchor, bam1_t *clip, int jSide, 
			std::string primaryContig, std::string clipCName,
			int lrIdx, std::ofstream & outbed){

    bool bViralAnchor = (jSide == JS_VIRUS);
    bool isLeftClip = (lrIdx % 2 == 0);
    bool isR1 = (lrIdx < 2);
    //Load Clip Alts
    std::vector<CXA> anchorMappings;
    std::vector<CXA> clipMappings;
    ParseReadXA(anchor,primaryContig,anchorMappings);
    ParseReadXA(clip,clipCName,clipMappings);
    
    //Process Primary
    
    std::string qname = bam_get_qname(anchor);
    qname += (isR1) ? "_1" : "_2";


    //Iterate over all pairs of anchor and clips
    //And note all uniq breakpoints this read supports
    std::unordered_set<std::string> uniqBPStrSet;
    for(int i = 0; i < anchorMappings.size(); i++){
	const CXA & anchorMap = anchorMappings[i];
        bool isViralAnchorXA = (VirusNameSet.count(anchorMap.chr));
        //Skip alt anchors which don't have the same virus status as the primary
        // anchor alignment
        if(bViralAnchor != isViralAnchorXA){
            continue;
        }
	for(int j = 0; j < clipMappings.size(); j++){
	    const CXA & clipMap = clipMappings[j];
	    hts_pos_t anchorPos = (isLeftClip) ? anchorMap.pos : anchorMap.endpos();
	    hts_pos_t clipPos = (isLeftClip) ? clipMap.endpos() : clipMap.pos;
	    std::array<char,2> strands = DetermineJunctionOrientation(bViralAnchor,
				    isLeftClip,anchorMap.bRev,clipMap.bRev,isR1);
            uniqBPStrSet.insert(ConstructCandidateString(anchorMap.chr,
                                                         anchorPos,
                                                         qname,
                                                         strands.front()));
            uniqBPStrSet.insert(ConstructCandidateString(clipMap.chr,
                                                         clipPos,
                                                         qname,
                                                         strands.back()));
	}
    }
    //Ouput each breakpoint
    for(const std::string & bpStr : uniqBPStrSet){
        outbed << bpStr << "\n";
    }
}

//Opens a given bam file and outputs all of the host/virus side junctions
//each clip or anchor supports
//Inputs - paths to both anchor and clip files the latter is only used for
//	    its header
//	 - a boolean indicating wether a clip or anchor file has been provided
//	 - an ofstream object to which to write
//Output - None, writes to outbed
void ProcessSplitReads(	std::string anchor_fname, std::string clip_fname,
			int jSide, std::ofstream & outbed){
    open_samFile_t* anchors_file = open_samFile(anchor_fname.c_str(),false,false);
    //Only needed for its header
    open_samFile_t* clips_file = open_samFile(clip_fname.c_str(),false,false);
    bam1_t* read = bam_init1();

    while (sam_read1(anchors_file->file, anchors_file->header, read) >= 0) {
        std::string qname = bam_get_qname(read);
	std::string cname = sam_hdr_tid2name(	anchors_file->header,
						read->core.tid);
	//Make sure there are any good clips for this anchor
	if(!GoodClipMap.count(qname)) continue;
	bool isR1 = (read->core.flag & BAM_FREAD1);
	int leftIdx = (isR1) ? GCT_R1L : GCT_R2L;
	int rightIdx = (isR1) ? GCT_R1R : GCT_R2R;
	//Iterate over both left and right clips
	for(int idx = leftIdx; idx <= rightIdx; idx++){
	    //Skip there isn't a matching read
	    if(!GoodClipMap[qname][idx]) continue;
	    //No good clip exists for this segment on this side;
	    std::string clip_cname = sam_hdr_tid2name(	clips_file->header,
						GoodClipMap[qname][idx]->core.tid);
	    ProcessSplitRead(	read,GoodClipMap[qname][idx],jSide,cname,
				clip_cname,idx,outbed);
	}
    }
    close_samFile(anchors_file);
    close_samFile(clips_file);
    bam_destroy1(read);
}

