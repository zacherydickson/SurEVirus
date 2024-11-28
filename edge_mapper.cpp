#include <array>
#include <iostream>
#include <htslib/sam.h>
#include <memory>
#include <unordered_map>
#include <unordered_set>

#include "config.h"
#include "utils.h"
#include "libs/ssw.h"
#include "libs/ssw_cpp.h"
#include "libs/cptl_stl.h"

//==== TYPE DECLARATIONS

//Object describing a read
struct Read_t {
    Read_t(std::string nm, bool bSplit, bool bViralR1) :
	name(nm), isSplit(bSplit), viralR1(bViralR1) {}
    std::string name;
    std::string hostSegment;
    std::string hostRC;
    std::string virusSegment;
    std::string virusRC;
    bool isSplit = false;
    bool viralR1 = false;
    int compare(const Read_t other) const {
	if(this->name != other.name) {
	    return (this->name < other.name) ? -1 : 1;
	}
	if(this->hostSegment.size() != other.hostSegment.size()){
	    return (this->hostSegment.size() < other.hostSegment.size())
		    ? -1 : 1;
	}
	if(this->virusSegment.size() != other.virusSegment.size()){
	    return (this->virusSegment.size() < other.virusSegment.size())
		    ? -1 : 1;
	}
	if(this->isSplit != other.isSplit){
	    return (other.isSplit) ? -1 : 1;
	}
	if(this->viralR1 != other.viralR1){
	    return (this->viralR1) ? -1 : 1;
	}
	return 0;
    }
    std::string to_string(bool seq = false) const {
	std::string str;
	if(!seq){
	    str =   name + '-' + std::to_string(hostSegment.size()) + '-' + 
		    std::to_string(virusSegment.size()) + '-' +
		    std::to_string(isSplit) + std::to_string(viralR1);
	} else {
	    if(hostSegment.size()){
		str += std::to_string('>') + name + '_' + std::to_string(1);
		str += std::to_string('\n') + hostSegment;
	    }
	    if(virusSegment.size()){
		str += std::to_string('>') + name + '_' + std::to_string(2);
		str += std::to_string('\n') + virusSegment;
	    }
	}
	return str;
    }
    const std::string & getSegment(bool bVirus, bool bRC) const {
	if(bVirus){
	   return (bRC) ? this->virusRC : this->virusSegment;
        } else {
	   return (bRC) ? this->hostRC : this->hostSegment;
        }
    }
};

typedef std::shared_ptr<Read_t> Read_pt;

struct Read_pt_EqFunctor {
    bool operator()(const Read_pt & a, const Read_pt & b) const {
	return a->compare(*b) == 0;
    }
};

struct Read_pt_HashFunctor {
    size_t operator()(const Read_pt & a) const {
	return std::hash<std::string>{}(a->to_string());
    }
};

typedef std::unordered_map<std::string,Read_pt> Name2ReadMap_t;

//Object describing a genomic location containing a
//candidate junction
struct Region_t {
    Region_t(kseq_t *seq, bool bVirus) : sequence(seq->seq.s), isVirus(bVirus){
	std::string name = seq->name.s;
	std::vector<std::string> fields = strsplit(name,',');
	this->chr = fields[0];
	this->left = std::stoul(fields[1]);
	this->right = std::stoul(fields[2]);
	this->strand = fields[3][0];
	for (auto & c: sequence) c = (char)toupper(c);
    }
    Region_t(const Region_t & other) :
	left(other.left), right(other.right), chr(other.chr), 
	strand(other.strand), sequence(other.sequence){}
    size_t left, right;
    std::string chr;
    char strand;
    std::string sequence;
    bool isVirus;
    int compare(const Region_t other){
	if(this->chr != other.chr){ 
	    return (this->chr < other.chr) ? -1 : 1;
	}
	if(this->strand != other.strand){
	    return (this->strand == '-') ? -1 : 1;
	}
	if(this->left != other.left){
	    return (this->left < other.left) ? -1 : 1;
	}
	if(this->right != other.right) {
	    return (this->right < other.right) ? -1 : 1;
	}
	return 0;
    }
    std::string to_string(bool seq=false) const {
	std::string str=    '>' + chr + ":" + strand +
			    std::to_string(left) + "-" +
			    std::to_string(right);
	if(seq) str += sequence;
	return str;
    }
};

typedef std::shared_ptr<Region_t> Region_pt;
typedef std::unordered_map<std::string,Region_pt> Name2RegionMap_t;

struct Region_pt_EqFunctor {
    bool operator()(const Region_pt & a, const Region_pt & b) const{
	return a->compare(*b) == 0;
    }
};

struct Region_pt_HashFunctor {
    size_t operator()(const Region_pt & a) const {
	return std::hash<std::string>{}(a->to_string());
    }
};



//Sets of shared pointers to reads and regions
typedef std::unordered_set<Read_pt,Read_pt_HashFunctor,Read_pt_EqFunctor>
	    ReadSet_t;
typedef std::unordered_set<Region_pt,Region_pt_HashFunctor,Region_pt_EqFunctor>
	    RegionSet_t;


typedef std::unordered_map< Read_pt,RegionSet_t,
			    Read_pt_HashFunctor,Read_pt_EqFunctor>
	    Read2RegionsMap_t;
typedef std::unordered_map< Region_pt, ReadSet_t,
			    Region_pt_HashFunctor,Region_pt_EqFunctor>
	    Region2ReadsMap_t;
typedef std::unordered_map< Read_pt, ReadSet_t,
			    Read_pt_HashFunctor,Read_pt_EqFunctor>
	    Read2ReadsMap_t;

//Object associating a pair of regions and the reads spanning the pair
struct Edge_t {
    Region_pt hostRegion;
    Region_pt virusRegion;
    ReadSet_t readSet;
    ReadSet_t uniqueReadSet;
    //Read2ReadsMap_t duplicatedReads;
    size_t hostOffset;
    size_t virusOffset;
    size_t nSplit = 0;
    Edge_t() : hostRegion(nullptr), virusRegion(nullptr), readSet() {}
    //Edge_t(const std::string & regStr, const std::string & readStr);
    Edge_t(Region_pt hostReg, Region_pt virusReg) :
	hostRegion(hostReg), virusRegion(virusReg), readSet() {}
    Edge_t(const Edge_t & other) :
	hostRegion(other.hostRegion), virusRegion(other.virusRegion),
	readSet(other.readSet), uniqueReadSet(other.uniqueReadSet),
	hostOffset(other.hostOffset), virusOffset(other.virusOffset), 
	nSplit(other.nSplit) {}
    public:
    bool addRead(const Read_pt & read){
	auto res = this->readSet.insert(read);
	if(res.second){
	    if(read->isSplit) nSplit++;
	    return true;
	}
	return false;
    }
    bool removeRead(const Read_pt & read){
	if(!this->readSet.erase(read)) return false;
	if(read->isSplit && nSplit) nSplit--;
	return true;
    }
    private:
    //void parseRegString(const std::string & regStr);
    //void parseReadString(const std::string & readStr);
};

typedef std::vector<Edge_t> EdgeVec_t;

//Object associating a query read with a subject region
struct SQPair_t {
    SQPair_t(const Region_pt & s, const Read_pt & q) : subject(s), query(q) {}
    Region_pt subject;
    Read_pt query;
    int compare(const SQPair_t & other) const {
	int res = this->subject->compare(*other.subject);
	if(res != 0) return res; 
	res = this->query->compare(*other.query);
	if(res != 0) return res;
	return 0;
    }
    std::string to_string() const {
	return	this->subject->to_string() + " vs " +
		this->query->to_string();
    }
};

struct SQPair_EqFunctor {
    bool operator()(const SQPair_t & a, const SQPair_t & b) const {
	return a.compare(b) == 0;
    }
};

struct SQPair_HashFunctor {
    size_t operator()(const SQPair_t &a) const {
	return std::hash<std::string>{}(a.to_string());
    }
};

//Mapping from a subject query pair to the resulting alignment
typedef std::unordered_map< SQPair_t,StripedSmithWaterman::Alignment,
			    SQPair_HashFunctor,SQPair_EqFunctor>
	    AlignmentMap_t;

typedef std::pair<size_t,size_t> PosPair_t;
struct PosPair_HashFunctor {
    size_t operator()(const PosPair_t & a) const {
	return std::hash<size_t>{}(a.first ^ a.second);
    }
};
typedef std::unordered_set<PosPair_t,PosPair_HashFunctor> PosPairSet_t;

struct JunctionInterval_t {
    size_t proximal, distal;
};

//==== GLOBAL VARIABLE DECLARATIONS

std::unordered_set<std::string> VirusNameSet;
config_t Config;
stats_t Stats;
bam_hdr_t* JointHeader;
//For an alignment to pass it must have a score of at least 30
StripedSmithWaterman::Filter AlnFilter;//(true,true,30,32767);
StripedSmithWaterman::Aligner Aligner(1,4,6,1,false);
int32_t AlnMaskLen;

static size_t MinimumReads = 4;
static size_t SplitBonus = 1;
static double MaxDiffRate = 0.06;

std::mutex Mtx;

//==== FUNCTION DECLARATIONS

void AlignRead(	int id, const Read_pt & read, const RegionSet_t & regSet,
		AlignmentMap_t & alnMap);
void AlignReads(const Read2RegionsMap_t &regMap,
		AlignmentMap_t & alnMap);
EdgeVec_t ConsensusSplitEdge(	int id, Edge_t & edge,
				const AlignmentMap_t & alnMap);
bool ConstructBamEntry(	const Read_pt & query, const Region_pt & subject,
			const Region_pt & mateSubject,
			const AlignmentMap_t & alnMap,
			bam1_t* entry);
breakpoint_t ConstructBreakpoint(const Region_pt & reg,size_t offset);
call_t ConstructCall(	int id, const Edge_t & edge,
			const AlignmentMap_t & alnMap);
JunctionInterval_t ConstructJIV(char strand, bool isVirus,
				const StripedSmithWaterman::Alignment & aln);
void DeduplicateEdge(Edge_t & edge,const AlignmentMap_t & alnMap);
size_t FillStringFromAlignment(	std::string & outseq,
				const std::string & inseq,
				size_t offset, size_t maxLen,
				const std::vector<uint32_t> & cigarVec);
void FilterEdgeVec(EdgeVec_t & edgeVec, const ReadSet_t * usedReads = nullptr);
void FilterHighInsertReads(Edge_t & edge, const AlignmentMap_t & alnMap);
void FilterSuspiciousReads(Edge_t & edge, const AlignmentMap_t & alnMap);
template<class T>
void FilterVector(  std::vector<T> & vec,
		    const std::unordered_set<size_t> & idxSet);
std::string GenerateConsensus(const std::vector<std::string> & rowVec,
				std::vector<size_t> & diffVec);
std::string GetAlignedSequence(	const Edge_t & edge, const Read_pt & read,
				const AlignmentMap_t & alnMap, size_t & nFill);
double GetEdgeScore(const Edge_t & edge,const AlignmentMap_t & alnMap,
		    const ReadSet_t & usedReads);
void IdentifyEdgeBreakpoints(Edge_t & edge, const AlignmentMap_t & alnMap);
void IdentifyEdgeSpecificReads(EdgeVec_t & edgeVec);
bool IsConsistent(const std::string & seq1, const std::string & seq2);
void LoadData(	const std::string & edgeFName,
		const std::string & regionsFName,
		const std::string & readsFName,
	      	Read2RegionsMap_t & read2regSetMap,
	      	Region2ReadsMap_t & reg2readSetMap,
	      	EdgeVec_t & edgeVec);
void LoadEdges(	std::string edgeFName,
	      	Read2RegionsMap_t & read2regSetMap, 
	      	Region2ReadsMap_t & reg2readSetMap,
	      	EdgeVec_t & edgeVec,
		const Name2ReadMap_t & readNameMap,
		const Name2RegionMap_t & regionNameMap);
void LoadReadSeq(   const std::string & readsFName,
		    Read2RegionsMap_t & read2regSetMap,
		    Name2ReadMap_t & nameMap);
void LoadRegionSeq( const std::string & regionsFName,
		    Region2ReadsMap_t & reg2readSetMap,
		    Name2RegionMap_t & nameMap);
void OutputEdge(const Edge_t & edge, const AlignmentMap_t & alnMap,
		const std::string & resFName, const std::string & readDir,
		ReadSet_t usedReads);
void OrderEdges(EdgeVec_t & edgeVec,const AlignmentMap_t & alnMap);
void OutputEdges(EdgeVec_t & edgeVec,const AlignmentMap_t & alnMap,
		 const std::string & resFName, const std::string & readDir);
bool PassesEffectiveReadCount(	const Edge_t & edge,
				const ReadSet_t * usedReads = nullptr);
void ProcessEdge(int id,Edge_t & edge, const AlignmentMap_t & alnMap);
void ProcessEdges(EdgeVec_t & edgeVec, const AlignmentMap_t & alnMap);
void RecursiveSplitEdge(Edge_t & edge, std::vector<Read_pt> & rowLabelVec,
			std::vector<std::string> & rowSeqVec,
			std::vector<size_t> & nFillVec,
			EdgeVec_t & newEdges);
void RemoveUnalignedReads(EdgeVec_t & edgeVec,const AlignmentMap_t & alnMap);
EdgeVec_t SplitEdges(EdgeVec_t & edgeVec, const AlignmentMap_t & alnMap);

//==== MAIN

//Parse Fasta Files for regions and reads as well as a file defining
//edges
//Regions sequences are expected to be the same strand as labelled
//Region names should be in the format CONTIG,OFF,END,STRAND(:+) 
//Reads are always on the reference strand regardless of fwd/rev
//Read names should be in the format READID_[12]
//
//After loading each read is mapped to all associated regions with
//poor alignments filtered away
//
//Edges are then processed, before each step checking if enough reads
//remain
//  Deduplication (based on alignment positions and cigars)
//  Split Edges By Consensus Sequence (Ex. R1->S1 R2->S2, R1->S1) -> 2 edges
//  Sum the score for each edge (conditional(from shared reads), and
//	unconditional (score from edge specific reads)
//  Sort edges on unconditional score, then by conditional score
//  Accept from best to worst until none pass anymore
int main(int argc, char* argv[]) {

    //## Parse Inputs
    std::string virus_ref_file_name = argv[1];
    std::string workdir = argv[2];
    std::string workspace = argv[3];
    //## Files to be used from the workdir
    std::string stats_file_name = workspace + "/stats.txt";
    std::string config_file_name = workdir + "/config.txt";
    std::string region_fasta_file_name = workdir + "/regions.fna";
    std::string read_fasta_file_name = workdir + "/edge_reads.fna";
    std::string edge_file_name = workdir + "/edges.tab";
    std::string bamFile = workspace + "/retained-pairs.namesorted.bam";
    //## Output Files
    std::string reg_file_name = workdir + "/results.txt";
    std::string reads_dir = workdir + "/readsx";
    //## Configuration
    LoadVirusNames(virus_ref_file_name,VirusNameSet);
    Config = parse_config(config_file_name);
    AlnMaskLen =  Config.read_len/2;
    Stats = parse_stats(stats_file_name);
    JointHeader = sam_hdr_read(sam_open(bamFile.c_str(),"r"));
    //## Raw Data
    Read2RegionsMap_t  read2regSetMap;
    Region2ReadsMap_t  reg2readSetMap;
    EdgeVec_t  edgeVec;
    LoadData(edge_file_name,region_fasta_file_name,read_fasta_file_name,
	     read2regSetMap,reg2readSetMap,edgeVec);
    //## Alignments 
    AlignmentMap_t alnMap;
    AlignReads(read2regSetMap,alnMap);
    RemoveUnalignedReads(edgeVec,alnMap);
    //## Edge Processing
    OrderEdges(edgeVec,alnMap);
    //## Output
    OutputEdges(edgeVec,alnMap,reg_file_name,reads_dir);
    //## Cleanup
    bam_hdr_destroy(JointHeader);
    fprintf(stderr,"edge_mapper Done\n");
}

//==== FUNCTION DEFINITIONS

//Aligns a read to all regions with which it is associated
//  filters the alignments based on a minimum score, and a minimum score
//  relative to the best alignment
//Inputs - an id, used by thead pool
//	 - a read object
//	 - a set of regions
//	 - a constant reference to a mapping of subject query pairs generated
//		alignment objects
//Output - None, modifies the alignment mapping
void AlignRead(	int id, const Read_pt & read, const RegionSet_t & regSet,
		AlignmentMap_t & alnMap)
{
    uint16_t bestScore = 0;
    //Iterate over regions and do the alignments
    std::vector<const SQPair_t *> sqPairVec;
    for( const Region_pt & reg : regSet){
        const std::string & query = read->getSegment(	reg->isVirus,
							reg->strand == '-');
	Mtx.lock();
        auto res = alnMap.emplace(std::make_pair(  SQPair_t(reg,read),
    				    StripedSmithWaterman::Alignment()));
	Mtx.unlock();
        if(!res.second) continue;
        StripedSmithWaterman::Alignment & aln = res.first->second;
        bool bPass = Aligner.Align( query.c_str(),reg->sequence.c_str(),
				    reg->sequence.length(),
				    AlnFilter, &(aln),AlnMaskLen);
	if(bPass) {
	    bPass = accept_alignment(aln, Config.min_sc_size);
	}
        if(!bPass) { //If the Alignment failed
	    Mtx.lock();
	    alnMap.erase(res.first);
	    Mtx.unlock();
	    continue;
        }
        sqPairVec.push_back(&(res.first->first));
        if(aln.sw_score > bestScore){
	   bestScore = aln.sw_score;
        }
    }
    double minScore = 0.75 * double(bestScore);
    //Iterate over sq pairs and erase any which are below threshold
    for( const SQPair_t * & pPair : sqPairVec){
	Mtx.lock();
        if(alnMap.at(*pPair).sw_score < minScore){
	    alnMap.erase(*pPair);
        }
	Mtx.unlock();
    }
}

//Iterate over each read and align it to all associated regions
//  Filters out low scoring alignments
//  Paralleizes on Reads
//Inputs - a reference to a mapping from reads to regions
//	 - a reference to a mapping from read-region pairs to alignments
//Output - none, Modifies the alignment map
void AlignReads(const Read2RegionsMap_t &regMap,
		AlignmentMap_t & alnMap) {
    fprintf(stderr,"Aligning reads ...\n");
    ctpl::thread_pool threadPool (Config.threads);
    std::vector<std::future<void>> futureVec;
    for(const auto & pair : regMap){
        std::future<void> future = threadPool.push( AlignRead,
        					    std::cref(pair.first),
        					    std::cref(pair.second),
        					    std::ref(alnMap));
        futureVec.push_back(std::move(future));
    }
    int pert = 0;
    size_t complete =0;
    for( auto & future : futureVec){
        future.get();
        complete++;
        double progress = complete / double(futureVec.size());
        if(1000.0 * progress > pert){
            pert = 1000 * progress;
            fprintf(stderr,"Progress: %0.1f%%\r",progress*100.0);
        }
    }
    //Single threaded
    //for(const auto & pair : regMap){
    //    AlignRead(0,pair.first,pair.second,alnMap);
    //}
    fprintf(stderr,"\nPassing Alignments: %lu\n",alnMap.size());
}

//Splits an edge into a number of edges for each unique consensus
//sequence of reads observed
//Inputs - an edge to process
//	 - an alignment map
//	 - a vector in which to store newly created edges
//Output - None, modifies the newEdges vector and edge object
EdgeVec_t ConsensusSplitEdge(	int id, Edge_t & edge,
				const AlignmentMap_t & alnMap)
{
    EdgeVec_t newEdges;
    //Build the Table of aligned sequences
    std::vector<Read_pt> rowLabelVec;
    std::vector<std::string> rowSeqVec;
    std::vector<size_t> nFillVec;
    //To track which rows are still be processed
    for(const Read_pt & read : edge.readSet){
	rowLabelVec.push_back(read);
	nFillVec.push_back(0);
	rowSeqVec.push_back(GetAlignedSequence(	edge,read,alnMap,
						nFillVec.back()));
    }
    RecursiveSplitEdge(edge,rowLabelVec,rowSeqVec,nFillVec,newEdges); 
    rowLabelVec.clear();
    rowSeqVec.clear();
    nFillVec.clear();
    return newEdges;
}



//Given a read region pair, and alignment info, construct a bam entry for
//the pair's alignment
//Inputs - a read-region pair
//	 - an alignment Map
//	 - a bam1_t pointer to store the result in
//Output - Boolean whether the construction was successful or not
//	 - Also modifes the entry object
bool ConstructBamEntry(	const Read_pt & query, const Region_pt & subject,
			const Region_pt & mateSubject,
			const AlignmentMap_t & alnMap,
			bam1_t* entry) {
    const StripedSmithWaterman::Alignment & aln =
	alnMap.at(SQPair_t(subject,query));
    const StripedSmithWaterman::Alignment & mateAln =
	alnMap.at(SQPair_t(mateSubject,query));
    bool bVirus = (VirusNameSet.count(subject->chr));
    bool bRev = (subject->strand == '-');
    uint16_t flag = BAM_FPAIRED;
    if(bRev) flag |= BAM_FREVERSE;
    if(mateSubject->strand == '-') flag |= BAM_FMREVERSE;
    flag |= (bVirus) ? BAM_FREAD2 : BAM_FREAD1;
    entry->core.qual = 255;
    entry->core.l_extranul = (4 - (query->name.length() % 4)) % 4;
    entry->core.l_qname = query->name.length() + entry->core.l_extranul;
    const std::string & qSeq = query->getSegment(bVirus,bRev);
    int l_qseq = qSeq.length();
    std::vector<char> qual(l_qseq,'<');
    int l_aux = 0;
    auto ql = bam_cigar2qlen(aln.cigar.size(),aln.cigar.data());
    std::vector<uint32_t> cigar = aln.cigar;
    if(l_qseq != ql) return false;
    bam_set1(	entry,query->name.length(),query->name.c_str(), flag,
		sam_hdr_name2tid(JointHeader,subject->chr.c_str()),
		subject->left + aln.ref_begin, 255, cigar.size(), cigar.data(),
		sam_hdr_name2tid(JointHeader,mateSubject->chr.c_str()),
		mateSubject->left + mateAln.ref_begin,
		0, l_qseq, qSeq.c_str(), qual.data(), l_aux);
    return true;
}

//Constructs a breakpoint from a known region
//Inputs - an offset describing where in the region the offset is
//	 - a region
//Output - a breakpoint_t object (see util.h)
breakpoint_t ConstructBreakpoint(const Region_pt & reg,size_t offset){
    bool bRev = (reg->strand == '+');
    int pos = (!bRev) ? reg->left + offset : reg->right - offset;
    return breakpoint_t(reg->chr,pos,pos,bRev);
}

//Given an edge and alignment information, calculates summary stats
//  hostPBS - the average score/aligned base on the host side
//  coverage - the half the total length of the host and virus sides as a
//		proportion of the max insert size
//Then constructs a call_t objects which can be output
//Inputs - an identifier for the output junction
//	 - an edge object
//	 - an alignment map
//Output - a call_t object (see utils.h)
call_t ConstructCall(	int id, const Edge_t & edge,
			const AlignmentMap_t & alnMap)
{
    size_t nReads = edge.readSet.size();
    breakpoint_t hostBP = ConstructBreakpoint(	edge.hostRegion,
						edge.hostOffset);
    breakpoint_t virusBP = ConstructBreakpoint(	edge.virusRegion,
						edge.virusOffset);
    double hostPBS = 0, virusPBS = 0;
    double hostCov = 0, virusCov = 0;
    size_t hostLeft = edge.hostRegion->sequence.length(), hostRight = 0;
    size_t virusLeft = edge.virusRegion->sequence.length(), virusRight = 0;
    int score = 0;
    for( const Read_pt & read : edge.readSet){
	SQPair_t hPair(edge.hostRegion,read);
	SQPair_t vPair(edge.virusRegion,read);
	const StripedSmithWaterman::Alignment & hAln = alnMap.at(hPair);
	const StripedSmithWaterman::Alignment & vAln = alnMap.at(vPair);
	if(hAln.ref_begin < hostLeft) hostLeft = hAln.ref_begin;
	if(vAln.ref_begin < virusLeft) virusLeft = vAln.ref_begin;
	if(hAln.ref_end > hostRight) hostRight = hAln.ref_end;
	if(vAln.ref_end > virusRight) virusRight = vAln.ref_end;
	double hLen = hAln.query_end - hAln.query_begin + 1;
	double vLen = vAln.query_end - vAln.query_begin + 1;
	hostPBS += double(hAln.sw_score) / hLen;
	virusPBS += double(vAln.sw_score) / vLen;
	score += hAln.sw_score + vAln.sw_score;
    }
    hostPBS /= double(nReads);
    virusPBS /= double(nReads);
    if(hostLeft <= hostRight) {
	hostCov = double(hostRight - hostLeft) / Stats.max_is;
    }
    if(virusLeft <= virusRight) {
	virusCov = double(virusRight - virusLeft) / Stats.max_is;
    }
    return call_t(id,hostBP,virusBP,nReads,nReads,edge.nSplit,0,0,
	    score,hostPBS,virusPBS,hostCov,virusCov);
}

//Given the an alignment from either/host or virus 
//Assigns the  reference begin/end offsets into juncton proximal/distal offsets
//  depending on the strand of the reference
// Host alignments are distal-proximal, unless reversed
// Viral alignments are proximal-distal, unless reversed
// Boils down to an XOR operation on the viralness and reversedness
//  if one is true proximal-distal, else distal-proximal
//Inputs - the strand of the reference
//	 - whether this alignment was against a viral reference
//	 - the alignment in question
//Output - a structure containing the junction distal and proximal positions
JunctionInterval_t ConstructJIV(char strand, bool isVirus,
				const StripedSmithWaterman::Alignment & aln)
{
    bool bRev = (strand == '-');
    JunctionInterval_t jIV;
    //!= is an XOR operation on boolean values
    if(bRev != isVirus){
	jIV.proximal = aln.ref_begin;
	jIV.distal = aln.ref_end;
    } else {
	jIV.distal = aln.ref_begin;
	jIV.proximal = aln.ref_end;
    }
    return jIV;
}

//Identifies Duplicate reads at an edge, suplicates are defined as having
//the same left map in the host and right map in the virus
//Inputs - an edge to process
//	 - an alignment map to inform the deduplication
//Output - None, modifies the given object
void DeduplicateEdge(Edge_t & edge,const AlignmentMap_t & alnMap) {
    PosPairSet_t outPosSet;
    std::vector<Read_pt> toRemoveVec;
    for(const Read_pt & read : edge.readSet){
	SQPair_t hPair(edge.hostRegion,read);
	SQPair_t vPair(edge.virusRegion,read);
	const StripedSmithWaterman::Alignment & hAln = alnMap.at(hPair);
	const StripedSmithWaterman::Alignment & vAln = alnMap.at(vPair);
	JunctionInterval_t hJIV = ConstructJIV(edge.hostRegion->strand,false,hAln);
	JunctionInterval_t vJIV = ConstructJIV(edge.virusRegion->strand,true,vAln);
	//If two reads have the same distal ends of their alignment, they are considered
	//duplicates
	PosPair_t pair = std::make_pair(hJIV.distal,vJIV.distal);
	auto res = outPosSet.insert(pair);
	if(!res.second){ // posPair not yet seen
	    toRemoveVec.push_back(read);
	}
    }
    for(const Read_pt & read : toRemoveVec){
	edge.removeRead(read);
    }
}

//Sets the characters of an output string to the appropriate characters
//from an input string according to a set of cigar operations
//Only up to maxpos in the output string will be altered
//The output is always in reference coordinates, as a result insertions
//are treated as alignments, that is they consume both the reference and
//query
//This may lead to the sequence being truncated at maxpos, which is fine
//The only situation this is a problem in is if a sequencing error exists
//at the end of a sequence with a long insertion *shrug*
//Inputs - a sequence to fill in
//	 - a sequence to fill from
//	 - where in the outseq to begin filling
//	 - where in the outsew to stop filling
//	 - a vector of bam formated cigar information
//Output - The number of defined characters filled in
size_t FillStringFromAlignment(	std::string & outseq,
				const std::string & inseq,
				size_t offset, size_t maxpos,
				const std::vector<uint32_t> & cigarVec) {
    size_t pos = offset;
    size_t qpos = 0;
    size_t nFill = 0;
    //Iterate over cigar operations and set characters in the outseq
    for(int i = 0; i < cigarVec.size(); i++){
	size_t opLen = bam_cigar_oplen(cigarVec[i]);
	char opchr = bam_cigar_opchr(cigarVec[i]);
	switch (opchr) {
	    case 'H':
	    case 'S':
		qpos += opLen;
		break;
	    case 'M':
	    case '=':
	    case 'X':
	    case 'I': //Insertions taken as is
		for(int j = 0;
		    j < opLen & pos < maxpos & qpos < inseq.size();
		    j++)
		{
		    outseq[pos++] = inseq[qpos++];
		    nFill++;
		}
		break;
	    case 'D':
		for(int j = 0;
		    j < opLen & pos < maxpos & qpos < inseq.size();
		    j++)
		{
		    outseq[pos++] = '-';
		    nFill++;
		}
		break;
	}
    }
    return nFill;
}

//Filters an edge Vector to only contain edges which have enough effective
//reads
//IMPORTANT: does not maintain element order, follow up with a sort!
//Inputs - an edge vector
//	 - a pointer to a set of used reads, may be null
//Output - none, modifies the given edge vector
void FilterEdgeVec(EdgeVec_t & edgeVec, const ReadSet_t * usedReads){
    size_t filtered = 0;
    for(auto it = edgeVec.begin(); it != edgeVec.end(); ){
        if(PassesEffectiveReadCount(*it,usedReads)){
            it++;
        } else {
	    //Delete the element by overwriting with last element,
	    //	then delete the last element
	    //Does not maintain element order, but is fast
	    *it = std::move(edgeVec.back());
	    edgeVec.pop_back();
            filtered++;
        }
    }
}


//Identifies read pairs with an apparent insert size which is too large
//and removes them
//Inputs - an edge to process
//	 - an alignment map 
//Output - None, modifies the edge object
void FilterHighInsertReads(Edge_t & edge, const AlignmentMap_t & alnMap){
    std::vector<Read_pt> toRemoveVec;
    for(const Read_pt & read : edge.readSet){
	SQPair_t hPair(edge.hostRegion,read);
	SQPair_t vPair(edge.virusRegion,read);
	const StripedSmithWaterman::Alignment & hAln = alnMap.at(hPair);
	const StripedSmithWaterman::Alignment & vAln = alnMap.at(vPair);
	JunctionInterval_t hJIV = ConstructJIV(edge.hostRegion->strand,false,hAln);
	JunctionInterval_t vJIV = ConstructJIV(edge.virusRegion->strand,true,vAln);
	size_t hIS = 1 + ((edge.hostOffset > hJIV.distal) ? (edge.hostOffset - hJIV.distal) :
			(hJIV.distal - edge.hostOffset));
	size_t vIS = 1 + ((edge.virusOffset > vJIV.distal) ? (edge.virusOffset - vJIV.distal) :
			(vJIV.distal - edge.virusOffset));
	size_t is = hIS + vIS;
	if(is > Stats.max_is){ // Insert size is too high
	    toRemoveVec.push_back(read);
	}
    }
    for(const Read_pt & read : toRemoveVec){
	edge.removeRead(read);
    }
}

//Remove reads that have suspicious alignments, alignments are considered
//suspicious if:
//  the aligned portion of the query or reference are low complexity
//  a split read's aligned position is too far from the breakpoint
//Inputs - an edge to process
//	 - an alignment map
//Output - None, modifies the edge object
void FilterSuspiciousReads(Edge_t & edge, const AlignmentMap_t & alnMap) {
    std::vector<Read_pt> toRemoveVec;
    for(const Read_pt & read : edge.readSet){
	std::array<Region_pt *,2> regArr = {&(edge.hostRegion),
					    &(edge.virusRegion)};
	for(const Region_pt * reg_p : regArr){
	    const StripedSmithWaterman::Alignment & aln =
		alnMap.at(SQPair_t(*reg_p,read));
	    bool bRev = ((*reg_p)->strand == '-');
	    const std::string & readSeq = read->getSegment( (*reg_p)->isVirus,
							    bRev);

	    bool qLC = is_low_complexity(readSeq.c_str(),
					aln.query_begin,aln.query_end);
	    bool rLC = is_low_complexity((*reg_p)->sequence.c_str(),
					aln.ref_begin,aln.ref_end);
	    if(qLC || rLC){
		toRemoveVec.push_back(read);
		continue;
	    }
	    uint32_t lClipLen = (bam_cigar_opchr(aln.cigar.front()) == 'S') ? 
				    bam_cigar_oplen(aln.cigar.front()) : 0;
	    uint32_t rClipLen = (bam_cigar_opchr(aln.cigar.back()) == 'S') ? 
				    bam_cigar_oplen(aln.cigar.back()) : 0;
	    //Only the matching clip is comparable to the breakpoint 
	    //	Left for virus, Right for host (based on all prior work
	    //	to make sure that's how things are arranged)
	    uint32_t clipLen = ((*reg_p)->isVirus) ? lClipLen : rClipLen;
	    uint32_t offset = ((*reg_p)->isVirus) ? edge.virusOffset :
						    edge.hostOffset;
	    //At this point breakpoints were defined by alignments
	    // :: no alignment to virus will start before the breakpoint
	    // :: no alignment to host will end after the breakpoint
	    uint32_t sMiss = (offset > aln.ref_begin) ?
			       offset - aln.ref_begin : aln.ref_begin - offset;
	    uint32_t eMiss = (offset > aln.ref_end) ?
			       offset - aln.ref_end : aln.ref_end - offset;
	    uint32_t miss = ((*reg_p)->isVirus) ? sMiss : eMiss;
	    //If the read is clipped enough, but starts/ends too far
	    //from the breakpoint it is wrongly clipped
	    if(clipLen > Config.max_sc_dist && miss > Config.max_sc_dist){
		toRemoveVec.push_back(read);
		continue;
	    }
	}
    }
    for(const Read_pt & read : toRemoveVec){
	edge.removeRead(read);
    }
}

//Generic function for filtering a vector to only a given set of indexes
//Inputs - a vector to process
//	 - a set of indexes
//Output - none, modifes the given vector
template<class T>
void FilterVector(  std::vector<T> & vec,
		    const std::unordered_set<size_t> & idxSet)
{
    size_t idx = 0;
    for( auto it = vec.begin(); it != vec.end(); idx++){
	if(!idxSet.count(idx)){
	    it = vec.erase(it);
	} else {
	    it++;
	}
    }
}



//Given a vecotr of strings (all assumed to be the same length) construct
//a majority rule consensus string
//Also record the number of differences from the generated consensus for
//each string
//Inputs - a vector of strings
//	 - a reference to a vector of size_t, will be scaled to rowVec
//	 Size, and witll contain the # of sites which differ from the
//	 conensus for each row
//Output - a majority rule consensus sequence
std::string GenerateConsensus(const std::vector<std::string> & rowVec,
				std::vector<size_t> & diffVec)
{
    
    if(!rowVec.size()) return std::string();
    if(diffVec.size()) diffVec.clear();
    std::string cons(rowVec.front().length(),'N');
    size_t nCol = rowVec[0].length();
    for(size_t col = 0; col < nCol; col++){
        std::unordered_map<char,size_t> nucCount;
        size_t max = 0;
        char best = 'N';
        //Iterate to determine consensus residue
        for(const std::string & seq : rowVec){
            char nuc = seq.at(col);
            size_t count = ++nucCount[nuc];
            if(count > max){
        	best = nuc;
        	max = count;
            }
        }
        cons[col] = best;
	diffVec.push_back(0);
        //Iterate to count differences from consensus
        for(size_t row = 0; row < rowVec.size(); row++){
            char nuc = rowVec.at(row).at(col);
            if(nuc != best && nuc != 'N'){
        	diffVec.back()++;
            }
        }
    }
    return cons;
}

//Uses the region information from an edge, and alignment information to
//generate the aligned sequence against a concatenation of both regions
//The Ooutput will have N's in positions with N's or no coverage
//Inputs - an edge, containing the host and virus region information
//	 - a read for which to contruct the sequence
//	 - an alignment map containing alignment information for the read
//	 - a reference to a size_t to store the number of filled characters
//Output - a string which represents the alignment to the host region
//	    concatenated to the alignmnet for the virus region
std::string GetAlignedSequence(	const Edge_t & edge, const Read_pt & read,
				const AlignmentMap_t & alnMap, size_t & nFill)
{
    const Region_pt & hReg = edge.hostRegion;
    const Region_pt & vReg = edge.virusRegion;
    SQPair_t hrPair(hReg,read);
    SQPair_t vrPair(vReg,read);
    const StripedSmithWaterman::Alignment & hostAln = alnMap.at(hrPair);
    const StripedSmithWaterman::Alignment & virusAln = alnMap.at(hrPair);
    size_t hostLen = hReg->sequence.length();
    size_t virusLen = vReg->sequence.length();
    size_t totalLen =	hostLen + virusLen;
    std::string outseq(totalLen,'N');
    //Determine which strand of was aligned to the region
    const std::string & hSeq = read->getSegment(false,hReg->strand == '-');
    const std::string & vSeq = read->getSegment(true,vReg->strand == '-');
    nFill = 0;
    //Fill in the host side of the alignment
    nFill += FillStringFromAlignment(	outseq,hSeq,hostAln.ref_begin,
					hostLen,hostAln.cigar);
    //Fill in the host side of the alignment
    nFill += FillStringFromAlignment(	outseq,vSeq,
					virusAln.ref_begin+hostLen,
					totalLen, virusAln.cigar);
    return outseq;
}


//Sums the scores on both the host and viral sides for all reads
//supporting the edges
//Edges with unique reads are given a boost in that all scores are
//multiplied by the proportion of their reads which are unique
//Inputs - an for which to calculate the score
//	 - an alignment map to get scores from
//	 - a set of reads which have been used in other junctions
//Output - A double value representing the edge's score
double GetEdgeScore(const Edge_t & edge,const AlignmentMap_t & alnMap,
		    const ReadSet_t & usedReads) {
    double score = 0;
    double numer = double(edge.uniqueReadSet.size() + 1);
    double denom = double(edge.readSet.size() + 2);
    double uniqueProp =  numer / denom;
			
    for(const Read_pt & read : edge.readSet){
	if(usedReads.count(read)) continue; 
	SQPair_t hPair(edge.hostRegion,read);
	SQPair_t vPair(edge.virusRegion,read);
	const StripedSmithWaterman::Alignment & hAln = alnMap.at(hPair);
	const StripedSmithWaterman::Alignment & vAln = alnMap.at(vPair);
	score += hAln.sw_score + vAln.sw_score;	
    }
    score *= uniqueProp;
    return score;
}


//Based on reads assigned to an edge, determine where within the host and
//viral regions the actual breakpoints are
//This goes through each read and finds the most extreme left and rightmost position
//The host breakpoint is the rightmost position in the host region, unless reversed
//The virus breakpoint is the leftmost position in the virus_region, unless reversed
//This position is only within the reference region, if the region was reversed remember to
//account for this during outputting
//Inputs - an edge to process
//	 - an alignment map
//Output - None, modifies the edge object
void IdentifyEdgeBreakpoints(Edge_t & edge, const AlignmentMap_t & alnMap){
    JunctionInterval_t hostIV = {0,edge.hostRegion->sequence.length()};
    JunctionInterval_t virusIV = {edge.virusRegion->sequence.length(),0};
    //Iterate over reads to find the extremes
    for( const Read_pt & read : edge.readSet){
	SQPair_t hPair(edge.hostRegion,read);
	SQPair_t vPair(edge.virusRegion,read);
	const StripedSmithWaterman::Alignment & hAln = alnMap.at(hPair);
	const StripedSmithWaterman::Alignment & vAln = alnMap.at(vPair);
	if(hAln.ref_end > hostIV.proximal) hostIV.proximal = hAln.ref_end;
	if(hAln.ref_begin < hostIV.distal) hostIV.distal = hAln.ref_begin;
	if(vAln.ref_begin < virusIV.proximal) virusIV.proximal = vAln.ref_begin;
	if(vAln.ref_end > virusIV.distal) virusIV.distal = vAln.ref_end;
    }
    if(edge.hostRegion->strand == '-'){
	std::swap(hostIV.proximal,hostIV.distal);
    }
    if(edge.virusRegion->strand == '-'){
	std::swap(virusIV.proximal,virusIV.distal);
    }
    edge.hostOffset = hostIV.proximal;
    edge.virusOffset = virusIV.proximal;
}

//Determines which if any reads assigned to an edge are only assigned to
//that edge
//Inputs - a vector if edges
//Output - None, modifies the edges
void IdentifyEdgeSpecificReads(EdgeVec_t & edgeVec) {
    std::unordered_map<Read_pt,size_t,Read_pt_HashFunctor,Read_pt_EqFunctor>
	edgeCount;
    for(const Edge_t & edge : edgeVec){
	for(const Read_pt & read : edge.readSet){
	    edgeCount[read]++;
	}
    }
    for(Edge_t & edge : edgeVec){
	for(const Read_pt & read : edge.readSet){
	    if(edgeCount[read] == 1){
		edge.uniqueReadSet.insert(read);
	    }
	}
    }
}

//Gicen two partial DNA sequences tests if the two have consistent
//sequences: that is they match at all non-N positions
//Inputs - two strings representing the two sequences
//Output - a boolean of whether they are consistent or not
bool IsConsistent(const std::string & seq1, const std::string & seq2){
    if(seq1.length() != seq2.length()) return false; 
    for(size_t i = 0; i < seq1.length(); i++){
	char c1 = seq1.at(i);
	char c2 = seq2.at(i);
	if(c1 != c2 && c1 != 'N' && c1 != 'N') return false;
    }
    return true;
}

//Wrapper function which handles all basic loading
//Inputs - strings representing the edges file, region sequence file,
//	    and bam file containing read sequences
//	 - references to objects to store the edges, and bidirectional
//	 mappings of reads and regions
//Output - None , modifies provided references
void LoadData(	const std::string & edgeFName,
		const std::string & regionsFName,
		const std::string & bamFName,
	      	Read2RegionsMap_t & read2regSetMap,
	      	Region2ReadsMap_t & reg2readSetMap,
	      	EdgeVec_t & edgeVec)
{
    Name2RegionMap_t regNameMap;
    Name2ReadMap_t readNameMap;
    LoadRegionSeq(regionsFName,reg2readSetMap,regNameMap);
    LoadReadSeq(bamFName,read2regSetMap,readNameMap);
    LoadEdges(	edgeFName,read2regSetMap,reg2readSetMap,edgeVec,
		readNameMap,regNameMap);
}

//Parses edge file to load the regions, reads, and their associations
//Inputs - a string representing the edge file name
//	 - a reference to mapping between reads and region sets
//	 - a reference to a mpaaing from regions to read sets
//	 - a reference to a vector of edges
//	 - a const reference to a mapping from read names to reads
//	 - a const reference to a mapping from region names to reads
//Output - None, modifies to provided references
void LoadEdges(	std::string edgeFName,
	      	Read2RegionsMap_t & read2regSetMap,
	      	Region2ReadsMap_t & reg2readSetMap,
	      	EdgeVec_t & edgeVec,
		const Name2ReadMap_t & readNameMap,
		const Name2RegionMap_t & regionNameMap)
{
    fprintf(stderr,"Loading Edges ...\n");
    std::ifstream in(edgeFName);
    if(!in.is_open()){
	fprintf(stderr,"[ERROR] Could not open %s for reading\n",
		edgeFName.c_str());
	throw std::invalid_argument("Failure to read Edge file");
    }
    std::string regStr,readStr;
    size_t nRead;
    std::string line;
    while(in >> regStr >> readStr >> nRead){
	std::vector<std::string> regionStringVec = strsplit(regStr,':');
	std::vector<std::string> readStringVec = strsplit(readStr,',');
        edgeVec.emplace_back(	regionNameMap.at(regionStringVec[0]),
				regionNameMap.at(regionStringVec[1]));
	Edge_t & edge = edgeVec.back();
	for(std::string & rName : readStringVec){
	    if(!rName.length()) continue;
	    char segment = rName[rName.length()-1];
	    if(	segment == 'H' || segment == 'V' ||
		segment == 'h' || segment == 'v')
	    {
		rName = rName.substr(0,rName.length()-2);
	    }
	    const Read_pt & read = readNameMap.at(rName);
	    edge.readSet.insert(read);
	    if(read->isSplit) edge.nSplit++;
	    //Update the read 2 region and region 2 read maps
	    read2regSetMap.at(read).insert(edge.hostRegion);
	    read2regSetMap.at(read).insert(edge.virusRegion);
	    reg2readSetMap.at(edge.hostRegion).insert(read);
	    reg2readSetMap.at(edge.virusRegion).insert(read);
	}
    }
    fprintf(stderr,"Loaded %lu Edges\n",edgeVec.size());
}

//Parses a fasta file containing reads and stores their sequences
//The header of the entries is parsed for a trailing _1 or _2 to
//seperate fwd and reverse reads
//Sequences are always read in as upper case
//Inputs - a string representing the bam file name
//	 - a reference to a mapping from reads to regions sets
//Output - None, modifies the key objects of the map
void LoadReadSeq(   const std::string & readsFName,
		    Read2RegionsMap_t & read2regSetMap,
		    Name2ReadMap_t & nameMap)
{
    fprintf(stderr,"Loading Reads ...\n");
    FILE* readsFasta = fopen(readsFName.c_str(),"r");
    kseq_t *seq = kseq_init(fileno(readsFasta));
    while(kseq_read(seq) >= 0){
	std::string name = seq->name.s;
	char segment = name[name.length() - 1];
	bool bSplit = true;
	bool bViralR1 = false;
	if( segment == 'H' || segment == 'V' ||
	    segment == 'h' || segment == 'v')
	{
	    if(segment == 'h' || segment == 'v') bViralR1 = true;
	    name = name.substr(0,name.length()-2);
	    bSplit = false;
	}
	if(!nameMap.count(name)){ //Create the read object if necessary
	    Read_pt read = std::make_shared<Read_t>(name,bSplit,bViralR1);
	    nameMap.insert(std::make_pair(name,read));
	}
	//Add the sequence to the read
	std::string sequence(seq->seq.s);
	for (auto & c: sequence) c = (char)toupper(c);
	switch (segment) {
	    case '1':
	    case '2':
		nameMap[name]->hostSegment = sequence;
		nameMap[name]->virusSegment = sequence;
		nameMap[name]->hostRC = get_seqrc(sequence);
		nameMap[name]->virusRC = get_seqrc(sequence);
		break;
	    case 'H':
	    case 'h':
		nameMap[name]->hostSegment = sequence;
		nameMap[name]->hostRC = get_seqrc(sequence);
		break;
	    case 'V':
	    case 'v':
		nameMap[name]->virusSegment = sequence;
		nameMap[name]->virusRC = get_seqrc(sequence);
		break;
	}
    }
    //Create the key-value pairs in the output map
    //Values are default constructed
    for( auto pair : nameMap){
	read2regSetMap[pair.second];
    }
    kseq_destroy(seq);
    fclose(readsFasta);
    fprintf(stderr,"Loaded %lu reads\n",read2regSetMap.size());
}

//Parses a fasta file containing region sequences and stores their
//sequences
//Sequences are always read in as upper case
//Inputs - a string representing the fasta file name
//	 - a reference to a mappinf from regions to read sets
//Output - None, modifes the key objects of the map
void LoadRegionSeq( const std::string & regionsFName,
		    Region2ReadsMap_t & reg2readSetMap,
		    Name2RegionMap_t & nameMap)
{
    fprintf(stderr,"Loading Regions ...\n");
    FILE* regionsFasta = fopen(regionsFName.c_str(),"r");
    kseq_t *seq = kseq_init(fileno(regionsFasta));
    while(kseq_read(seq) >= 0){
	std::string name = seq->name.s;
	size_t pos = name.find(':',1);
	if(pos == std::string::npos) continue;
	seq->name.l = pos; //set the current end of the name earlier
	seq->name.s[pos] = '\0';
	std::string contig = strsplit(std::string(seq->name.s),',')[0];
	bool bVirus = (VirusNameSet.count(contig));
	Region_pt reg = std::make_shared<Region_t>(seq,bVirus);
	nameMap.insert(std::make_pair(std::string(seq->name.s),reg));
	auto res = reg2readSetMap.emplace(std::make_pair(reg,ReadSet_t()));
	if(!res.second){
	    fprintf(stderr,"[WARNING] Duplicate Region Sequence ignored\n");
	}
    }
    kseq_destroy(seq);
    fclose(regionsFasta);
    fprintf(stderr,"Loaded %lu regions\n",reg2readSetMap.size());
}

//Process Edges and puts them into an order from most likely to be real to
//least
//Processing includes:
//  Spliting Edges which appear to have multiple
//  consensus sequences
//  Finding the actual breakpoints
//  Deduplicating
//  Removing High insert size reads
//  Identifying which if any reads are unique to the edge
//Inputs - a reference to an edge vector
//	 - a const reference to an alignment map
//Output - None, modifies the edge vector
void OrderEdges(EdgeVec_t & edgeVec,const AlignmentMap_t & alnMap) {
    fprintf(stderr,"Ordering Edges ...\n");
    EdgeVec_t newEdges = SplitEdges(edgeVec,alnMap); 
    //Eliminate Edges with low read counts
    FilterEdgeVec(edgeVec);
    //Add any new edges back in (these are already filtered)
    edgeVec.insert(edgeVec.end(),newEdges.begin(),newEdges.end());
    //Process all edges
    ProcessEdges(edgeVec,alnMap);
    //Remove insufficiently supported Edges
    //Find the edges which are unique to a particular edge
    IdentifyEdgeSpecificReads(edgeVec);
    ReadSet_t usedReads; //Required by GetEdgeScore, but not acually needed yet
    std::sort(edgeVec.begin(), edgeVec.end(),
	    [&alnMap,&usedReads](Edge_t & a, Edge_t & b){
		//Sort in descending order (we'll process from the back)
		return (GetEdgeScore(a,alnMap,usedReads) <
			GetEdgeScore(b,alnMap,usedReads));
	    });
    fprintf(stderr,"Ordered %lu edges\n",edgeVec.size());
}

void OutputEdge(int id, const Edge_t & edge, const AlignmentMap_t & alnMap,
		std::ofstream & out, const std::string & readDir,
		ReadSet_t usedReads)
{
    //Output the call
    call_t call = ConstructCall(id, edge,alnMap);
    out << call.to_string() << "\n"; 
    //Output the reads
    samFile* writer = open_bam_writer(	readDir,std::to_string(id)+".bam",
					JointHeader);
    bam1_t* entry = bam_init1();
    for(const Read_pt & read : edge.readSet){
	if(usedReads.count(read)) continue;
	bool bCons = ConstructBamEntry( read,edge.hostRegion,
	    			    edge.virusRegion,alnMap,
	    			    entry);
	if(!bCons) throw std::runtime_error("Cigar failure");
	int ok = sam_write1(writer,JointHeader,entry);
	if(ok < 0) throw std::runtime_error("Failed to write to " +
					    std::string(writer->fn));
	ConstructBamEntry(	read,edge.virusRegion,edge.hostRegion,alnMap,
	    		entry);
	ok = sam_write1(writer,JointHeader,entry);
	if(ok < 0) throw std::runtime_error("Failed to write to " +
					    std::string(writer->fn));
    }

    bam_destroy1(entry);
    sam_close(writer);
}

//Proceeds from high confidence edges to low confidence edges, ensuring
//each read is used exactly once
//  The edges may be reordered as 
//Inputs - a reference to a vector of edges sorted from most to
//	    least confident
//	 - a const reference to an alignment map
//	 - a string representing the results file
//	 - a string representing the reads directory 
//Output - None, prints to outfile
void OutputEdges(EdgeVec_t & edgeVec,const AlignmentMap_t & alnMap,
		 const std::string & resFName, const std::string & readDir)
{
    fprintf(stderr,"Outputting Edges ensuring reads support only one edge...\n");
    ReadSet_t usedReads;
    std::ofstream out(resFName);
    int nextJunctionID = 0;
    int start = edgeVec.size();
    int pert = 0;
    while(edgeVec.size()) {
	for(const Read_pt & read : edgeVec.back().readSet){
	    usedReads.insert(read);
	}
	OutputEdge( nextJunctionID++,edgeVec.back(),alnMap,out,readDir,
		    usedReads);
	edgeVec.pop_back();
	FilterEdgeVec(edgeVec,&usedReads);
	std::sort(edgeVec.begin(), edgeVec.end(),
	    [&alnMap,&usedReads](Edge_t & a, Edge_t & b){
		return (GetEdgeScore(a,alnMap,usedReads) <
			GetEdgeScore(b,alnMap,usedReads));
	    });
	double progress = (start - edgeVec.size()) / double(start);
    	    if(1000.0 * progress > pert){
    	        pert = 1000 * progress;
    	        fprintf(stderr,"Progress: %0.1f%%\r",progress*100.0);
    	    }
    	}
    fprintf(stderr,"\nOutput %d Edges...\n",nextJunctionID);
}

//Given an edge reports wheteher it has enough effective reads
//Inputs - an Edge
//Output - boolean wether it has enough reads
bool PassesEffectiveReadCount(	const Edge_t & edge,
				const ReadSet_t * usedReads)
{
    size_t count = edge.readSet.size() + ((edge.nSplit) ? SplitBonus : 0);
    if(usedReads){
	count = 0;
	bool bSplit = false;
	for(const Read_pt & read : edge.readSet){
	    if(!usedReads->count(read)){
		count++;
		if(read->isSplit) bSplit = true;
	    }
	}
	if(bSplit) count += SplitBonus;
    }
    return (count >= MinimumReads);
}

///Performs all filtering steps on the edge
//  Identifying break point locations
//  Deduplicating reads
//  Removing High insert size reads
//Inputs - an id, used by thread_pool
//	 - an edge
//	 - an alignment map
void ProcessEdge(int id,Edge_t & edge, const AlignmentMap_t & alnMap){
    IdentifyEdgeBreakpoints(edge,alnMap);
    DeduplicateEdge(edge,alnMap);
    FilterHighInsertReads(edge,alnMap);
    FilterSuspiciousReads(edge,alnMap);
}

void ProcessEdges(EdgeVec_t & edgeVec, const AlignmentMap_t & alnMap){
    fprintf(stderr,"Processing %ld Edges ...\n",edgeVec.size());
    ctpl::thread_pool threadPool (Config.threads);
    std::vector<std::future<void>> futureVec;
    for( Edge_t & edge : edgeVec){
	auto future = threadPool.push(	ProcessEdge,std::ref(edge),
					std::cref(alnMap));
	futureVec.push_back(std::move(future));
    }
    int pert = 0;
    size_t complete = 0;
    for (auto & future : futureVec){
	future.get();
	complete++;
	double progress = complete / double(futureVec.size());
	if(1000.0 * progress > pert){
	    pert = 1000 * progress;
	    fprintf(stderr,"Progress: %0.1f%%\r",progress*100.0);
	}
    }
    FilterEdgeVec(edgeVec);
    fprintf(stderr,"\nProcessed and retained %ld Edges\n",edgeVec.size());
}

//Recursivly processes prepared data describing the sequences of an edge
//First a consensus sequence is generated for the edge
//then all reads which are too different from the consensus (adjusted for
//the length of the read) are identified
//These reads are removed fromt he parent edge and moved to a child edge
//Additionally any reads from the parent edge which are completely
//consistent with any of the discarded reads are also included
//Overall this allows the potential for multiple alleles of junctions
//Any child edges with too fiew reads are ignored
//Continues until no valid child edge is made
//Inputs - an edge to process
//	 - a vector of the reads in the edge
//	 - a vector of the aligned sequences of the reads
//	 - a vector of the aligned length of the reads
//	 - a reference to a vector of edges in which to store new edges
//Output - None, modifies all inputs
void RecursiveSplitEdge(Edge_t & edge, std::vector<Read_pt> & rowLabelVec,
			std::vector<std::string> & rowSeqVec,
			std::vector<size_t> & nFillVec,
			EdgeVec_t & newEdges)
{
    size_t nRowIn = rowLabelVec.size();
    std::vector<size_t> diffCount;
    std::string consensus = GenerateConsensus(rowSeqVec,diffCount);
    Edge_t * newEdge_p = nullptr;
    std::unordered_set<size_t> roiSet;
    for(size_t a = 0; a < rowSeqVec.size(); a++){
        //Calculate the # of diffs per defined site
        double diffRate = double(diffCount[a]) / double(nFillVec[a]);
        if(diffRate < MaxDiffRate) continue;
        //Remove this read from the parent edge
        edge.removeRead(rowLabelVec[a]);
        //Build the new edge if necessary
        if(!newEdge_p){
            newEdges.emplace_back(edge.hostRegion,edge.virusRegion);
            newEdge_p = & newEdges.back();
        }
        //Add this read to the new Edge
        newEdge_p->addRead(rowLabelVec[a]);
        roiSet.insert(a);
        //Any reads consistent with this read will be included
        for(size_t b = a + 1; b < rowSeqVec.size(); b++){
            if(!IsConsistent(rowSeqVec[a],rowSeqVec[b])) continue;
            newEdge_p->addRead(rowLabelVec[b]);
            roiSet.insert(b);
        }
    }
    //We are done if no new edge was created
    if(!newEdge_p) return;
    //We are also done if the new edge is too small
    if(!PassesEffectiveReadCount(*newEdge_p)){
	newEdges.pop_back();
	newEdge_p = nullptr;
	return;
    }
    //Reduce the vectors to only the rows of interest for the new edge
    FilterVector(rowLabelVec,roiSet); 
    FilterVector(rowSeqVec,roiSet); 
    FilterVector(nFillVec,roiSet); 
    //Prevent infinite recursion by requiring that the recursion stops if
    //the next round isn't smaller
    if(rowLabelVec.size() >= nRowIn) return;
    RecursiveSplitEdge(*newEdge_p,rowLabelVec,rowSeqVec,nFillVec,newEdges);
}

//Some reads may have failed during alignment, remove them from the edges
//Inputs - a vector of edges to modify
//	 - an alignment map to check
//Output - None, modifes the edge vector
void RemoveUnalignedReads(EdgeVec_t & edgeVec,const AlignmentMap_t & alnMap){
    fprintf(stderr,"Removing Unaligned reads from %ld edges...\n",edgeVec.size());
    for( Edge_t & edge : edgeVec){
	std::vector<Read_pt> toRemoveVec;
	for( const Read_pt & read : edge.readSet){
	    //Check if both the host and virus alignments passed
	    if(	!alnMap.count(SQPair_t(edge.hostRegion,read)) ||
		!alnMap.count(SQPair_t(edge.virusRegion,read)))
	    {
		toRemoveVec.push_back(read);
	    }
	}
	for( const Read_pt & read : toRemoveVec){
	    edge.removeRead(read);
	}
    }
    FilterEdgeVec(edgeVec);
    fprintf(stderr,"Edges remaining: %ld\n",edgeVec.size());
}

//Given a vector of edges, in parallel splits each into edges for
//each consensus sequence present
//Inputs - a vector of edges
//	 - an alignment map
//Output - a vector of new edges, also modifies the edges in the input vecto
EdgeVec_t SplitEdges(EdgeVec_t & edgeVec, const AlignmentMap_t & alnMap){
    fprintf(stderr,"Splitting Edges based on consensus sequences ...\n");
    ctpl::thread_pool threadPool (Config.threads);
    std::vector<std::future<EdgeVec_t>> futureVec;
    for( Edge_t & edge : edgeVec){
        auto future = threadPool.push(	ConsensusSplitEdge,std::ref(edge),
        				std::cref(alnMap));
        futureVec.push_back(std::move(future));
    }
    EdgeVec_t newEdges;
    for(auto & future : futureVec){
        EdgeVec_t localNew = future.get();
        newEdges.insert(newEdges.end(),localNew.begin(),localNew.end());
    }
    /*//Single Threaded version
    EdgeVec_t newEdges;
    for( Edge_t & edge : edgeVec){
        EdgeVec_t localNew = ConsensusSplitEdge(1,edge,alnMap);
        newEdges.insert(newEdges.end(),localNew.begin(),localNew.end());
    }*/
    fprintf(stderr,"Identified %lu new Edges ...\n", newEdges.size());
    return newEdges;
}







