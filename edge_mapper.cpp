#include <iostream>
#include <htslib/sam.h>
#include <memory>
#include <unordered_map>
#include <unordered_set>

#include "config.h"
#include "utils.h"
#include "libs/ssw_cpp.h"
#include "libs/cptl_stl.h"

//==== TYPE DECLARATIONS

//Object describing a read
struct Read_t {
    Read_t(std::string nm) : name(nm) {}
    std::string name;
    std::string hostSegment;
    std::string virusSegment;
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
	return 0;
    }
    std::string to_string(bool seq = false) const {
	std::string str;
	if(!seq){
	    str =   name + std::to_string(hostSegment.size()) +
		    std::to_string(virusSegment.size());
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
	std::string str=    std::to_string('>') + chr + ":" + strand +
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
    double score;
    Edge_t() : hostRegion(nullptr), virusRegion(nullptr), readSet() {}
    //Edge_t(const std::string & regStr, const std::string & readStr);
    Edge_t(Region_pt hostReg, Region_pt virusReg) :
	hostRegion(hostReg), virusRegion(virusReg), readSet() {}
    Edge_t(const Edge_t & other) :
	hostRegion(other.hostRegion), virusRegion(other.virusRegion),
	readSet(other.readSet) {}
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


//==== GLOBAL VARIABLE DECLARATIONS

std::unordered_set<std::string> VirusNameSet;
config_t Config;
stats_t Stats;
//For an alignment to pass it must have a score of at least 30
static StripedSmithWaterman::Filter AlnFilter(true,true,30,32767);
static StripedSmithWaterman::Aligner Aligner(1,4,6,1,false);
int32_t AlnMaskLen;

std::mutex Mtx;

//==== FUNCTION DECLARATIONS

void AlignRead(	int id, const Read_pt & read, const RegionSet_t & regSet,
		AlignmentMap_t & alnMap);
void AlignReads(const Read2RegionsMap_t &regMap,
		AlignmentMap_t & alnMap);
void ConsensusSplitEdge(Edge_t & edge,const AlignmentMap_t & alnMap,
			EdgeVec_t & newEdges);
void DeduplicateEdge(Edge_t & edge,const AlignmentMap_t & alnMap);
void FilterHighInsertReads(Edge_t & edge, const AlignmentMap_t & alnMap);
double GetEdgeScore(const Edge_t & edge,const AlignmentMap_t & alnMap,
		    const ReadSet_t & usedReads, bool bUniqueOnly);
void IdentifyEdgeBreakpoints(Edge_t & edge, const AlignmentMap_t & alnMap);
void IdentifyEdgeSpecificReads(EdgeVec_t & edgeVec);
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
void OrderEdges(EdgeVec_t & edgeVec,const AlignmentMap_t & alnMap);
void OutputEdges(const EdgeVec_t & edgeVec,const AlignmentMap_t & alnMap,
		 const std::string & resFName, const std::string & readDir);

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
    //#Parse Inputs
    std::string virus_ref_file_name = argv[1];
    std::string workdir = argv[2];
    std::string workspace = argv[3];

    //## Files to be used from the workdir
    std::string stats_file_name = workspace + "/stats.txt";
    std::string config_file_name = workdir + "/config.txt";
    std::string region_fasta_file_name = workdir + "/regions.fna";
    std::string read_fasta_file_name = workdir + "/edge_reads.fna";
    std::string edge_file_name = workdir + "edges.tab";
    //## Output Files
    std::string reg_file_name = workdir + "/results.remapped.txt";
    std::string reads_dir = workdir + "/readsx";

    LoadVirusNames(virus_ref_file_name,VirusNameSet);
    Config = parse_config(config_file_name);
    AlnMaskLen =  Config.read_len/2;
    Stats = parse_stats(stats_file_name);

    
    Read2RegionsMap_t  read2regSetMap;
    Region2ReadsMap_t  reg2readSetMap;
    EdgeVec_t  edgeVec;
    LoadData(edge_file_name,region_fasta_file_name,read_fasta_file_name,
	     read2regSetMap,reg2readSetMap,edgeVec);
    
    AlignmentMap_t alnMap;
    AlignReads(read2regSetMap,alnMap);
    
    OrderEdges(edgeVec,alnMap);
    OutputEdges(edgeVec,alnMap,reg_file_name,reads_dir);
}

//==== FUNCTION DEFINITIONS

////Single-Theaded version see other function
//void AlignRead(const Read_pt & read, const RegionSet_t & regSet, AlignmentMap_t & alnMap){
//    AlignRead(0,read,regSet,alnMap);
//}

//Aligns a read to all regions with which it is associated
//  filters the alignments based on a minimum score, and a minimum score
//  relative to the best alignment
//Inputs - an id, used by thead pool
//	 - a read object
//	 - a set of regions
//	 - a constant reference to a mapping of subject query pairs generated
//		alignment objects
//Output - None, modifies the alignment mapping
void AlignRead(int id, const Read_pt & read, const RegionSet_t & regSet, AlignmentMap_t & alnMap){
    std::string hostRC = get_seqrc(read->hostSegment);
    std::string virusRC = get_seqrc(read->virusSegment);
    uint16_t bestScore = 0;
    //Iterate over regions and do the alignments
    std::vector<const SQPair_t *> sqPairVec;
    for( const Region_pt & reg : regSet){
        std::string * query;
        if(reg->isVirus){
	   query = (reg->strand == '+') ?   &(read->virusSegment) :
					    &(virusRC);
        } else {
	   query = (reg->strand == '+') ?   &(read->hostSegment) :
					    &(hostRC);
        }
	Mtx.lock();
        auto res = alnMap.emplace(std::make_pair(  SQPair_t(reg,read),
    				    StripedSmithWaterman::Alignment()));
	Mtx.unlock();
        if(!res.second) continue;
        StripedSmithWaterman::Alignment & aln = res.first->second;
        bool bPass = Aligner.Align(	query->c_str(),AlnFilter,
    				&(aln),AlnMaskLen);
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
    double minScore = 0.75 * bestScore;
    //Iterate over sq pairs and erase any which are below threshold
    for( const SQPair_t * & pPair : sqPairVec){
        if(alnMap.at(*pPair).sw_score < minScore){
	    Mtx.lock();
	    alnMap.erase(*pPair);
	    Mtx.unlock();
        }
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
    ctpl::thread_pool threadPool (Config.threads);
    std::vector<std::future<void>> futureVec;
    for(const auto & pair : regMap){
	std::future<void> future = threadPool.push( AlignRead,
						    std::cref(pair.first),
						    std::cref(pair.second),
						    std::ref(alnMap));
	futureVec.push_back(std::move(future));
    }
    for( auto & future : futureVec){
	future.get();
    }
}

//Splits an edge into a number of edges for each unique consensus
//sequence of reads observed
//Inputs - an edge to process
//	 - an alignment map
//	 - a vector in which to store newly created edges
//Output - None, modifies the newEdges vector and edge object
void ConsensusSplitEdge(Edge_t & edge,const AlignmentMap_t & alnMap,
			EdgeVec_t & newEdges) {
    //TODO:
}

//Identifies Duplicate reads at an edge 
//Inputs - an edge to process
//	 - an alignment map to inform the deduplication
//Output - None, modifies the given object
void DeduplicateEdge(Edge_t & edge,const AlignmentMap_t & alnMap) {
    //TODO:
}

//Identifies read pairs with an apparent insert size which is too large
//and removes them
//Inputs - an edge to process
//	 - an alignment map 
//Output - None, modifies the edge object
void FilterHighInsertReads(Edge_t & edge, const AlignmentMap_t & alnMap){
    //TODO:
}

//Sums the scores on both the host and viral sides for all reads
//supporting the edges
//  may limit to reads which are unique to the edge
//Inputs - an for which to calculate the score
//	 - an alignment map to get scores from
//	 - a set of reads which have been used in other junctions
//	 - a boolean on whether to use only reads which are unique to the
//	    junction
//Output - A double value representing the edge's score
double GetEdgeScore(const Edge_t & edge,const AlignmentMap_t & alnMap,
		    const ReadSet_t & usedReads, bool bUniqueOnly) {
    double score = 0;
    //TODO:
    return score;
}


//Based on reads assigned to an edge, determine where within the host and
//viral regions the actual breakpoints are
//Inputs - an edge to process
//	 - an alignment map
//Output - None, modifies the edge object
void IdentifyEdgeBreakpoints(Edge_t & edge, const AlignmentMap_t & alnMap){
    //TODO:
}

//Determines which if any reads assigned to an edge are only assigned to
//that edge
//Inputs - a vector if edges
//Output - None, modifies the edges
void IdentifyEdgeSpecificReads(EdgeVec_t & edgeVec) {
    //TODO:
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
    std::ifstream in(edgeFName);
    std::string regStr,readStr;
    size_t nRead;
    while(in >> regStr >> readStr >> nRead){
	std::vector<std::string> regionStringVec = strsplit(regStr,':');
	std::vector<std::string> readStringVec = strsplit(readStr,',');
        edgeVec.emplace_back(	regionNameMap.at(regionStringVec[1]),
				regionNameMap.at(regionStringVec[2]));
	for(std::string & rName : readStringVec){
	    rName = rName.substr(0,rName.length()-2);
	    edgeVec.back().readSet.insert(readNameMap.at(rName));
	}
    }
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
    FILE* readsFasta = fopen(readsFName.c_str(),"r");
    kseq_t *seq = kseq_init(fileno(readsFasta));
    while(kseq_read(seq) >= 0){
	std::string name = seq->name.s;
	char segment = name[name.length() - 1];
	if(segment == 'H' || segment == 'V'){
	    name = name.substr(0,name.length()-2);
	}
	if(!nameMap.count(name)){ //Create the read object if necessary
	    Read_pt read = std::make_shared<Read_t>(name);
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
		break;
	    case 'H':
		nameMap[name]->hostSegment = sequence;
		break;
	    case 'V':
		nameMap[name]->virusSegment = sequence;
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
	    fprintf(stderr,"[WARNING] Duplicate Region Sequence ignored");
	}
    }
    kseq_destroy(seq);
    fclose(regionsFasta);
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
    EdgeVec_t newEdges;
    for( Edge_t & edge : edgeVec){
	ConsensusSplitEdge(edge,alnMap,newEdges);
    }
    edgeVec.insert(edgeVec.end(),newEdges.begin(),newEdges.end());
    for( Edge_t & edge : edgeVec){
	IdentifyEdgeBreakpoints(edge,alnMap);
	DeduplicateEdge(edge,alnMap);
	FilterHighInsertReads(edge,alnMap);
    }
    IdentifyEdgeSpecificReads(edgeVec);
    ReadSet_t usedReads; //Required by GetEdgeScore, but not acually needed yet
    std::sort(edgeVec.begin(), edgeVec.end(),
	    [&alnMap,&usedReads](Edge_t & a, Edge_t & b){
		return (GetEdgeScore(a,alnMap,usedReads,true) >
			GetEdgeScore(b,alnMap,usedReads,true));
	    });
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
    //TODO:
}
