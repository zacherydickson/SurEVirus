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
    std::string r1;
    std::string r2;
    int compare(const Read_t other) const {
	if(this->name != other.name) {
	    return (this->name < other.name) ? -1 : 1;
	}
	if(this->r1.size() != other.r1.size()){
	    return (this->r1.size() < other.r1.size()) ? -1 : 1;
	}
	if(this->r2.size() != other.r2.size()){
	    return (this->r2.size() < other.r2.size()) ? -1 : 1;
	}
	return 0;
    }
    std::string to_string(bool seq = false) const {
	std::string str;
	if(!seq){
	    str =   name + std::to_string(r1.size()) +
		    std::to_string(r2.size());
	} else {
	    if(r1.size()){
		str += std::to_string('>') + name + '_' + std::to_string(1);
		str += std::to_string('\n') + r1;
	    }
	    if(r2.size()){
		str += std::to_string('>') + name + '_' + std::to_string(2);
		str += std::to_string('\n') + r2;
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

//Object describing a genomic location containing a
//candidate junction
struct Region_t {
    Region_t(kseq_t *seq) : sequence(seq->seq.s){
	std::string name = seq->name.s;
	std::vector<std::string> fields = strsplit(name,',');
	this->chr = fields[0];
	this->left = std::stoul(fields[1]);
	this->right = std::stoul(fields[2]);
	this->strand = fields[3][0];
    }
    Region_t(const Region_t & other) :
	left(other.left), right(other.right), chr(other.chr), 
	strand(other.strand), sequence(other.sequence){}
    size_t left, right;
    std::string chr;
    char strand;
    std::string sequence;
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

//Object associating a pair of regions and the reads spanning the pair
struct Edge_t {
    Region_pt hostRegion;
    Region_pt virusRegion;
    ReadSet_t readSet;
    double score;
    Edge_t() : hostRegion(nullptr), virusRegion(nullptr), readSet() {}
    //Edge_t(const std::string & regStr, const std::string & readStr);
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

config_t Config;
stats_t Stats;
//For an alignment to pass it must have a score of at least 30
static StripedSmithWaterman::Filter AlnFilter(true,true,30,32767);
static StripedSmithWaterman::Aligner Aligner(1,4,6,1,false);

//==== FUNCTION DECLARATIONS


void AlignReads(const Read2RegionsMap_t &regMap,
		AlignmentMap_t & alnMap);
void LoadData(	const std::string & edgeFName,
		const std::string & regionsFName,
		const std::string & readsFName,
	      	Read2RegionsMap_t & read2regSetMap,
	      	Region2ReadsMap_t & reg2readSetMap,
	      	EdgeVec_t & edgeVec);
void LoadEdges(	std::string edgeFName,
	      	Read2RegionsMap_t & read2regSetMap,
	      	Region2ReadsMap_t & reg2readSetMap,
	      	EdgeVec_t & edgeVec);
void LoadReadSeq(   const std::string & readsFName,
		    Read2RegionsMap_t & read2regSetMap);
void LoadRegionSeq( const std::string & regionsFName,
		    Region2ReadsMap_t & reg2readSetMap);

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
    std::string workdir = argv[1];
    std::string workspace = argv[2];

    //## Files to be used from the workdir
    std::string stats_file_name = workspace + "/stats.txt";
    std::string config_file_name = workdir + "/config.txt";
    std::string region_fasta_file_name = workdir + "/regions.fna";
    std::string read_fasta_file_name = workdir + "/edge_reads.fna";
    std::string edge_file_name = workdir + "edges.tab";
    //## Output Files
    std::string reg_file_name = workdir + "/results.remapped.txt";
    std::string reads_dir = workdir + "/readsx";

    Config = parse_config(config_file_name);
    Stats = parse_stats(stats_file_name);

    
    Read2RegionsMap_t  read2regSetMap;
    Region2ReadsMap_t  reg2readSetMap;
    EdgeVec_t  edgeVec;
    LoadData(edge_file_name,region_fasta_file_name,read_fasta_file_name,
	     read2regSetMap,reg2readSetMap,edgeVec);
    
    AlignmentMap_t alnMap;
    AlignReads(read2regSetMap,alnMap);
    //Need to figure out how best to orient sequences for alignment 
   
    //Identify and remove PCR duplicates
    //Restrict to only reads which agree with the consensus
    //Score Edges and process from best to worst score
}

//==== FUNCTION DEFINITIONS

//Iterate over each read and align it to all associated regions
//  Filters out low scoring alignments
void AlignReads(const Read2RegionsMap_t &regMap,
		AlignmentMap_t & alnMap) {
    //Need to align all reads to associated regions
    //Accept alignments with scores higher than 30
    //Reject alignments whose score is lower than 75% of best score
    // both the host and viral sides
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
    LoadRegionSeq(regionsFName,reg2readSetMap);
    LoadReadSeq(bamFName,read2regSetMap);
    LoadEdges(edgeFName,read2regSetMap,reg2readSetMap,edgeVec);
}

//Parses edge file to load the regions, reads, and their associations
//Inputs - a string representing the edge file name
//	 - a reference to mapping between reads and region sets
//	 - a reference to a mpaaing from regions to read sets
//	 - a reference to a vector of edges
//Output - None, modifies to provided references
void LoadEdges(	std::string edgeFName,
	      	Read2RegionsMap_t & read2regSetMap,
	      	Region2ReadsMap_t & reg2readSetMap,
	      	EdgeVec_t & edgeVec)
{
    std::ifstream in(edgeFName);
    std::string regStr,readStr;
    size_t nRead;
    while(in >> regStr >> readStr >> nRead){
	//TODO: Need a connection from read and region names to the
	//objects, so that the edge can be properly formed
        edgeVec.emplace_back(regStr,readStr);
        Edge_t & edge = edgeVec.back();
        //Iterate over the host and virus regions
        for(Region_pt * curRegion : {&(edge.hostRegion),&(edge.virusRegion)}){
	    //Iterate over reads
            for(const Read_pt & read : edge.readSet){
        	//Add the read to the region map
        	reg2readSetMap.at(*curRegion).insert(read);
        	//Add the region to the read map
        	read2regSetMap.at(read).insert(*curRegion);
            }
        }
    }
}

//Parses a fasta file containing reads and stores their sequences
//The header of the entries is parsed for a trailing _1 or _2 to
//seperate fwd and reverse reads
//Inputs - a string representing the bam file name
//	 - a reference to a mapping from reads to regions sets
//Output - None, modifies the key objects of the map
void LoadReadSeq(   const std::string & readsFName,
		    Read2RegionsMap_t & read2regSetMap)
{
    FILE* readsFasta = fopen(readsFName.c_str(),"r");
    kseq_t *seq = kseq_init(fileno(readsFasta));
    std::unordered_map<std::string,Read_pt> nameMap;
    while(kseq_read(seq) >= 0){
	std::string name = seq->name.s;
	char segment = name[name.length() - 1];
	name = name.substr(0,name.length()-2);
	if(!nameMap.count(name)){ //Create the read object if necessary
	    Read_pt read = std::make_shared<Read_t>(name);
	    nameMap.insert(std::make_pair(name,read));
	}
	//Add the sequence to the read
	switch (segment) {
	    case '1':
		nameMap[name]->r1 = std::string(seq->seq.s);
		break;
	    case '2':
		nameMap[name]->r2 = std::string(seq->seq.s);
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
//Inputs - a string representing the fasta file name
//	 - a reference to a mappinf from regions to read sets
//Output - None, modifes the key objects of the map
void LoadRegionSeq( const std::string & regionsFName,
		    Region2ReadsMap_t & reg2readSetMap)
{
    FILE* regionsFasta = fopen(regionsFName.c_str(),"r");
    kseq_t *seq = kseq_init(fileno(regionsFasta));
    while(kseq_read(seq) >= 0){
	std::string name = seq->name.s;
	size_t pos = name.find(':',1);
	if(pos == std::string::npos) continue;
	seq->name.l = pos; //set the current end of the name earlier
	seq->name.s[pos] = '\0';
	Region_pt reg = std::make_shared<Region_t>(seq);
	auto res = reg2readSetMap.emplace(std::make_pair(reg,ReadSet_t()));
	if(!res.second){
	    fprintf(stderr,"[WARNING] Duplicate Region Sequence ignored");
	}
    }
    kseq_destroy(seq);
    fclose(regionsFasta);
}

////Constructor for Edge_t structure from region and read strings
//Edge_t::Edge_t(const std::string & regStr, const std::string & readStr) :
//    Edge_t()
//{
//    this->parseRegString(regStr);
//    this->parseReadString(readStr);
//}
//
////Internal Edge_t function for parsing a region string
////Populates the host and virus regions members
//void Edge_t::parseRegString(const std::string & regStr) {
//    //TODO
//}
//
////Internal Edge_t function for parsing region string
////Populates the readSet member
//void Edge_t::parseReadString(const std::string & readStr) {
//    //TODO
//}
