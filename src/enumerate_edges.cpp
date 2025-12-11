#include <algorithm>
#include <iostream>
#include <iterator>
#include <fstream>
#include <queue>
#include <string>
#include "str_utils.h"
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include "utils.h"
#include "edge_utils.h"

// ==== TYPE DEFINITIONS

typedef std::unordered_map<std::string,std::vector<std::string>> RegVecByReadNameMap_t;
typedef std::unordered_map<std::string,std::vector<std::string>> ReadVecByRegMap_t;

// ==== FUNCTION DECLARATIONS
EdgeList_t BuildEdges(const RegVecByReadNameMap_t & vRegVecByRead,
                     const RegVecByReadNameMap_t & hRegVecByRead,
                     const Name2RegionMap_t & regMap);
bool CheckForEdgeOverlap(const Edge_t & a, const Edge_t & b);
bool CheckForEdgeContainment(const Edge_t & a, const Edge_t & b);
void PrintUsage();
EdgeList_t::iterator ProcessBlock(EdgeList_t & edgeList, EdgeList_t::iterator first,
                                 EdgeList_t::iterator last);
std::string Region2String(const Region_t & reg);
ReadVecByRegMap_t LoadRegionReadAssociations(const std::string & fname);
void LoadRegions(   const std::string & fname,
                    const std::unordered_set<std::string> & vNameSet,
                    const ReadVecByRegMap_t & readByReg,
                    RegVecByReadNameMap_t & vRegVecByRead,
                    RegVecByReadNameMap_t & hRegVecByRead,
                    Name2RegionMap_t & regMap);
void SubsumeSubEdges(EdgeList_t & edgeList);

// === GLOBAL VARIABLES

const size_t MinReads = 3;

struct {
    bool operator()(const Edge_t & a, const Edge_t & b ) const {
        int hostCmp = a.hostRegion->compare(*b.hostRegion);
        if(hostCmp != 0) return (hostCmp == -1);
        int virusCmp = a.virusRegion->compare(*b.virusRegion);
        return (virusCmp == -1);
    }
} EdgePositionLess;

struct {
    bool operator()(const Edge_t & a, const Edge_t & b ) const {
        return (a.readSet.size() > b.readSet.size()); 
    }
} EdgeSupportGreater;

// ==== MAIN

int main(int argc, const char* argv[]) {
    //Handle Input
    if(argc < 3){
	PrintUsage();
	return 1;
    }
    std::string vNameFileName = argv[1];
    std::string regFileName = argv[2];
    std::string readFileName = argv[3];
    //LoadVirus Names
    std::unordered_set<std::string> vNameSet;
    LoadVirusNames(vNameFileName,vNameSet);
    //Load the region - read associations
    ReadVecByRegMap_t readByReg = LoadRegionReadAssociations(readFileName);
    //Object for Recording the loaded regions with values split apart for comparison purposes
    Name2RegionMap_t regMap;
    //Load the regions
    RegVecByReadNameMap_t vRegVecByRead;
    RegVecByReadNameMap_t hRegVecByRead;
    LoadRegions(regFileName,vNameSet,readByReg,vRegVecByRead,hRegVecByRead,regMap);

    //Iterate over reads associated with viral regions
    //std::unordered_map<std::string,std::vector<std::string>> readVecByEdge;
    EdgeList_t edgeList = BuildEdges(vRegVecByRead,hRegVecByRead,regMap);

    //Check for any edges wholly contained within other edges
    SubsumeSubEdges(edgeList); 

    //Sort the edges by descending number of reads
    edgeList.sort(EdgeSupportGreater);
//#std::sort(edgeList.begin(),edgeList.end(),EdgeSupportGreater);

    fprintf(stderr,"Filtering Edges ...\n");
    size_t nEdges = 0;
    //Iterate over the list and stop when you get to an edge with insufficient support
    for(auto it = edgeList.begin(); it != edgeList.end() && it->readSet.size() >= MinReads; it++){
        nEdges++;
        std::string edgeStr = Region2String(*(it->hostRegion)) + ':' +
                              Region2String(*(it->virusRegion));
        std::cout << edgeStr << '\t';
        auto readIt = it->readSet.begin();
        std::cout << (*readIt)->name;
        while(++readIt != it->readSet.end()){
            std::cout << ',' << (*readIt)->name;
        }
        std::cout << '\t' << it->readSet.size() << '\n';
    }
    fprintf(stderr,"Enumerated %zu edges with sufficent support ...\n",nEdges);

    return 0;
}

// ==== FUNCTION DEFINITIONS

//Iterates over reads assocaiated with any viral region and checks if they are
// associated with a host region. If so, then generates all combinations of host
// and virus regions that this read is associated with. For each unique combination
// an edge is created and the read is added to each each edge it supports
//Inputs - a mapping from read names to vectors of viral region names
//       - a mapping from read names to vectors of host region names
//Output - a vector of Edge_t objects
EdgeList_t BuildEdges(const RegVecByReadNameMap_t & vRegVecByRead,
                     const RegVecByReadNameMap_t & hRegVecByRead,
                     const Name2RegionMap_t & regMap)
{
    fprintf(stderr,"Building Edges ...\n");
    std::unordered_map<std::string,Edge_t*> edgePtrByName;
    EdgeList_t edgeList;
    for(const auto & pair : vRegVecByRead){
	const std::string & rID = pair.first;
	//Iterate over human regions associated with this read
	if(!hRegVecByRead.count(rID)) continue;
        Read_pt read = std::make_shared<Read_t>(rID,false,false); 
	for(const std::string & hReg : hRegVecByRead.at(rID)){
	    for(const std::string & vReg : pair.second){
		std::string edgeStr = hReg + ':' + vReg;
                if(!edgePtrByName.count(edgeStr)){
                    edgeList.emplace_back(regMap.at(hReg),regMap.at(vReg));
                    edgePtrByName[edgeStr] = &(edgeList.back());
                }
                edgePtrByName[edgeStr]->addRead(read);
	    }
	}
    }
    fprintf(stderr,"Built %zu edges\n",edgePtrByName.size());
    return edgeList;
}

//Given two edges, checks if both the host regions and virus regions overlap
//For an overlap check the order doesn't matter so we flip the regions so that 
//  a is always less than b which
bool CheckForEdgeOverlap(const Edge_t & a, const Edge_t & b) {
    const Region_pt & hostA = (a.hostRegion->compare(*b.hostRegion) == 1) ?
                              b.hostRegion : a.hostRegion;
    const Region_pt & hostB = (a.hostRegion->compare(*b.hostRegion) == 1) ?
                              a.hostRegion : b.hostRegion;
    const Region_pt & virusA = (a.virusRegion->compare(*b.virusRegion) == 1) ?
                               b.virusRegion : a.virusRegion;
    const Region_pt & virusB = (a.virusRegion->compare(*b.virusRegion) == 1) ?
                               a.virusRegion : b.virusRegion;
    //Check if the strand and chromosomes match
    if(hostA->chr + hostA->strand != hostB->chr + hostB->strand ||
       virusA->chr + virusA->strand != virusB->chr + virusB->strand)
    {
        return false;
    }
    //If both B's start at or before the end of Both A's they must overlap
    if(hostA->seqRight >= hostB->seqLeft && 
       virusA->seqRight >= virusB->seqLeft)
    {
        return true;
    }
    return false;
}

//Given two edges checks if both the host and virus regions of B are subregions
//  of their respective counterparts in A
bool CheckForEdgeContainment(const Edge_t & a, const Edge_t & b) {
    const Region_pt & hostA = a.hostRegion;
    const Region_pt & hostB = b.hostRegion;
    const Region_pt & virusA = a.virusRegion;
    const Region_pt & virusB = b.virusRegion;
    //Check if the strand and chromosomes match
    if(hostA->chr + hostA->strand != hostB->chr + hostB->strand ||
       virusA->chr + virusA->strand != virusB->chr + virusB->strand)
    {
        return false;
    }
    if( hostA->seqLeft <= hostB->seqLeft &&
        hostA->seqRight >= hostB->seqRight &&
        virusA->seqLeft <= virusB->seqLeft &&
        virusA->seqRight >= virusB->seqRight)
    {
        return true;
    }
    return false;
}

//Open a file of which details which reads are associated with a given region
//Inputs - a path to a file with 2 columns, read name and regionID
//Output - a ReadVecByReagMap_t object
ReadVecByRegMap_t LoadRegionReadAssociations(const std::string & fname) {
    fprintf(stderr,"Loading read-region associations from %s ...\n", fname.c_str());
    std::ifstream readFIn(fname.c_str());
    std::string readName, regionID;
    ReadVecByRegMap_t readByReg;
    size_t nAssoc = 0;
    while (readFIn >> readName >> regionID){
        readByReg[regionID].push_back(readName);
        nAssoc++;
    }
    fprintf(stderr,"Loaded %zu associations across %zu regions\n",
            nAssoc,readByReg.size());
    return readByReg;
}

//Given a file of regions, the set of references which are viral, and a mapping
// between regions and reads, loads the regions associated with viral and host
//Notes: fname is a path to a file with 6 columns:
// the name of the read, the offset and end within the region,
// the name of the region, the score of the mapping (unused), and the strand
void LoadRegions(   const std::string & fname,
                    const std::unordered_set<std::string> & vNameSet,
                    const ReadVecByRegMap_t & readByReg,
                    RegVecByReadNameMap_t & vRegVecByRead,
                    RegVecByReadNameMap_t & hRegVecByRead,
                    Name2RegionMap_t & regMap)
{
    fprintf(stderr,"Loading regions from %s ...\n", fname.c_str());
    std::string rname, score, strand, regionID;
    size_t off, end;
    std::ifstream regFIn(fname.c_str());
    while (regFIn >> rname >> off >> end >> regionID >> score >> strand){
        bool isViral = vNameSet.count(rname);
	std::string regStr = rname + ',' + std::to_string(off) + ',' + std::to_string(end) + ',' + strand;
        std::string coordStr = std::to_string(off) + '-' + std::to_string(end);
        regMap[regStr] = std::make_shared<Region_t>(regStr,"",isViral,coordStr);
	//for(const std::string & rID : strsplit(readListStr,',')){
        for(const std::string & rID : readByReg.at(regionID)){
	    if(isViral){
		vRegVecByRead[rID].push_back(regStr);
	    } else {
		hRegVecByRead[rID].push_back(regStr);
	    }
	}
    }
    fprintf(stderr,"Loaded %zu regions associated with %zu viral Reads and %zu host Reads\n",
            regMap.size(), vRegVecByRead.size(), hRegVecByRead.size());
}

//Iterator over the upper triangle of comparisons and subsume any sub edges
EdgeList_t::iterator ProcessBlock(EdgeList_t & edgeList, EdgeList_t::iterator first,
                                 EdgeList_t::iterator last)
{
    for(; first != last && std::next(first,1) != last; first++){
        for(auto it = std::next(first); it != last;){
            if(CheckForEdgeContainment(*first,*it)){
                //Add the subedge's reads to this edge
                for(const Read_pt & read : it->readSet){
                    first->addRead(read);
                }
                //The erase will invalidate last so we need to be able to
                // reconstruct it by advancing from the new it
                //first shouldn't be invalidated as it is before it
                //auto distance = std::distance(it,last);
                it = edgeList.erase(it);
                //last = std::next(it,distance-1);
                //distance = std::distance(it,last);
            } else {
                it++;
            }
        }
    }
    return last;
}

//Prints the Usage Message
//Inputs - None
//Outputs - None, writes to stderr
void PrintUsage() {
    std::cerr	<< "===Description\n"
	<< "\tGive a fasta file of viral seqs, a file containing region\n"
	<< " candidates, and a file connecting readIDs to regions,\n"
        << " outputs a set of edges between viral and non-viral regions\n"
	<< "===Usage\n"
	<< "\tenumerate_edges virus.fna region_candidates.bed read_regionMap.tab\n"
	<< "===ARGUMENTS\n"
	<< "\tvirus PATH\tPath to a fasta file containing viral seqs\n"
	<< "\tregion_candidates PATH\tPath to a tab delim file describing\n"
	<< "\t regions; BED formatted with a comma sep list of reads in the"
	<< "\t name field\n"
        << "\tread_regionMap.tab PATH\tPath to a tab delim file with two columns:\n"
        << "\t readID and regionID; any ID may appear multiple times, but each\n"
        << "\t combination will be unique\n"
	<< "===Output\n"
	<< "\tOutput is tab delim with columns of edge, readList, and nReads\n"
	<< "\tThe output will be sorted in descending order of nReads\n"
	;
}


std::string Region2String(const Region_t & reg){
    return reg.chr + ',' + std::to_string(reg.left) + ',' +
           std::to_string(reg.right) + ',' + reg.strand;
}

//Identifies edges for which both regions are enclosed inside of the regions of
// some other edge, the reads associated with these subedges are added to the
// edge with large regions and the subedge is deleted
//Inputs - a reference to a vector of edges
//Output - None, modifies the edge vector
void SubsumeSubEdges(EdgeList_t & edgeList){
    fprintf(stderr,"Subsuming edges...\n");
    //First the Edges must be sorted by:
    //   hostChr+strand -> hostStart -> hostEnd -> virusChr+strand -> virusStart -> virusEnd
    edgeList.sort(EdgePositionLess);
//#std::sort(edgeList.begin(),edgeList.end(),EdgePositionLess);
    //Next we iterate over the edges whenever we encounter the end of a block of overlapping edges we process
    //that block
    auto blockStart = edgeList.begin();
    for(auto it = std::next(edgeList.begin()); it != edgeList.end(); it++){
        if(!CheckForEdgeOverlap(*std::prev(it,1),*it)){
            it = ProcessBlock(edgeList,blockStart,it);
            blockStart=it;
        }
    }
    ProcessBlock(edgeList,blockStart,edgeList.end());
    fprintf(stderr,"After subsumption, %zu edges remain...\n",edgeList.size());
}
