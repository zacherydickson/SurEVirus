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

// ==== FUNCTION DECLARATIONS
EdgeVec_t BuildEdges(const RegVecByReadNameMap_t & vRegVecByRead,
                     const RegVecByReadNameMap_t & hRegVecByRead,
                     const Name2RegionMap_t & regMap);
bool CheckForEdgeOverlap(const Edge_t & a, const Edge_t & b);
bool CheckForEdgeContainment(const Edge_t & a, const Edge_t & b);
void PrintUsage();
EdgeVec_t::iterator ProcessBlock(EdgeVec_t & edgeVec, EdgeVec_t::iterator first,
                                 EdgeVec_t::iterator last);
std::string Region2String(const Region_t & reg);
void SubsumeSubEdges(EdgeVec_t & edgeVec);

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
    std::ifstream readFIn(readFileName.c_str());
    std::string readName, regionID;
    std::unordered_map<std::string,std::vector<std::string>> readByReg;
    while (readFIn >> readName >> regionID){
        readByReg[regionID].push_back(readName);
    }
    //Load the regions
    std::ifstream regFIn(regFileName.c_str());
    RegVecByReadNameMap_t vRegVecByRead;
    RegVecByReadNameMap_t hRegVecByRead;
    //Record the loaded regions with values split appart for comparison purposes
    Name2RegionMap_t regMap;
    std::string rname, score, strand;
    size_t off, end;
    while (regFIn >> rname >> off >> end >> regionID >> score >> strand){
        bool isViral = vNameSet.count(rname);
	std::string regStr = rname + ',' + std::to_string(off) + ',' + std::to_string(end) + ',' + strand;
        std::string coordStr = std::to_string(off) + '-' + std::to_string(end);
        regMap[regStr] = std::make_shared<Region_t>(regStr,"",isViral,coordStr);
	//for(const std::string & rID : strsplit(readListStr,',')){
        for(const std::string & rID : readByReg[regionID]){
	    if(isViral){
		vRegVecByRead[rID].push_back(regStr);
	    } else {
		hRegVecByRead[rID].push_back(regStr);
	    }
	}
    }

    //Iterate over reads associated with viral regions
    //std::unordered_map<std::string,std::vector<std::string>> readVecByEdge;
    EdgeVec_t edgeVec = BuildEdges(vRegVecByRead,hRegVecByRead,regMap);

    SubsumeSubEdges(edgeVec); 

    //Sort the edges by descending number of reads
    std::sort(edgeVec.begin(),edgeVec.end(),EdgeSupportGreater);

    //Iterate over the vector and stop when you get to an edge with insufficient support
    for(auto it = edgeVec.begin(); it != edgeVec.end() && it->readSet.size() >= MinReads; it++){
        std::string edgeStr = Region2String(*(it->hostRegion)) + ',' +
                              Region2String(*(it->virusRegion));
        std::cout << edgeStr << '\t';
        bool bFirst = true;
        auto readIt = it->readSet.begin();
        std::cout << (*readIt)->name;
        while(++readIt != it->readSet.end()){
            std::cout << ',' << (*readIt)->name;
        }
        std::cout << '\t' << it->readSet.size() << '\n';
    }

    ////Define a comparator for sorting edges
    //struct compfunctor {
    //    const std::unordered_map<std::string,std::vector<std::string>> * prVbE;
    //    compfunctor(std::unordered_map<std::string,std::vector<std::string>> * p)
    //        : prVbE(p) {}
    //    bool operator() (const std::string & a,const std::string & b) const {
    //        return this->prVbE->at(a).size() < this->prVbE->at(b).size();
    //    }
    //};
    ////Sort the edges
    //std::priority_queue<std::string,std::vector<std::string>,compfunctor> pq((compfunctor(&readVecByEdge)));
    //for(const auto & pair : readVecByEdge){
    //    if(pair.second.size() >= MinReads)
    //        pq.push(pair.first);
    //}
    ////Output the edges in sorted order
    //for(;!pq.empty(); pq.pop()){
    //    const std::string & edge = pq.top();
    //    const std::vector<std::string> & rIDVec = readVecByEdge[edge];
    //    std::string s = rIDVec.front();
    //    for(int i = 2; i < rIDVec.size(); i++){
    //        s += ',' + rIDVec[i];
    //    }
    //    std::cout << edge << '\t' << s << '\t' << rIDVec.size() << "\n";
    //}
    //for(const auto & pair : readVecByEdge){
    //    if(pair.second.size() < MinReads) continue;
    //    
    //}
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
EdgeVec_t BuildEdges(const RegVecByReadNameMap_t & vRegVecByRead,
                     const RegVecByReadNameMap_t & hRegVecByRead,
                     const Name2RegionMap_t & regMap)
{
    std::unordered_map<std::string,size_t> edgeIdxByName;
    EdgeVec_t edgeVec;
    for(const auto & pair : vRegVecByRead){
	const std::string & rID = pair.first;
	//Iterate over human regions associated with this read
	if(!hRegVecByRead.count(rID)) continue;
        Read_pt read = std::make_shared<Read_t>(rID,false,false); 
	for(const std::string & hReg : hRegVecByRead.at(rID)){
	    for(const std::string & vReg : pair.second){
		std::string edgeStr = hReg + ':' + vReg;
		//readVecByEdge[edgeStr].push_back(rID);
                if(!edgeIdxByName.count(edgeStr)){
                    edgeIdxByName[edgeStr] = edgeVec.size();
                    edgeVec.emplace_back(regMap.at(hReg),regMap.at(vReg));
                }
                size_t edgeIdx = edgeIdxByName[edgeStr];
                edgeVec[edgeIdx].addRead(read);
	    }
	}
    }
    return edgeVec;
}

//Given two edges, checks if both the host regions and virus regions overlap
//For an overlap check the order doesn't matter so we flip the regions so that 
//  a is always less than b which
bool CheckForEdgeOverlap(const Edge_t & a, const Edge_t & b) {
    const Region_pt & hostA = (a.hostRegion->compare(*b.hostRegion) == 1) ?
                              b.hostRegion : a.hostRegion;
    const Region_pt & hostB = (a.hostRegion->compare(*b.hostRegion) == 1) ?
                              b.hostRegion : a.hostRegion;
    const Region_pt & virusA = (a.virusRegion->compare(*b.virusRegion) == 1) ?
                               b.virusRegion : a.virusRegion;
    const Region_pt & virusB = (a.virusRegion->compare(*b.virusRegion) == 1) ?
                               b.virusRegion : a.virusRegion;
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

//Iterator over the upper triangle of comparisons and subsume any sub edges
EdgeVec_t::iterator ProcessBlock(EdgeVec_t & edgeVec, EdgeVec_t::iterator first,
                                 EdgeVec_t::iterator last)
{
    for(; std::next(first,1) != last; first++){
        for(auto it = std::next(first); it != last;){
            if(CheckForEdgeContainment(*first,*it)){
                //Add the subedge's reads to this edge
                for(const Read_pt & read : it->readSet){
                    first->addRead(read);
                }
                //The erase will invalidate last so we need to be able to
                // reconstruct it by advancing from the new it
                //first shouldn't be invalidated as it is before it
                auto distance = std::distance(it,last);
                it = edgeVec.erase(it);
                last = std::next(it,distance-1);
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
void SubsumeSubEdges(EdgeVec_t & edgeVec){
    //First the Edges must be sorted by:
    //   hostChr+strand -> hostStart -> hostEnd -> virusChr+strand -> virusStart -> virusEnd
    std::sort(edgeVec.begin(),edgeVec.end(),EdgePositionLess);
    //Next we iterate over the edges whenever we encounter the end of a block of overlapping edges we process
    //that block
    auto blockStart = edgeVec.begin();
    for(auto it = edgeVec.begin()+1; it != edgeVec.end(); it++){
        if(!CheckForEdgeOverlap(*std::prev(it,1),*it)){
            it = ProcessBlock(edgeVec,blockStart,it);
            blockStart=it;
        }
    }
    ProcessBlock(edgeVec,blockStart,edgeVec.end());
}
