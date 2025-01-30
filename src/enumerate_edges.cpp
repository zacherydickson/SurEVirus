#include <iostream>
#include <fstream>
#include <queue>
#include <string>
#include "str_utils.h"
#include <unordered_map>
#include <unordered_set>
//#include "utils.h"

// ==== TYPE DEFINITIONS


// ==== FUNCTION DECLARATIONS
void PrintUsage();

// ==== MAIN

int main(int argc, const char* argv[]) {
    //Handle Input
    if(argc < 2){
	PrintUsage();
	return 1;
    }
    std::string vNameFileName = argv[1];
    std::string regFileName = argv[2];
    //LoadVirus Names
    std::unordered_set<std::string> vNameSet;
    LoadVirusNames(vNameFileName,vNameSet);
    //Load the reagion - read associations
    std::ifstream regFIn(regFileName.c_str());
    std::unordered_map<std::string,std::vector<std::string>> vRegVecByRead;
    std::unordered_map<std::string,std::vector<std::string>> hRegVecByRead;
    std::string rname, off, end, readListStr, score, strand;
    while (regFIn >> rname >> off >> end >> score >> strand){
	std::string regStr = rname + ',' + off + ',' + end + ',' + strand;
	for(const std::string & rID : strsplit(readListStr,',')){
	    if(vNameSet.count(rname)){
		vRegVecByRead[rID].push_back(regStr);
	    } else {
		hRegVecByRead[rID].push_back(regStr);
	    }
	}
    }
    //Iterate over reads associated with viral regions
    std::unordered_map<std::string,std::vector<std::string>> readVecByEdge;
    for(const auto & pair : vRegVecByRead){
	const std::string & rID = pair.first;
	//Iterate over human regions associated with this read
	if(!hRegVecByRead.count(rID)) continue;
	for(const std::string & hReg : hRegVecByRead.at(rID)){
	    for(const std::string & vReg : pair.second){
		std::string edge = hReg + ':' + vReg;
		readVecByEdge[edge].push_back(rID);
	    }
	}
    }
    //Define a comparator for sorting edges
    struct compfunctor {
	const std::unordered_map<std::string,std::vector<std::string>> * prVbE;
	compfunctor(std::unordered_map<std::string,std::vector<std::string>> * p)
	    : prVbE(p) {}
	bool operator() (const std::string & a,const std::string & b) const {
	    return this->prVbE->at(a).size() > this->prVbE->at(b).size();
	}
    };
    //Sort the edges
    std::priority_queue<std::string,std::vector<std::string>,compfunctor> pq((compfunctor(&readVecByEdge)));
    for(const auto & pair : readVecByEdge){
	if(pair.second.size() > 3)
	    pq.push(pair.first);
    }
    //Output the edges in sorted order
    for(;!pq.empty(); pq.pop()){
	const std::string & edge = pq.top();
	const std::vector<std::string> & rIDVec = readVecByEdge[edge];
	std::string s = rIDVec.front();
	for(int i = 2; i < rIDVec.size(); i++){
	    s += rIDVec[i];
	}
	std::cout << edge << '\t' << s << '\t' << rIDVec.size();
    }
    for(const auto & pair : readVecByEdge){
	if(pair.second.size() < 3) continue;
	
    }
    return 0;
}

// ==== FUNCTION DEFINITIONS

//Prints the Usage Message
//Inputs - None
//Outputs - None, writes to stderr
void PrintUsage() {
    std::cerr	<< "===Description\n"
	<< "\tGive a list file of virus names, and file containing region\n"
	<< " candidates, outputs a set of edges between viral and non-viral\n"
	<< " regions\n"
	<< "===Usage\n"
	<< "\tenumerate_edges virusNames.list region_candidates.bed\n"
	<< "===ARGUMENTS\n"
	<< "\tvirusNames PATH\tPath to a fasta file containing viral seqs\n"
	<< "\tregion_candidates PATH\tPath to a tab delim file describing\n"
	<< "\t regions; BED formatted with a comma sep list of reads in the"
	<< "\t name field\n"
	<< "===Output\n"
	<< "\tOutput is tab delim with columns of edge, readList, and nReads\n"
	<< "\tThe output will be sorted in descending order of nReads\n"
	;
}
