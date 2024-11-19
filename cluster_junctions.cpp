#include <algorithm>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "config.h"
#include "utils.h"

//==== TYPE DECLARATIONS

struct jRegLabel_t {
    std::string chr;
    char strand;
    size_t pos;
    int compare(const jRegLabel_t & other, int dist = 0) const {
	if(this->chr != other.chr) {
	    return (this->chr < other.chr) ? -1 : 1;
	}
	if(this->strand != other.strand) {
	    return (this->strand == '+') ? 1 : -1;
	}
	if(this->pos > other.pos && this->pos - other.pos > dist) return 1;
	if(this->pos < other.pos && other.pos - this->pos > dist) return -1;
	return 0;
    }
};

struct jRegLabel_EqFunctor {
    bool operator()(jRegLabel_t a, jRegLabel_t b) const {
	return (a.compare(b) == 0);
    }
};

struct jRegLabel_HashFunctor {
    size_t operator()(jRegLabel_t a) const {
	return std::hash<std::string>{}(a.chr + a.strand + std::to_string(a.pos));
    }
};

struct junctionRegion_t {
    junctionRegion_t(size_t pos)
	: left(pos), right(pos), hasSplit(false) {}
    size_t left;
    size_t right;
    bool hasSplit;
    std::unordered_set<std::string> QNameSet;
};


typedef std::vector<jRegLabel_t> jRegLabelVector_t;
typedef std::unordered_map<jRegLabel_t,size_t,jRegLabel_HashFunctor,jRegLabel_EqFunctor>
	    jRegLabelCount_t;
typedef std::unordered_map<jRegLabel_t,junctionRegion_t,jRegLabel_HashFunctor,jRegLabel_EqFunctor>
	    jRegMap_t;

//==== GLOBAL VARIABLE DECLARATIONS

static int MinimumReads = 4;
static int SplitBonus = 1;
int MaxInsertSize;
int ReadLength;
std::unordered_set<std::string> VirusNameSet;

//==== FUNCTION DECLARATIONS

void OrderJunctions(const std::string fname, jRegLabelCount_t & labelCount,
		    jRegLabelVector_t & labelVec);
void ClusterRegions(const std::string fname, const jRegLabelCount_t & labelCount,
		    const jRegLabelVector_t & labelVec, jRegMap_t & regionMap);
void FilterRegions(const std::string ofname, const jRegMap_t & regionMap);

//==== MAIN

//Passes over a candidate junction file twice
//  The first time counts the instances of each junction and puts them into order from most
//  to least prevalent
//  Second pass assigns each read to the most prevalent region in range
//  Finally the regions are filtered by the number of assigned reads and printed
int main(int argc, char* argv[]) {
    //#Parse Inputs
    std::string virus_ref_fname = argv[1];
    std::string workdir = argv[2];
    std::string workspace = argv[3];

    //## Files to be used from the workdir
    std::string stats_file_name = workspace + "/stats.txt";
    std::string candidate_file_name = workdir + "/junction-candidates.bed";
    std::string config_file_name = workdir + "/config.txt";
    //## Output Files
    std::string reg_file_name = workdir + "/region-candidates.bed";

    LoadVirusNames(virus_ref_fname,VirusNameSet);

    MaxInsertSize = parse_stats(stats_file_name).max_is;
    ReadLength = parse_config(config_file_name).read_len;

    jRegLabelVector_t labelVec;
    jRegLabelCount_t labelCount;
    OrderJunctions(candidate_file_name,labelCount,labelVec);
    jRegMap_t regionMap;
    ClusterRegions(candidate_file_name,labelCount,labelVec,regionMap);

    FilterRegions(reg_file_name,regionMap);
}

//==== FUNCTION DEFINITIONS

//Parsed the candidate junctions and counts the instances of each junction
//Then outputs a vector of junction labels ordered from most to least numerous
//Inputs - a string representing the file containing candidate junctions
//	 - a reference to a mapping of counts for labels
//	 - a reference to a vector of labels
//Output - None, modifies the label count map and vector
void OrderJunctions(const std::string fname, jRegLabelCount_t & labelCount,
		    jRegLabelVector_t & labelVec) {
    std::ifstream in(fname);
    std::string chr,qname;
    size_t off, end;
    char score, strand;
    ;
    while (in >> chr >> off >> end >> qname >> score >> strand){
	jRegLabel_t label = {chr,strand,off};
	if(!labelCount.count(label)) {
	    labelCount[label] = 0;
	    labelVec.push_back(label);
	}
	labelCount[label]++;
    }
    //Sort the vector of labels in descending order by count
    //sort's comparator is true if less (i.e earlier in sorted order),
    //	so we must return true if greater
    std::sort(	labelVec.begin(),labelVec.end(),
		[&labelCount](jRegLabel_t & a, jRegLabel_t & b){
		    return (labelCount[a] > labelCount[b]);
		});
}

//Parses candidate junctions and assigns each to the most prevalent region within the max
//insert size; the range for paired reads is max insert size, the range for split reads is
//read length
//Inputs - a string representing the file containing candidate junctions
//	 - a reference to an ordered vector of region labels
//	 - a reference to a regionMap to stor results in
//	 - also used global max insert size and read length
//Output - None, modifes the regionMap
void ClusterRegions(const std::string fname, const jRegLabelCount_t & labelCount,
		    const jRegLabelVector_t & labelVec, jRegMap_t & regionMap){
     std::ifstream in(fname);

     std::string chr, qname;
     size_t off, end;
     char score, strand;
     //Iterate over candidate junctions
     while (in >> chr >> off >> end >> qname >> score >> strand){
	jRegLabel_t aLabel = {chr,strand,off};
	bool bSplit = (qname[qname.length()-2] == '_');
	int range = (bSplit) ? ReadLength : MaxInsertSize;
	bool bFound = false;
	//Search over labels from most to least prevalent
	for (int i = 0; i < labelVec.size() && !bFound; i++){
	    //start at i: test if in range
	    //if i passes, then we can keep checking for ties
	    //add the qname to every tied region
	    int j = i;
	    do {
		const jRegLabel_t bLabel = labelVec[j];
		if(aLabel.compare(bLabel,range) == 0){
		    bFound = true;
		    if(!regionMap.count(bLabel)){
			regionMap[bLabel] = junctionRegion_t(bLabel.pos);
		    }
		    //Extend the bounds of the region
		    if(off < regionMap[bLabel].left){
			regionMap[bLabel].left = off;
		    } else if(off > regionMap[bLabel].right){
			regionMap[bLabel].right = off;
		    }
		    //Indicate if the region has split reads if necessary
		    if(bSplit) regionMap[bLabel].hasSplit = true;
		    regionMap[bLabel].QNameSet.insert(qname);
		}
	    } while(bFound && j+1 < labelVec.size() &&
	    	labelCount.at(labelVec[j+1]) == labelCount.at(labelVec[i]));
	}
     }
}

//Outputs regions which have enough reads assigned, enough is defined as the minimum reads
//minus the split bonus if any split reads are present
//After this check, removes reads which now map to only host or only virus
//Then the final set of junctions is written
//Inputs - a string representing the file to which to write
//	 - a mapping of region labels to regions
//	 - also uses global min reads and split bonus, and viral names
//Output - None, writes to ofname
void FilterRegions(const std::string ofname, const jRegMap_t & regionMap){
    if(MinimumReads == SplitBonus){/*TODO:TRASHME*/}
    //TODO: Write ME
}
