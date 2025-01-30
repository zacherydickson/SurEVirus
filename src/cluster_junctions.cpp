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
    bool operator()(const jRegLabel_t & a, const jRegLabel_t & b) const {
	return (a.compare(b) == 0);
    }
};

struct jRegLabel_HashFunctor {
    size_t operator()(const jRegLabel_t & a) const {
	return std::hash<std::string>{}(a.chr + a.strand + std::to_string(a.pos));
    }
};

struct junctionRegion_t {
    junctionRegion_t() : left(0), right(0), nSplit(0) {}
    junctionRegion_t(size_t pos)
	: left(pos), right(pos), nSplit(0) {}
    junctionRegion_t(const junctionRegion_t & other)
	:   left(other.left), right(other.right), nSplit(other.nSplit),
	    QNameSet(other.QNameSet) {}
    size_t left;
    size_t right;
    size_t nSplit;
    std::unordered_set<std::string> QNameSet;
};


typedef std::vector<jRegLabel_t> jRegLabelVector_t;
typedef std::unordered_map< jRegLabel_t,size_t,jRegLabel_HashFunctor,
			    jRegLabel_EqFunctor>
	    jRegLabelCount_t;
typedef std::unordered_map< jRegLabel_t,bool,jRegLabel_HashFunctor,
			    jRegLabel_EqFunctor>
	    jRegSplitStatus_t;
typedef std::unordered_map< jRegLabel_t,junctionRegion_t,jRegLabel_HashFunctor,
			    jRegLabel_EqFunctor>
	    jRegMap_t;
typedef std::unordered_set< jRegLabel_t, jRegLabel_HashFunctor,
			    jRegLabel_EqFunctor>
	    jRegLabelSet_t;
typedef std::unordered_map< jRegLabel_t, jRegLabelSet_t,
			    jRegLabel_HashFunctor,jRegLabel_EqFunctor> 
	    BestJRegSetMap_t;

//==== GLOBAL VARIABLE DECLARATIONS

static int MinimumReads = 4;
static int SplitBonus = 1;
int MaxInsertSize;
int ReadLength;
std::unordered_set<std::string> VirusNameSet;

//==== FUNCTION DECLARATIONS

void OrderJunctions(const std::string fname, jRegLabelCount_t & labelCount,
		    jRegLabelVector_t & labelVec);
void IdentifyBestJunctions( const std::string fname,BestJRegSetMap_t & bestSetMap);
void ClusterRegions(const std::string fname, const BestJRegSetMap_t & bestSetMap,
		    jRegMap_t & regionMap);
void FilterRegions(jRegMap_t & regionMap);
void OutputRegions(std::string fname, const jRegMap_t & regionMap);

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

    BestJRegSetMap_t regionAssignments;
    IdentifyBestJunctions(candidate_file_name,regionAssignments);
    jRegMap_t regionMap;
    ClusterRegions(candidate_file_name,regionAssignments,regionMap);
    FilterRegions(regionMap);
    OutputRegions(reg_file_name,regionMap);
    fprintf(stderr,"Done - cluster_junctions\n");
}

//==== FUNCTION DEFINITIONS

//Parsed the candidate junctions and counts the instances of each junction
//Then outputs a vector of junction labels ordered from most to least numerous
//As a side effect the split status for the labels is recorded
//Inputs - a string representing the file containing candidate junctions
//	 - a reference to a mapping of counts for labels
//	 - a reference to a vector of labels
//	 - a reference to a mapping of split status for labels
//Output - None, modifies the label count map and vector and split status map
void OrderJunctions(const std::string fname, jRegLabelCount_t & labelCount,
		    jRegLabelVector_t & labelVec,
		    jRegSplitStatus_t & labelFromSplitOnly)
{
    fprintf(stderr,"Ordering Junctions ...\n");
    std::ifstream in(fname);
    std::string chr,qname;
    size_t off, end;
    char score, strand;
    ;
    while (in >> chr >> off >> end >> qname >> score >> strand){
	jRegLabel_t label = {chr,strand,off};
	bool bSplit = (qname[qname.length()-2] == '_');
	if(!labelCount.count(label)) {
	    labelCount[label] = 0;
	    labelVec.push_back(label);
	    labelFromSplitOnly[label] = true;
	}
	labelCount[label]++;
	if(!bSplit) labelFromSplitOnly[label] = false;
    }
    //Sort the vector of labels in descending order by count
    //sort's comparator is true if less (i.e earlier in sorted order),
    //	so we must return true if greater
    std::sort(	labelVec.begin(),labelVec.end(),
		[&labelCount](jRegLabel_t & a, jRegLabel_t & b){
		    return (labelCount[a] > labelCount[b]);
		});

    fprintf(stderr,"Ordered %lu unique Junctions\n",labelVec.size());
}

//First loads all unique junctions and sorts them by order of prevalance
//Then assigns each junction to the best junction within range
//Inputs - a string representing the file containing candidate junctions
//	 - a reference to a map from label to best label set
//Output - None, modifies the best set map
void IdentifyBestJunctions( const std::string fname,
			    BestJRegSetMap_t & bestSetMap) {

    fprintf(stderr,"Assigning junctions ...\n");
    jRegLabelVector_t labelVec;
    jRegLabelCount_t labelCount;
    jRegSplitStatus_t labelFromSplitOnly;
    OrderJunctions(fname,labelCount,labelVec,labelFromSplitOnly);

    //Iterate over junctions from most prevalent to least
    double perc = 0;
    for(size_t i = 0; i < labelVec.size(); ){
	size_t first = i; 
        //Iterate over all junctions with the same score as this one
	do {
	    jRegLabel_t & aLabel = labelVec[i];
	    size_t aRange = (labelFromSplitOnly[aLabel]) ?	ReadLength : 
								MaxInsertSize;
	    size_t left = (aLabel.pos > aRange) ? aLabel.pos - aRange : 0;
	    size_t right = aLabel.pos + aRange;
            //Iterate over all junctions in range of this junction
	    for(size_t pos = left; pos <= right; pos++){
		jRegLabel_t bLabel = {aLabel.chr, aLabel.strand, pos};
		//Check if the potential junction exists
		if(!labelCount.count(bLabel)) continue;
		//Make sure the other junction is mutually in range
		size_t bRange = (labelFromSplitOnly[bLabel]) ?	ReadLength : 
								MaxInsertSize;
		size_t dist = (pos < aLabel.pos) ?  aLabel.pos - pos :
						    pos - aLabel.pos;	
		if(dist > bRange) continue;
		//Check if the bLabel has already been assigned
		if(!bestSetMap.count(bLabel)){
		    bestSetMap[bLabel] = jRegLabelSet_t();
		    bestSetMap[bLabel].insert(aLabel);
		} else if(bestSetMap[bLabel].count(labelVec.at(first))){
		    //bLabel has been assigned to a label with the same
		    //score as aLabel
		    bestSetMap[bLabel].insert(aLabel);
		}
	    }
	    i++; 
	    while( (i * 100.0) / double(labelVec.size()) > perc){
	        perc += 1;
	        fprintf(stderr,"Progress %lu of %lu (~%0.0f%%)\r",i,labelVec.size(),perc);
	    }
	} while(i < labelVec.size() && 
		labelCount[labelVec[i]] == labelCount[labelVec[first]]);
    }
    fprintf(stderr,"\nJunctions Assigned\n");
}


//Parses candidate junctions and assigns each to the most prevalent region within the max
//insert size; the range for paired reads is max insert size, the range for split reads is
//read length
//Inputs - a string representing the file containing candidate junctions
//	 - a reference to an ordered vector of region labels
//	 - a reference to a regionMap to stor results in
//	 - also used global max insert size and read length
//Output - None, modifes the regionMap
void ClusterRegions(const std::string fname, const BestJRegSetMap_t & bestSetMap,
		    jRegMap_t & regionMap) {
     fprintf(stderr,"Clustering Regions ...\n");
     std::ifstream in(fname);
     std::string chr, qname;
     size_t off, end;
     char score, strand;
     //Iterate over candidate junctions
     while (in >> chr >> off >> end >> qname >> score >> strand){
	 jRegLabel_t aLabel = {chr,strand,off};
	 bool bSplit = (qname[qname.length()-2] == '_');
	 for( const auto & bLabel : bestSetMap.at(aLabel) ){
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
	    if(bSplit) regionMap[bLabel].nSplit++;
	    regionMap[bLabel].QNameSet.insert(qname);
	 }
     }
     fprintf(stderr,"Selected %lu unique regions\n",regionMap.size());
}

//Outputs regions which have enough reads assigned,
//  enough is defined as the minimum reads
//  minus the split bonus if any split reads are present
//After this check, removes reads which now map to only host or only virus
//This pair of filters is repeated until no filtering occurs
//Then the final set of junctions is written
//Inputs - a mapping of region labels to regions
//	 - also uses global min reads and split bonus, and viral names
//Output - None, modifes the regionMap
void FilterRegions(jRegMap_t & regionMap){
    //FUTURE:
    //Current implementation does a lot of unecessary filtering
    //One could map each qname to a region and vice versa so that the
    //on each iteration only regions which may have changed are checked
    //Cost - extra complexity and memory
    bool bFilter;
    size_t filterRound= 0;
    do {
	fprintf(stderr,"Filtering Regions - Round %lu ...\n",filterRound++);
	bFilter=false;
	//First Pass Remove regions with too few reads
	//Containers tracking the number of host/viral regions a particular
    	//qname is present in After the first filter pass
    	std::unordered_map<std::string,size_t> hostCount;
    	std::unordered_map<std::string,size_t> viralCount;
    	for( auto it = regionMap.begin(); it != regionMap.end();){
    	    const jRegLabel_t & label = it->first;
    	    const junctionRegion_t & reg = it->second;
    	    size_t effectiveReads = reg.QNameSet.size();
    	    if(reg.nSplit) effectiveReads += SplitBonus;
    	    if(effectiveReads < MinimumReads){ //Fails filter remove
		bFilter=true;
		it = regionMap.erase(it);
    	    } else { // Keep the region
		std::unordered_map<std::string,size_t> * pCountObj =
		    (VirusNameSet.count(label.chr)) ? &viralCount : &hostCount;
		//Increment the host/virus reg count for the qname
		for(const std::string & qname : reg.QNameSet){
		    if(!pCountObj->count(qname)){
			(*pCountObj)[qname] = 0;
		    }
		    (*pCountObj)[qname]++;
		}
    	        it++;
    	    }
    	}
	//If no regions were removed, second pass is unecessary
	if(!bFilter) continue;
	fprintf(stderr,"\tAfter 1st Pass: %lu Regions remain\n",regionMap.size());
	//Second Pass - Remove qnames which now map to only host or 
	//If no qnames get removed the next iteration won't remove any regions
	bFilter=false; 
	for( auto & pair : regionMap){
	    for (   auto it = pair.second.QNameSet.begin();
		    it != pair.second.QNameSet.end(); )
	    {
		//If the qname is associated with both a host and virus
		//region, it may stay
		if(hostCount[*it] && viralCount[*it]){
		    it++;
		} else { // erase the non-junction qname
		    bFilter = true;
		    //Update the split read count for the region if
		    //necessary
		    bool bSplit = ((*it)[it->length()-2] == '_');
		    if(bSplit){
			//The weird syntax is to prevent underflow in a
			// case which shouldn't happen
			pair.second.nSplit += (pair.second.nSplit) ? -1 : 0;
		    }
		    it = pair.second.QNameSet.erase(it);
		}
	    }
	}
    } while(bFilter);

    fprintf(stderr,"Filtered Down to %lu Regions\n",regionMap.size());
}

//Given a map of regions, outputs to a given file
//Inputs - a string represnting the output file name
//	 - a reference to a region map containing regions to output
//Output - None, writes to the povided file
void OutputRegions(std::string fname, const jRegMap_t & regionMap){
    fprintf(stderr,"Printing Regions ...\n");
    std::ofstream output(fname);
    for(auto & pair : regionMap){
    	const jRegLabel_t & label = pair.first;
    	const junctionRegion_t & reg = pair.second;
	std::string qnameStr = *(reg.QNameSet.begin());
	for( const std::string & qname : reg.QNameSet){
	    if(qname == qnameStr) continue;
	    qnameStr = qnameStr + "," + qname;
	}
	output	<< label.chr << "\t" << reg.left << "\t" << reg.right + 1 << "\t"
		<< qnameStr << "\t.\t" << label.strand << "\n";
    }
}

