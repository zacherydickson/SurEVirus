#ifndef SURVERYOR_EDGE_UTILS_H
#define SURVERYOR_EDGE_UTILS_H

#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <htslib/kseq.h>
#include <vector>

#include "utils.h"
#include <ssw.h>
#include <ssw_cpp.h>

struct Read_t;

typedef std::shared_ptr<Read_t> Read_pt;

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
    Read_pt mate = nullptr;
    int compare(const Read_t & other) const {
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
    Region_t(kseq_t *seq, bool bVirus,const std::string & coordStr) : 
        sequence(seq->seq.s), isVirus(bVirus)
    {
        std::string name = seq->name.s;
        std::vector<std::string> fields = strsplit(name,',');
        std::vector<std::string> coords = strsplit(coordStr,'-');
        this->chr = fields[0];
        this->left = std::stoul(fields[1]);
        this->right = std::stoul(fields[2]);
        this->seqLeft = std::stoul(coords[0]);
        this->seqRight = std::stoul(coords[1]);
        this->strand = fields[3][0];
        for (auto & c: sequence) c = (char)toupper(c);
    }
    Region_t(const Region_t & other) :
        left(other.left), right(other.right),
        seqLeft(other.seqLeft), seqRight(other.seqRight),
        chr(other.chr), strand(other.strand), sequence(other.sequence) {}
    size_t left, right; //These are labels defining the region
    //These are the actual genomic corrdinates of the sequence associated with this region
    size_t seqLeft, seqRight; 
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

struct Edge_t; 

typedef std::vector<Edge_t> EdgeVec_t;

//Object associating a query read with a subject region
struct SQPair_t {
    SQPair_t(const Region_pt & s, const Read_pt & q) : subject(s), query(q) {}
    SQPair_t(const SQPair_t & other) :
        subject(other.subject), query(other.query) {}
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
        return        this->subject->to_string() + " vs " +
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
    double lastScore = -1;
    Edge_t() :        hostRegion(nullptr), virusRegion(nullptr), readSet(),
                uniqueReadSet(), hostOffset(0), virusOffset(0) {}
    //Edge_t(const std::string & regStr, const std::string & readStr);
    Edge_t(Region_pt hostReg, Region_pt virusReg) :
            hostRegion(hostReg), virusRegion(virusReg), readSet(),
            uniqueReadSet(), hostOffset(0), virusOffset(0) {}
    Edge_t(const Edge_t & other) :
        hostRegion(other.hostRegion), virusRegion(other.virusRegion),
        readSet(other.readSet), uniqueReadSet(other.uniqueReadSet),
        hostOffset(other.hostOffset), virusOffset(other.virusOffset), 
        nSplit(other.nSplit) {}
    public:
    bool addRead(const Read_pt & read){
        auto res = this->readSet.insert(read);
        if(res.second){
            this->lastScore = -1;
            if(read->isSplit) nSplit++;
            return true;
        }
        return false;
    }
    double cachedScore(   const AlignmentMap_t & alnMap,
                    const ReadSet_t & usedReads)
    {
        if(this->lastScore == -1)
            this->lastScore = this->score(alnMap,usedReads);
        return this->lastScore;
    }
    bool removeRead(const Read_pt & read){
        //if(this->readSet.empty()) return false;
        if(!this->readSet.erase(read)) {
            return false;
        }
        if(read->isSplit && nSplit) nSplit--;
        this->lastScore = -1;
        return true;
    }
    double score(   const AlignmentMap_t & alnMap,
                    const ReadSet_t & usedReads,
                    bool useCache = false) const
    {
        double my_score = 0;
        for(const Read_pt & read : this->readSet){
           if(usedReads.count(read)) continue; 
               SQPair_t hPair(this->hostRegion,read);
               SQPair_t vPair(this->virusRegion,read);
               const StripedSmithWaterman::Alignment & hAln =
               alnMap.at(hPair);
               const StripedSmithWaterman::Alignment & vAln =
               alnMap.at(vPair);
               my_score += hAln.sw_score + vAln.sw_score;        
        }
        return my_score;
    }
    private:
    //void parseRegString(const std::string & regStr);
    //void parseReadString(const std::string & readStr);
};


#endif //SURVERYOR_EDGE_UTILS_H



