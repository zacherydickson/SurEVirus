#ifndef SURVEYOR_BRANCHED_QUEUE_H
#define SURVEYOR_BRANCHED_QUEUE_H

#include <forward_list>
#include <memory>

#include "edge_utils.h"

struct BranchedEdgeQueueNode_t;

class CBranchedEdgeQueue;

typedef std::unique_ptr<BranchedEdgeQueueNode_t> BranchedEdgeQueueNode_pt;
typedef std::forward_list<BranchedEdgeQueueNode_pt> BranchedEdgeQueueNodeList_t;
typedef std::unique_ptr<CBranchedEdgeQueue> CBranchedEdgeQueue_p;

struct BranchedEdgeQueueNode_t {
    public:
        friend CBranchedEdgeQueue;
    protected:
    //Members
    CBranchedEdgeQueue_p    childQueue;
    ReadSet_t               descendentReads;
    Edge_t                  edge;
    ReadSet_t               originalReads;
    //Con-/Destruction
    public:
    BranchedEdgeQueueNode_t(const Edge_t & e)
        : childQueue(nullptr), descendentReads(e.readSet),edge(e),
          originalReads(e.readSet) {}
    protected:
    //Methods
    void addNode(   BranchedEdgeQueueNode_pt node,
                    const AlignmentMap_t & alnMap,
                    const ReadSet_t & used,
                    bool bNewElement);
    void addReads(ReadSet_t reads);
    size_t queueSize();
    void removeReads(const ReadSet_t & reads,bool bRecurse=false);
    double score(const AlignmentMap_t & alnMap, const ReadSet_t & used){
        return this->edge.cachedScore(alnMap,used);
    }
    
};



//Data structure for efficiently outputting edges such that any read
//supports only one edge
//The Elements of the Queue are sorted by score
//Any element which shares reads with a higher score element is 
//  placed in the higher scoring subqueue, and those shared reads are
//  deleted from the sub-edge
//When an element is removed, its subqueues are added to the queue
//Note: behaves best when edges are added to the queue in descending order
//of score
//Data is maintained in a forward_list where the last element always has the
//best score
class CBranchedEdgeQueue {
    //Friendship
    public:
        friend BranchedEdgeQueueNode_t;
    //Members
    public:
        const AlignmentMap_t * pAlnMap;
        const ReadSet_t * pUsedReads;
    protected:
        BranchedEdgeQueueNodeList_t data;
        size_t nEdges;
        size_t nElements;
        bool validNElementCache;
    //Con-/Destruction
    public:
        CBranchedEdgeQueue(AlignmentMap_t const * pMap, ReadSet_t const * pSet)
            : pAlnMap(pMap), pUsedReads(pSet), data(), nEdges(0), nElements(0),
              validNElementCache(false)
        {}
        ~CBranchedEdgeQueue() {}
    //Accessors
        bool empty() const {return this->data.empty();}
        size_t queueSize() {
            if(validNElementCache){
                return this->nElements;
            }
            nElements = this->nEdges;
            for(auto & node : this->data){
                nElements += node->queueSize();
            }
            validNElementCache = true;
            return this->nElements;
        }
        size_t size() const {return this->nEdges;}
        const Edge_t & top() const {
            if(this->data.empty()){
                throw std::logic_error( "Attempt to call top on empty"
                                        " CBranchEdgeQueue");
            }
            return this->data.front()->edge;
        }
    //Methods
    public:
        void addEdge(const Edge_t & edge) {
            BranchedEdgeQueueNode_pt node =
                std::make_unique<BranchedEdgeQueueNode_t>(edge);
            this->addNode(std::move(node),true);
        }
        //Remove the first element, return any unused reads 
        // to its descendents
        //Then remove the top element
        void cannabalize() { this->removeTop(true); }
        //Just remove the top element
        void pop() { this->removeTop(false); }
    protected:
        //Return the change in the number of elements in the queue
        void addNode(BranchedEdgeQueueNode_pt node,bool bNewElement)
        {
            if(bNewElement) {
                this->validNElementCache = false;
            }
            const AlignmentMap_t & alnMap = *(this->pAlnMap);
            const ReadSet_t & usedReads = *(this->pUsedReads);
            if(this->data.empty()){
                data.push_front(std::move(node));
                this->nEdges++;
                return;
            }
            //Set the inertion and merge points to values
            //corresponding to:
            //        insert at front, do not merge
            BranchedEdgeQueueNodeList_t::iterator insIt = 
                this->data.before_begin();
            BranchedEdgeQueueNodeList_t::iterator mergeIt = 
                this->data.end();
            BranchedEdgeQueueNodeList_t::iterator preMergeIt = 
                this->data.before_begin();
            //Locate the merge and insertion points
            //int idx,ins,merge;
            //idx=0;
            //ins = -1;
            //merge = this->nEdges+1;
            for(auto it = this->data.begin(),
                    pit=this->data.before_begin(); 
                it != this->data.end() && mergeIt == this->data.end();
                it++,pit++)//,idx++)
            {
                const BranchedEdgeQueueNode_pt & curNode = *it;
                //Check If the current Node has a higher score
                //  The last one for which this is true is the
                //  insert after for
                if( curNode->score(alnMap,usedReads) >
                    node->score(alnMap,usedReads))
                { //Track that for later
                    insIt = it;
                    //ins = idx;
                }
                //
                for(const Read_pt & read : node->edge.readSet){
                    if(curNode->descendentReads.count(read) ||
                       (read->mate &&
                        curNode->descendentReads.count(read->mate)))
                    {
                        preMergeIt = pit;
                        mergeIt = it;
                        //merge = idx;
                        break;
                    }
                }
            }
            //The Case where This Node has no shared reads
            //        with any node already in the queue
            //Simply put it in place
            if(mergeIt == this->data.end()){
                this->data.insert_after(insIt,std::move(node));
                this->nEdges++;
                return;
            }
            BranchedEdgeQueueNode_pt & mergeNode = *mergeIt;
            //The Case where the node shares reads with some
            //other node in the queue
            //Favouring nodes already in the queue as it
            //is less expensive
            if( mergeNode->score(alnMap,usedReads) >=
                node->score(alnMap,usedReads))
            {
                //This node is lower priority, add it to the
                //mergePoint
                //Add this node's reads to the mergeNodes's 
                //  set of reads
                mergeNode->descendentReads.insert(
                        node->descendentReads.begin(),
                        node->descendentReads.end());
                //Remove the merge Node's reads from 
                //  this node's edge's readSet
                node->removeReads(mergeNode->edge.readSet,false);
                //Don't bother adding a node that is now empty
                mergeNode->addNode( std::move(node),alnMap,usedReads,
                                    bNewElement);
                return;
            }
            //Final Case, this node is better so we have to 
            // pull up the mergeNode, put this node in place
            // and add the merge Node to this Node's queue
            node->descendentReads.insert(
                mergeNode->descendentReads.begin(),
                mergeNode->descendentReads.end());
            mergeNode->removeReads(node->edge.readSet,true);
            //Invalidate the entry in this queue by pulling up the merge
            //Node and adding to this node's queue
            //If this returns false the node to add became empty and was
            //removed
            node->addNode(std::move(mergeNode),alnMap,usedReads,bNewElement);
            //Remove the invalidated entry
            this->data.erase_after(preMergeIt);
            //Put this node into this queue at the correct spot
            this->data.insert_after(insIt,std::move(node));
            return;
        }
        void removeReads(const ReadSet_t & set){
            for(auto & node : this->data){
                node->removeReads(set,true);
            }
        }
        void removeTop(bool cannabalize = false);
};

void BranchedEdgeQueueNode_t::addNode(
            BranchedEdgeQueueNode_pt node,
            const AlignmentMap_t & alnMap,
            const ReadSet_t & used,
            bool bNewElement)
{
    //Provided an empty node: do not add it
    // but attempt to recover its children
    if(!node->edge.readSet.size()){
        if(node->childQueue){
            for(BranchedEdgeQueueNode_pt & cNode: node->childQueue->data){
                this->addNode(std::move(cNode),alnMap,used,false);
            }
        }
        return;
    }
    if(!this->childQueue){
        this->childQueue = std::make_unique<CBranchedEdgeQueue>(
                &alnMap,&used);
    }
    this->childQueue->addNode(std::move(node),bNewElement);
}

void BranchedEdgeQueueNode_t::addReads(ReadSet_t reads){
    //Add any reads which belong to this node to it's edge's readSet
    //  and remove them from the list of reads to add
    for(auto it = reads.begin(); it != reads.end();){
        const Read_pt & read = *it;
        if(!this->descendentReads.count(read)){
            it = reads.erase(it);
        } else if(this->originalReads.count(read)){
            this->edge.addRead(read);
            if(read->mate && this->originalReads.count(read->mate)){
                this->edge.addRead(read->mate);
            }
            it = reads.erase(it);
        } else if(read->mate && this->originalReads.count(read->mate)){
            this->edge.addRead(read->mate);
            it = reads.erase(it);
        } else{
            it++;
        }
    }
    if(this->childQueue){
        for(auto & cNode : this->childQueue->data){
            cNode->addReads(reads);
        }
    }
}
    
size_t BranchedEdgeQueueNode_t::queueSize() {
    if(this->childQueue)
        return this->childQueue->queueSize();
    return 0;
}

void BranchedEdgeQueueNode_t::removeReads(
            const ReadSet_t & reads,
            bool bRecurse)
{
    for(const Read_pt & read : reads){
        this->edge.removeRead(read);
        if(read->mate)
            this->edge.removeRead(read->mate);
        //Won't remove from descendedReads
        //        Expensive operation, for no benefit
    }
    if(bRecurse && this->childQueue){
        this->childQueue->removeReads(reads);
    }
}

void CBranchedEdgeQueue::removeTop(bool cannabalize){
    if(this->data.empty()){
        throw std::logic_error( "Attempt to remove top element on empty"
                                " CBranchEdgeQueue");
    }
    //Pull off the first element
    BranchedEdgeQueueNode_pt node = std::move(this->data.front());
    //Erase it from the vector
    this->data.pop_front();
    this->nEdges--;
    this->nElements--;
    //Add the nodes' children to this queue
    if(node->childQueue){
        BranchedEdgeQueueNodeList_t & childVec =
            node->childQueue->data;
            for(auto it = childVec.begin();
                it != childVec.end(); it++)
            {
                if(cannabalize)
                    (*it)->addReads(node->edge.readSet);
                this->addNode(std::move(*it),false);
            }
    }
}



#endif //SURVEYOR_BRANCHED_QUEUE_H




