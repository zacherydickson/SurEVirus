#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <htslib/hts.h>
#include <htslib/sam.h>

#include "sam_utils.h"
#include "config.h"
#include "libs/cptl_stl.h"
#include "utils.h"

enum JunctionSide_t {JS_HOST, JS_VIRUS};
typedef std::unordered_map<std::string, std::pair<char, char> > GoodClipMap_t;
typedef std::unordered_map<std::string, bool> R1AnchorOrientationMap_t;

std::unordered_set<std::string> VirusNameSet;
std::mutex mtx;

//void categorize(int id, std::string contig, std::string bam_fname,
//                std::vector<bam1_t*>* reads, std::vector<bam1_t*>* anchor_reads,
//                std::unordered_map<std::string, std::pair<char, char> >& good_clips) {
//
//    open_samFile_t* bam_file_open = open_samFile(bam_fname.data());
//
//    samFile* bam_file = bam_file_open->file;
//    hts_idx_t* idx = bam_file_open->idx;
//    bam_hdr_t* header = bam_file_open->header;
//
//    bam1_t* read = bam_init1();
//
//    hts_itr_t* iter = sam_itr_querys(idx, header, contig.c_str());
//    while (sam_itr_next(bam_file, iter, read) >= 0) {
//        if (is_unmapped(read)) continue;
//
//        if (good_clips.count(bam_get_qname(read))) {
//            std::pair<char, char> good_clip_flag = good_clips[bam_get_qname(read)];
//            bool is_first_read = read->core.flag & BAM_FREAD1;
//            mtx.lock();
//            if ((is_first_read && good_clip_flag.first) || (!is_first_read && good_clip_flag.second)) {
//                uint8_t dir = is_first_read ? good_clip_flag.first : good_clip_flag.second;
//                bam1_t* copy = bam_dup1(read);
//                bam_aux_append(copy, "CD", 'A', 1, &dir);
//                anchor_reads->push_back(copy);
//            } else {
//                reads->push_back(bam_dup1(read));
//		//TODO: Determine if this can be skipped
//            }
//            mtx.unlock();
//        } else if (is_dc_pair(read) && !is_poly_ACGT(read, true)) {
//            std::string target_name = contig;
//            std::string mate_target_name = header->target_name[read->core.mtid];
//	    //TODO: test more extensivly to exclude pairs which both
//	    //map to either host or virus in  alt alignments
//            if (VirusNameSet.count(target_name) != VirusNameSet.count(mate_target_name)) {
//                mtx.lock();
//                reads->push_back(bam_dup1(read));
//                mtx.unlock();
//            }
//        }
//    }
//
//    close_samFile(bam_file_open);
//    bam_destroy1(read);
//    bam_itr_destroy(iter);
//}
//
//std::string print_fq(bam1_t* r) {
//    std::stringstream ss;
//    ss << "@" << bam_get_qname(r) << "\n";
//    std::string seq = get_sequence(r);
//    if (bam_is_rev(r)) get_rc(seq);
//    ss << seq << "\n";
//    ss << "+" << "\n";
//    for (int i = 0; i < r->core.l_qseq; i++) {
//        int idx = bam_is_rev(r) ? r->core.l_qseq-i-1 : i;
//        ss << (char) (bam_get_qual(r)[idx] + 33);
//    }
//    ss << std::endl;
//    return ss.str();
//}

//Opens a file (presumed to be a viral fasta reference) and
//  stores all accessions in the global VirusNameSet variable
//Inputs - a cstring representing the file name
//Output - none, modifies global vars
void LoadVirusNames(char* file){
    FILE* virus_ref_fasta = fopen(file, "r");
    kseq_t *seq = kseq_init(fileno(virus_ref_fasta));
    while (kseq_read(seq) >= 0) {
        VirusNameSet.insert(seq->name.s);
    }
    kseq_destroy(seq);
    fclose(virus_ref_fasta);
}


//Parses the configuration file in the workdir an extracts the number
//of threads
//Inputs - a string representing the workdir
//Output - an integer representing the number of threads
int GetConfigThreads(std::string workdir){
    return parse_config(workdir + "/config.txt").threads;
}

//given a read and some information about its clip status, return
//either the 5`-most position in the reference of the read
//or the 3`-most position
//Inputs - 0x1 - isClip, 0x2 - isLeft, 0x4 - isRev
//Output - Either zero or the reference length of the alignment
bool IsLeftOfJunction(uint8_t flag){
    if(flag & 0x1){ // The read is a clip
	if(flag & 0x2){ // The clip is a left clip
	    return true;
	}
    } else if(!(flag & 0x4)) { // A paired read which is not reversed
	return true;
    }
    return false;
}

//Note Junctions fall into one of 4 situations:
//  (A) H+V+, (B) H+V-, (C) H-V+, and (D) H-V-
//  tehnically there are 4 more where the order of Host and virus and the strands are flipped but these are equivalent (i.e. V+H+ === H-V-)
//  A clip which is reverse complemented must come from either B or C
//  B is the the case where the clip is to the left of the junction
//  and human, or to the right and viral
//Inputs - 3 boolean values: isRev, is Viral, isLeft
//Output - the strand of the junction side (either + or -)
char DetermineJunctionStrand (bool isRev, bool isViral, bool isLeft){
    char strand = '+';
    if(isRev){ // Case B or Case C
	if(isViral != isLeft){
	    return '-';
	}
    } else { // Case A or Case D
	if(!isViral & isLeft
    }
    return '+';
}

//Given a read (assumed to be good), outputs all bed entries
//corresponding to the junction side the read supports

//Inputs - an ofstream to which to write
//	 - an bam1_t object to process
//	 - a string representing the contig name
//	 - a byte containing two boolean values
//	    0x1 - isClip, 0x2 - isLeft
//	 - (optional) the read identifier, by default the read's identier is taken
//Output - none, writes to the ofstream
void OutputBEDEntries(	std::ofstream & outbed, const bam1_t* read,
			std::string cname, uint8_t clipflag,
			std::string qname=std::string()
		     )
{
    //Get the core information
    if(qname.empty()){
        qname = bam_get_qname(read);
    }

    //Get the junction position (either the left or right side of the
    //alignment depending on paired (rev vs fwd) vs clip (right vs left)
    hts_pos_t pos = read->core.pos;
    //The shift moves the bit from the 16s place of BAM_FREVERSE
    //to the 4's place of flag
    bool bLeft = IsLeftOfJunction(clipflag | 
				  ((read->core.flag & BAM_FREVERSE) >> 2));
    if(bLeft){
	pos = bam_endpos(read);
    }

    char strand = '+';
    if(read->core.flag & BAM_FREVERSE) {
    } else 
    }

    //Output Core Entry
    outbed  << cname << '\t'	<< pos  << '\t' << pos +1 << '\t'
	    << qname << "\t.\t" << strand   << '\n';

    //Check if alternate positions exist
    uint8_t *xa = bam_aux_get(read,"XA");
    if(!xa) return; //No need to continue if they don't
		    
    std::string xaListStr = bam_aux2Z(xa);
    size_t spos = 0, prev = 0;
    while((spos = xaListStr.find(';',spos+1)) != std::string::npos){
	std::string xaStr = xaListStr.substr(prev,pos=prev);
	prev=pos;
	CXA xaObj(xaStr);
	hts_pos_t xpos = xaObj.pos;
	bool xbLeft = IsLeftOfJunction(clipflag | (xaObj.bRev ? 0x8 : 0));
	if(xbLeft){
	    xpos = xaObj.endpos();
	}
	char xstrand = (xaObj.bRev) ? '-' : '+';
	outbed  << xaObj.chr	<< '\t'	<< xpos  << '\t' << xpos + 1 << '\t'
		<< qname	<< "\t.\t" << xstrand   << '\n';
    }

}

//Opens a bam file containing mapped clips.
//It is assumed that the bam file only contains primary mappped clips
//(no secondary/supplementary/unmapped)
//All reads in these files are assumed to define good clips
//Input  - a string reperesenting a file name
//	 - a reference to a GoodClipMap_t object to store the results
//Output - none, modifies the provided object
void LoadGoodClips(std::string fname,GoodClipMap_t & good_clips){
    open_samFile_t* clips_file = open_samFile(fname.c_str(), true);
    bam1_t* read = bam_init1();

    while (sam_read1(clips_file->file, clips_file->header, read) >= 0) {
        std::string clip_name = bam_get_qname(read);
        std::string qname = clip_name.substr(0, clip_name.length()-4);
        char dir = clip_name[clip_name.length()-3];
        if (clip_name[clip_name.length()-1] == '1') {
            good_clips[qname].first = dir;
        } else {
            good_clips[qname].second = dir;
        }
    }
    close_samFile(clips_file);
    bam_destroy1(read);
}

void LoadAnchorOrientation(std::string fname, R1AnchorOrientationMap_t & anchor_orient){
    open_samFile_t* anchors_file = open_samFile(fname.c_str(), true);
    bam1_t* read = bam_init1();


    while (sam_read1(clips_file->file, clips_file->header, read) >= 0) {
    }

    close_samFile(anchors_file);
    bam_destroy1(read);
}

//Process BAM files from the working environment to extract all
//candidate junctions based on BWA alignment
//Inputs - A path to the viral reference in fasta format
//	 - A path to the working directory
//	 - A path to the bam workspace
//Outputs - A bed6 formatted file detailing potential junction positions
//Each entry describes a side of a junction, whether host or virus can
//be determined by checing the contig id in column 1
//The indexes describe the interval in which the junction could be
//  For split reads this is a single position (Defined by the split)
//  For chimeric reads this ranges up to the maxIS from the junction
//  side of the mapped portion of the read (this is limited by 
//The name column is id of the fragment supporting the candidate
//  a single fragment can support multiple candidates
//The score column is unused 
//The strand column
int main(int argc, char* argv[]) {

    LoadVirusNames(argv[1]);

    std::string workdir = argv[2];
    std::string workspace = argv[3];
    //Files to be used from the workspace
    std::string bam_fname = workspace + "/retained-pairs.namesorted.bam";
    //Note: virus-clips are on the host side of a junction
    //	and host-clips are on the virus side of a junction
    std::string clip_bam_fnames[2];
    clip_bam_fnames[JS_HOST] = workspace + "/virus-clips.cs.bam";
    clip_bam_fnames[JS_VIRUS] = workspace + "/host-clips.cs.bam";
    //Host anchors = host-side, virus anchors = virus side
    std::string anchor_bam_fnames[2];
    anchor_bam_fnames[JS_HOST] = workspace + "/host-anchors.bam";
    anchor_bam_fnames[JS_HOST] = workspace + "/virus-anchors.bam";
    //Output file
    std::string bed_fname = workdir + "/junction-candidates.bed";
    std::ofstream outbed(bed_fname);

    //Set up the thread pool
    int nThread = GetConfigThreads(workdir);
    ctpl::thread_pool thread_pool(nThread);

    //Load the ids and directions of clips which map properly
    GoodClipMap_t good_clips;
    for (int clipSide = JS_HOST; clipSide <= JS_VIRUS; clipSide++){
	//TODO: Have this pull out candidates on the clip side
	LoadGoodClips(clip_bam_fnames[clipSide],good_clips);
    }

    R1AnchorOrientationMap_t

    //TODO: Parse the anchors for their supported junctions
    //TODO: Parse the paired reads for their supported junctions
    //	Note, have to filter out the good clips


    //std::vector<bam1_t*> virus_side_reads, host_side_reads;
    //std::vector<bam1_t*> virus_anchors, host_anchors;

    //std::vector<std::future<void> > futures;
    //for (int i = 0; i < header->n_targets; i++) {
    //    bool is_virus = VirusNameSet.count(header->target_name[i]);
    //    std::future<void> future = thread_pool.push(categorize, header->target_name[i], bam_fname,
    //                                                is_virus ? &virus_side_reads : &host_side_reads,
    //                                                is_virus ? &virus_anchors : &host_anchors, good_clips);
    //    futures.push_back(std::move(future));
    //}
    //thread_pool.stop(true);
    //for (int i = 0; i < futures.size(); i++) {
    //    try {
    //        futures[i].get();
    //    } catch (char const* s) {
    //        std::cout << s << std::endl;
    //    }
    //}

    //std::ofstream virus_fq(workspace + "/virus-side.fq"), host_fq(workspace + "/host-side.fq");
    //for (bam1_t* r : host_side_reads) {
    //    host_fq << print_fq(r);
    //    bam_destroy1(r);
    //}
    //for (bam1_t* r : virus_side_reads) {
    //    virus_fq << print_fq(r);
    //    bam_destroy1(r);
    //}

    //samFile* virus_anchor_writer = open_bam_writer(workspace, "virus-anchors.bam", header);
    //for (bam1_t* anchor : virus_anchors) {
    //    int ok = sam_write1(virus_anchor_writer, header, anchor);
    //    if (ok < 0) throw "Failed to write to " + std::string(virus_anchor_writer->fn);
    //    bam_destroy1(anchor);
    //}
    //sam_close(virus_anchor_writer);

    //samFile* host_anchor_writer = open_bam_writer(workspace, "host-anchors.bam", header);
    //for (bam1_t* anchor : host_anchors) {
    //    int ok = sam_write1(host_anchor_writer, header, anchor);
    //    if (ok < 0) throw "Failed to write to " + std::string(host_anchor_writer->fn);
    //    bam_destroy1(anchor);
    //}
    //sam_close(host_anchor_writer);

}
