#ifndef SURVEYOR_SAM_UTILS_H
#define SURVEYOR_SAM_UTILS_H

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <queue>
#include <cmath>
#include <cstring>
#include <string>
#include "htslib/sam.h"
#include "config.h"
#include "str_utils.h"

extern const int MIN_CLIP_LEN;
extern const int MAX_READ_SUPPORTED;

class CXA {
    public:
    std::string chr;
    uint64_t pos;
    uint32_t * cigar;
    size_t nCigar;
    uint32_t nm;
    bool bRev;
    CXA(std::string xaStr) {
	size_t field = 0;
	this->nCigar=1;
	this->cigar = (uint32_t*) malloc(sizeof(uint32_t));
	for(std::string & fieldStr : strsplit(xaStr,',')){
	    switch(field){
		case 0:
		    this->chr = fieldStr;
		    break;
		case 1:
		    this->bRev = (fieldStr.front() == '-');
		    try{
			this->pos = std::stoull(fieldStr.substr(1,
						fieldStr.length()-1));
		    } catch (std::invalid_argument & e){
			throw std::invalid_argument("Bad XA pos (" +
						    fieldStr
						    + ") from " + xaStr + "\n" +
						    e.what());
		    }
		    break;
		case 2:
		    if(sam_parse_cigar(	fieldStr.c_str(),nullptr,
					&(this->cigar),&(this->nCigar)) == -1)
			throw std::invalid_argument("parse XA cigar(" + fieldStr + ") failure for " + xaStr);
		    break;
		case 3:
		    this->nm = std::stoul(fieldStr);
		    break;
	    }
	    field++;
	}
    }
    CXA(const CXA & other) {
	this->chr = other.chr;
	this->pos = other.pos;
	this->nCigar = other.nCigar;
	this->cigar = (uint32_t*) malloc(this->nCigar*sizeof(uint32_t));
	memcpy(this->cigar,other.cigar,this->nCigar*sizeof(uint32_t));
	this->nm = other.nm;
	this->bRev = other.bRev;
    }
    ~CXA() {if(this->cigar) {free(this->cigar); this->cigar=nullptr;}}
    size_t endpos () const {
	return this->pos + bam_cigar2rlen(this->nCigar,this->cigar); 
    }
    uint8_t clipSide(){
        //0b00 = 0 - Not clipped
        //0b01 = 1 - Left Clipped
        //0b10 = 2 - Right Clipped
        //0b11 = 3 - Double Clipped
        uint8_t clipFlag = 0;
        //Get the first and last cigar operations
        char op[2] = {
            bam_cigar_opchr(this->cigar[0]),
            bam_cigar_opchr(this->cigar[this->nCigar-1])
        };
        //Determine the clip state
        for(int side = 0; side < 2; side++){
            if(op[side] == 'S' || op[side] == 'H')
                clipFlag |= (1 << side);
        }
        return clipFlag;
    }
};

int get_mate_endpos(bam1_t* r);

bool is_unmapped(bam1_t* r) {
    return r->core.flag & BAM_FUNMAP;
}
bool is_mate_unmapped(bam1_t* r) {
    return r->core.flag & BAM_FMUNMAP;
}

bool is_dup(bam1_t* r) {
    return r->core.flag & BAM_FDUP;
}

bool is_primary(bam1_t* r) {
    return !(r->core.flag & BAM_FSECONDARY) && !(r->core.flag & BAM_FSUPPLEMENTARY);
}

bool is_samechr(bam1_t* r) {
    return r->core.tid == r->core.mtid;
}

bool is_dc_pair(bam1_t* r) {
    return !is_samechr(r) || std::abs(r->core.isize) > 100000;
}

bool is_valid(bam1_t* r) {
    return !is_unmapped(r) && !is_mate_unmapped(r) && is_primary(r);
}

int get_left_clip_len(bam1_t* r, bool include_HC = false) {
    uint32_t* cigar = bam_get_cigar(r);
    return bam_cigar_opchr(cigar[0]) == 'S' || (include_HC && bam_cigar_opchr(cigar[0]) == 'H') ?
            bam_cigar_oplen(cigar[0]) : 0;
}

bool is_left_clipped(bam1_t* r, bool include_HC = false) {
    return get_left_clip_len(r, include_HC) >= MIN_CLIP_LEN;
}

bool no_fullAln_alt(bam1_t* r) {
    uint8_t *xa = bam_aux_get(r,"XA");
    if(!xa) return true;
    std::string xaListStr = bam_aux2Z(xa);
    for(std::string & xaStr : strsplit(xaListStr,';')){
        CXA xaObj(xaStr);
        bool allAln = true;
        for(size_t i = 0; i < xaObj.nCigar && allAln; i++){
            char opChr = bam_cigar_opchr(xaObj.cigar[i]);
            if(opChr != 'M' && opChr != '=' && opChr != 'X'){
        	allAln = false;
            }
        }
        if(allAln){
            return false;
        }
    }
    return true;
}

int get_right_clip_len(bam1_t *r, bool include_HC = false) {
    uint32_t* cigar = bam_get_cigar(r);
    return bam_cigar_opchr(cigar[r->core.n_cigar-1]) == 'S' || (include_HC && bam_cigar_opchr(cigar[r->core.n_cigar-1]) == 'H') ?
            bam_cigar_oplen(cigar[r->core.n_cigar-1]) : 0;
}

bool is_right_clipped(bam1_t* r, bool include_HC = false) {
    return get_right_clip_len(r, include_HC) >= MIN_CLIP_LEN;
}

std::string get_cigar_code(const uint32_t* cigar, int n_cigar) {
    std::stringstream ss;
    for (int i = 0; i < n_cigar; i++) {
        ss << bam_cigar_oplen(cigar[i]) << bam_cigar_opchr(cigar[i]);
    }
    return ss.str();
}
std::string get_cigar_code(bam1_t* r) {
    return get_cigar_code(bam_get_cigar(r), r->core.n_cigar);
}

int get_mate_endpos(bam1_t* r) {
    uint8_t *mcs = bam_aux_get(r, "MC");
    if (mcs == NULL) return r->core.mpos; // if no MC, return mpos

    char* mc = bam_aux2Z(mcs);
    int i = 0, mclen = strlen(mc);

    int len = 0, pos = r->core.mpos;
    while (i < mclen) {
        if (mc[i] >= '0' && mc[i] <= '9') {
            len = (len*10) + (mc[i]-'0');
        } else {
            if (mc[i] != 'I' && mc[i] != 'S') {
                pos += len;
            }
            len = 0;
        }
        i++;
    }
    return pos-1;
}


char get_base(const uint8_t* seq, int i) {
    char nucl2chr[16];
    nucl2chr[1] = 'A'; nucl2chr[2] = 'C'; nucl2chr[4] = 'G'; nucl2chr[8] = 'T'; nucl2chr[15] = 'N';
    return nucl2chr[bam_seqi(seq, i)];
}

int get_avg_qual(bam1_t* r, int start, int end) {
    if (start >= end) return 0;
    int q = 0;
    for (int i = start; i < end; i++) {
        q += bam_get_qual(r)[i];
    }
    return q/(end-start);
}
int get_avg_qual(bam1_t* r, bool aligned_only = true) {
    int start = aligned_only ? get_left_clip_len(r) : 0;
    int end = r->core.l_qseq - (aligned_only ? get_right_clip_len(r) : 0);
    return get_avg_qual(r, start, end);
}

bool is_poly_ACGT(bam1_t* r, bool aligned_only = false) {
    int a = 0, c = 0, g = 0, t = 0;
    int start = aligned_only ? get_left_clip_len(r) : 0;
    int end = r->core.l_qseq - (aligned_only ? get_right_clip_len(r) : 0);
    const uint8_t* bam_seq = bam_get_seq(r);
    for (int i = start; i < end; i++) {
        char base = get_base(bam_seq, i);
        if (base == 'A') a++;
        else if (base == 'C') c++;
        else if (base == 'G') g++;
        else if (base == 'T') t++;
    }

    int maxc = std::max(std::max(a,c), std::max(g,t));
    return double(maxc)/(end-start) >= 0.8;
}

void reverse(std::string& s) {
    int len = s.length();
    for (int i = 0; i < len/2; i++) {
        std::swap(s[i], s[len-i-1]);
    }
}

void get_rc(std::string& read) {
    int len = read.length();
    reverse(read);
    for (int i = 0; i < len; i++) {
        if (read[i] == 'A') read[i] = 'T';
        else if (read[i] == 'C') read[i] = 'G';
        else if (read[i] == 'G') read[i] = 'C';
        else if (read[i] == 'T') read[i] = 'A';
        else read[i] = 'N';
    }
}

void reverse(char* read, int len) {
    for (int i = 0; i < len/2; i++) {
        std::swap(read[i], read[len-i-1]);
    }
}

void get_rc(char* read, int len) {
    reverse(read, len);
    for (int i = 0; i < len; i++) {
        if (read[i] == 'A') read[i] = 'T';
        else if (read[i] == 'C') read[i] = 'G';
        else if (read[i] == 'G') read[i] = 'C';
        else if (read[i] == 'T') read[i] = 'A';
        else read[i] = 'N';
    }
}

std::string get_sequence(bam1_t* r, bool original_seq = false) { // if original_seq == true, return the sequence found in fastx file
    char seq[MAX_READ_SUPPORTED];
    const uint8_t* bam_seq = bam_get_seq(r);
    for (int i = 0; i < r->core.l_qseq; i++) {
        seq[i] = get_base(bam_seq, i);
    }
    seq[r->core.l_qseq] = '\0';
    if (original_seq && bam_is_rev(r)) get_rc(seq, r->core.l_qseq);
    return std::string(seq);
}

#define bam1_seq_seti(s, i, c) ( (s)[(i)>>1] = ((s)[(i)>>1] & 0xf<<(((i)&1)<<2)) | (c)<<((~(i)&1)<<2) )

void rc_sequence(bam1_t* read) {
    std::string seq = get_sequence(read);
    get_rc(seq);
    uint8_t* s = bam_get_seq(read);
    for (int i = 0; i < read->core.l_qseq; ++i){
        bam1_seq_seti(s, i, seq_nt16_table[(int)seq[i]]);
    }
}
void set_to_forward(bam1_t* read) {
    if (bam_is_rev(read)) {
        rc_sequence(read);
    }
    read->core.flag &= ~BAM_FREVERSE;
}
void set_to_reverse(bam1_t* read) {
    if (!bam_is_rev(read)) {
        rc_sequence(read);
    }
    read->core.flag |= BAM_FREVERSE;
}
bam1_t* rc_read(bam1_t* read) {
    if (bam_is_rev(read)) set_to_forward(read);
    else set_to_reverse(read);
    return read;
}

std::string get_qualities(bam1_t* r, bool original_quals = false) {
    char quals[MAX_READ_SUPPORTED];
    char* q = (char*) bam_get_qual(r);
    for (int i = 0; i < r->core.l_qseq; i++) {
        quals[i] = q[i]+33;
    } quals[r->core.l_qseq] = '\0';
    if (original_quals && bam_is_rev(r)) reverse(quals, r->core.l_qseq);
    return quals;
}

//Requires that the start and end positions be in the current sequences orientation
bool is_low_complexity(const char* seq, int start, int end) {

    int count[256] = {};

    char chr2nucl[256];
    chr2nucl[(int)'A'] = 1;
    chr2nucl[(int)'C'] = 2;
    chr2nucl[(int)'G'] = 4;
    chr2nucl[(int)'T'] = 8;
    chr2nucl[(int)'N'] = 15;
    chr2nucl[(int)'a'] = 1;
    chr2nucl[(int)'c'] = 2;
    chr2nucl[(int)'g'] = 4;
    chr2nucl[(int)'t'] = 8;
    chr2nucl[(int)'n'] = 15;

    for (int i = start+1; i < end; i++) {
        int twobases = (chr2nucl[(int)seq[i-1]] << 4) | chr2nucl[(int)seq[i]];
        count[twobases]++;
    }

    uint8_t top1 = 0, top2 = 0;
    for (uint8_t i = 1; i < 16; i*=2) {
        for (uint8_t j = 1; j < 16; j*=2) {
            uint8_t c = (i << 4) | j;
            if (count[c] > count[top1]) {
                top2 = top1;
                top1 = c;
            } else if (count[c] > count[top2]) {
                top2 = c;
            }
        }
    }

    bool is_lc = count[top1] + count[top2] >= (end-start+1)*0.75;
    return is_lc;
}
bool is_low_complexity(char* seq, bool rc, int start, int end) {

    int count[256] = {};

    if (rc) get_rc(seq, strlen(seq)); // what's the use of RC when determining low-complexity?

    char chr2nucl[256];
    chr2nucl[(int)'A'] = 1;
    chr2nucl[(int)'C'] = 2;
    chr2nucl[(int)'G'] = 4;
    chr2nucl[(int)'T'] = 8;
    chr2nucl[(int)'N'] = 15;
    chr2nucl[(int)'a'] = 1;
    chr2nucl[(int)'c'] = 2;
    chr2nucl[(int)'g'] = 4;
    chr2nucl[(int)'t'] = 8;
    chr2nucl[(int)'n'] = 15;

    for (int i = start+1; i < end; i++) {
        int twobases = (chr2nucl[(int)seq[i-1]] << 4) | chr2nucl[(int)seq[i]];
        count[twobases]++;
    }

    uint8_t top1 = 0, top2 = 0;
    for (uint8_t i = 1; i < 16; i*=2) {
        for (uint8_t j = 1; j < 16; j*=2) {
            uint8_t c = (i << 4) | j;
            if (count[c] > count[top1]) {
                top2 = top1;
                top1 = c;
            } else if (count[c] > count[top2]) {
                top2 = c;
            }
        }
    }

    bool is_lc = count[top1] + count[top2] >= (end-start+1)*0.75;
    return is_lc;
}
bool is_low_complexity(bam1_t* read, bool rc, int start, int end) {
    std::string seq = get_sequence(read);
    return is_low_complexity((char*) seq.data(), rc, start, end);
}

std::string get_left_clip(bam1_t* read) {
    std::string seq = get_sequence(read);
    return seq.substr(0, get_left_clip_len(read));
}
std::string get_right_clip(bam1_t* read) {
    std::string seq = get_sequence(read);
    return seq.substr(seq.length()- get_right_clip_len(read)-1, get_right_clip_len(read));
}

void del_aux(bam1_t* read, const char* tag) {
    uint8_t* b = bam_aux_get(read, tag);
    if (b == NULL) return;
    bam_aux_del(read, b);
}

std::pair<int, const uint32_t*> cigar_str_to_array(std::string& cigar_str) {
    std::vector<uint32_t> opv;

    int pos = 0, prev = 0;
    std::string bam_ops = BAM_CIGAR_STR;
    while ((pos = cigar_str.find_first_of(bam_ops, prev)) != std::string::npos) {
        opv.push_back(bam_cigar_gen(std::stoi(cigar_str.substr(prev, pos-prev)), bam_ops.find(cigar_str[pos])));
        prev = pos+1;
    }

    uint32_t* opa = new uint32_t[opv.size()];
    std::copy(opv.begin(), opv.end(), opa);
    return std::make_pair(opv.size(), opa);
}
std::string cigar_array_to_str(int cigar_len, const uint32_t* cigar) {
    std::stringstream ss;
    for (int i = 0; i < cigar_len; i++) {
        ss << bam_cigar_oplen(cigar[i]) << bam_cigar_opchr(cigar[i]);
    }
    return ss.str();
}

struct open_samFile_t {
    samFile* file;
    bam_hdr_t* header;
    hts_idx_t* idx;

    open_samFile_t() : file(NULL), header(NULL), idx(NULL) {}
};

open_samFile_t* open_samFile(const char* fname, bool index_file = false, bool open_index = true) {
    open_samFile_t* sam_file = new open_samFile_t;
    sam_file->file = sam_open(fname, "r");
    if (sam_file->file == NULL) {
        throw "Could not open " + std::string(fname);
    }

    if (index_file) {
        int code = sam_index_build(fname, 0);
        if (code != 0) {
            throw "Cannot index " + std::string(fname);
        }
    }

    if (open_index) {
        sam_file->idx = sam_index_load(sam_file->file, sam_file->file->fn);
        if (sam_file->idx == NULL) {
            throw "Unable to open index for " + std::string(fname);
        }
    }

    sam_file->header = sam_hdr_read(sam_file->file);
    if (sam_file->header == NULL) {
        throw "Unable to open header for " + std::string(fname);
    }

    return sam_file;
}

void close_samFile(open_samFile_t* f) {
    hts_idx_destroy(f->idx);
    bam_hdr_destroy(f->header);
    sam_close(f->file);
    delete f;
}

samFile* open_bam_writer(std::string dir, std::string name, bam_hdr_t* header, bool sam = false) {
    std::string filename = dir + "/" + name;
    samFile* writer = sam_open(filename.c_str(), sam ? "w" : "wb");
    if (sam_hdr_write(writer, header) != 0) {
        throw "Could not write file " + filename;
    }
    return writer;
}

void load_reads(std::string bam_fname, std::vector<bam1_t*>& reads, bool load_duplicates = true) {
    bam1_t* read = bam_init1();
    open_samFile_t* host_file = open_samFile(bam_fname.c_str(), false, false);
    while (sam_read1(host_file->file, host_file->header, read) >= 0) {
        if (!load_duplicates && is_dup(read)) continue;
        reads.push_back(bam_dup1(read));
    }
    close_samFile(host_file);
    bam_destroy1(read);
}


template<typename T>
inline T max(T a, T b, T c) { return std::max(std::max(a,b), c); }

template<typename T>
inline T max(T a, T b, T c, T d) { return std::max(std::max(a,b), std::max(c,d)); }

template<typename T>
inline T max(T a, T b, T c, T d, T e) { return max(std::max(a,b), std::max(c,d), e); }

template<typename T>
inline T min(T a, T b, T c) { return std::min(std::min(a,b), c); }

template<typename T>
inline T min(T a, T b, T c, T d) { return std::min(std::min(a,b), std::min(c,d)); }

int compareBamByPos(bam1_t* & a, bam1_t* & b){
    if(a->core.mtid != b->core.mtid){
        return (a->core.mtid < b->core.mtid) ? -1 : 1;
    } else if(a->core.pos != b->core.pos){
        return (a->core.pos < b->core.pos) ? -1 : 1;
    }
    return 0;
}

#endif //SURVEYOR_SAM_UTILS_H
