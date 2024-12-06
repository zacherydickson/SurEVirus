#ifndef SURVEYOR_CLUSTER_H
#define SURVEYOR_CLUSTER_H

#include <algorithm>
#include <htslib/kseq.h>
#include <htslib/sam.h>
#include <iostream>
#include <sstream>
#include <map>
#include <unordered_set>
#include <unistd.h>
#include <vector>

#include "config.h"
#include "sam_utils.h"
#include "libs/ssw.h"
#include "libs/ssw_cpp.h"

static const size_t MinimumAlignmentLength = 30;


KSEQ_INIT(int, read)

struct breakpoint_t {
    std::string chr;
    bool rev;
    int start, end;

    int pos() {return rev ? start : end;}

    breakpoint_t() {}
    explicit breakpoint_t(const std::string& str) {
        char contig[1001], strand;
        sscanf(str.data(), "%1000[^:]:%c:%d:%d", contig, &strand, &start, &end);
        this->chr = contig;
        rev = (strand == '-');
    }
    breakpoint_t(std::string chr, int start, int end, bool rev) : chr(chr), rev(rev), start(start), end(end) {}

    bool operator == (breakpoint_t& other) {
        return chr == other.chr && rev == other.rev && pos() == other.pos();
    }

    std::string to_string() {
        std::stringstream ssout;
        ssout << chr << ':' << (rev ? '-' : '+') << ":" << start << ":" << end;
        return ssout.str();
    }

    std::string to_human_string() {
        std::stringstream ssout;
        ssout << chr << ':' << (rev ? '-' : '+') << pos();
        return ssout.str();
    }
};

struct call_t {
    int id;
    breakpoint_t host_bp, virus_bp;
    int reads, good_pairs, split_reads, reads_w_dups, unique_reads_w_dups, score;
    double host_pbs, virus_pbs;
    double host_cov, virus_cov;
    bool removed = false;
    int paired_with = -1;

    explicit call_t(std::string& line) {
        std::stringstream ssin(line);
        std::string host_bp_str, virus_bp_str;
        ssin >> id >> host_bp_str >> virus_bp_str >> reads >> good_pairs >> split_reads >> score >> host_pbs >> virus_pbs
             >> reads_w_dups >> unique_reads_w_dups >> host_cov >> virus_cov;
        host_bp = breakpoint_t(host_bp_str);
        virus_bp = breakpoint_t(virus_bp_str);
    }

    call_t(int id, breakpoint_t& host_bp, breakpoint_t& virus_bp, int reads, int good_pairs, int split_reads,
        int reads_w_dups, int uniquer_reads_w_dups, int score, double host_pbs, double virus_pbs,
        double host_cov, double virus_cov) :
    id(id), host_bp(host_bp), virus_bp(virus_bp), reads(reads), good_pairs(good_pairs), split_reads(split_reads),
    reads_w_dups(reads_w_dups), unique_reads_w_dups(uniquer_reads_w_dups), score(score),
    host_pbs(host_pbs), virus_pbs(virus_pbs),
    host_cov(host_cov), virus_cov(virus_cov) {}

    double coverage() { return (host_cov + virus_cov)/2; }

    bool is_paired() { return paired_with >= 0; }

    std::string to_string() {
        char buffer[4096];
        sprintf(buffer, "%d %s %s %d %d %d %d %lf %lf %d %d %lf %lf", id,
                host_bp.to_string().c_str(), virus_bp.to_string().c_str(), reads, good_pairs, split_reads, score,
                host_pbs, virus_pbs, reads_w_dups, unique_reads_w_dups, host_cov, virus_cov);
        return buffer;
    }

    std::string to_human_string() {
        char buffer[10000];
        sprintf(buffer, "ID=%d %s %s SUPPORTING_PAIRS=%d SPLIT_READS=%d HOST_PBS=%lf COVERAGE=%lf",
                id, host_bp.to_human_string().c_str(), virus_bp.to_human_string().c_str(), good_pairs, split_reads, host_pbs, coverage());
        std::string sout = buffer;
        if (is_paired()) sout += " PAIRED WITH ID=" + std::to_string(paired_with);
        return sout;
    }
};

int pair_dist(call_t& c1, call_t& c2) {
    if (c1.host_bp.chr == c2.host_bp.chr && c1.host_bp.rev != c2.host_bp.rev && c1.virus_bp.rev != c2.virus_bp.rev) {
        int rev_pos, fwd_pos;
        if (c1.host_bp.rev) {
            rev_pos = c1.host_bp.pos();
            fwd_pos = c2.host_bp.pos();
        } else {
            rev_pos = c2.host_bp.pos();
            fwd_pos = c1.host_bp.pos();
        }

        if (rev_pos-fwd_pos >= -50 && rev_pos-fwd_pos <= 1000) {
            return rev_pos-fwd_pos;
        }
    }
    return INT32_MAX;
}

char _cigar_int_to_op(uint32_t c) {
    char op = cigar_int_to_op(c);
    return (op != 'X' && op != '=') ? op : 'M';
}

std::string alignment_cigar_to_bam_cigar(const std::vector<uint32_t> & cigar) {
    std::stringstream ss;
    char op = ' '; int len = 0;
    for (uint32_t c : cigar) {
        if (op != _cigar_int_to_op(c)) {
            if (op != ' ') ss << len << op;
            op = _cigar_int_to_op(c);
            len = cigar_int_to_len(c);
        } else {
            len += cigar_int_to_len(c);
        }
    }
    ss << len << op;
    return ss.str();
}

bool accept_alignment(const StripedSmithWaterman::Alignment & alignment, int min_sc_size) {
    bool long_enough = alignment.ref_end-alignment.ref_begin+1 >= MinimumAlignmentLength;
    bool score_enough = alignment.sw_score >= 30;
    uint32_t c0 = alignment.cigar[0], cl = alignment.cigar[alignment.cigar.size()-1];
    bool left_clipped = cigar_int_to_op(c0) == 'S' && cigar_int_to_len(c0) >= min_sc_size;
    bool right_clipped = cigar_int_to_op(cl) == 'S' && cigar_int_to_len(cl) >= min_sc_size;
    bool clipped_both_sides = left_clipped && right_clipped;
    return long_enough && score_enough && !clipped_both_sides;
}

//Opens a file (presumed to be a viral fasta reference) and
//  stores all accessions in the provided unordered set
//Inputs - a cstring representing the file name
//	 - a reference to a vector of strings in which to to store results
//Output - none, modifies provided set
void LoadVirusNames(const std::string & file,std::unordered_set<std::string> & virusNameSet){
    FILE* virus_ref_fasta = fopen(file.c_str(), "r");
    kseq_t *seq = kseq_init(fileno(virus_ref_fasta));
    while (kseq_read(seq) >= 0) {
        virusNameSet.insert(seq->name.s);
    }
    kseq_destroy(seq);
    fclose(virus_ref_fasta);
}

//From Arafat Hasan: Answer to stack Overflow Question 14265581
std::vector<std::string> strsplit (const std::string &s, char delim) {
    std::vector<std::string> result;
    std::stringstream ss (s);
    std::string item;

    while (getline (ss, item, delim)) {
        result.push_back (item);
    }

    return result;
}

std::string get_seqrc(const std::string & seq) {
    std::string rc(seq);
    std::reverse(rc.begin(),rc.end());
    for( char & c : rc){
	char r = 'N';
	switch(c) {
	    case 'A': r = 'T'; break;
	    case 'C': r = 'G'; break;
	    case 'T': r = 'A'; break;
	    case 'G': r = 'C'; break;
	}
	c=r;
    }
    return rc;
}

#endif //SURVEYOR_CLUSTER_H
