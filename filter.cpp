#include <iostream>
#include <vector>
#include <cstring>
#include <unistd.h>
#include <htslib/kseq.h>

KSEQ_INIT(int, read)

#include "config.h"
#include "libs/ssw.h"
#include "libs/ssw_cpp.h"


void get_rc(std::string& read) {
    int len = read.length();
    for (int i = 0; i < len/2; i++) {
        std::swap(read[i], read[len-i-1]);
    }
    for (int i = 0; i < len; i++) {
        if (toupper(read[i]) == 'A') read[i] = 'T';
        else if (toupper(read[i]) == 'C') read[i] = 'G';
        else if (toupper(read[i]) == 'G') read[i] = 'C';
        else if (toupper(read[i]) == 'T') read[i] = 'A';
        else read[i] = 'N';
    }
}

struct breakpoint_t {
    std::string chr;
    bool rev;
    int start, end;

    int pos() {return rev ? start : end;}

    breakpoint_t() {}
    breakpoint_t(std::string& str) {
        char chr[1000], strand;
        sscanf(str.data(), "%[^:]:%c:%d:%d", chr, &strand, &start, &end);
        this->chr = chr;
        rev = (strand == '-');
    }

    bool operator == (breakpoint_t& other) {
        return chr == other.chr && rev == other.rev && pos() == other.pos();
    }

    std::string to_string() {
        std::stringstream ssout;
        ssout << chr << ':' << (rev ? '-' : '+') << pos();
        return ssout.str();
    }
};

struct call_t {
    int id;
    breakpoint_t host_bp, virus_bp;
    int reads, good_reads, split_reads, reads_w_dups, unique_reads_w_dups, score;
    double host_pbs, virus_pbs;
    double host_cov, virus_cov;
    bool removed = false;
    int paired_with = -1;

    call_t(std::string& line) {
        std::stringstream ssin(line);
        std::string host_bp_str, virus_bp_str;
        ssin >> id >> host_bp_str >> virus_bp_str >> reads >> good_reads >> split_reads >> score >> host_pbs >> virus_pbs
            >> reads_w_dups >> unique_reads_w_dups >> host_cov >> virus_cov;
        host_bp = breakpoint_t(host_bp_str);
        virus_bp = breakpoint_t(virus_bp_str);
    }

    double coverage() { return (host_cov + virus_cov)/2; }

    bool is_paired() { return paired_with >= 0; }

    std::string to_string() {
        char buffer[10000];
        sprintf(buffer, "ID=%d %s %s GOOD_READS=%d TOT_READS=%d SPLIT_READS=%d HOST_PBS=%lf COVERAGE=%lf",
                id, host_bp.to_string().c_str(), virus_bp.to_string().c_str(), good_reads, reads, split_reads, host_pbs, coverage());
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

void print_calls(std::vector<call_t>& calls) {
    for (int i = 0; i < calls.size(); i++) {
        call_t call = calls[i];
        if (!call.removed) {
            std::cout << call.to_string() << std::endl;
        }
    }
}

StripedSmithWaterman::Alignment align(StripedSmithWaterman::Aligner& aligner, std::string& ref, std::string& query) {
    std::string rc_query = query;
    get_rc(rc_query);

    StripedSmithWaterman::Alignment h_aln, h_aln_rc;
    StripedSmithWaterman::Filter filter;
    aligner.Align(query.c_str(), ref.c_str(), ref.length(), filter, &h_aln, 0);
    aligner.Align(rc_query.c_str(), ref.c_str(), ref.length(), filter, &h_aln_rc, 0);

    if (h_aln.sw_score > h_aln_rc.sw_score) return h_aln;
    else return h_aln_rc;
}

int main(int argc, char* argv[]) {

    std::string workdir = argv[1];

    config_t config = parse_config(workdir + "/config.txt");

    std::ifstream fin(workdir + "/results.txt");

    std::vector<std::string> args;
    bool print_rejected = false, accept_clipped = false;
    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--print-rejected") == 0) {
            print_rejected = true;
        } else if (strcmp(argv[i], "--accept-all-clipped") == 0) {
            accept_clipped = true;
        } else {
            args.push_back(argv[i]);
        }
    }

    double min_reads = 4;
    if (args.size() >= 2) min_reads = std::stoi(args[1]);

    double min_host_pbs = 0.8;
    if (args.size() >= 3) min_host_pbs = std::stod(args[2]);

    double min_cov_for_accept = 0.5;
    if (args.size() >= 4) min_cov_for_accept = std::stod(args[3]);


    // read host and virus sequences
    std::string host_bp_seqs_fname = workdir + "/host_bp_seqs.fa", virus_bp_seqs_fname = workdir + "/virus_bp_seqs.fa";
    FILE* host_bp_fasta = fopen(host_bp_seqs_fname.c_str(), "r"), *virus_bp_fasta = fopen(virus_bp_seqs_fname.c_str(), "r");
    kseq_t* hseq = kseq_init(fileno(host_bp_fasta)), *vseq = kseq_init(fileno(virus_bp_fasta));
    std::vector<std::pair<std::string, std::string> > host_virus_bp_seqs;
    while (kseq_read(hseq) >= 0 && kseq_read(vseq) >= 0) {
        host_virus_bp_seqs.push_back({hseq->seq.s, vseq->seq.s});
    }


    std::vector<call_t> accepted_calls, rejected_calls;
    std::string line;
    while (std::getline(fin, line)) {
        call_t call(line);

        bool accept = true;
        if (call.host_pbs < min_host_pbs) accept = false;
        if (call.good_reads < min_reads-1 || (call.good_reads == min_reads-1 && call.split_reads == 0)) accept = false;
        if (host_virus_bp_seqs[call.id].first.length() < config.read_len/2 || host_virus_bp_seqs[call.id].second.length() < config.read_len/2) {
            accept = false;
        }

        if (call.coverage() >= min_cov_for_accept && call.good_reads > 2) accept = true;
        if (accept_clipped && call.split_reads > 0) accept = true;

        if (accept) {
            accepted_calls.push_back(call);
        } else {
            rejected_calls.push_back(call);
        }
    }


    // find pairs
    for (int i = 0; i < accepted_calls.size(); i++) {
        call_t a_call1 = accepted_calls[i];

        int pair_j = -1, min_dist = INT32_MAX;
        for (int j = i+1; j < accepted_calls.size(); j++) {
            call_t a_call2 = accepted_calls[j];
            if (a_call2.is_paired()) continue;

            int dist = pair_dist(a_call1, a_call2);
            if (abs(dist) < min_dist) {
                pair_j = j;
                min_dist = abs(dist);
            }
        }

        if (pair_j != -1) {
            accepted_calls[i].paired_with = accepted_calls[pair_j].id;
            accepted_calls[pair_j].paired_with = accepted_calls[i].id;
        }
    }

    for (int i = 0; i < accepted_calls.size(); i++) {
        call_t a_call = accepted_calls[i];
        if (a_call.is_paired()) continue;

        int pair_j = -1, min_dist = INT32_MAX;
        for (int j = 0; j < rejected_calls.size(); j++) {
            if (rejected_calls[j].is_paired()) continue;

            call_t r_call = rejected_calls[j];
            int dist = pair_dist(a_call, r_call);
            if (abs(dist) < min_dist) {
                pair_j = j;
                min_dist = abs(dist);
            }
        }

        if (pair_j != -1) {
            a_call.paired_with = rejected_calls[pair_j].id;
            rejected_calls[pair_j].paired_with = a_call.id;
            accepted_calls.push_back(rejected_calls[pair_j]);
        }
    }


    // find duplicates
    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
    for (int j = 0; j < accepted_calls.size(); j++) {
        std::string hseq_j = host_virus_bp_seqs[accepted_calls[j].id].first;
        std::string vseq_j = host_virus_bp_seqs[accepted_calls[j].id].second;

        for (int i = 0; i < j; i++) {
            std::string hseq_i = host_virus_bp_seqs[accepted_calls[i].id].first;
            std::string vseq_i = host_virus_bp_seqs[accepted_calls[i].id].second;

            std::string h_ref = hseq_i.length() >= hseq_j.length() ? hseq_i : hseq_j;
            std::string h_query = hseq_i.length() >= hseq_j.length() ? hseq_j : hseq_i;
            std::string v_ref = vseq_i.length() >= vseq_j.length() ? vseq_i : vseq_j;
            std::string v_query = vseq_i.length() >= vseq_j.length() ? vseq_j : vseq_i;

            StripedSmithWaterman::Alignment h_aln = align(aligner, h_ref, h_query), v_aln = align(aligner, v_ref, v_query);
            bool same_hseq = h_aln.query_end-h_aln.query_begin >= h_query.length()*0.8;
            bool same_vseq = v_aln.query_end-v_aln.query_begin >= v_query.length()*0.8;

            if (same_hseq && same_vseq) {
//                if (hseq_j.length() > hseq_i.length()) {
//                    std::cerr << "SWAPPING " << accepted_calls[i].id << " " << accepted_calls[j].id << std::endl;
//                    accepted_calls[i].host_bp = accepted_calls[j].host_bp;
//                    host_virus_bp_seqs[i].first = host_virus_bp_seqs[j].first;
//                }
//                if (vseq_j.length() > vseq_i.length()) {
//                    accepted_calls[i].virus_bp = accepted_calls[j].virus_bp;
//                    host_virus_bp_seqs[i].second = host_virus_bp_seqs[j].second;
//                }
                accepted_calls[j].removed = true;
                break;
            }
        }
    }


    if (print_rejected) {
        print_calls(rejected_calls);
    } else {
        print_calls(accepted_calls);
    }
}