#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <htslib/sam.h>

#include "config.h"
#include "sam_utils.h"
#include <ssw.h>
#include <ssw_cpp.h>
#include "utils.h"

void print_calls(std::vector<call_t>& calls) {
    for (int i = 0; i < calls.size(); i++) {
        call_t call = calls[i];
        if (!call.removed) {
            std::cout << call.to_human_string() << std::endl;
        }
    }
}

StripedSmithWaterman::Alignment align(StripedSmithWaterman::Aligner& aligner, std::string& ref, std::string& query, bool align_rc) {
    std::string rc_query = query;
    get_rc(rc_query);

    StripedSmithWaterman::Alignment aln;
    StripedSmithWaterman::Filter filter;
    if (align_rc) {
        aligner.Align(rc_query.c_str(), ref.c_str(), ref.length(), filter, &aln, 0);
    } else {
        aligner.Align(query.c_str(), ref.c_str(), ref.length(), filter, &aln, 0);
    }
    return aln;
}

int main(int argc, char* argv[]) {

    std::string workdir = argv[1];

    config_t config = parse_config(workdir + "/config.txt");

    std::vector<std::string> args;
    bool print_rejected = false, accept_clipped = false, remapped = false;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--print-rejected") == 0) {
            print_rejected = true;
        } else if (strcmp(argv[i], "--accept-all-clipped") == 0) {
            accept_clipped = true;
        } else if (strcmp(argv[i], "--remapped") == 0) {
            remapped = true;
        } else {
            args.push_back(argv[i]);
        }
    }

    std::ifstream fin(workdir + (remapped ? "/results.remapped.txt" : "/results.txt"));

    double min_pairs = 4;
    if (args.size() >= 2) min_pairs = std::stoi(args[1]);

    double min_host_pbs = 0.8;
    if (args.size() >= 3) min_host_pbs = std::stod(args[2]);

    double min_cov_for_accept = 0.5;
    if (args.size() >= 4) min_cov_for_accept = std::stod(args[3]);


    // read host and virus sequences
    std::string host_bp_seqs_fname = workdir + "/host_bp_seqs.fa", virus_bp_seqs_fname = workdir + "/virus_bp_seqs.fa";
    FILE* host_bp_fasta = fopen(host_bp_seqs_fname.c_str(), "r"), *virus_bp_fasta = fopen(virus_bp_seqs_fname.c_str(), "r");
    kseq_t* hseq = kseq_init(fileno(host_bp_fasta)), *vseq = kseq_init(fileno(virus_bp_fasta));
    std::vector<std::pair<std::string, std::string> > host_virus_bp_seqs;
    std::vector<std::pair<int, int> > host_virus_bp_Ns;
    while (kseq_read(hseq) >= 0 && kseq_read(vseq) >= 0) {
        host_virus_bp_seqs.push_back({hseq->seq.s, vseq->seq.s});
        host_virus_bp_Ns.push_back({0, 0});
    }

    std::string host_bp_seqs_masks_fname = workdir + "/host_bp_seqs.masked.bed",
                virus_bp_seqs_masks_fname = workdir + "/virus_bp_seqs.masked.bed";
    std::ifstream host_bp_masks_fin(host_bp_seqs_masks_fname), virus_bp_masks_fin(virus_bp_seqs_masks_fname);
    std::vector<std::pair<int, int> > host_virus_bp_seqs_masked;
    std::string id;
    int start, end;
    while (host_bp_masks_fin >> id >> start >> end) {
        int n = std::stoi(id.substr(0, id.length()-2));
        host_virus_bp_Ns[n].first = end-start+1;
    }
    while (virus_bp_masks_fin >> id >> start >> end) {
        int n = std::stoi(id.substr(0, id.length()-2));
        host_virus_bp_Ns[n].second = end-start+1;
    }

    std::vector<call_t> accepted_calls, rejected_calls;
    std::string line;
    while (std::getline(fin, line)) {
        call_t call(line);

        bool accept = true;
        if (call.host_pbs < min_host_pbs) accept = false;
        if (call.good_pairs < min_pairs-1 || (call.good_pairs == min_pairs-1 && call.split_reads == 0)) accept = false;

        std::string host_seq = host_virus_bp_seqs[call.id].first;
        std::string virus_seq = host_virus_bp_seqs[call.id].second;
        if (host_seq.length() < config.read_len/2 || virus_seq.length() < config.read_len/2) {
            accept = false;
        }

        if (call.coverage() >= min_cov_for_accept && call.good_pairs > 2) accept = true;
        if (accept_clipped && call.split_reads > 0) accept = true;

        int host_Ns = host_virus_bp_Ns[call.id].first;
        int virus_Ns = host_virus_bp_Ns[call.id].second;
        if (double(host_Ns)/host_seq.length() > 0.80 || double(virus_Ns)/virus_seq.length() > 0.80) {
            accept = false;
        }

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

            // TODO: if the two calls are on opposite strand, align rc, otherwise, no
            bool h_rc = accepted_calls[i].host_bp.rev != accepted_calls[j].host_bp.rev;
            bool v_rc = accepted_calls[i].virus_bp.rev != accepted_calls[j].virus_bp.rev;
            StripedSmithWaterman::Alignment h_aln = align(aligner, h_ref, h_query, h_rc),
                                            v_aln = align(aligner, v_ref, v_query, v_rc);
            bool same_hseq = h_aln.query_end-h_aln.query_begin >= h_query.length()*0.8;
            bool same_vseq = v_aln.query_end-v_aln.query_begin >= v_query.length()*0.8;

            if (same_hseq && same_vseq) {
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
