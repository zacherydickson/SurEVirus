#include <iostream>
#include <unordered_set>
#include <mutex>
#include <bitset>
#include <unistd.h>
#include <zlib.h>
#include <htslib/kseq.h>
#include <cptl_stl.h>

#define USE_BITSET


KSEQ_INIT(gzFile, gzread)

#include "sam_utils.h"
#include "config.h"

typedef unsigned long long ull;

static const int KMER_LEN = 18;
static const int KMER_STR_LEN = KMER_LEN + 1;
const int NUMBER_OF_SEGS = 6;

const int KMER_BITS = KMER_LEN * 2;

const int SEG_LEN = KMER_LEN/NUMBER_OF_SEGS;
const int SEG_BITS = SEG_LEN * 2;

const ull KMER_MASK = (1ll << KMER_BITS)-1;

const int MASKED_KMER_LEN = KMER_LEN - SEG_LEN;
const int MASKED_KMER_BITS = MASKED_KMER_LEN * 2;
const ull MASKED_KMER_MASK = (1ll << MASKED_KMER_BITS)-1;

const int MAX_READS_TO_PROCESS = 10000;


ull nucl_bm[256] = { 0 };
char bm_nucl[4] = { 'A', 'C', 'G', 'T' };
char nucl2chr[16];

#ifdef USE_BITSET
std::bitset<(1ll << MASKED_KMER_BITS)> segs[NUMBER_OF_SEGS];
#else
std::unordered_set<ull>  segs[NUMBER_OF_SEGS];
#endif

inline bool check(ull masked_kmer, int seg_n) {
#ifdef USE_BITSET
    return segs[seg_n].test(masked_kmer);
#else
    return segs[seg_n].count(masked_kmer);
#endif
}

inline void insert(ull masked_kmer, int seg_n) {
#ifdef USE_BITSET
    segs[seg_n].set(masked_kmer);
#else
    segs[seg_n].insert(masked_kmer);
#endif
}

std::mutex mtx, mtx_out;
config_t config;
std::ofstream RetainedFQStream1, RetainedFQStream2;

std::string print(ull kmer, int len) {
    char s[KMER_STR_LEN];
    s[len] = '\0';
    while (len > 0) {
        len--;
        s[len] = bm_nucl[kmer%4];
        kmer /= 4;
    }
    return s;
}

inline ull mask(ull kmer, int seg_n) {
    ull first_n_segs_mask = (1ll << SEG_BITS*seg_n)-1;
    ull first_n_segs = kmer & first_n_segs_mask;
    ull remaining_segs_shifted = (kmer >> SEG_BITS) & ~first_n_segs_mask;
    return (first_n_segs | remaining_segs_shifted) & MASKED_KMER_MASK;
}

inline bool valid_kmer(ull kmer, int len) {
    int count[256];
    count[int('A')] = count[int('C')] = count[int('G')] = count[int('T')] = 0;
    for (int i = 0; i < len; i++) {
        count[int(bm_nucl[kmer%4])]++;
        kmer /= 4;
    }

    // filter poly-(ACGT)
    int max_freq = std::max(std::max(count[int('A')], count[int('C')]), std::max(count[int('G')], count[int('T')]));
    if (max_freq >= len-2) return false;
    return true;
}


void index_seq(char* seq, size_t len) {
    ull kmer = 0;
    for (int i = 0; i < len; i++) {
        ull nv = nucl_bm[int(seq[i])];
        kmer = ((kmer << 2) | nv) & KMER_MASK;

        if (i >= KMER_LEN-1) {
            for (int j = 0; j < NUMBER_OF_SEGS; j++) {
                ull seg = mask(kmer, j);

                if (valid_kmer(seg, MASKED_KMER_LEN)) {
                    insert(seg, j);
                }
            }
        }
    }
}

char* kstring_to_cstr(kstring_t kstring) {
    char* cstr = (char*) malloc(kstring.l+1);
    strncpy(cstr, kstring.s, kstring.l);
    cstr[kstring.l] = '\0';
    return cstr;
}

struct read_t {
    char* name;
    char* seq;
    char* qual;

    read_t(kseq_t* kseq) {
        name = kstring_to_cstr(kseq->name);
        seq = kstring_to_cstr(kseq->seq);
        qual = kstring_to_cstr(kseq->qual);
    }
    read_t(const read_t & other){
        size_t nameSize = std::strlen(other.name)+1;
        this->name = (char*) std::malloc(sizeof(char) * nameSize);
        std::strncpy(this->name,other.name,nameSize);
        size_t seqSize = std::strlen(other.seq)+1;
        this->seq = (char*) std::malloc(sizeof(char) * seqSize);
        std::strncpy(this->seq,other.seq,seqSize);
        size_t qualSize = std::strlen(other.qual)+1;
        this->qual = (char*) std::malloc(sizeof(char) * qualSize);
        std::strncpy(this->qual,other.qual,qualSize);
    }
    ~read_t(){
        free(name);
        free(seq);
        free(qual);
    }
};
typedef std::pair<read_t, read_t> ReadPair_t;
typedef std::vector<ReadPair_t> ReadBlock_t;

bool is_virus_read(read_t read) {
    ull kmer = 0;
    int hit = 0, len = 0;

    const int seq1_len = strlen(read.seq);
    for (int i = 0; i < seq1_len; i++) {
        ull nv = nucl_bm[int(read.seq[i])];
        kmer = ((kmer << 2) | nv) & KMER_MASK;
        len++;

        if (len >= KMER_LEN) {
            print(kmer, KMER_LEN);
            for (int j = 0; j < NUMBER_OF_SEGS; j++) {
                ull seg = mask(kmer, j);
                if (check(seg, j)) {
                    hit++;
                    break;
                }
            }
        }

        if (hit >= 2) break;
    }

    return hit >= 2;
}

void outputReadPair(const ReadPair_t & rp){
    RetainedFQStream1 << "@" << rp.first.name << std::endl;
    RetainedFQStream1 << rp.first.seq << std::endl;
    RetainedFQStream1 << "+" << rp.first.name << std::endl;
    RetainedFQStream1 << rp.first.qual << std::endl;

    RetainedFQStream2 << "@" << rp.second.name << std::endl;
    RetainedFQStream2 << rp.second.seq << std::endl;
    RetainedFQStream2 << "+" << rp.second.name << std::endl;
    RetainedFQStream2 << rp.second.qual << std::endl;
}

//Processes a block of read pairs and identifies those which may be viral
//Input - an id (for use by thread_pool)
//      - a constant reference to a block of read pairs to process
//      - a reference to a block of read pairs to which the identified reads may be added
//Output: None, Modifies the to_write read block
void isolate(int id, const ReadBlock_t & read_pairs, ReadBlock_t & to_write) {
    for (const ReadPair_t& rp : read_pairs) {
        if (is_virus_read(rp.first) || is_virus_read(rp.second)) {
            to_write.emplace_back(rp);
        }
    }
}


int main(int argc, char* argv[]) {

    nucl_bm[int('A')] = 0;
    nucl_bm[int('C')] = 1;
    nucl_bm[int('G')] = 2;
    nucl_bm[int('T')] = 3;
    nucl_bm[int('N')] = 0;

    std::string fq1_fname = argv[1];
    std::string fq2_fname = argv[2];
    std::string host_ref = argv[3];
    std::string virus_ref = argv[4];
    std::string workdir = argv[5];
    std::string workspace = argv[6];

    gzFile fastaf = gzopen(virus_ref.c_str(), "r");
    kseq_t* seq = kseq_init(fastaf);
    while (kseq_read(seq) >= 0) {
        for (int i = 0; i < seq->seq.l; i++) {
            seq->seq.s[i] = toupper(seq->seq.s[i]);
        }
        index_seq(seq->seq.s, seq->seq.l);
        get_rc(seq->seq.s, seq->seq.l);
        index_seq(seq->seq.s, seq->seq.l);
    }
    kseq_destroy(seq);
    gzclose(fastaf);

    nucl2chr[1] = 'A'; nucl2chr[2] = 'C'; nucl2chr[4] = 'G'; nucl2chr[8] = 'T'; nucl2chr[15] = 'N';

    config = parse_config(workdir + "/config.txt");

    gzFile fq1f = gzopen(fq1_fname.c_str(), "r");
    gzFile fq2f = gzopen(fq2_fname.c_str(), "r");
    kseq_t* seq1 = kseq_init(fq1f);
    kseq_t* seq2 = kseq_init(fq2f);
    RetainedFQStream1.open(workspace + "/retained-pairs_1.fq");
    RetainedFQStream2.open(workspace + "/retained-pairs_2.fq");

    ctpl::thread_pool thread_pool(config.threads);
    //Multithreading is achieved by maintaining two global queues one containing reads to process
    // The other containing processed reads
    //An empty read block is created in the Output queue
    //A Block of reads is read from the input fastq files and stored in the Input queue
    //A job of processing the block is passed off to the thread pool
    //This continues until all reads have been read in and passed off
    //Processing waits for the blocks in order, frees the Input block, and outputs then frees the Output block
    //This ensures reads are output in the same order they were input
    std::vector<std::future<void> > futures;
    std::queue<ReadBlock_t> toWrite;
    std::queue<ReadBlock_t> toProcess;
    while(kseq_read(seq1) >= 0 && kseq_read(seq2) >= 0){
        toWrite.emplace();
        toProcess.emplace();
        ReadBlock_t & block = toProcess.back();
        block.push_back({read_t(seq1),read_t(seq2)});
        for (int i = 1; i < MAX_READS_TO_PROCESS && kseq_read(seq1) >= 0 && kseq_read(seq2) >= 0; i++) {
            read_t r1(seq1), r2(seq2);
            block.push_back({r1, r2});
        }
        std::future<void> future = thread_pool.push(isolate, std::cref(block), std::ref(toWrite.back()));
        futures.push_back(std::move(future));
    }
    for (int i = 0; i < futures.size(); i++) {
        try {
            futures[i].get();
            //mtx.lock();
            toProcess.pop();
            //mtx.unlock();
            ReadBlock_t block = toWrite.front();
            for (const ReadPair_t & rp : block){
                outputReadPair(rp);
            }
           //mtx_out.lock();
           toWrite.pop();
           //mtx_out.unlock();
        } catch (char const* s) {
            std::cout << s << std::endl;
        }
    }

    kseq_destroy(seq1);
    kseq_destroy(seq2);
    gzclose(fq1f);
    gzclose(fq2f);
    RetainedFQStream1.close();
    RetainedFQStream2.close();
}
