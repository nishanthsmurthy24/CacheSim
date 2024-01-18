// Author: Priyank Faldu
// This code was submitted to
// THE 2ND CACHE REPLACEMENT CHAMPIONSHIP, 2017
// URL: http://crc2.ece.tamu.edu
// CRC2 2017 Submission #6
// Title: Reuse-Aware Management for Last-Level Caches

#ifndef _LEEWAY_NRU_H_
#define _LEEWAY_NRU_H_
#include <cstdlib>
#include <cassert>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <utility>
#include <iomanip>
#include <random>
#include <cstdint>
#include <cmath>

#define MAXIMUM_CORES 4
#define MAXIMUM_SETS MAXIMUM_CORES*2048

typedef enum {
    FOLLOWER_SET = 0,
    ROP_SET = 1,
    BOP_SET = 2,
    MAX_SETS = 3,
} SetRole; 

enum {
    ROP = ROP_SET-1,
    BOP = BOP_SET-1,
    REUSE_TYPES = MAX_SETS-1
};

struct SetType {
    uint8_t owner;
    SetRole type;
};


// This is used to identify whether a given set is a sampler set or a follower set
// It need not cost any storage as it will be hardwired in the logic
SetType* set_types;
int32_t follower[MAXIMUM_CORES];

// Registers
uint32_t random_counter[MAXIMUM_CORES][REUSE_TYPES] = {0}; // 9 bit saturating counter 
int32_t misses[MAXIMUM_CORES][MAXIMUM_CORES][REUSE_TYPES]; // 16 bit saturating counter
int32_t accesses[MAXIMUM_CORES][MAXIMUM_CORES][REUSE_TYPES]; // 16 bit saturating counter
uint32_t last_cycle_count = 0; // 27-bit saturating counter

// Constants
const int32_t max_nru_bits = 3; // 2-bit NRU
int32_t max_conservative_conf = 7;
int32_t max_aggressive_conf = 2;
int stat_threshold = 0;
uint32_t bypass_th[REUSE_TYPES] = {0};
uint32_t interval = 100000000;
const uint64_t prefetch_sig = 0xfeedfeed;
const uint64_t writeback_sig =0xbeefbeef;
int32_t num_sample_sets = 64;
uint32_t llc_sets;
int32_t num_cores;
const int32_t llc_assoc = 16;
uint32_t ldpt_sets;
uint32_t ldpt_assoc;
uint32_t config_number;

// Supporting stats which are not part of the implementation
uint64_t prediction_stat[REUSE_TYPES+1][max_nru_bits+2];
uint64_t bypass_applied = 0, bypass_skipped = 0;
uint64_t update_applied = 0, update_skipped = 0;
uint64_t invalid_blocks = 0;
uint64_t sampler_bp[REUSE_TYPES] = {0};
uint64_t sampler_update[REUSE_TYPES] = {0};
uint64_t sampler_read[REUSE_TYPES] = {0};

string policyStr(int p) {
    switch(p) {
        case ROP:
            return "Reuse-Oriented-Policy";
            break;
        case BOP:
            return "Bypass-Oriented-Policy";
            break;
        default:
            assert(0);
            return "Unknown Policy";
            break;
    }
    assert(0);
    return "Unknown Policy";
}

void CheckPolicyStatus() {
    if ( ( current_core_cycle[0] - last_cycle_count) >= interval ) {
        last_cycle_count = current_core_cycle[0];

        cout << endl;
        for ( int i = 0 ; i < num_cores; i++ ) {
            if ( i == 0 ) {
                cout << setw(15) << "ROP/BOP/DIFF";
            }
            for ( int j = 0 ; j < num_cores; j++ ) {
                if ( i == 0 ) {
                    cout << setw(28) << "CORE " << j << "\t";
                }
            }
            if ( i == 0  ) {
                cout << endl;
            }
            cout << setw(15) << i;
            for ( int j = 0 ; j < num_cores; j++ ) {
                cout << setw(9) << misses[i][j][ROP] << "/" << setw(9) << misses[i][j][BOP] << "/" << setw(9) << misses[i][j][ROP] - misses[i][j][BOP] << "\t";
            }
            cout << endl;
        }
        cout << endl;


        cout << endl;
        for ( int i = 0 ; i < num_cores; i++ ) {
            if ( i == 0 ) {
                cout << setw(15) << "ROP/BOP/DIFF";
            }
            for ( int j = 0 ; j < num_cores; j++ ) {
                if ( i == 0 ) {
                    cout << setw(28) << "CORE " << j << "\t";
                }
            }
            if ( i == 0  ) {
                cout << endl;
            }
            cout << setw(15) << i;
            for ( int j = 0 ; j < num_cores; j++ ) {
                cout << setw(9) << accesses[i][j][ROP] << "/" << setw(9) << accesses[i][j][BOP] << "/" << setw(9) << accesses[i][j][ROP] - accesses[i][j][BOP] << "\t";
            }
            cout << endl;
        }
        cout << endl;


        cout << endl;
        for ( int i = 0 ; i < num_cores; i++ ) {
            if ( i == 0 ) {
                cout << setw(15) << "ROP/BOP/DIFF";
            }
            for ( int j = 0 ; j < num_cores; j++ ) {
                if ( i == 0 ) {
                    cout << setw(28) << "CORE " << j << "\t";
                }
            }
            if ( i == 0  ) {
                cout << endl;
            }
            cout << setw(15) << i;
            for ( int j = 0 ; j < num_cores; j++ ) {
                cout << setw(9) << misses[i][j][ROP] * 100.0 / accesses[i][j][ROP] << "/" << setw(9) << misses[i][j][BOP] * 100.0 / accesses[i][j][BOP] << "/" << setw(9) << (misses[i][j][ROP] * 100.0 / accesses[i][j][ROP]) - (misses[i][j][BOP] * 100.0 / accesses[i][j][BOP]) << "\t";
            }
            cout << endl;
        }
        cout << endl;

        for ( int i = 0 ; i < num_cores; i++ ) {
            cout << "Thread: " << i << " Old Policy: " << policyStr(follower[i]);
            int64_t total[REUSE_TYPES];
            for ( int l = 0 ; l < REUSE_TYPES ; l++ ) {
                total[l] = 0;
            }
            for ( int j = 0 ; j < num_cores; j++ ) {
                total[ROP] += misses[i][j][ROP];
                total[BOP] += misses[i][j][BOP];
            }
            if ( total[ROP] < total[BOP] ) {
                //if ( misses[i][i][ROP] < misses[i][i][BOP] )
                follower[i] = ROP;
            } else {
                follower[i] = BOP;
            }
            cout << " New Policy: " << policyStr(follower[i]) << endl;
            cout << " total[ROP]: " << total[ROP] << " total[BOP]:" << total[BOP] << " diff: " << int64_t(total[ROP] - total[BOP]) << endl;
        }
        for ( int i = 0 ; i < num_cores; i++ ) {
            for ( int j = 0 ; j < num_cores; j++ ) {
                misses[i][j][ROP] = 0;
                misses[i][j][BOP] = 0;
                accesses[i][j][ROP] = 0;
                accesses[i][j][BOP] = 0;
            }
        }
    }
}

struct LeewayBlock {
    uint64_t signature; // 22-bits
    int32_t nru_bits; // 2-bits {0, 1, 2, 3}
    int32_t current_ld; // 2-bits; {-1, 0, 1, 2}; -1 if no hits; hit on 2 or 3 is forced to 2
    int32_t predicted_ld; /// 2-bits; only possible values {-1, 0, 1, 2}  if live distance is 2 or 3, it is forced to 2
    int32_t inst_type; // 1-bit; In multicore setting, it is used to identify which policy to update
    // Total 29 bits per block if sampler block
    // 4-bit otherwise
    LeewayBlock() : signature(0), nru_bits(max_nru_bits), current_ld(-1), predicted_ld(max_nru_bits), inst_type(REUSE_TYPES) {}
};

// This function assigns roles to different sets
// Every set is either a follower set or leader set for some core 
void RandomSetSelection(SetType* sets, int32_t total_sets, int32_t num_threads, int32_t num_types, int32_t sets_to_sample) {
    cout << "Selecting " << sets_to_sample << " sets for " << num_threads << " threads for " << num_types << " types" << " from " << total_sets << " sets " << endl; 
    assert ( num_types * num_threads * sets_to_sample <= total_sets);
    for ( int t = 0 ; t < num_threads ; t++ ) {
        follower[t] = ROP;
    }
    for ( int i = 0 ; i < total_sets ; i ++ ) {
        sets[i].type = FOLLOWER_SET;
        sets[i].owner = num_threads;
    }
    int total_sets_to_sample = num_types * num_threads * sets_to_sample;
    assert(total_sets_to_sample <= total_sets);
    std::default_random_engine generator(15485863);
    std::uniform_int_distribution<int> distribution(0,total_sets-1);
    for (int next_type = 1; next_type <= num_types ; next_type++) {
        for (int thread_t = 0; thread_t < num_threads ; thread_t++) {
            int32_t remaining_sets_to_sample = sets_to_sample;
            while (remaining_sets_to_sample > 0 ) {
                int local_i = distribution(generator);
                if ( sets[local_i].type == FOLLOWER_SET ) {
                    sets[local_i].type = (SetRole)next_type;
                    sets[local_i].owner = thread_t;
                    remaining_sets_to_sample--;
                } else {
                    continue;
                }
            }
        }
    }

    cout << setw(10) << "Thread"
        << setw(10) << "Type"
        << setw(10) << "Count"
        << endl;
    for ( int t = 0 ; t <= num_threads ; t++ ) {
        for ( int type = 0 ; type <= num_types ; type++ ) {
            int count = 0;
            for ( int s = 0 ; s < total_sets ; s++ ) {
                if ( sets[s].type == type && sets[s].owner == t) {
                    count++;
                }
            }
            cout << setw(10) << t
                << setw(10) << type
                << setw(10) << count
                << endl;
        }
    }
}

struct LDPT_ENTRY {
    int32_t stable_ld; // 2-bit
    int32_t variance_conf; // 3-bit
    int32_t variance_dir; // 1-bit
    LDPT_ENTRY() : stable_ld(0), variance_conf(0), variance_dir(0) {}
};

struct Meta {
    uint32_t lru; // 2-bit
    uint32_t tag; // 13-bit, less for multi-core
    bool valid; // 1-bit
    LDPT_ENTRY* ldpt_entry; // 12-bit
    // total 28-bits
};
class LDPT {
    Meta** meta; // NUM_CORE * 512 * 4 * 28-bits
    const uint32_t instances;
    const uint32_t sets;
    const uint32_t ways;
    uint64_t hits;
    uint64_t misses;
    int32_t* inc_threshold;
    int32_t* dec_threshold;
    uint32_t set_bits;
    uint32_t tag_bits;
    uint32_t address_bits;
    uint32_t tag_mask;
    uint32_t set_mask;
    uint32_t address_mask;
    public:
    LDPT(uint32_t _instances, uint32_t _sets, uint32_t _ways) : instances(_instances), sets(_sets), ways(_ways), hits(0), misses(0) {
        inc_threshold = new int32_t[instances];
        dec_threshold = new int32_t[instances];
        assert(inc_threshold);
        assert(dec_threshold);
        for ( uint32_t i = 0 ; i < instances ; ++i ) {
            inc_threshold[i] = 0;
            dec_threshold[i] = 0;
        }
        meta = new Meta*[sets];
        assert(meta);
        for ( uint32_t j = 0 ; j < sets ; j++ ) {
            meta[j] = new Meta[ways];
            assert(meta[j]);
            for ( uint32_t k = 0 ; k < ways ; k++ ) {
                meta[j][k].lru = k;
                meta[j][k].tag = 0;
                meta[j][k].valid = false;
                meta[j][k].ldpt_entry = new LDPT_ENTRY[instances];
                assert(meta[j][k].ldpt_entry);
                for ( uint32_t i = 0 ; i < instances; i++ ) {
                    init(meta[j][k].ldpt_entry[i]);
                }
            }
        }
        set_bits = log2(sets);
        set_mask = (1 << set_bits)-1;
        address_bits = 22;
        address_mask = 0;
        for ( unsigned int i = 0 ; i < address_bits; i++ ) { 
            address_mask <<= 1;
            address_mask = address_mask | 1;
        }
        tag_bits = address_bits - set_bits;
        tag_mask = (1 << tag_bits)-1;

        cout << "LDPT sets: " << sets << " assoc:" << ways << " instances: " << instances << endl;
        cout << "address_bits: " << address_bits
            << " set_bits: " << set_bits
            << " tag_bits: " << tag_bits << endl;
        cout << hex << "address_mask: " << address_mask
            << " set_mask: " << set_mask
            << " tag_mask: " << tag_mask << dec << endl;

        get_storage_state();

#if 0
        for ( uint32_t i = 0 ; i < 10000 ;  ) {
            cout << hex << setw(10) << mask_address(i) << setw(10) << get_set(mask_address(i)) << setw(10) << get_tag(mask_address(i)) << dec << endl;
            i++;
            cout << hex << setw(10) << mask_address(i) << setw(10) << get_set(mask_address(i)) << setw(10) << get_tag(mask_address(i)) << dec << endl;
            i += sets -1;
        }
#endif
    }
    ~LDPT() {
        for ( unsigned int i = 0 ; i < sets ; i++ ) {
            for ( unsigned int j = 0 ; j < ways ; j++ ) {
                for ( unsigned int k = 0 ; k < instances ; j++ ) {
                    delete[] meta[i][j].ldpt_entry;
                }
                delete[] meta[i];
            }
            delete[] meta;
        }
    }
    uint32_t mask_address(uint32_t addr)  {
        return ( addr ) & ( address_mask);
    }
    uint32_t get_set(uint32_t addr) {
        return  ( addr ) & ( set_mask);
    }
    uint32_t get_tag(uint32_t addr) {
        return (addr >> set_bits);
    }
    void setThreshold(int32_t inst, int32_t inc, int32_t dec) {
        inc_threshold[inst] = inc;
        dec_threshold[inst] = dec;
        cout << " Setting threshold for instance : " << inst << " inc: " << inc_threshold[inst] << " dec: " << dec_threshold[inst] << endl;
    }
    void init(LDPT_ENTRY& e) {
        e.stable_ld = max_nru_bits;
        e.variance_conf = 0;
    }
    void updateEntry(uint64_t sig, int32_t current_ld, uint32_t inst) {
        if ( inst >= instances ) {
            cout << " " << inst << " " << instances << endl;
            assert(0);
        }
        LDPT_ENTRY* entry = findEntry(sig, inst); 
        if ( entry->stable_ld == current_ld ) {
            //entry->conf--;
            entry->variance_conf = 0;
        } else {
            if ( current_ld > entry->stable_ld ) {
                if ( entry->variance_dir == 1 ) {
                    entry->variance_dir = 0;
                    entry->variance_conf = 1;
                } else {
                    entry->variance_conf++;
                }
                if ( entry->variance_conf >= inc_threshold[inst] ) {
                    entry->variance_conf = 0;
                    entry->stable_ld = current_ld;
                }
            } else {
                if ( entry->variance_dir == 0 ) {
                    entry->variance_dir = 1;
                    entry->variance_conf = 1;
                } else {
                    entry->variance_conf++;
                }
                if ( entry->variance_conf >= dec_threshold[inst] ) {
                    entry->variance_conf = 0;
                    entry->stable_ld = current_ld;
                }
            }
        }
    }
    LDPT_ENTRY* findEntry(uint64_t sig, uint32_t inst) {
        if ( inst >= instances ) {
            cout << " " << inst << " " << instances << endl;
            assert(0);
        }
        return &(meta[0][0].ldpt_entry[0]);
        uint32_t addr = mask_address(sig);
        uint32_t set_index = get_set(addr);
        uint32_t tag = get_tag(addr);
        uint32_t lru_index = -1;
        assert(set_index < sets);
        assert(set_index > 0);
        for ( uint32_t i = 0 ; i < ways ; i++ ) {
            if ( meta[set_index][i].tag == tag && meta[set_index][i].valid) {
                hits++;
                for ( uint32_t k = 0 ; k < ways ; k++ ) {
                    if ( meta[set_index][k].lru < meta[set_index][i].lru ) {
                        meta[set_index][k].lru++;
                    }
                }
                assert(meta[set_index][i].lru < ways);
                assert(meta[set_index][i].lru >= 0);
                meta[set_index][i].lru = 0;
                return &(meta[set_index][i].ldpt_entry[inst]);
            }
            if ( meta[set_index][i].lru == (ways-1) ) {
                assert(lru_index == (static_cast<uint32_t>(-1)));
                lru_index = i;
            }
        }

        assert(lru_index >= 0);
        assert(lru_index < ways);
        for ( uint32_t i = 0 ; i < instances ; i++ ) {
            init(meta[set_index][lru_index].ldpt_entry[i]);
        }
        for ( uint32_t i = 0 ; i < ways; i++ ) {
            if ( meta[set_index][i].lru < meta[set_index][lru_index].lru ) {
                meta[set_index][i].lru++;
            }
        }
        meta[set_index][lru_index].lru = 0;
        meta[set_index][lru_index].tag = tag;
        meta[set_index][lru_index].valid = true;
        misses++;
        return &(meta[set_index][lru_index].ldpt_entry[inst]);
    }

    int32_t getMyInstanceId(int32_t thread_index, int32_t set_index) {
        // Identify policy for a given set
        if ( (thread_index < 0) || (thread_index >= num_cores) ) {
            cout << thread_index << " " << num_cores << endl;
            assert(0);
        }
        if ( set_index < 0 || set_index >= (int32_t)llc_sets ) {
            cout << set_index << " " << llc_sets << endl;
            assert(0);
        }
        if ( set_types[set_index].owner == thread_index && (set_types[set_index].type == ROP_SET || set_types[set_index].type == BOP_SET) ) {
            return set_types[set_index].type-1;
        } else {
            return follower[thread_index];
        }
    }
    void display_stats() {
        cout << " ldpt hits: " << setw(12) << hits
            << " misses: " << setw(12) <<  misses
            << " hit-rate: " << setw(12) << fixed << setprecision(4) << hits * 100.0 / (hits + misses) << "%" << endl;
    }

    void get_storage_state() {
#if 0
        cout << setw(80) << "Number of LDPT entries: " << setw(20) << sets * ways << endl;
        cout << setw(80) << "Fields of LDPT" << endl; 
        int total_bits = 0;
        cout << setw(80) << "lru: " << setw(20) << log2(ways) << endl;
        total_bits += log2(ways);
        cout << setw(80) << "tag: " << setw(20) << tag_bits << endl;
        total_bits += tag_bits;
        cout << setw(80) << "valid: " << setw(20) << 1 << endl;
        total_bits += 1;
        int ldpt_entry  = 0;
        cout << setw(80) << "ldpt_entry[stable_ld]: " << setw(20) << log2(max_nru_bits+1) + 1<< endl;
        ldpt_entry += log2(max_nru_bits+1) + 1;
        cout << setw(80) << "ldpt_entry[variance_conf]: " << setw(20) << log2(max_conservative_conf+1) << endl;
        ldpt_entry += log2(max_conservative_conf+1);
        cout << setw(80) << "ldpt_entry[variance_dir]: " << setw(20) << 1 << endl;
        ldpt_entry += 1;
        cout << setw(80) << "LDPT entry size" << setw(20) << total_bits + ( ldpt_entry * instances)  << endl; 

        cout << setw(80) << "Total LDPT size ( bit) = " << setw(20) << (total_bits + (ldpt_entry * instances) ) * sets * ways << endl;
        cout << setw(80) << "Total LDPT size (kbytes) = " << setw(20) << (total_bits + (ldpt_entry * instances) ) * sets * ways  / ( 1024.0 * 8) << endl;

#endif
    }
};

// Following functions are used from source code of 
// Dead Block Replacement and Bypass with a Sampling Predictor by 
// D. Jimenez (University of Texas at San Antonio, USA)
// URL: http://www.jilp.org/jwac-1/online/code/001_jiminez.tgz
uint64_t mix (uint64_t a, uint64_t b, uint64_t c) {
    a=a-b;  a=a-c;  a=a^(c >> 13);
    b=b-c;  b=b-a;  b=b^(a << 8);
    c=c-a;  c=c-b;  c=c^(b >> 13);
    return c;
}

uint64_t f1 (uint64_t x) {
    return mix (0xfeedface, 0xdeadb10c, x);
}

uint64_t f2 (uint64_t x) {
    return mix (0xc001d00d, 0xfade2b1c, x);
}

uint64_t fi (uint64_t x, int i) {
    return f1 (x) + (f2 (x) >> i);
}

unsigned int get_pc_trace (uint64_t PC, int32_t thread) {
    unsigned int trace_bits = 22;
    uint64_t trace_mask = ( 1 << trace_bits ) -1;
    uint32_t t = ((PC & 0xFFFFFFFFFF) != 0);
    uint64_t trace = PC;
    uint64_t x = fi (trace ^ (t << 2), thread);
    uint32_t hash = x & trace_mask;
    return hash;
}

#endif
