#include <bits/stdc++.h>
#include <stdlib.h>
#include <pthread.h>
#include <stdint.h>
#include <gmp.h>
#define printf gmp_printf

using namespace std;

const int DATABITS = 1 << 20;

const double LOG_SECONDS = 1.0;

// coefficient with respect to the table size (2^27)
const double BATCH_COEF = 1.0/2;
const int JUMP_AFTER_BATCHES = 10; // >= 1

#ifdef P40
    // 40 bits
    // answer: 0x4e2dea1f3c
    const char *P_HEX = "ffffffffa9";
    const char *FNAME = "test.bin";

    const int TAB1 = 1 << 10;
    const int TAB2 = 1 << 8;
#endif

#ifdef P64
    // 64 bits
    // answer: 0x90644c931a3fba5
    const char *P_HEX = "ffffffffffffffc5";
    const char *FNAME = "0.bin";

    const int TAB1 = 1 << 14;
    const int TAB2 = 1 << 8;
#endif

#ifdef P74
    // 74 bits
    // answer: 0x384f17db02976dcf63d
    const char *P_HEX = "3ffffffffffffffffdd";
    const char *FNAME = "1.bin";

    const int TAB1 = 1 << 14;
    const int TAB2 = 1 << 8;
#endif

#ifdef P84
    // 84 bits
    // answer: ???
    const char *P_HEX = "fffffffffffffffffffdd";
    const char *FNAME = "2.bin";

    const int TAB1 = 1 << 14;
    const int TAB2 = 1 << 8;
#endif



bitset<DATABITS> PRF_output;
mpz_t mP;

// separate tables to improve locality
// (we access info rarely)
vector<uint64_t> table_keys;
vector<uint64_t> table_info;

uint64_t pack_info(int d, int off, int rev) {
    return (uint64_t(d) << 40) | (uint64_t(off) << 20) | uint64_t(rev);
}
tuple<int, int, int> unpack_info(uint64_t v) {
    int d = (v >> 40) & 0xfffff;
    int off = (v >> 20) & 0xfffff;
    int rev = (v >> 0) & 0xfffff;
    return {d, off, rev};
}


// 0 if is quadratic residue,
// 1 otherwise
uint8_t LegendreF2(mpz_t ma) {
    return (1 - mpz_jacobi(ma, mP)) >> 1;
}
uint8_t LegendreF2(uint64_t a) {
    mpz_t ma;
    mpz_init(ma);
    mpz_set_ui(ma, a);
    auto res = LegendreF2(ma);
    mpz_clear(ma);
    return res;
}
uint8_t LegendreF2precomp[TAB2+5] = {};

int FOUND_SOLUTION = 0;
void* worker(void * vdata);

int main(int argc, char *argv[]) {
    mpz_init(mP);
    mpz_set_str(mP, P_HEX, 16);

    FILE *fd = fopen(FNAME, "r");
    assert(fd);
    uint8_t data[DATABITS/8];
    fread(data, sizeof(data), 1, fd);
    fclose(fd);

    for (int i = 0; i < DATABITS/8; i++) {
        for (int j = 0; j < 8; j++) {
            uint64_t bit = (data[i] >> (7 - j)) & 1;
            PRF_output[(i<<3)|j] = bit^1;
        }
    }

    printf("=========================\n");
    printf("Prime p = 0x%Zx\n", mP);
    printf("TAB1: %d TAB2: %d\n", TAB1, TAB2);
    printf("NTHREADS: %d\n", NTHREADS);
    printf("BATCH_COEF: %.3lf\n", BATCH_COEF);
    printf("JUMP_AFTER_BATCHES: %d\n", JUMP_AFTER_BATCHES);
    printf("=========================\n");
    fflush(stdout);

    printf("\n-------------------------\n");
    printf("Step 1:\n");
    printf("-------------------------\n");
    vector<pair<uint64_t, uint64_t>> table;
    for (int b = 1; b <= TAB1; b++) {
        if (64 * b > DATABITS) break;

        uint64_t bflag = LegendreF2(b);
        for(int a = 0; a < b; a++) {
            uint64_t v = 0;
            uint64_t vr = 0;
            int j = a;
            for(int i = 0; i < 64; i++) {
                v = (v << 1) | (bflag ^ PRF_output[j]);
                vr = (vr >> 1) | ((bflag ^ PRF_output[j]) << 63);
                j += b;
            }
            int rev = (v > vr);
            uint64_t info = pack_info(b, a, rev);
            table.push_back({rev ? vr : v, info});
        }
    }
    for(int b = 1; b <= TAB2; b++) {
        LegendreF2precomp[b] = LegendreF2(b);
    }

    printf("table size: %lu (%.3f GB + %.3f GB extra)\n", table.size(), table.size()*8 / 1e9, table.size()*8 / 1e9);
    sort(table.begin(), table.end());
    for(auto kv: table) {
        table_keys.push_back(kv.first);
        table_info.push_back(kv.second);
    }
    printf("sorted\n");
    fflush(stdout);

    FOUND_SOLUTION = 0;

    printf("\n-------------------------\n");
    printf("Step 2:\n");
    printf("-------------------------\n");
    pthread_t threads[NTHREADS];
    int thread_id[NTHREADS];
    for (int i = 0; i < NTHREADS; i++) {
        thread_id[i] = i;
        assert(!pthread_create(&threads[i], NULL, worker, thread_id+i));
    }
    for (int i = 0; i < NTHREADS; i++) {
        pthread_join(threads[i], NULL);
    }
    return 0;
}


// annoying time funcs
void sleep_ms(long msec) {
    struct timespec ts;
    ts.tv_sec = msec / 1000;
    ts.tv_nsec = (msec % 1000) * 1000000;
    nanosleep(&ts, &ts);
}

chrono::system_clock::time_point time_now() {
    return chrono::high_resolution_clock::now();
}
chrono::system_clock::time_point t_start, t_log;

template<typename T>
double elapsed(T start) {
    chrono::duration<double, micro> ts = time_now() - start;
    return ts.count();
}

// batch membership (sorted arrays)
template<typename T>
pair<T, T> batch_contains(T qi, T qend, T pi, T pend) {
    while (pi != pend && qi != qend) {
        uint64_t vp = *pi;
        uint64_t vq = *qi;
        if (vp == vq)
            return {qi, pi};
        else if (vp < vq)
            pi++;
        else
            qi++;
    }
    return {qend, pend};
}

void fast_sort(vector<uint64_t> &vec);
void fast_sort(uint64_t *begin, uint64_t *end);
void recover(uint64_t val, mpz_t mGuess, uint64_t table_index, uint64_t itr_last, int d);

static thread_local uint64_t total = 1;
static thread_local uint64_t total_batches = 0;

void * worker(void * vdata) {
    int *data = (int*)vdata;
    int thread_id = *data;

    // desynchronize threads so they don't all access memory at the same time
    srand(rand() ^ time(0));
    srand(rand() ^ (uint64_t)elapsed(time_now()));
    srand(rand() ^ thread_id);
    // if (thread_id)
    //     sleep_ms(rand() % 10000);

    gmp_randstate_t mRand;
    gmp_randinit_default(mRand);
    gmp_randseed_ui(mRand, rand());

    mpz_t mGuess;
    mpz_init(mGuess);
    mpz_urandomm(mGuess, mRand, mP);

    printf("thread %02d guess %Zx\n", thread_id, mGuess);

    pair<uint64_t, uint64_t> seqs[TAB2+1][TAB2+1];

    vector<uint64_t> batch;
    t_start = time_now();
    t_log = time_now();
    for(uint64_t itr = 0; ; itr++, total++) {
        uint64_t jsbit = LegendreF2(mGuess);

        for(int b = 1; b <= TAB2; b++) {
            int a = itr % b;
            uint64_t & v = seqs[b][a].first;
            v = (v << 1) | (jsbit ^ LegendreF2precomp[b]);
            uint64_t & vr = seqs[b][a].second;
            vr = (vr >> 1) | ((jsbit ^ LegendreF2precomp[b]) << 63);

            // consider only full sequences
            // (ignore initial zero-filled parts)
            if (itr >= b * 64 + a) {
                batch.push_back( (v < vr) ? v : vr );
            }
            else {
                batch.push_back(0); // to keep easy index recovery
            }
        }

        mpz_add_ui(mGuess, mGuess, 1);

        // batch processing etc.

        if (batch.size() >= BATCH_COEF * table_keys.size()) {
            if (FOUND_SOLUTION) return 0;
            total_batches++;

            int print_log = 0;
            if (elapsed(t_log) > LOG_SECONDS * 1e6) {
                t_log = time_now();
                print_log = 1;
            }
            if (print_log) printf("thread %02d guess %Zx : start batch %.3f mil tab1 %.3f mil buf\n", thread_id, mGuess, table_keys.size()/1e6, batch.size()/1e6);

            vector<uint64_t> batch_orig = batch;
            fast_sort(batch);

            auto pairpos = batch_contains(table_keys.begin(), table_keys.end(), batch.begin(), batch.end());
            auto table_index = pairpos.first - table_keys.begin();
            auto batch_index = pairpos.second - batch.begin();
            if (batch_index < batch.size()) {
                uint64_t val = batch[batch_index];
                assert(table_keys[table_index] == val);

                printf("match after 0x%lx iterations : value 0x%lx\n", total, val);
                for(uint64_t i = 0; i < batch_orig.size(); i++) {
                    if (batch_orig[i] == val) {
                        uint64_t itr_last = batch_orig.size() / TAB2 - i / TAB2;
                        printf("itr offset %lu : itr_last 0x%lx mod %ld (i %lu)\n", i, itr_last, i % TAB2, i);
                        int b = i % TAB2 + 1;
                        recover(val, mGuess, table_index, itr_last, b);
                        if (FOUND_SOLUTION) return 0;
                    }
                }
            }
            batch.clear();

            double ms_per_call = elapsed(t_start) / ((double)total * TAB2);
            double ms_per_itr = elapsed(t_start) / ((double)total);
            if (print_log) printf("thread %02d guess %Zx : finish batch : 0x%lx iterations, %.7lf seconds per itr, %.7lf microseconds per check\n\n",
                thread_id, mGuess, total, ms_per_itr, ms_per_call);
            fflush(stdout);

            // re-randomize? jump
            if (total_batches && (total_batches % JUMP_AFTER_BATCHES) == 0) {
                mpz_urandomm(mGuess, mRand, mP);
                printf("thread %02d randomize : guess %Zx\n", thread_id, mGuess);
                itr = 0;
                itr--;
            }
        }
    }
    return NULL;
}

void recover(uint64_t val, mpz_t mGuess, uint64_t table_index, uint64_t itr_last, int b) {
    // get L0 position
    mpz_t base, coef, ptr;
    mpz_init(base);
    mpz_init(coef);
    mpz_init(ptr);

    mpz_t key;
    mpz_init(key);

    mpz_t chk;
    mpz_init(chk);

    mpz_set_ui(coef, b);
    mpz_invert(coef, coef, mP);

    mpz_sub_ui(base, mGuess, itr_last);
    mpz_mul(base, base, coef);

    mpz_sub_ui(base, base, 63);
    mpz_mod(base, base, mP);

    int rev;
    uint64_t v = 0, vr = 0;
    for (rev = 0; rev < 2; rev ++) {
        mpz_set(ptr, base);
        for (int i = 0; i < 64; i++) {
            uint64_t jsbit = LegendreF2(ptr);
            v = (v << 1) | jsbit;
            vr = (vr >> 1) | (jsbit << 63);
            mpz_add_ui(ptr, ptr, 1);
        }
        if (v == val) break;
        mpz_sub(base, mP, base);
        mpz_sub_ui(base, base, 63);
    }
    assert(v == val);
    printf("match L0: b=%d base=0x%Zx rev=%d\n", b, base, rev);

    // get Lk position
    int b2, a2, rev2;
    tie(b2, a2, rev2) = unpack_info(table_info[table_index]);

    mpz_mul_ui(key, base, b2);
    mpz_sub_ui(key, key, a2);
    mpz_mod(key, key, mP);

    if (rev2) {
        printf("match rev Lk: b=%d a=%d\n", b2, a2);

        mpz_sub(key, mP, key);
        mpz_sub_ui(key, key, 63 * b2 + a2 * 2);
    }
    else {
        printf("match     Lk: b=%d a=%d\n", b2, a2);
    }

    // check more bits
    uint64_t v2, v2r;
    mpz_set(chk, key);
    int i;
    for(i = 0; i < 256; i++) {
        uint8_t jsbit = LegendreF2(chk);
        if (PRF_output[i] != jsbit)
            break;
        mpz_add_ui(chk, chk, 1);
    }
    if (i == 256) {
        printf("\n=========================\n");
        printf("# final key: 0x%Zx\n", key);
        printf("=========================\n\n");
        FOUND_SOLUTION = 1;
        system("date >> SOLUTION");
        system("hostname >> SOLUTION");
        FILE *fd = fopen("SOLUTION", "a");
        gmp_fprintf(fd, "iterations %lu : key 0x%Zx\n", total, key);
        gmp_fprintf(fd, "time: %lf seconds\n\n", elapsed(t_start) / 1e6);
        fclose(fd);
        return;
    }
    else {
        printf("false positive: %d extra bits\n", i);
    }
    fflush(stdout);
}


static thread_local vector<uint64_t> vtemp;

void lsb_radix_sort(uint64_t *begin, uint64_t *end);

void fast_sort(uint64_t *begin, uint64_t *end) {
    uint64_t sz = end - begin;
    if (sz > vtemp.size())
        vtemp.resize(sz);
    return lsb_radix_sort(begin, end);
}
void fast_sort(vector<uint64_t> &vec) {
    uint64_t *begin = vec.data();
    uint64_t *end = begin + vec.size();
    return fast_sort(begin, end);
}

/*
Based on the answer by Andreas Kaseorg
https://www.quora.com/What-is-the-most-efficient-way-to-sort-a-million-32-bit-integers
*/
void lsb_radix_sort(uint64_t *begin, uint64_t *end) {
    uint64_t *begin1 = vtemp.data();
    uint64_t *end1 = begin1 + (end - begin);
    if (0) {
        for (uint64_t shift = 0; shift < 64; shift += 8) {
            size_t count[0x100] = {};
            for (uint64_t *p = begin; p != end; p++)
                count[(*p >> shift) & 0xFF]++;
            uint64_t *bucket[0x100], *q = begin1;
            for (int i = 0; i < 0x100; q += count[i++])
                bucket[i] = q;
            for (uint64_t *p = begin; p != end; p++)
                *bucket[(*p >> shift) & 0xFF]++ = *p;
            swap(begin, begin1);
            swap(end, end1);
        }
    }
    else {
        // minor optimization

        // sort only 32 MSBits
        // should be almost-sorted
        for (uint64_t shift = 32; shift < 64; shift += 8) {
            size_t count[0x100] = {};
            for (uint64_t *p = begin; p != end; p++)
                count[(*p >> shift) & 0xFF]++;
            uint64_t *bucket[0x100], *q = begin1;
            for (int i = 0; i < 0x100; q += count[i++])
                bucket[i] = q;
            for (uint64_t *p = begin; p != end; p++)
                *bucket[(*p >> shift) & 0xFF]++ = *p;
            swap(begin, begin1);
            swap(end, end1);
        }
        // fix few discrepancies using Bubble sort
        uint64_t *p = begin;
        uint64_t *endm1 = --end;
        while(p < endm1) {
            while (p[0] > p[1]) {
                swap(p[0], p[1]);
                if (p > begin)
                    p--;
                else
                    break;
            }
            p++;
        }
    }
}

