#include <errno.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <wordexp.h>
#include "shimmer.h"
#include "khash.h"
#include "kvec.h"
#include "kalloc.h"

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

extern char *optarg;
extern int optind, opterr, optopt;

#define handle_error(msg) \
	do { perror(msg); exit(EXIT_FAILURE); } while (0)

#define MMER_COUNT_LOWER_BOUND 2
#define MMER_COUNT_UPPER_BOUND 72
#define ORIGINAL 0
#define REVERSED 1
#define READ_END_FUZZINESS 48
#define LOCAL_OVERLAP_UPPERBOUND 72
#define BESTN 2
#define OVERLAP 0
#define CONTAINMENT 1

typedef struct {
	mm128_v * mmers;
	void * mmer0_map;
	void * rlmap;
	void * mcmap;} py_mmer_t;

KHASH_MAP_INIT_INT64(RPAIR, uint8_t);

uint8_t rmap[] = {
	0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
	16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
	32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
	48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
	64,  84,  66,  71,  68,  69,  70,  67,  72,  73,  74,  75,  76,  77,  78,  79,
	80,  81,  82,  83,  65,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,
	96, 116,  98, 103, 100, 101, 102,  99, 104, 105, 106, 107, 108, 109, 110, 111,
	112, 113, 114, 115,  97, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127,
	128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
	144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
	160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
	176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
	192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
	208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
	224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
	240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255
};



void reverse_complement(char * seq, size_t len) {
	size_t p, rp;
	char tmp;
	p = 0;
	for (;;) {
		rp = len - p - 1;
		if ((rp - p)  <= 1) break;
		tmp = seq[rp];
		seq[rp] = rmap[(uint8_t)seq[p]];
		seq[p] = rmap[(uint8_t)tmp];
		p ++;
	}
}

int rp128_comp(const void * a, const void * b) {
	rp128_t * a0 = (rp128_t *) a;
	rp128_t * b0 = (rp128_t *) b;
	return ((a0->y0 & 0xFFFFFFFF) >> 1) < ((b0->y0 & 0xFFFFFFFF) >> 1);
}


void build_shimmer_map4py(py_mmer_t * py_mmer,        
		char * seqdb_prefix, 
		char * shimmer_prefix, 
		uint32_t mychunk, 
		uint32_t total_chunk,
		uint32_t lowerbound,
		uint32_t upperbound) {

	char mmc_file_path[8192];
	char mmer_file_path[8192]; 
	char seq_idx_file_path[8192];
	char seqdb_file_path[8291];

	wordexp_t p; 
	char **mmc_fns; 
	char **shimmer_fns;

	mm128_v mmers_;
	mm_count_v mmc;
	
	khash_t(RLEN) *rlmap_; 
	khash_t(MMC) *mcmap_ = kh_init(MMC);
	khash_t(MMER0) * mmer0_map_;

	assert(total_chunk > 0);
	assert(mychunk > 0 && mychunk <= total_chunk);
	
	if (seqdb_prefix == NULL) {
		seqdb_prefix = (char *) calloc(8192, 1);
		snprintf( seqdb_prefix, 8191, "seq_dataset" );
	}

	if (shimmer_prefix == NULL) {
		seqdb_prefix = (char *) calloc(8192, 1);
		snprintf( shimmer_prefix, 8191, "shimmer-L2" );
	}

	

	int written;
	written = snprintf(seq_idx_file_path, sizeof(seq_idx_file_path), "%s.idx", seqdb_prefix);
	assert(written < sizeof(seq_idx_file_path));
	fprintf(stderr, "using index file: %s\n", seq_idx_file_path);

	rlmap_ = get_read_length_map(seq_idx_file_path);
	
	written = snprintf(seqdb_file_path, sizeof(seqdb_file_path), "%s.seqdb", seqdb_prefix);
	assert(written < sizeof(seqdb_file_path));
	fprintf(stderr, "using seqdb file: %s\n", seqdb_file_path);
	
	py_mmer->mmers = malloc(sizeof(mm128_v));
	py_mmer->mmers->n=0;
	py_mmer->mmers->m=0;
	py_mmer->mmers->a=0;

	written = snprintf(mmer_file_path, sizeof(mmer_file_path), "%s-[0-9]*-of-[0-9]*.dat", shimmer_prefix);
	assert(written < sizeof(mmer_file_path));
	wordexp(mmer_file_path, &p, 0);
	shimmer_fns = p.we_wordv;
	for (int i = 0; i < p.we_wordc; i++) {
		fprintf(stderr, "useing shimmer data file: %s\n", shimmer_fns[i]);
		mmers_ = read_mmlist(shimmer_fns[i]);
		append_mmlist(py_mmer->mmers, &mmers_);
		kv_destroy(mmers_);
	}
	wordfree(&p);	


	written = snprintf(mmc_file_path, sizeof(mmc_file_path), "%s-MC-[0-9]*-of-[0-9]*.dat", shimmer_prefix);
	assert(written < sizeof(mmc_file_path));
	wordexp(mmc_file_path, &p, 0);
	mmc_fns = p.we_wordv;
	for (int i = 0; i < p.we_wordc; i++) {
		fprintf(stderr, "using shimmer count file: %s\n", mmc_fns[i]);
		mmc = read_mm_count(mmc_fns[i]);
		aggregate_mm_count(mcmap_, &mmc);
		kv_destroy(mmc);
	}

	wordfree(&p);	

	mmer0_map_ = kh_init(MMER0);

	build_map(
			py_mmer->mmers,
			mmer0_map_,
			rlmap_,
			mcmap_,
			mychunk,
			total_chunk,
			lowerbound,
			upperbound);
	py_mmer->mmer0_map = (void *) mmer0_map_;
	py_mmer->rlmap = (void *) rlmap_;
	py_mmer->mcmap = (void *) mcmap_;
}

void process_overlaps(py_mmer_t * py_mmer) {

	khash_t(MMER0) * mmer0_map_ = (khash_t(MMER0) *) py_mmer->mmer0_map; 
	//khash_t(RLEN) * rlmap_ = (khash_t(RLEN) *) rlmap_;
	//khash_t(MMC) * mcmap_ = (khash_t(MMC) *) mcmap_; 

	rp128_v * rpv;
	uint64_t mhash0, mhash1;

	khash_t(MMER1) * mmer1_map;
    uint32_t count = 0;

	for (khiter_t __i = kh_begin(mmer0_map_); __i != kh_end(mmer0_map_); ++__i) {
		if (!kh_exist(mmer0_map_,__i)) continue;
		mhash0 = kh_key(mmer0_map_, __i);
	    mhash0 >>= 8;	
		printf("%09ld\n", mhash0);
		mmer1_map = kh_val(mmer0_map_, __i);
		for (khiter_t __j = kh_begin(mmer1_map); __j != kh_end(mmer1_map); ++__j) {
		    if (!kh_exist(mmer1_map,__j)) continue;
			mhash1 = kh_key(mmer1_map, __j);
			mhash1 >>= 8;
			rpv = kh_val(mmer1_map, __j);
			//if (rpv->n <= 2 || rpv->n > LOCAL_OVERLAP_UPPERBOUND) continue;
			//qsort(rpv->a, rpv->n, sizeof(rp128_t), rp128_comp);
		    printf("X %lu %lu %lu\n", mhash0, mhash1, rpv->n);
			count++;
		}
	}
	printf("X2: %u\n", count);
}
