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

int mp128_comp(const void * a, const void * b) {
	mp128_t * a0 = (mp128_t *) a;
	mp128_t * b0 = (mp128_t *) b;
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


uint32_t get_mmer_count( py_mmer_t * py_mmer, uint64_t mhash ) {
	khash_t(MMC) *mcmap_ = (khash_t(MMC) *) py_mmer->mcmap;
	khiter_t k = kh_get(MMC, mcmap_, mhash);
	if (k != kh_end(mcmap_)) {
		return kh_val(mcmap_, k);
	} else {
		return 0;
	}

}


typedef struct { uint64_t x0, x1, y0, y1; uint8_t direction;} mp256_t;
typedef struct { size_t n, m; mp256_t *a; } mp256_v;

mp256_v *get_shimmer_hits(py_mmer_t * py_mmer, uint64_t mhash0, uint32_t span) {

	khash_t(MMER0) * mmer0_map_ = (khash_t(MMER0) *) py_mmer->mmer0_map; 
	//khash_t(RLEN) * rlmap_ = (khash_t(RLEN) *) rlmap_;
	//khash_t(MMC) * mcmap_ = (khash_t(MMC) *) mcmap_; 

	mp128_v * mpv;
	mp256_t mp256;
	mp256_v * mpv_out;
	uint64_t mhash1;
	khiter_t k;

	mpv_out = calloc(sizeof(mp256_v), 1);

	khash_t(MMER1) * mmer1_map;
    mhash0 <<= 8;
	mhash0 |= span;
	mp256.x0 = mhash0;

	k = kh_get(MMER0, mmer0_map_, mhash0);
	if ( k == kh_end(mmer0_map_) ) {
		return mpv_out;
	}
	mmer1_map = kh_val(mmer0_map_, k);
	for (khiter_t __j = kh_begin(mmer1_map); __j != kh_end(mmer1_map); ++__j) {

		if (!kh_exist(mmer1_map,__j)) continue;
		mhash1 = kh_key(mmer1_map, __j);
		mp256.x1 = mhash1;
		mhash1 >>= 8;
		mpv = kh_val(mmer1_map, __j);
		qsort(mpv->a, mpv->n, sizeof(mp128_t), mp128_comp);
		for (size_t __k0 = 0; __k0 < (mpv->n); __k0++) {  
			mp256.y0 = mpv->a[__k0].y0;
			mp256.y1 = mpv->a[__k0].y1;
			mp256.direction = mpv->a[__k0].direction;

			kv_push(mp256_t, NULL, *mpv_out, mp256); 
		}
	}
	return mpv_out;
}
