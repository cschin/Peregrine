#include <errno.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include "shimmer.h"
#include "khash.h"
#include "kvec.h"
#include "kalloc.h"
#include "ksort.h"

#define LOWERBOUND 2
#define UPPERBOUND 30

#define sort_key_128x(a) ((a).x)
KRADIX_SORT_INIT(128x, mm128_t, sort_key_128x, 8);

#define sort_key_256x(a) ( ( (((uint64_t) (a).x0) << 16) & 0xFFFFFFFF00000000) |  ( (uint64_t) ((a).x1 >> 24 )))
KRADIX_SORT_INIT(256x, mm256_t, sort_key_256x, 8);

void main() {
	char mmc_file_path[] = "../test/test/hmmer-L2-MC-01-of-01.dat";
	char mmer_file_path[] = "../test/test/hmmer-L2-01-of-01.dat";
	char seq_idx_file_path[] = "../test/seq_dataset.idx";

	mm128_v mmers;
	mm_count_v mmc;
	mm256_v mmpairs = {0, 0, 0};
	mm256_t mmpair;

    uint64_t mhash;
    mm128_t mmer0, mmer1;
	uint32_t rid;
	uint32_t strand;
    uint32_t pos, rpos;
	uint32_t span;
	uint32_t mcount = 0;
	int32_t absent;
	rp128_t rp;
	rp128_v *rpv;
	khiter_t k;
	khash_t(RLEN) *rlmap; 
	khash_t(MMC) *mcmap = kh_init(MMC);
    
	khash_t(MMER1) * mmer1_map;
	khash_t(MMER0) * mmer0_map = kh_init(MMER0);

	mmers = read_mmlist(mmer_file_path);
	mmc = read_mm_count(mmc_file_path);
	aggregate_mm_count(mcmap, &mmc);
	rlmap = get_read_length_map(seq_idx_file_path);

	size_t s=0;
	while (1) {
		mmer0 = mmers.a[s];
		mhash = mmer0.x >> 8;
		k = kh_get(MMC, mcmap, mhash);
		assert(k != kh_end(mcmap));
        mcount = kh_val(mcmap, k); 
		if (mcount >= LOWERBOUND && mcount < UPPERBOUND) break; 
	}
	for( size_t i=s+1; i < mmers.n; i++ ){
		mmer1 = mmers.a[i];
		mhash = mmer1.x >> 8;
		k = kh_get(MMC, mcmap, mhash);
		assert(k != kh_end(mcmap));
        mcount = kh_val(mcmap, k); 
		if (mcount < LOWERBOUND ||  mcount > UPPERBOUND) continue; 

	    if ( (mmer0.y >> 32) == (mmer1.y >> 32) ) {  // the pairs are in the same read
			k = kh_put(MMER0, mmer0_map, mmer0.x, &absent);
			if (absent) kh_value(mmer0_map, k) = kh_init(MMER1);

			mmer1_map = kh_value(mmer0_map, k);
			k = kh_put(MMER1, mmer1_map, mmer1.x, &absent);
			if (absent) {
			    rpv = kmalloc(NULL, sizeof(rp128_v));
				rpv->n = 0;
				rpv->m = 0;
				rpv->a = NULL;
				kh_value(mmer1_map, k) = rpv;
			} else {
				rpv = kh_value(mmer1_map, k);
			}
			rp.y0 = mmer0.y;
			rp.y1 = mmer1.y;
			kv_push(rp128_t, NULL, *rpv, rp);
			
			// reverse direction
			k = kh_put(MMER0, mmer0_map, mmer1.x, &absent);
			if (absent) kh_value(mmer0_map, k) = kh_init(MMER1);

			mmer1_map = kh_value(mmer0_map, k);
			k = kh_put(MMER1, mmer1_map, mmer0.x, &absent);
			if (absent) {
			    rpv = kmalloc(NULL, sizeof(rp128_v));
				rpv->n = 0;
				rpv->m = 0;
				rpv->a = NULL;
				kh_value(mmer1_map, k) = rpv;
			} else {
				rpv = kh_value(mmer1_map, k);
			}
			rp.y0 = mmer1.y;
			span = mmer1.x & 0xFF;
		    rid = rp.y0 >> 32;
            pos = ((rp.y0 & 0xFFFFFFFF) >> 1) + 1;
			k = kh_get(RLEN, rlmap, rid);
			assert(k != kh_end(rlmap));
			rpos = kh_val(rlmap, k) - pos + span -1;
			rp.y0 = ((rp.y0 & 0xFFFFFFFF00000001) | (rpos << 1 )) ^ 0x1; // ^0x1 -> switch strand

			rp.y1 = mmer0.y;
			span = mmer0.x & 0xFF;
		    rid = rp.y1 >> 32;
            pos = ((rp.y1 & 0xFFFFFFFF) >> 1) + 1;
			k = kh_get(RLEN, rlmap, rid);
			assert(k != kh_end(rlmap));
			rpos = kh_val(rlmap, k) - pos + span - 1;
			rp.y1 = ((rp.y1 & 0xFFFFFFFF00000001) | (rpos << 1 )) ^ 0x1; // ^0x1 -> switch strand

			kv_push(rp128_t, NULL, *rpv, rp);
		}	
		mmer0 = mmer1;
	}

	khint_t __i;
	khint_t __j;

	uint64_t mhash0, mhash1;
	uint32_t rid0, rid1;
	uint32_t strand0, strand1;
	uint32_t pos0, pos1;
	uint32_t rpos0, rpos1;
	uint32_t span0, span1;
	uint32_t mcount0, mcount1;
	uint64_t y0, y1;
	for (__i = kh_begin(mmer0_map); __i != kh_end(mmer0_map); ++__i) {
		if (!kh_exist(mmer0_map,__i)) continue;
		mhash0 = kh_key(mmer0_map, __i);
		span0 = mhash0 & 0xFF;
	    mhash0 >>= 8;	
		mmer1_map = kh_val(mmer0_map, __i);
		for (__j = kh_begin(mmer1_map); __j != kh_end(mmer1_map); ++__j) {
		    if (!kh_exist(mmer1_map,__j)) continue;
			mhash1 = kh_key(mmer1_map, __j);
			span1 = mhash1 & 0xFF;
			mhash1 >>= 8;
			rpv = kh_val(mmer1_map, __j);
			for (size_t __k=0; __k < rpv->n; __k++) {
				y0 = rpv->a[__k].y0;
				y1 = rpv->a[__k].y1;
				rid0 = (uint32_t) (y0 >> 32);
				rid1 = (uint32_t) (y1 >> 32);
			    strand0 = (uint32_t) (y0 & 0x1);
			    strand1 = (uint32_t) (y1 & 0x1);
				pos0 = (uint32_t) ((y0 & 0xFFFFFFFF) >> 1) + 1;
				pos1 = (uint32_t) ((y1 & 0xFFFFFFFF) >> 1) + 1;
				
				k = kh_get(RLEN, rlmap, rid0);
				assert(k != kh_end(rlmap));
				rpos0 = kh_val(rlmap, k) - pos0 + span0;

				k = kh_get(RLEN, rlmap, rid1);
				assert(k != kh_end(rlmap));
				rpos1 = kh_val(rlmap, k) - pos1 + span1;

				k = kh_get(MMC, mcmap, mhash0);
				assert(k != kh_end(mcmap));
				mcount0 = kh_val(mcmap, k);

				k = kh_get(MMC, mcmap, mhash1);
				assert(k != kh_end(mcmap));
				mcount1 = kh_val(mcmap, k);

				printf("%lu %u %lu %u %u %u %u %u %u %u %u %u\n", mhash0, strand0, mhash1, strand1,
						rid0, pos0, rpos0, rid1, pos1, rpos1, mcount0, mcount1);
			}
		}
	}
	// TODO: clean up memory
}

