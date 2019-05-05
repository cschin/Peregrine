#include <errno.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <wordexp.h>
#include "shimmer.h"
#include <time.h>
#include "khash.h"
#include "kvec.h"
#include "kalloc.h"

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>


shmr_aln_v * shmr_aln( 
        mm128_v *mmers0, 
        mm128_v *mmers1, 
        uint8_t direction, 
        double maxdiff) {
    /* generate a list of co-aligned mimimizer from two 
     * minimizer lists: mv1 and mv2
     */
    uint64_t mhash;
    mm128_t mmer0, mmer1;
    khash_t(MMIDX) *mmidx_map = kh_init(MMIDX);
    shmr_aln_v *alns;
	khiter_t k;
	int32_t absent;

    mm_idx_v * idx_tmp;

    alns = calloc(sizeof(shmr_aln_v), 1);

	mm_idx_t s=0;
    /* build a hasmap from mhash to the index of the minimizer array */
	for(;;) {
		if (s >= mmers0->n) break;
		mmer0 = mmers0->a[s];
		mhash = mmer0.x >> 8;

        k = kh_put(MMIDX, mmidx_map, mhash, &absent);
        if (absent) {
            idx_tmp = calloc(sizeof(mm_idx_v), 1);
            kv_push(mm_idx_t, 0, *idx_tmp, s);
            kh_val(mmidx_map, k) = idx_tmp;
        } else {
            k = kh_get(MMIDX, mmidx_map, mhash);
            assert( k != kh_end(mmidx_map));
            idx_tmp = kh_val(mmidx_map, k);
            kv_push(mm_idx_t, 0, *idx_tmp, s);
        }
		s++;
	}

    /* loop through 2nd shimmer list to build alginements */
	mm_idx_t ss=0;
	for(;;) {
		if (ss >= mmers1->n) break;
        if (direction == 1) {  // reversed
            s =  mmers1->n - ss; 
        } else {
            s = ss;
        }
		mmer1 = mmers1->a[s];
		mhash = mmer1.x >> 8;
        k = kh_get(MMIDX, mmidx_map, mhash);
        if (k == kh_end(mmidx_map)) {
            ss++;
            continue;
        }
        idx_tmp = kh_val(mmidx_map, k);
        for (uint32_t i=0; i < idx_tmp -> n; i++) {
            mmer0 = mmers0->a[idx_tmp->a[i]];
            int64_t delta0, delta1;
            int64_t mm_dist;
            uint32_t grouped = 0;
            if (direction == 1) {
                delta0 = abs(mmer_pos(&mmer0) + mmer_pos(&mmer1));
            } else {
                delta0 = abs(mmer_pos(&mmer0) - mmer_pos(&mmer1));
            }
            for (uint32_t aln_idx = 0; aln_idx < alns->n; aln_idx ++ ){
                mm128_t m0, m1;
                shmr_aln_t * aln; 
                aln = alns->a + aln_idx;
                m0 = aln->m0.a[aln->m0.n-1];
                m1 = aln->m1.a[aln->m1.n-1];

                if (direction == 1) {
                    delta1 = abs(mmer_pos(&m0) + mmer_pos(&m1));
                } else {
                    delta1 = abs(mmer_pos(&m0) - mmer_pos(&m1));
                }
                mm_dist = abs(mmer_pos(&mmer0) - mmer_pos(&m0));
                // should we group the new miminiizer pair to the min-diff one?
                if ( (double) abs(delta0 - delta1) / (double) (mm_dist) < maxdiff ) {
                    kv_push(mm128_t, 0, aln->m0, mmer0);
                    kv_push(mm128_t, 0, aln->m1, mmer1);
                    kv_push(mm_idx_t, 0, aln->idx0, idx_tmp->a[i]);
                    kv_push(mm_idx_t, 0, aln->idx1, s); 
                    grouped = 1;

                    break;
                }
            }
            if (grouped == 0) {
                shmr_aln_t *aln;
                aln = calloc(sizeof(shmr_aln_t),1);
                kv_push(mm128_t, 0, aln->m0, mmer0);
                kv_push(mm128_t, 0, aln->m1, mmer1);
                kv_push(mm_idx_t, 0, aln->idx0, idx_tmp->a[i]);
                kv_push(mm_idx_t, 0, aln->idx1, s); 
                kv_push(shmr_aln_t, 0, *alns, *aln);
            }
        }
        ss++;
    }

	for (khiter_t __i = kh_begin(mmidx_map); __i != kh_end(mmidx_map); ++__i) {
		if (!kh_exist(mmidx_map,__i)) continue;
		idx_tmp = kh_val(mmidx_map, __i);
        kv_destroy(*idx_tmp);
	}
    kh_destroy(MMIDX, mmidx_map);
    return alns;
}

void free_shmr_alns(shmr_aln_v * alns) {
    for (uint32_t aln_idx = 0; aln_idx < alns->n; aln_idx ++ ){
        for (uint32_t idx =0; idx < alns[aln_idx].n; idx ++) {
            kv_destroy(alns[aln_idx].a[idx].m0);
            kv_destroy(alns[aln_idx].a[idx].m1);
            kv_destroy(alns[aln_idx].a[idx].idx0);
            kv_destroy(alns[aln_idx].a[idx].idx1);
        }
            kv_destroy(alns[aln_idx]);
    }
}
