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

#define LOWERBOUND 2
#define UPPERBOUND 30
#define REVERSED 1

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
		seq[rp] = rmap[seq[p]];
		seq[p] = rmap[tmp];
		p ++;
	}
}

int rp128_comp(const void * a, const void * b) {
	rp128_t * a0 = (rp128_t *) a;
	rp128_t * b0 = (rp128_t *) b;
	return ((a0->y0 & 0xFFFFFFFF) >> 1) < ((b0->y0 & 0xFFFFFFFF) >> 1);
}

void main() {
	char mmc_file_path[] = "../test/test/hmmer-L2-MC-01-of-01.dat";
	char mmer_file_path[] = "../test/test/hmmer-L2-01-of-01.dat";
	char seq_idx_file_path[] = "../test/test/seq_dataset.idx";
	char seq_db_file_path[] = "../test/test/seq_dataset.seqdb";

	FILE * seqdb;
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
	seqdb = fopen(seq_db_file_path, "rb");

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
			rp.direction = 0;
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
			rpos = kh_val(rlmap, k).len - pos + span -1;
			rp.y0 = ((rp.y0 & 0xFFFFFFFF00000001) | (rpos << 1 )) ^ 0x1; // ^0x1 -> switch strand

			rp.y1 = mmer0.y;
			span = mmer0.x & 0xFF;
		    rid = rp.y1 >> 32;
            pos = ((rp.y1 & 0xFFFFFFFF) >> 1) + 1;
			k = kh_get(RLEN, rlmap, rid);
			assert(k != kh_end(rlmap));
			rpos = kh_val(rlmap, k).len - pos + span - 1;
			rp.y1 = ((rp.y1 & 0xFFFFFFFF00000001) | (rpos << 1 )) ^ 0x1; // ^0x1 -> switch strand
            rp.direction = REVERSED;
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
	uint32_t rlen; 
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
			if (rpv->n <= 1) continue;
			qsort(rpv->a, rpv->n, sizeof(rp128_t), rp128_comp);

			char *seq;
			char *pseq = NULL;
			uint32_t ppos = 0;
			uint32_t prid = 0;
			uint32_t pstrand = 0;
			uint32_t plen = 0;
			for (size_t __k=0; __k < rpv->n; __k++) {
				y0 = rpv->a[__k].y0;
				y1 = rpv->a[__k].y1;
				rid0 = (uint32_t) (y0 >> 32);
				rid1 = (uint32_t) (y1 >> 32);
				assert( rid0 == rid1 );
			    strand0 = (uint32_t) (y0 & 0x1);
			    strand1 = (uint32_t) (y1 & 0x1);
				pos0 = (uint32_t) ((y0 & 0xFFFFFFFF) >> 1) + 1;
				pos1 = (uint32_t) ((y1 & 0xFFFFFFFF) >> 1) + 1;
			    
				k = kh_get(RLEN, rlmap, rid0);
			    rlen = kh_val(rlmap, k).len;	
				assert(k != kh_end(rlmap));
				rpos0 = rlen - pos0 + span0;
				rpos1 = rlen - pos1 + span1;

				k = kh_get(MMC, mcmap, mhash0);
				assert(k != kh_end(mcmap));
				mcount0 = kh_val(mcmap, k);

				k = kh_get(MMC, mcmap, mhash1);
				assert(k != kh_end(mcmap));
				mcount1 = kh_val(mcmap, k);

	   		    alignment * aln;
				seq = get_read_seq(seqdb, rid0, rlmap);
				strand = rpv->a[__k].direction;
				if (strand == REVERSED) {
					reverse_complement(seq, rlen);
				}
				//printf("%014lX %u %014lX %u %u %09u %u %u %09u %u %u %u %u\n", 
				//		mhash0, strand0, mhash1, strand1, rpv->a[__k].direction,
				//		rid0, pos0, rpos0, rid1, pos1, rpos1, mcount0, mcount1);
				if (pseq == NULL) {
					pseq = seq;
					ppos = pos0;
					prid = rid0;
					pstrand = strand;
					plen = rlen;
				} else {
					//printf("seq1:%s\nseq2:%s\n\n", pseq+ppos-pos0, seq);
					aln = align(pseq+ppos-pos0, strlen(pseq+ppos-pos0), seq, strlen(seq), 100);
					printf("%d %d %d %d %d %d\n", 
							aln->aln_str_size, aln->dist,
							aln->aln_q_s, aln->aln_q_e, 
							aln->aln_t_s, aln->aln_t_e);

					uint32_t a_bgn, a_end, b_bgn, b_end;
					if (pstrand == 0) {
						a_bgn = ppos-pos0;
						a_end = ppos-pos0 + aln->aln_q_e - aln->aln_q_s;
					} else {
						a_bgn = plen - (ppos - pos0) - (aln->aln_q_e - aln->aln_q_s);
						a_end = plen - (ppos - pos0);
					}
					if (strand == 0) {
						b_bgn = aln->aln_t_s;
						b_end = aln->aln_t_e;
					} else {
						b_bgn = rlen - aln->aln_t_e;
						b_end = rlen - aln->aln_t_s;;
					}
					if ( aln->aln_t_e > 500 && aln->aln_q_e > 500) {
						printf("%09d %09d %d %f %u %u %u %u %u %u %u %u %s\n", 
								prid, rid0, 0, 99.0,
								pstrand, a_bgn, a_end, plen,
								strand, b_bgn, b_end, rlen, "overlap");	   
					}
				    kfree(NULL, pseq);
					pseq = seq;
					ppos = pos0;
					prid = rid0;
					pstrand = rpv->a[__k].direction;
					plen = rlen;
				    free_alignment(aln);
				}
			}
			kfree(NULL, pseq);
		}
	}
	// TODO: clean up memory
}

