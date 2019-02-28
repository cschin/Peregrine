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

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

#define handle_error(msg) \
	do { perror(msg); exit(EXIT_FAILURE); } while (0)

#define LOWERBOUND 2
#define UPPERBOUND 30
#define REVERSED 1
#define READENDFUZZINESS 48

KHASH_SET_INIT_INT64(RPAIR);

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

void build_map(
		mm128_v * mmers,
		khash_t(MMER0) * mmer0_map, 
		khash_t(RLEN) *rlmap, 
		khash_t(MMC) *mcmap) {

    uint64_t mhash;
    mm128_t mmer0, mmer1;
	rp128_v * rpv;
	uint32_t rid;
    uint32_t pos, rpos;
	uint32_t span;
	uint32_t mcount = 0;
	int32_t absent;
	rp128_t rp;
	khiter_t k;
	khash_t(MMER1) * mmer1_map;
	size_t s=0;

	while (1) {
		mmer0 = mmers->a[s];
		mhash = mmer0.x >> 8;
		k = kh_get(MMC, mcmap, mhash);
		assert(k != kh_end(mcmap));
        mcount = kh_val(mcmap, k); 
		if (mcount >= LOWERBOUND && mcount < UPPERBOUND) break; 
		s++;
	}
	for( size_t i=s+1; i < mmers->n; i++ ){
		mmer1 = mmers->a[i];
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
				rpv->n = 0; rpv->m = 0; rpv->a = NULL;
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
				rpv->n = 0; rpv->m = 0; rpv->a = NULL;
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

}


void pring_overlaps(
		rp128_v * rpv,
		khash_t(RLEN) * rlmap,
		khash_t(RPAIR) * rid_pairs,
		char * seq_p){
	uint64_t y0;
	uint64_t ridp;
	uint32_t rid0;
	uint32_t pos0;
	uint32_t rlen; 
	uint32_t strand;
	khiter_t k;
	char *seq;
	char *pseq = NULL;
	uint32_t ppos = 0;
	uint32_t prid = 0;
	uint32_t pstrand = 0;
	uint32_t plen = 0;
	int32_t absent;
	bool contains;
	for (size_t __k=0; __k < rpv->n; __k++) {
		y0 = rpv->a[__k].y0;
		rid0 = (uint32_t) (y0 >> 32);
		pos0 = (uint32_t) ((y0 & 0xFFFFFFFF) >> 1) + 1;

		k = kh_get(RLEN, rlmap, rid0);
		rlen = kh_val(rlmap, k).len;	
		assert(k != kh_end(rlmap));
		alignment * aln;
		seq = get_read_seq_mmap(seq_p, rid0, rlmap);
		strand = rpv->a[__k].direction;

		if (strand == REVERSED) {
			reverse_complement(seq, rlen);
		}
		if (pseq == NULL) {
			pseq = seq;
			ppos = pos0;
			prid = rid0;
			pstrand = strand;
			plen = rlen;
		} else {
			assert(ppos-pos0 >= 0 );

			ridp = prid < rid0 ? ((uint64_t) prid) << 32 | ((uint64_t) rid0) : ((uint64_t) rid0) << 32 | ((uint64_t) prid);
			k = kh_put(RPAIR, rid_pairs, ridp, &absent);

			if (absent) {
				uint32_t pslen = plen-ppos+pos0;
				uint32_t slen = rlen;
				aln = align(pseq+ppos-pos0, pslen, seq, slen, 100);

				seq_coor_t q_bgn, q_end, t_bgn, t_end;
				q_bgn = aln->aln_q_s; q_end = aln->aln_q_e; t_bgn = aln->aln_t_s; t_end = aln->aln_t_e;


				if ((q_bgn < READENDFUZZINESS &&
							t_bgn < READENDFUZZINESS &&
							(abs(pslen - q_end) < READENDFUZZINESS ||
							 abs(slen - t_end) < READENDFUZZINESS)) &&  
						q_end > 500 && t_end > 500) {
				    //printf("%d %d %d %d\n", q_bgn, q_end, t_bgn, t_end);
					//printf("%s\n%s\n", pseq+ppos-pos0, seq); 
					seq_coor_t a_bgn, a_end, b_bgn, b_end;

					double err_est;
					err_est	= 100.0 - 100.0 * (double) (aln->dist) / (double) (aln->aln_str_size);

					q_bgn -= t_bgn;
					t_bgn = 0;

					char type[128];
					if ( slen < pslen ) {
						t_end = slen; 
						q_end = slen + (q_end - t_end);
						strcpy(type ,"contains");
						contains = true;
					} else {
						t_end = pslen  - (q_end - t_end); 
						q_end = pslen;
						strcpy(type ,"overlap");
						contains = false;
					}
					if (pstrand == 0) {
						a_bgn = (seq_coor_t) (ppos-pos0) + q_bgn;
						a_end = (seq_coor_t) (ppos-pos0) + q_end;
					} else {
						q_bgn -= t_bgn;
						t_bgn = 0;
						a_bgn = (seq_coor_t) plen - (seq_coor_t) (ppos - pos0) - q_end;
						a_end = (seq_coor_t) plen - (seq_coor_t) (ppos - pos0) - q_bgn;
					}
					if (strand == 0) {
						b_bgn = t_bgn;
						b_end = t_end;
					} else {
						b_bgn = (seq_coor_t) rlen - t_end;
						b_end = (seq_coor_t) rlen - t_bgn;
					}
					printf("%09d %09d %d %0.1f %u %d %d %u %u %d %d %u %s\n", 
							prid, rid0, -aln->aln_str_size, err_est,
							0, a_bgn, a_end, plen,
							pstrand == 0 ? strand : 1-strand, b_bgn, b_end, rlen, type);

				}
				free_alignment(aln);
			}
			if (contains == false) {
				kfree(NULL, pseq);
				pseq = seq;
				ppos = pos0;
				prid = rid0;
				pstrand = strand;
				plen = rlen;
			}
		}
	}
	kfree(NULL, pseq);
}

void process_overlaps(
		khash_t(MMER0) * mmer0_map, 
		khash_t(RLEN) *rlmap, 
		khash_t(MMC) *mcmap) {

	int fd;
	struct stat sb;
	char * seq_p;
	rp128_v * rpv;
	uint64_t mhash0, mhash1;

	char seq_db_file_path[] = "../test/test/seq_dataset.seqdb";
	khash_t(MMER1) * mmer1_map;

	fd = open(seq_db_file_path, O_RDONLY);
	if (fd == -1)
		handle_error("open");
	
	if (fstat(fd, &sb) == -1)           /* To obtain file size */
		handle_error("fstat");

	seq_p = mmap((caddr_t)0, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);

	khash_t(RPAIR) * rid_pairs = kh_init(RPAIR);
	for (khiter_t __i = kh_begin(mmer0_map); __i != kh_end(mmer0_map); ++__i) {
		if (!kh_exist(mmer0_map,__i)) continue;
		mhash0 = kh_key(mmer0_map, __i);
	    mhash0 >>= 8;	
		mmer1_map = kh_val(mmer0_map, __i);
		for (khiter_t __j = kh_begin(mmer1_map); __j != kh_end(mmer1_map); ++__j) {
		    if (!kh_exist(mmer1_map,__j)) continue;
			mhash1 = kh_key(mmer1_map, __j);
			mhash1 >>= 8;
			rpv = kh_val(mmer1_map, __j);
			if (rpv->n <= 1) continue;
			qsort(rpv->a, rpv->n, sizeof(rp128_t), rp128_comp);
			pring_overlaps(rpv, rlmap, rid_pairs, seq_p);	   
		}
	}
	kh_destroy(RPAIR, rid_pairs);
	munmap(seq_p, sb.st_size);
}

int main() {
	char mmc_file_path[] = "../test/test/hmmer-L2-MC-01-of-01.dat";
	char mmer_file_path[] = "../test/test/hmmer-L2-01-of-01.dat";
	char seq_idx_file_path[] = "../test/test/seq_dataset.idx";

	mm128_v mmers;
	mm_count_v mmc;

	khash_t(RLEN) *rlmap; 
	khash_t(MMC) *mcmap = kh_init(MMC);
    
	khash_t(MMER0) * mmer0_map = kh_init(MMER0);

	mmers = read_mmlist(mmer_file_path);
	mmc = read_mm_count(mmc_file_path);
	aggregate_mm_count(mcmap, &mmc);
	rlmap = get_read_length_map(seq_idx_file_path);

    build_map( &mmers, mmer0_map, rlmap, mcmap);
	process_overlaps(mmer0_map, rlmap, mcmap);
	// TODO: clean up memory
}

