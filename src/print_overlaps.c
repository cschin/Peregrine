#include <errno.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
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
#define BESTN 1

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
		khash_t(MMC) *mcmap,
		uint32_t chunk,
		uint32_t total_chunk) {

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

	for(;;) {
		mmer0 = mmers->a[s];
		mhash = mmer0.x >> 8;
		k = kh_get(MMC, mcmap, mhash);
		assert(k != kh_end(mcmap));
        mcount = kh_val(mcmap, k); 
		if (mcount >= MMER_COUNT_LOWER_BOUND && mcount < MMER_COUNT_UPPER_BOUND) break;
		s++;
	}
	for( size_t i=s+1; i < mmers->n; i++ ){
		mmer1 = mmers->a[i];
		mhash = mmer1.x >> 8;
		k = kh_get(MMC, mcmap, mhash);
		assert(k != kh_end(mcmap));
        mcount = kh_val(mcmap, k); 
		if (mcount < MMER_COUNT_LOWER_BOUND || mcount > MMER_COUNT_UPPER_BOUND) continue;

	    if ((mmer0.y >> 32) == (mmer1.y >> 32)) {  // the pairs are in the same read

			// don't use two minimers that are too close to each other
			//if (((mmer1.y >> 1) & 0xFFFFFFF) - ((mmer0.y >> 1) & 0xFFFFFFF) < 100) {
			//	mmer0 = mmer1;
			//	continue;
			//}

			if ((mmer0.x >> 8) % total_chunk == chunk % total_chunk) {
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

				rp.direction = ORIGINAL;
				kv_push(rp128_t, NULL, *rpv, rp);
				// printf("Y %lu %lu %09u\n", mmer0.x >> 8, mmer1.x >> 8, (uint32_t) (rp.y0 >> 32));
			}
			
			// reverse direction
			if ((mmer1.x >> 8) % total_chunk == chunk % total_chunk) {
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
				// printf("Y %lu %lu %09u\n", mmer1.x >> 8, mmer0.x >> 8, rid);
			}
		}	
		mmer0 = mmer1;
	}

}


void print_overlaps(
		rp128_v * rpv,
		khash_t(RLEN) * rlmap,
		khash_t(RPAIR) * rid_pairs,
		char * seq_p){

	uint64_t ridp;
	uint64_t y0, y1;
	uint32_t rid0, pos0, rlen0, strand0;
	uint32_t rid1, pos1, rlen1, strand1;
	khiter_t k;
	char *seq0 = NULL;
	char *seq1 = NULL;
	int32_t absent;
	for (size_t __k0=0; __k0 < rpv->n; __k0++) {
		y0 = rpv->a[__k0].y0;
		rid0 = (uint32_t) (y0 >> 32);
		pos0 = (uint32_t) ((y0 & 0xFFFFFFFF) >> 1) + 1;

		k = kh_get(RLEN, rlmap, rid0);
		rlen0 = kh_val(rlmap, k).len;	
		assert(k != kh_end(rlmap));
		seq0 = get_read_seq_mmap(seq_p, rid0, rlmap);
		strand0 = rpv->a[__k0].direction;

		if (strand0 == REVERSED) {
			reverse_complement(seq0, rlen0);
		}

		size_t overlap_count = 0;
		for (size_t __k1=1; (__k0+__k1 < rpv->n) && (overlap_count <= BESTN); __k1++ ) {
			y0 = rpv->a[__k0+__k1].y0;
			rid1 = (uint32_t) (y0 >> 32);
			
			ridp = rid0 < rid1 ? ( ((uint64_t) rid0) << 32) | ((uint64_t) rid1) : (((uint64_t) rid1) << 32) | ((uint64_t) rid0);
			k = kh_get(RPAIR, rid_pairs, ridp);
			if (rid0 == rid1) continue;
			if (k != kh_end(rid_pairs) || rid0 == rid1) {
				overlap_count += 1;
				continue;
			}

			pos1 = (uint32_t) ((y0 & 0xFFFFFFFF) >> 1) + 1;
			assert(pos0-pos1 >= 0 );
			k = kh_get(RLEN, rlmap, rid1);
			rlen1 = kh_val(rlmap, k).len;	
			assert(k != kh_end(rlmap));
			seq1 = get_read_seq_mmap(seq_p, rid1, rlmap);
			strand1 = rpv->a[__k0+__k1].direction;

			//printf("X1: %lu %lu %lu %lu %09u %09u\n", rpv->n, __k0, __k1, overlap_count, rid0, rid1);
			if (strand1 == REVERSED) {
				reverse_complement(seq1, rlen1);
			}

			//printf("%09d %s\n%09d %s\n",rid0, seq0+pos0-pos1,rid1, seq1); 
			alignment * aln;
			uint32_t slen0 = rlen0 - pos0 + pos1;
			uint32_t slen1 = rlen1;
			aln = align(seq0 + pos0 - pos1, slen0, seq1, slen1, 100);
			seq_coor_t q_bgn, q_end, t_bgn, t_end;
			q_bgn = aln->aln_q_s; q_end = aln->aln_q_e; 
			t_bgn = aln->aln_t_s; t_end = aln->aln_t_e;
			//printf("X2: %u %u %d %d %d %d\n", pos0, pos1, q_bgn, q_end, t_bgn, t_end);

			if ((q_bgn < READ_END_FUZZINESS &&
						t_bgn < READ_END_FUZZINESS &&
						(abs(slen0 - q_end) < READ_END_FUZZINESS ||
						 abs(slen1 - t_end) < READ_END_FUZZINESS)) &&
					q_end > 500 && t_end > 500) {
				//printf("X2: %u %u %d %d %d %d\n", pos0, pos1, q_bgn, q_end, t_bgn, t_end);
				// printf("%s\n%s\n", seq0+pos0-pos1, seq1); 


				seq_coor_t a_bgn, a_end, b_bgn, b_end;
				double err_est;
				err_est	= 100.0 - 100.0 * (double) (aln->dist) / (double) (aln->aln_str_size);

				q_bgn -= t_bgn;
				t_bgn = 0;

				char ovlp_type[128];
				if ( slen1 < slen0 ) {
					t_end = slen1; 
					q_end = slen1 + (q_end - t_end);
					strcpy(ovlp_type ,"contains");
				} else {
					t_end = slen0  - (q_end - t_end); 
					q_end = slen0;
					strcpy(ovlp_type ,"overlap");
					overlap_count ++;
				}
				if (strand0 == ORIGINAL) {
					a_bgn = (seq_coor_t) (pos0-pos1) + q_bgn;
					a_end = (seq_coor_t) (pos0-pos1) + q_end;
				} else {
					q_bgn -= t_bgn;
					t_bgn = 0;
					a_bgn = (seq_coor_t) rlen0 - (seq_coor_t) (pos0 - pos1) - q_end;
					a_end = (seq_coor_t) rlen0 - (seq_coor_t) (pos0 - pos1) - q_bgn;
				}
				if (strand1 == ORIGINAL) {
					b_bgn = t_bgn;
					b_end = t_end;
				} else {
					b_bgn = (seq_coor_t) rlen1 - t_end;
					b_end = (seq_coor_t) rlen1 - t_bgn;
				}
				kh_put(RPAIR, rid_pairs, ridp, &absent);
				//assert(absent == 1);
				printf("%09d %09d %d %0.1f %u %d %d %u %u %d %d %u %s\n", 
						rid0, rid1, -aln->aln_str_size, err_est,
						ORIGINAL, a_bgn, a_end, rlen0,
						strand0 == ORIGINAL ? strand1 : 1-strand1, b_bgn, b_end, rlen1, ovlp_type);
			}
			free_alignment(aln);
			kfree(NULL, seq1);
		}
		kfree(NULL, seq0);
	}
}

void process_overlaps(char * seq_db_file_path,
		khash_t(MMER0) * mmer0_map, 
		khash_t(RLEN) *rlmap, 
		khash_t(MMC) *mcmap) {

	int fd;
	struct stat sb;
	char * seq_p;
	rp128_v * rpv;
	uint64_t mhash0, mhash1;

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
			if (rpv->n <= 1 || rpv->n > LOCAL_OVERLAP_UPPERBOUND) continue;
			qsort(rpv->a, rpv->n, sizeof(rp128_t), rp128_comp);
		    //printf("X %lu %lu %lu\n", mhash0, mhash1, rpv->n);
			print_overlaps(rpv, rlmap, rid_pairs, seq_p);	   
		}
	}
	kh_destroy(RPAIR, rid_pairs);
	munmap(seq_p, sb.st_size);
}

int main(int argc, char *argv[]) {
	char * seqdb_prefix = NULL;
	char * shimmer_prefix = NULL;

	char mmc_file_path[8192];
	char mmer_file_path[8192]; 
	char seq_idx_file_path[8192];
	char seqdb_file_path[8291];
        int c;	
	uint32_t total_chunk, mychunk;

	mm128_v mmers = {0, 0, 0};
	mm128_v mmers_;
	mm_count_v mmc;

	khash_t(RLEN) *rlmap; 
	khash_t(MMC) *mcmap = kh_init(MMC);
    
	khash_t(MMER0) * mmer0_map;
	khash_t(MMER1) * mmer1_map;

	rp128_v * rpv;
	
	opterr = 0;

	while ((c = getopt(argc, argv, "p:l:t:c:")) != -1) {
		switch (c) {
			case 'p':
				seqdb_prefix = optarg;
				break;
			case 'l':
				shimmer_prefix = optarg;
				break;
			case 't':
				total_chunk = atoi(optarg);
				break;
			case 'c':
				mychunk = atoi(optarg);
				break;
			case '?':
				if (optopt == 'p') {
					fprintf (stderr, "Option -%c not specified, using 'seq_dataset' as the sequence db prefix\n", optopt);
				}
				if (optopt == 'l') {
					fprintf(stderr, "Option -%c not specified, using 'shimmer-L2' as the L2 index prefix\n", optopt);
				}
				return 1;
			default:
				abort();
		}
	}

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

	rlmap = get_read_length_map(seq_idx_file_path);
	
	written = snprintf(seqdb_file_path, sizeof(seqdb_file_path), "%s.seqdb", seqdb_prefix);
	assert(written < sizeof(seqdb_file_path));
	fprintf(stderr, "using seqdb file: %s\n", seqdb_file_path);
	
        for (int c=1; c<=total_chunk; c++) {	
		int written;
		written = snprintf(mmer_file_path, sizeof(mmer_file_path), "%s-%02d-of-%02d.dat", shimmer_prefix, c, total_chunk);
		assert(written < sizeof(mmer_file_path));
		fprintf(stderr, "useing data file: %s\n", mmer_file_path);
		
		written = snprintf(mmc_file_path, sizeof(mmc_file_path), "%s-MC-%02d-of-%02d.dat", shimmer_prefix, c, total_chunk);
		assert(written < sizeof(mmc_file_path));
		fprintf(stderr, "output data file: %s\n", mmc_file_path);

		mmers_ = read_mmlist(mmer_file_path);
		append_mmlist(&mmers, &mmers_);
		
		mmc = read_mm_count(mmc_file_path);
		aggregate_mm_count(mcmap, &mmc);
		kv_destroy(mmc);
		kv_destroy(mmers_);
	}

	mmer0_map = kh_init(MMER0);
	build_map(&mmers, mmer0_map, rlmap, mcmap, mychunk, total_chunk);
	process_overlaps(seqdb_file_path, mmer0_map, rlmap, mcmap);
	
	for (khiter_t __i = kh_begin(mmer0_map); __i != kh_end(mmer0_map); ++__i) {
		if (!kh_exist(mmer0_map,__i)) continue;
		mmer1_map = kh_val(mmer0_map, __i);
		for (khiter_t __j = kh_begin(mmer1_map); __j != kh_end(mmer1_map); ++__j) {
			if (!kh_exist(mmer1_map,__j)) continue;
			rpv = kh_val(mmer1_map, __j);
			kv_destroy(*rpv);
		}
		kh_destroy(MMER1, mmer1_map);
	}
	
	kh_destroy(MMER0, mmer0_map);
	kh_destroy(MMC, mcmap);
	kh_destroy(RLEN, rlmap);
	// TODO: clean up memory
}

