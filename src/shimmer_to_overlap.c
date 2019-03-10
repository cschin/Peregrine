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

extern char *optarg;
extern int optind, opterr, optopt;

#define handle_error(msg) \
	do { perror(msg); exit(EXIT_FAILURE); } while (0)

#define MMER_COUNT_LOWER_BOUND 2
#define MMER_COUNT_UPPER_BOUND 240
#ifndef ORIGINAL
#define ORIGINAL 0
#endif
#ifndef REVERSED
#define REVERSED 1
#endif
#define READ_END_FUZZINESS 48
#define LOCAL_OVERLAP_UPPERBOUND 120
#define BESTN 4
#define OVERLAP 0
#define CONTAINMENT 1
#define ALNBANDSIZE 100


KHASH_MAP_INIT_INT64(RPAIR, uint8_t);


int mp128_comp(const void * a, const void * b) {
	mp128_t * a0 = (mp128_t *) a;
	mp128_t * b0 = (mp128_t *) b;
	return ((a0->y0 & 0xFFFFFFFF) >> 1) < ((b0->y0 & 0xFFFFFFFF) >> 1);
}


void shimmer_to_overlap(
		mp128_v * mpv,
		khash_t(RLEN) * rlmap,
		khash_t(RPAIR) * rid_pairs,
		char * seq_p){

	uint64_t ridp;
	uint64_t y0;
	uint32_t rid0, pos0, rlen0, strand0;
	uint32_t rid1, pos1, rlen1, strand1;
	uint32_t right_ext = 0;
	khiter_t k;
	char *seq0 = NULL;
	char *seq1 = NULL;
	int32_t absent;
	uint8_t * contained;

	contained = (uint8_t *) calloc(mpv->n, sizeof(uint8_t)); // use calloc to set element to zero

	// clock_t time_begin = clock();
	// clock_t time_end;
	for (size_t __k0 = (mpv->n)-1; __k0 > 0; __k0--) {  // note: k0 is an unsigned type
		y0 = mpv->a[__k0-1].y0;
		rid0 = (uint32_t) (y0 >> 32);
		pos0 = (uint32_t) ((y0 & 0xFFFFFFFF) >> 1) + 1;
		k = kh_get(RLEN, rlmap, rid0);
		assert(k != kh_end(rlmap));
		rlen0 = kh_val(rlmap, k).len;	
		strand0 = mpv->a[__k0-1].direction;
		seq0 = get_read_seq_mmap(seq_p, rid0, rlmap, strand0);

		if (right_ext == 0) {
			right_ext = rlen0;
		}

		size_t overlap_count = 0;
		for (size_t __k1=1; (__k0+__k1-1 < mpv->n) && (overlap_count < BESTN); __k1++ ) {

			if ( contained[__k0+__k1-1] == 1 ) continue;

			y0 = mpv->a[__k0+__k1-1].y0;
			rid1 = (uint32_t) (y0 >> 32);
			
			if (rid0 == rid1) continue;
			// time_end = clock();
			// printf("X0: %lu %lu %lu %lu %09u %09u %lu\n", mpv->n, __k0, __k1, overlap_count, rid0, rid1, time_end-time_begin);
			// time_begin = time_end;
			ridp = rid0 < rid1 ? ( ((uint64_t) rid0) << 32) | ((uint64_t) rid1) : (((uint64_t) rid1) << 32) | ((uint64_t) rid0);
			k = kh_get(RPAIR, rid_pairs, ridp);
			if (k != kh_end(rid_pairs)) {
				if ( kh_val(rid_pairs, k) == OVERLAP ) overlap_count += 1;
				continue;
			}

			pos1 = (uint32_t) ((y0 & 0xFFFFFFFF) >> 1) + 1;
			assert(pos0-pos1 >= 0 );
			k = kh_get(RLEN, rlmap, rid1);
			assert(k != kh_end(rlmap));
			rlen1 = kh_val(rlmap, k).len;
			strand1 = mpv->a[__k0+__k1-1].direction;
			seq1 = get_read_seq_mmap(seq_p, rid1, rlmap, strand1);

			//printf("X1: %lu %lu %lu %lu %09u %09u\n", mpv->n, __k0, __k1, overlap_count, rid0, rid1);

			//printf("%09d %s\n%09d %s\n",rid0, seq0+pos0-pos1,rid1, seq1); 
			uint32_t slen0 = rlen0 - pos0 + pos1;
			uint32_t slen1 = rlen1;
			alignment_t * aln;
			aln = align(seq0 + pos0 - pos1, slen0, seq1, slen1, ALNBANDSIZE);
			seq_coor_t q_bgn, q_end, t_bgn, t_end;
			q_bgn = aln->q_bgn; q_end = aln->q_end; 
			t_bgn = aln->t_bgn; t_end = aln->t_end;
			//printf("X2: %u %u %d %d %d %d\n", pos0, pos1, q_bgn, q_end, t_bgn, t_end);

			if ((q_bgn < READ_END_FUZZINESS &&
						t_bgn < READ_END_FUZZINESS &&
						(abs(slen0 - q_end) < READ_END_FUZZINESS ||
						 abs(slen1 - t_end) < READ_END_FUZZINESS)) &&
					q_end > 500 && t_end > 500) {
				// printf("X3: %u %u %d %d %d %d\n", pos0, pos1, q_bgn, q_end, t_bgn, t_end);
				// printf("%s\n%s\n", seq0+pos0-pos1, seq1); 


				seq_coor_t a_bgn, a_end, b_bgn, b_end;
				double err_est;
				err_est	= 100.0 - 100.0 * (double) (aln->dist) / (double) (aln->astr_size);

				q_bgn -= t_bgn;
				t_bgn = 0;

				char ovlp_type[128];
				if ( slen1 < slen0 + READ_END_FUZZINESS) {
					t_end = slen1; 
					q_end = slen1 + (q_end - t_end);
					strcpy(ovlp_type ,"contains");
					k = kh_put(RPAIR, rid_pairs, ridp, &absent);
					kh_val(rid_pairs, k) = CONTAINMENT;
					contained[__k0+__k1-1] = 1;
				} else {
					t_end = slen0  - (q_end - t_end); 
					q_end = slen0;
					strcpy(ovlp_type ,"overlap");
					overlap_count ++;
					k = kh_put(RPAIR, rid_pairs, ridp, &absent);
					kh_val(rid_pairs, k) = OVERLAP;
				}
				if (strand0 == ORIGINAL) {
					a_bgn = (seq_coor_t) (pos0-pos1) + q_bgn;
					a_end = (seq_coor_t) (pos0-pos1) + q_end;
					a_bgn = a_bgn < 0 ? 0 : a_bgn;              //this ad-hoc fix, read should be stiched by alignment
					a_end = a_end >= rlen0 ? rlen0 : a_end;
				} else {
					q_bgn -= t_bgn;
					t_bgn = 0;
					a_bgn = (seq_coor_t) rlen0 - (seq_coor_t) (pos0 - pos1) - q_end;
					a_end = (seq_coor_t) rlen0 - (seq_coor_t) (pos0 - pos1) - q_bgn;
					a_bgn = a_bgn < 0 ? 0 : a_bgn;              //this ad-hoc fix
					a_end = a_end >= rlen0 ? rlen0 : a_end;
				}
				if (strand1 == ORIGINAL) {
					b_bgn = t_bgn;
					b_end = t_end;
					b_bgn = b_bgn < 0 ? 0 : b_bgn;              //this ad-hoc fix
					b_end = b_end >= rlen1 ? rlen1 : b_end;
				} else {
					b_bgn = (seq_coor_t) rlen1 - t_end;
					b_end = (seq_coor_t) rlen1 - t_bgn;
					b_bgn = b_bgn < 0 ? 0 : b_bgn;              //this ad-hoc fix
					b_end = b_end >= rlen1 ? rlen1 : b_end;
				}
				assert(absent == 1);
				fprintf(stdout,"%09d %09d %d %0.1f %u %d %d %u %u %d %d %u %s\n",
						rid0, rid1, -(aln->astr_size), err_est,
						ORIGINAL, a_bgn, a_end, rlen0,
						(strand0 == ORIGINAL ? strand1 : 1-strand1), b_bgn, b_end, rlen1,
						ovlp_type);
			}
			free_alignment(aln);
			kfree(NULL, seq1);
		}
		kfree(NULL, seq0);
	}
	free(contained);
}

void process_overlaps(char * seq_db_file_path,
		khash_t(MMER0) * mmer0_map, 
		khash_t(RLEN) *rlmap, 
		khash_t(MMC) *mcmap) {

	int fd;
	struct stat sb;
	char * seq_p;
	mp128_v * mpv;
	uint64_t mhash0, mhash1;

	khash_t(MMER1) * mmer1_map;

	fd = open(seq_db_file_path, O_RDONLY);
	if (fd == -1)
		handle_error("open");
	
	if (fstat(fd, &sb) == -1)           /* To obtain file size */
		handle_error("fstat");

	seq_p = mmap((caddr_t)0, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);

	khash_t(RPAIR) * rid_pairs = kh_init(RPAIR);
	uint32_t iter_count = 0;
	clock_t begin = clock();
	clock_t end;
	for (khiter_t __i = kh_begin(mmer0_map); __i != kh_end(mmer0_map); ++__i) {
		if (!kh_exist(mmer0_map,__i)) continue;
		mhash0 = kh_key(mmer0_map, __i);
	    mhash0 >>= 8;	
		mmer1_map = kh_val(mmer0_map, __i);
		for (khiter_t __j = kh_begin(mmer1_map); __j != kh_end(mmer1_map); ++__j) {
		    if (!kh_exist(mmer1_map,__j)) continue;
			mhash1 = kh_key(mmer1_map, __j);
			mhash1 >>= 8;
			mpv = kh_val(mmer1_map, __j);
			if (mpv->n <= 2 || mpv->n > LOCAL_OVERLAP_UPPERBOUND) continue;
			qsort(mpv->a, mpv->n, sizeof(mp128_t), mp128_comp);
			shimmer_to_overlap(mpv, rlmap, rid_pairs, seq_p);
			iter_count++;
			if (iter_count % 10000 == 0) {
				end = clock();
				fprintf(stderr, "iter_count: %d, %f s for 10000 candidates\n", iter_count, (double)(end - begin) / CLOCKS_PER_SEC);
				begin = end;
			}
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
	uint32_t total_chunk = 1, mychunk = 1;

	wordexp_t p; 
	char **mmc_fns; 
	char **shimmer_fns;

	mm128_v mmers = {0, 0, 0};
	mm128_v mmers_;
	mm_count_v mmc;

	khash_t(RLEN) *rlmap; 
	khash_t(MMC) *mcmap = kh_init(MMC);
    
	khash_t(MMER0) * mmer0_map;
	khash_t(MMER1) * mmer1_map;

	mp128_v * mpv;
	
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
	

	written = snprintf(mmer_file_path, sizeof(mmer_file_path), "%s-[0-9]*-of-[0-9]*.dat", shimmer_prefix);
	assert(written < sizeof(mmer_file_path));
	wordexp(mmer_file_path, &p, 0);
	shimmer_fns = p.we_wordv;
	for (int i = 0; i < p.we_wordc; i++) {
		fprintf(stderr, "useing shimmer data file: %s\n", shimmer_fns[i]);
		mmers_ = read_mmlist(shimmer_fns[i]);
		append_mmlist(&mmers, &mmers_);
		kv_destroy(mmers_);
	}
	wordfree(&p);	

	char buffer[32768];

	setvbuf(stdout, buffer, _IOFBF, sizeof(buffer));


	written = snprintf(mmc_file_path, sizeof(mmc_file_path), "%s-MC-[0-9]*-of-[0-9]*.dat", shimmer_prefix);
	assert(written < sizeof(mmc_file_path));
	wordexp(mmc_file_path, &p, 0);
	mmc_fns = p.we_wordv;
	for (int i = 0; i < p.we_wordc; i++) {
		fprintf(stderr, "using shimmer count file: %s\n", mmc_fns[i]);
		mmc = read_mm_count(mmc_fns[i]);
		aggregate_mm_count(mcmap, &mmc);
		kv_destroy(mmc);
	}

	wordfree(&p);	
	


	mmer0_map = kh_init(MMER0);
	build_map(&mmers, mmer0_map, rlmap, mcmap, mychunk, total_chunk, MMER_COUNT_LOWER_BOUND, MMER_COUNT_UPPER_BOUND);
	process_overlaps(seqdb_file_path, mmer0_map, rlmap, mcmap);
	
	for (khiter_t __i = kh_begin(mmer0_map); __i != kh_end(mmer0_map); ++__i) {
		if (!kh_exist(mmer0_map,__i)) continue;
		mmer1_map = kh_val(mmer0_map, __i);
		for (khiter_t __j = kh_begin(mmer1_map); __j != kh_end(mmer1_map); ++__j) {
			if (!kh_exist(mmer1_map,__j)) continue;
			mpv = kh_val(mmer1_map, __j);
			kv_destroy(*mpv);
		}
		kh_destroy(MMER1, mmer1_map);
	}
	
	kh_destroy(MMER0, mmer0_map);
	kh_destroy(MMC, mcmap);
	kh_destroy(RLEN, rlmap);
	// TODO: clean up memory
}

