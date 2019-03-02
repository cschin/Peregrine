#include <assert.h>
#include <errno.h>
#include <zlib.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include "shimmer.h"
#include "kseq.h"
#include "kvec.h"
#include "khash.h"
#include "kseq.h"

#define REDUCTION_FACTOR 6

KSEQ_INIT(gzFile, gzread);

extern char *optarg;
extern int optind, opterr, optopt;

void write_mc_count_mm128(char *fn, mm128_v * shimmer) {
	mm_count_v mmc = {0,0,0};
	khash_t(MMC) *mcmap = kh_init(MMC);
	mm_count(shimmer, mcmap, &mmc);
	write_mm_count(fn, &mmc);
	kv_destroy(mmc);
	kh_destroy(MMC, mcmap);
}

int main(int argc, char *argv[])
{
	gzFile fp;
	FILE *seq_dataset_file;
	kseq_t *seq;
	uint32_t rid;
	char *seq_dataset_path = NULL;
	char *out_prefix = NULL;
	char *index_path = NULL;
	char shimmer_output_path[8192];
	char seq_file_path[8192];
	int l, c, written;
	int total_chunk = 1;
	int mychunk = 1;
	int reduction_factor = REDUCTION_FACTOR;
	int number_layers = 2;
	int seq_file_counter;
	int output_L0_mc = 1;
	mm128_v shimmerL0 = {0,0,0};
	mm128_v shimmerL1 = {0,0,0};
	mm128_v shimmerL2 = {0,0,0};
	seq_data_v seq_data = {0, 0, 0};

    khash_t(RID) *hmap;
    khash_t(RLEN) *rlmap=kh_init(RLEN);
	opterr = 0;

	while ((c = getopt(argc, argv, "d:i:o:t:c:l:r:m:")) != -1) {
		switch (c) {
			case 'd':
				seq_dataset_path = optarg;
				break;
			case 'i':
				index_path = optarg;
				break;
			case 'o':
				out_prefix = optarg;
				break;
			case 't':
				total_chunk = atoi(optarg);
				break;
			case 'c':
				mychunk = atoi(optarg);
				break;
			case 'r':
				reduction_factor = atoi(optarg);
				break;
			case 'l':
				number_layers = atoi(optarg);
				break;
			case 'm':
				output_L0_mc = atoi(optarg);
				break;
			case '?':
				if (optopt == 'd') {
					fprintf (stderr, "Option -%c not specified, using 'seq_dataset.lst' as the input sequence dataset file list\n", optopt);
				}
				if (optopt == 'i') {
					fprintf(stderr, "Option -%c not specified, using 'seq_dataset.idx' as the input index path\n", optopt);
                }
				else if (optopt == 'o') {
					fprintf (stderr, "Option -%c not specified, using 'shimmer' as the output prefix\n", optopt);
				}
				return 1;
			default:
				abort();
		}
	}

	assert(total_chunk > 0);
	assert(mychunk > 0 && mychunk <= total_chunk);
	assert(reduction_factor < 256);
	
	fprintf(stderr, "reduction factor= %d\n", reduction_factor);

	if (seq_dataset_path == NULL) {
		seq_dataset_path = (char *) calloc(8192, 1);
		snprintf( seq_dataset_path, 8191, "seq_dataset.lst" );
	}

	if (index_path == NULL) {
		index_path = (char *) calloc(8192, 1);
		snprintf( index_path, 8191, "seq_dataset.idx" );
	}

	if (out_prefix == NULL) {
		out_prefix = (char *) calloc(8192, 1);
		snprintf( out_prefix, 8191, "shimmer" );
	}


	seq_dataset_file = fopen(seq_dataset_path, "r");
	fprintf(stderr, "input sequqnece dataset file: '%s'\n", seq_dataset_path);
	if (!seq_dataset_file) {
		fprintf(stderr, "file '%s' open error: %s\n", seq_dataset_path, strerror(errno));
		exit(1);
	}

	fprintf(stderr, "using sequqnece index file: '%s'\n", index_path);
    hmap = build_read_index(index_path, &seq_data, rlmap);
	
	seq_file_counter = 0;
	while( fscanf(seq_dataset_file, "%s", seq_file_path) != EOF ) {
		seq_file_counter++;
		if ( (seq_file_counter % total_chunk) != (mychunk % total_chunk)) continue;
		fprintf(stderr, "prcess file #%d: %s in chunk %d of %d\n", seq_file_counter, seq_file_path, mychunk, total_chunk);
		fp = gzopen(seq_file_path, "r");
		if (!fp) {
			fprintf(stderr, "file '%s' open error: %s\n", seq_file_path, strerror(errno));
			exit(1);
		}
		seq = kseq_init(fp);

		while ((l = kseq_read(seq)) >= 0) {
			int is_missing;
			khiter_t k;
			if (seq->seq.l < 2000) continue;
			k = kh_get(RID, hmap, seq->name.s);
			is_missing = (k == kh_end(hmap));
			rid = kh_value(hmap, k);
			if (is_missing) {
			    fprintf(stderr, "WARNNING: cannot map read:%s to an internal unique identifier. The read will be ignored. \n", seq->name.s);
				continue;
			};
			mm_sketch(NULL, seq->seq.s, seq->seq.l, 80, 16, rid, 0, &shimmerL0);
		}
		gzclose(fp);
	}

	written = snprintf(shimmer_output_path, sizeof shimmer_output_path, "%s-L0-%02d-of-%02d.dat", out_prefix, mychunk, total_chunk);
	assert(written < sizeof(shimmer_output_path));
	fprintf(stderr, "output data file: %s\n", shimmer_output_path);
	write_mmlist(shimmer_output_path, &shimmerL0);

    /* temporary disable this as it is not used for now	
	mm128_v shimmerE5 = {0,0,0};
	mm128_v shimmerE3 = {0,0,0};

    mm_end_filter(&shimmerL0, &shimmerE5, &shimmerE3, rlmap, 250);
	written = snprintf(shimmer_output_path, sizeof shimmer_output_path, "%s-E5-%02d-of-%02d.dat", out_prefix, mychunk, total_chunk);
	assert(written < sizeof(shimmer_output_path));
	printf("output data file: %s\n", shimmer_output_path);
	write_mmlist(shimmer_output_path, &shimmerE5);
	kv_destroy(shimmerE5);

	written = snprintf(shimmer_output_path, sizeof shimmer_output_path, "%s-E3-%02d-of-%02d.dat", out_prefix, mychunk, total_chunk);
	assert(written < sizeof(shimmer_output_path));
	printf("output data file: %s\n", shimmer_output_path);
	write_mmlist(shimmer_output_path, &shimmerE3);
	kv_destroy(shimmerE3);
    */

	if (output_L0_mc == 1) {
		written = snprintf(shimmer_output_path, sizeof shimmer_output_path, "%s-L0-MC-%02d-of-%02d.dat", out_prefix, mychunk, total_chunk);
		assert(written < sizeof(shimmer_output_path));
		printf("output data file: %s\n", shimmer_output_path);
		write_mc_count_mm128(shimmer_output_path, &shimmerL0);
	}

    mm_reduce(&shimmerL0, &shimmerL1, reduction_factor);
	kv_destroy(shimmerL0);
	if (number_layers == 1) {
		written = snprintf(shimmer_output_path, sizeof shimmer_output_path, "%s-L1-%02d-of-%02d.dat", out_prefix, mychunk, total_chunk);
		assert(written < sizeof(shimmer_output_path));
		printf("output data file: %s\n", shimmer_output_path);
		write_mmlist(shimmer_output_path, &shimmerL1);
		
		written = snprintf(shimmer_output_path, sizeof shimmer_output_path, "%s-L1-MC-%02d-of-%02d.dat", out_prefix, mychunk, total_chunk);
		assert(written < sizeof(shimmer_output_path));
		printf("output data file: %s\n", shimmer_output_path);
		write_mc_count_mm128(shimmer_output_path, &shimmerL1);
	}
	else if (number_layers > 1) {
		mm_reduce(&shimmerL1, &shimmerL2, reduction_factor);
	    kv_destroy(shimmerL1);
		written = snprintf(shimmer_output_path, sizeof shimmer_output_path, "%s-L2-%02d-of-%02d.dat", out_prefix, mychunk, total_chunk);
		assert(written < sizeof(shimmer_output_path));
		fprintf(stderr, "output data file: %s\n", shimmer_output_path);
		write_mmlist(shimmer_output_path, &shimmerL2);
		
		written = snprintf(shimmer_output_path, sizeof shimmer_output_path, "%s-L2-MC-%02d-of-%02d.dat", out_prefix, mychunk, total_chunk);
		assert(written < sizeof(shimmer_output_path));
		printf("output data file: %s\n", shimmer_output_path);
		write_mc_count_mm128(shimmer_output_path, &shimmerL2);

	    kv_destroy(shimmerL2);
	}

	kh_destroy(RID, hmap);
	kh_destroy(RLEN, rlmap);
	for (size_t _i=0; _i < seq_data.n; _i++) {
		kfree(NULL, seq_data.a[_i].name);
	}
	kv_destroy(seq_data);

	fclose(seq_dataset_file);
	if (!seq_dataset_path) free(seq_dataset_path);
	if (!out_prefix) free(out_prefix);
	return 0;
}
