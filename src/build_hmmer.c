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


KSEQ_INIT(gzFile, gzread);
KHASH_MAP_INIT_STR(RID, uint32_t);

extern char *optarg;
extern int optind, opterr, optopt;

khash_t(RID) * build_read_index(char *fpath, seq_data_v *seq_data, khash_t(RLEN) *rlmap ) {
	khash_t(RID) *hmap = kh_init(RID);
	int absent;
	khiter_t k;
	FILE *index_file;
	index_file = fopen(fpath, "r");
	uint32_t i, rid, rlen;
	char name_buf[256];
	char * name;
	while (fscanf(index_file, "%u %255s %u", &rid, name_buf, &rlen) != EOF) {
		kv_resize(seq_data_t, NULL, *seq_data, seq_data->n + 8192);
		name = kmalloc(NULL, 256 * sizeof(char));
		strncpy(name, name_buf, 256);
		i = seq_data->n;
		seq_data->a[i].name = name;
		seq_data->a[i].rid = rid;
		seq_data->n ++;
		k=kh_put(RID, hmap, seq_data->a[i].name, &absent);
		if (absent) kh_value(hmap, k) = rid;
		k=kh_put(RLEN, rlmap, rid, &absent);
		kh_value(rlmap, k) = rlen;
		/* 
        { 
			khint_t __i;		
			for (__i = kh_begin(hmap); __i != kh_end(hmap); ++__i) {		
				if (!kh_exist(hmap,__i)) continue;						
				printf("testx:%s %u\n", kh_key(hmap,__i), kh_val(hmap,__i));
			} 
		}
	    */	
	}
	fclose(index_file);
	return hmap;
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
	char hmmer_output_path[8192];
	char seq_file_path[8192];
	int l, c, written;
	int total_chunk = 1;
	int mychunk = 1;
	int reduction_factor = 8;
	int number_layers = 2;
	int seq_file_counter;
	mm128_v hmmerL0 = {0,0,0};
	mm128_v hmmerL1 = {0,0,0};
	mm128_v hmmerL2 = {0,0,0};
	seq_data_v seq_data = {0, 0, 0};

    khash_t(RID) *hmap;
    khash_t(RLEN) *rlmap=kh_init(RLEN);
	opterr = 0;

	while ((c = getopt(argc, argv, "i:o:t:c:l:r:")) != -1) {
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
			case '?':
				if (optopt == 'd') {
					fprintf (stderr, "Option -%c not specified, using 'seq_dataset.lst' as the input sequence dataset file list\n", optopt);
				}
				if (optopt == 'i') {
					fprintf(stderr, "Option -%c not specified, using 'seq_dataset.idx' as the input index path\n", optopt);
                }
				else if (optopt == 'o') {
					fprintf (stderr, "Option -%c not specified, using 'hmmer' as the output prefix\n", optopt);
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
		snprintf( out_prefix, 8191, "hmmer" );
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
			    fprintf(stderr, "WARNNING: cannot map read:%s to an internal unique identifier. The read will be ignore. The read will be ignored. \n", seq->name.s);
				continue;
			};
			mm_sketch(NULL, seq->seq.s, seq->seq.l, 80, 16, rid, 0, &hmmerL0);
		}
		gzclose(fp);
	}

	written = snprintf(hmmer_output_path, sizeof hmmer_output_path, "%s-L0-%02d-of-%02d.dat", out_prefix, mychunk, total_chunk);
	assert(written < sizeof(hmmer_output_path));
	fprintf(stderr, "output data file: %s\n", hmmer_output_path);
	write_mmlist(hmmer_output_path, &hmmerL0);

	
	mm128_v hmmerE5 = {0,0,0};
	mm128_v hmmerE3 = {0,0,0};

    mm_end_filter(&hmmerL0, &hmmerE5, &hmmerE3, rlmap, 250);
	written = snprintf(hmmer_output_path, sizeof hmmer_output_path, "%s-E5-%02d-of-%02d.dat", out_prefix, mychunk, total_chunk);
	assert(written < sizeof(hmmer_output_path));
	printf("output data file: %s\n", hmmer_output_path);
	write_mmlist(hmmer_output_path, &hmmerE5);
	kv_destroy(hmmerE5);

	written = snprintf(hmmer_output_path, sizeof hmmer_output_path, "%s-E3-%02d-of-%02d.dat", out_prefix, mychunk, total_chunk);
	assert(written < sizeof(hmmer_output_path));
	printf("output data file: %s\n", hmmer_output_path);
	write_mmlist(hmmer_output_path, &hmmerE3);
	kv_destroy(hmmerE3);

	mm_count_v mmc = {0,0,0};
	mm_count(&hmmerL0, &mmc);
	written = snprintf(hmmer_output_path, sizeof hmmer_output_path, "%s-MC-%02d-of-%02d.dat", out_prefix, mychunk, total_chunk);
	assert(written < sizeof(hmmer_output_path));
	printf("output data file: %s\n", hmmer_output_path);
	write_mm_count(hmmer_output_path, &mmc);
	kv_destroy(mmc);

    mm_reduce(&hmmerL0, &hmmerL1, reduction_factor);
	kv_destroy(hmmerL0);
	if (number_layers == 1) {
		written = snprintf(hmmer_output_path, sizeof hmmer_output_path, "%s-L1-%02d-of-%02d.dat", out_prefix, mychunk, total_chunk);
		assert(written < sizeof(hmmer_output_path));
		printf("output data file: %s\n", hmmer_output_path);
		write_mmlist(hmmer_output_path, &hmmerL1);
	}
	else if (number_layers > 1) {
		mm_reduce(&hmmerL1, &hmmerL2, reduction_factor);
	    kv_destroy(hmmerL1);
		written = snprintf(hmmer_output_path, sizeof hmmer_output_path, "%s-L2-%02d-of-%02d.dat", out_prefix, mychunk, total_chunk);
		assert(written < sizeof(hmmer_output_path));
		fprintf(stderr, "output data file: %s\n", hmmer_output_path);
		write_mmlist(hmmer_output_path, &hmmerL2);
	    kv_destroy(hmmerL2);
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
