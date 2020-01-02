#include <assert.h>
#include <errno.h>
#include <zlib.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "shimmer.h"
#include "kseq.h"
#include "kvec.h"
#include "khash.h"

#define handle_error(msg) \
	do { perror(msg); exit(EXIT_FAILURE); } while (0)

#define REDUCTION_FACTOR 6
#define DEFAULT_WINDOW_SIZE 80
#define DEFAULT_KMER_SIZE 16

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
	char * seqdb_prefix = NULL;
	char * shimmer_prefix = NULL;

	char seq_idx_file_path[8192];
	char seqdb_file_path[8291];
	char shimmer_output_path[8192];
	
	FILE * seq_index_file;
	uint32_t rid;

	int c, written;
	int total_chunk = 1;
	int mychunk = 1;
	int reduction_factor = REDUCTION_FACTOR;
	int number_layers = 2;
	int output_L0 = 1;
	int window_size = DEFAULT_WINDOW_SIZE;
	int kmer_size = DEFAULT_KMER_SIZE;
	mm128_v shimmerL0 = {0,0,0};
	mm128_v shimmerL1 = {0,0,0};
	mm128_v shimmerL2 = {0,0,0};
	seq_data_v seq_data = {0, 0, 0};

    khash_t(RLEN) *rlmap=kh_init(RLEN);
	opterr = 0;

	while ((c = getopt(argc, argv, "p:o:t:c:l:r:m:w:k:")) != -1) {
		switch (c) {
			case 'p':
				seqdb_prefix = optarg;
				break;
			case 'o':
				shimmer_prefix = optarg;
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
				output_L0 = atoi(optarg);
				break;
			case 'w':
				window_size = atoi(optarg);
				break;
			case 'k':
				kmer_size = atoi(optarg);
				break;
			case '?':
				if (optopt == 'p') {
					fprintf (stderr, "Option -%c not specified, using 'seq_dataset' as the sequence db prefix\n", optopt);
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
	assert(window_size >= 24 && kmer_size >= 12 && window_size > kmer_size);
	
	fprintf(stderr, "reduction factor= %d\n", reduction_factor);

	if (seqdb_prefix == NULL) {
		seqdb_prefix = (char *) calloc(8192, 1);
		snprintf( seqdb_prefix, 8191, "seq_dataset" );
	}

	if (shimmer_prefix == NULL) {
		shimmer_prefix = (char *) calloc(8192, 1);
		snprintf( shimmer_prefix, 8191, "shimmer" );
	}

	written = snprintf(seq_idx_file_path, sizeof(seq_idx_file_path), "%s.idx", seqdb_prefix);
	assert(written < sizeof(seq_idx_file_path));
	fprintf(stderr, "using index file: %s\n", seq_idx_file_path);

	rlmap = get_read_length_map(seq_idx_file_path);
	
	written = snprintf(seqdb_file_path, sizeof(seqdb_file_path), "%s.seqdb", seqdb_prefix);
	assert(written < sizeof(seqdb_file_path));
	fprintf(stderr, "using seqdb file: %s\n", seqdb_file_path);
	
    int fd;	
	struct stat sb;
	uint8_t * seq_p;
	fd = open(seqdb_file_path, O_RDONLY);
	if (fd == -1)
		handle_error("open");
	
	if (fstat(fd, &sb) == -1)           /* To obtain file size */
		handle_error("fstat");

	seq_p = (uint8_t *)  mmap((void *)0, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);

	seq_index_file = fopen(seq_idx_file_path, "r");
	char name_buf[512];
	uint32_t rlen;
	size_t offset;
	while (fscanf(seq_index_file, "%u %255s %u %lu", &rid, name_buf, &rlen, &offset) != EOF) {
		if ( (rid % total_chunk) != (mychunk % total_chunk)) continue;
		char * seq = malloc(rlen+1); 
		decode_biseq(seq_p + offset, seq, rlen, 0);
		seq[rlen] = '\0'; 
		mm_sketch(NULL, seq, rlen, window_size, kmer_size, rid, 0, &shimmerL0);
		free(seq);
	}

	if (output_L0 == 1) {
        written = snprintf(shimmer_output_path, sizeof shimmer_output_path, "%s-L0-%02d-of-%02d.dat", shimmer_prefix, mychunk, total_chunk);
        assert(written < sizeof(shimmer_output_path));
        fprintf(stderr, "output data file: %s\n", shimmer_output_path);
        write_mmlist(shimmer_output_path, &shimmerL0);

        /* temporary disable this as it is not used for now	
        mm128_v shimmerE5 = {0,0,0};
        mm128_v shimmerE3 = {0,0,0};

        mm_end_filter(&shimmerL0, &shimmerE5, &shimmerE3, rlmap, 250);
        written = snprintf(shimmer_output_path, sizeof shimmer_output_path, "%s-E5-%02d-of-%02d.dat", shimmer_prefix, mychunk, total_chunk);
        assert(written < sizeof(shimmer_output_path));
        printf("output data file: %s\n", shimmer_output_path);
        write_mmlist(shimmer_output_path, &shimmerE5);
        kv_destroy(shimmerE5);

        written = snprintf(shimmer_output_path, sizeof shimmer_output_path, "%s-E3-%02d-of-%02d.dat", shimmer_prefix, mychunk, total_chunk);
        assert(written < sizeof(shimmer_output_path));
        printf("output data file: %s\n", shimmer_output_path);
        write_mmlist(shimmer_output_path, &shimmerE3);
        kv_destroy(shimmerE3);
    */

		written = snprintf(shimmer_output_path, sizeof shimmer_output_path, "%s-L0-MC-%02d-of-%02d.dat", shimmer_prefix, mychunk, total_chunk);
		assert(written < sizeof(shimmer_output_path));
		printf("output data file: %s\n", shimmer_output_path);
		write_mc_count_mm128(shimmer_output_path, &shimmerL0);
	}

    mm_reduce(&shimmerL0, &shimmerL1, reduction_factor);
	kv_destroy(shimmerL0);
	if (number_layers == 1) {
		written = snprintf(shimmer_output_path, sizeof shimmer_output_path, "%s-L1-%02d-of-%02d.dat", shimmer_prefix, mychunk, total_chunk);
		assert(written < sizeof(shimmer_output_path));
		printf("output data file: %s\n", shimmer_output_path);
		write_mmlist(shimmer_output_path, &shimmerL1);
		
		written = snprintf(shimmer_output_path, sizeof shimmer_output_path, "%s-L1-MC-%02d-of-%02d.dat", shimmer_prefix, mychunk, total_chunk);
		assert(written < sizeof(shimmer_output_path));
		printf("output data file: %s\n", shimmer_output_path);
		write_mc_count_mm128(shimmer_output_path, &shimmerL1);
	}
	else if (number_layers > 1) {
		mm_reduce(&shimmerL1, &shimmerL2, reduction_factor);
	    kv_destroy(shimmerL1);
		written = snprintf(shimmer_output_path, sizeof shimmer_output_path, "%s-L2-%02d-of-%02d.dat", shimmer_prefix, mychunk, total_chunk);
		assert(written < sizeof(shimmer_output_path));
		fprintf(stderr, "output data file: %s\n", shimmer_output_path);
		write_mmlist(shimmer_output_path, &shimmerL2);
		
		written = snprintf(shimmer_output_path, sizeof shimmer_output_path, "%s-L2-MC-%02d-of-%02d.dat", shimmer_prefix, mychunk, total_chunk);
		assert(written < sizeof(shimmer_output_path));
		printf("output data file: %s\n", shimmer_output_path);
		write_mc_count_mm128(shimmer_output_path, &shimmerL2);

	    kv_destroy(shimmerL2);
	}

	kh_destroy(RLEN, rlmap);
	for (size_t _i=0; _i < seq_data.n; _i++) {
		kfree(NULL, seq_data.a[_i].name);
	}
	kv_destroy(seq_data);

	munmap(seq_p, sb.st_size);
	if (!seqdb_prefix) free(seqdb_prefix);
	if (!shimmer_prefix) free(shimmer_prefix);
	return 0;
}
