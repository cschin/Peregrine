#ifndef SHIMMER_H
#define SHIMMER_H

#include <assert.h>
#include <stdint.h>
#include <unistd.h>
#include "khash.h"

#ifdef __cplusplus
extern "C" {
#endif

void encode_4bit_bidirection(uint8_t *, char *, size_t ); 

void decode_4bit_bidirection(uint8_t *, char *, size_t, uint8_t);

void reverse_complement(char *, size_t);

typedef struct { uint64_t x, y; } mm128_t;
typedef struct { size_t n, m; mm128_t *a; } mm128_v;
void mm_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p);
void mm_reduce(mm128_v *, mm128_v *,  uint8_t);
uint32_t ks_ksmall_uint32_t(size_t n, uint32_t arr[], size_t kk);

void write_mmlist(char *, mm128_v *); 
mm128_v read_mmlist(char *); 
void append_mmlist(mm128_v *, mm128_v *);

typedef struct { char * name; uint32_t rid;} seq_data_t; 
typedef struct { size_t n, m; seq_data_t * a; } seq_data_v;

typedef struct { uint32_t len; size_t offset; } rl_t;
KHASH_MAP_INIT_STR(RID, uint32_t);
KHASH_MAP_INIT_INT(RLEN, rl_t);
khash_t(RID) * build_read_index(char *, seq_data_v *, khash_t(RLEN) * );
khash_t(RLEN) * get_read_length_map(char *);

void mm_end_filter(mm128_v *, mm128_v *, mm128_v *, khash_t(RLEN) *, uint32_t); 

KHASH_MAP_INIT_INT64(MMC, uint32_t);
typedef struct { uint64_t mer; uint32_t count; } mm_count_t;
typedef struct { size_t n, m; mm_count_t *a; } mm_count_v;
void mm_count(mm128_v *, khash_t(MMC) *, mm_count_v *); 
void write_mm_count(char *, mm_count_v *); 
void mm_count_to_vec(khash_t(MMC) *, mm_count_v *);
mm_count_v read_mm_count(char *fn);

void aggregate_mm_count(khash_t(MMC) *,  mm_count_v *); 

typedef struct { uint64_t x0, x1, y0, y1; } mm256_t;
typedef struct { size_t n, m; mm256_t *a; } mm256_v;

typedef struct { uint64_t y0, y1; uint8_t direction;} mp128_t;
typedef struct { size_t n, m; mp128_t *a; } mp128_v;
KHASH_MAP_INIT_INT64(MMER1, mp128_v *);
typedef khash_t(MMER1) * mmert1_p_t;
KHASH_MAP_INIT_INT64(MMER0, mmert1_p_t);

void build_map( mm128_v *, khash_t(MMER0) *, khash_t(RLEN) *, khash_t(MMC) *, 
		uint32_t, uint32_t, uint32_t, uint32_t); 

char * get_read_seq(FILE *, uint32_t, khash_t(RLEN) *);
char * get_read_seq_mmap(char *, uint32_t, khash_t(RLEN) *, uint8_t);

// For DWalign
typedef int32_t seq_coor_t;

typedef struct {
	seq_coor_t aln_str_size, dist ;
	seq_coor_t aln_q_s, aln_q_e;
	seq_coor_t aln_t_s, aln_t_e;
	char * q_aln_str; char * t_aln_str;

} alignment;

typedef struct {
	seq_coor_t pre_k;
	seq_coor_t x1, y1;
	seq_coor_t x2, y2;
} d_path_data;

typedef struct {
	seq_coor_t d, k;
	seq_coor_t pre_k;
	seq_coor_t x1, y1;
	seq_coor_t x2, y2;
} d_path_data2;

typedef struct {
	seq_coor_t x;
	seq_coor_t y;
} path_point;

typedef struct {
    seq_coor_t s1, e1;
    seq_coor_t s2, e2;
    long int score;
} aln_range;

alignment * align(char *, seq_coor_t,
                  char *, seq_coor_t,
                  seq_coor_t); 

void free_alignment(alignment *);

#ifdef __cplusplus
}
#endif

#endif
