#ifndef SHIMMER_H
#define SHIMMER_H

#include <assert.h>
#include "khash.h"
#include "ksort.h"


#ifdef __cplusplus
extern "C" {
#endif

typedef struct { uint64_t x, y; } mm128_t;
typedef struct { size_t n, m; mm128_t *a; } mm128_v;
void mm_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p);
void mm_reduce(mm128_v *, mm128_v *,  uint8_t);
uint32_t ks_ksmall_uint32_t(size_t n, uint32_t arr[], size_t kk);

void write_mmlist(char *, mm128_v *); 
mm128_v read_mmlist(char *); 

typedef struct { char * name; uint32_t rid;} seq_data_t; 
typedef struct { size_t n, m; seq_data_t * a; } seq_data_v;

KHASH_MAP_INIT_STR(RID, uint32_t);
KHASH_MAP_INIT_INT(RLEN, uint32_t);
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

#define sort_key_128x(a) ((a).x)
KRADIX_SORT_INIT(128x, mm128_t, sort_key_128x, 8);

#define sort_key_256x(a) ( ( (((uint64_t) (a).x0) << 16) & 0xFFFFFFFF00000000) |  ( (uint64_t) ((a).x1 >> 24 )))
KRADIX_SORT_INIT(256x, mm256_t, sort_key_256x, 8);

#define mm256x_lt(a, b) ( 1 ? ((a).x0 < (b).x0) :  ((a).x1 < (b).x1)   )
KSORT_INIT(256x, mm256_t, mm256x_lt);
#ifdef __cplusplus
}
#endif

#endif
