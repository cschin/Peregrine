#ifndef SHIMMER_H
#define SHIMMER_H

#include <assert.h>
#include "khash.h"


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

KHASH_MAP_INIT_INT(RLEN, uint32_t);
void mm_end_filter(mm128_v *, mm128_v *, mm128_v *, khash_t(RLEN) *, uint32_t); 

#ifdef __cplusplus
}
#endif

#endif
