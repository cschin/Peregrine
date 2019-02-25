#ifndef SHIMMER_H
#define SHIMMER_H

#include <assert.h>


#ifdef __cplusplus
extern "C" {
#endif

typedef struct { uint64_t x, y; } mm128_t;
typedef struct { size_t n, m; mm128_t *a; } mm128_v;
void mm_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p);
void mm_select(mm128_v *, mm128_v *,  uint8_t);
uint32_t ks_ksmall_uint32_t(size_t n, uint32_t arr[], size_t kk);

void write_mmlist(char *, mm128_v *); 
mm128_v read_mmlist(char *); 

#ifdef __cplusplus
}
#endif

#endif
