#ifndef SHIMMER_H
#define SHIMMER_H

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>
#include "khash.h"
#include "kvec.h"

#define ORIGINAL 0
#define REVERSED 1

#ifdef __cplusplus
extern "C" {
#endif

void encode_biseq(uint8_t *, char *, size_t);

void decode_biseq(uint8_t *, char *, size_t, uint8_t);

void reverse_complement(char *, size_t);

typedef struct {
  uint64_t x, y;
} mm128_t;
typedef struct {
  size_t n, m;
  mm128_t *a;
} mm128_v;
void mm_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid,
               int is_hpc, mm128_v *p);
void mm_reduce(mm128_v *, mm128_v *, uint8_t);
uint32_t ks_ksmall_uint32_t(size_t n, uint32_t arr[], size_t kk);

void write_mmlist(char *, mm128_v *);
mm128_v read_mmlist(char *);
void append_mmlist(mm128_v *, mm128_v *);

typedef struct {
  char *name;
  uint32_t rid;
} seq_data_t;
typedef struct {
  size_t n, m;
  seq_data_t *a;
} seq_data_v;

typedef struct {
  uint32_t len;
  size_t offset;
} rl_t;
KHASH_MAP_INIT_STR(RID, uint32_t);
KHASH_MAP_INIT_INT(RLEN, rl_t);
khash_t(RID) * build_read_index(char *, seq_data_v *, khash_t(RLEN) *);
khash_t(RLEN) * get_read_length_map(char *);

void mm_end_filter(mm128_v *, mm128_v *, mm128_v *, khash_t(RLEN) *, uint32_t);

KHASH_MAP_INIT_INT64(MMC, uint32_t);
typedef struct {
  uint64_t mer;
  uint32_t count;
} mm_count_t;
typedef struct {
  size_t n, m;
  mm_count_t *a;
} mm_count_v;
void mm_count(mm128_v *, khash_t(MMC) *, mm_count_v *);
void write_mm_count(char *, mm_count_v *);
void mm_count_to_vec(khash_t(MMC) *, mm_count_v *);
mm_count_v read_mm_count(char *fn);

void aggregate_mm_count(khash_t(MMC) *, mm_count_v *);

typedef struct {
  uint64_t y0, y1;
  uint8_t direction;
} mp128_t;
typedef struct {
  size_t n, m;
  mp128_t *a;
} mp128_v;
KHASH_MAP_INIT_INT64(MMER1, mp128_v *);
typedef khash_t(MMER1) * mmert1_p_t;
KHASH_MAP_INIT_INT64(MMER0, mmert1_p_t);

void build_map(mm128_v *, khash_t(MMER0) *, khash_t(RLEN) *, khash_t(MMC) *,
               uint32_t, uint32_t, uint32_t, uint32_t);

char *get_read_seq(FILE *, uint32_t, khash_t(RLEN) *);
uint8_t *get_read_seq_mmap_ptr(uint8_t *, uint32_t, khash_t(RLEN) *);

// For DWmatch
typedef int32_t seq_coor_t;

typedef struct {
  seq_coor_t m_size, dist;
  seq_coor_t q_bgn, q_end;
  seq_coor_t t_bgn, t_end;
  seq_coor_t t_m_end, q_m_end;
} ovlp_match_t;

typedef struct {
  uint64_t y0, y1;
  uint32_t rl0, rl1;
  uint8_t strand0, strand1;
  uint8_t ovlp_type;
  ovlp_match_t match;
} ovlp_t;

typedef struct {
  seq_coor_t s1, e1;
  seq_coor_t s2, e2;
  long int score;
} match_range;

ovlp_match_t *ovlp_match(uint8_t *, seq_coor_t, uint8_t, uint8_t *, seq_coor_t,
                         uint8_t, seq_coor_t);

void free_ovlp_match(ovlp_match_t *);

typedef struct {
  uint64_t x0, x1, y0, y1;
  uint8_t direction;
} mp256_t;
typedef struct {
  size_t n, m;
  mp256_t *a;
} mp256_v;

typedef struct {
  mm128_v *mmers;
  void *mmer0_map;
  void *rlmap;
  void *mcmap;
  void *ridmm;
} py_mmer_t;

// for shmr_align
typedef uint32_t mm_idx_t;
typedef kvec_t(mm_idx_t) mm_idx_v;

typedef struct {
  mm_idx_v idx0;
  mm_idx_v idx1;
} shmr_aln_t;

typedef kvec_t(shmr_aln_t) shmr_aln_v;

KHASH_MAP_INIT_INT64(MMIDX, mm_idx_v *);
shmr_aln_v *shmr_aln(mm128_v *, mm128_v *, uint8_t, uint32_t, uint32_t,
                     uint32_t);

void free_shmr_alns(shmr_aln_v *);

KHASH_MAP_INIT_INT(RIDMM, mm128_v *);
void get_ridmm(khash_t(RIDMM) *, mm128_v *);
uint32_t mmer_pos(mm128_t *);

#ifdef __cplusplus
}
#endif

#endif
