#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "kalloc.h"
#include "khash.h"
#include "kvec.h"
#include "shimmer.h"

#ifndef ORIGINAL
#define ORIGINAL 0
#endif
#ifndef REVERSED
#define REVERSED 1
#endif

uint8_t fourbit_map_f[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 4,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 2, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

uint8_t fourbit_map_r[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 4, 0, 0, 0, 2,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 8, 0, 4, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

void encode_biseq(uint8_t *target, char *seq, size_t len) {
  size_t p, rp;

  for (p = 0, rp = len - 1; p < len; p++, rp--) {
    *(target + p) = ((fourbit_map_r[(uint8_t) * (seq + rp)] << 4) |
                     fourbit_map_f[(uint8_t) * (seq + p)]);
  }
};

char bits_to_base[] = {'N', 'A', 'C', 'N', 'G', 'N', 'N', 'N',
                       'T', 'N', 'N', 'N', 'N', 'N', 'N', 'N'};

void decode_biseq(uint8_t *src, char *seq, size_t len, uint8_t strand) {
  size_t p;
  for (p = 0; p < len; p++) {
    *(seq + p) = strand == ORIGINAL ? bits_to_base[(*(src + p)) & 0x0F]
                                    : bits_to_base[(*(src + p)) >> 4];
  }
};

uint8_t rmap[] = {
    0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,
    15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,
    30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,
    45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,
    60,  61,  62,  63,  64,  84,  66,  71,  68,  69,  70,  67,  72,  73,  74,
    75,  76,  77,  78,  79,  80,  81,  82,  83,  65,  85,  86,  87,  88,  89,
    90,  91,  92,  93,  94,  95,  96,  116, 98,  103, 100, 101, 102, 99,  104,
    105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 97,  117, 118, 119,
    120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134,
    135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149,
    150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164,
    165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179,
    180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194,
    195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209,
    210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224,
    225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
    240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254,
    255};

void reverse_complement(char *seq, size_t len) {
  size_t p, rp;
  char tmp;
  p = 0;
  for (;;) {
    rp = len - p - 1;
    if ((rp - p) <= 1) break;
    tmp = seq[rp];
    seq[rp] = rmap[(uint8_t)seq[p]];
    seq[p] = rmap[(uint8_t)tmp];
    p++;
  }
}

void write_mmlist(char *fn, mm128_v *p) {
  FILE *out_data;
  out_data = fopen(fn, "wb");
  if (!out_data) {
    fprintf(stderr, "file '%s' open error: %s\n", fn, strerror(errno));
    exit(1);
  }
  fwrite(&(kv_size(*p)), sizeof(size_t), 1, out_data);
  fwrite(p->a, sizeof(mm128_t), kv_size(*p), out_data);
  fclose(out_data);
};

mm128_v read_mmlist(char *fn) {
  mm128_v p = {0, 0, 0};
  FILE *in_data;
  in_data = fopen(fn, "rb");
  if (!in_data) {
    fprintf(stderr, "file '%s' open error: %s\n", fn, strerror(errno));
    exit(1);
  }
  fread(&(kv_size(p)), sizeof(size_t), 1, in_data);
  kv_resize(mm128_t, NULL, p, kv_size(p));
  fread(p.a, sizeof(mm128_t), kv_size(p), in_data);
  fclose(in_data);
  return p;
};

void append_mmlist(mm128_v *target, mm128_v *source) {
  kv_resize(mm128_t, NULL, *target, kv_size(*target) + kv_size(*source));
  memcpy(target->a + target->n, source->a, sizeof(mm128_t) * source->n);
  target->n += source->n;
}

void mm_count(mm128_v *p, khash_t(MMC) * mcmap, mm_count_v *cp) {
  uint32_t idx;
  khiter_t k;
  mm128_t mmer;
  uint64_t mhash;
  int32_t absent;

  for (idx = 0; idx < p->n; idx++) {
    mmer = p->a[idx];
    mhash = mmer.x >> 8;
    k = kh_put(MMC, mcmap, mhash, &absent);
    if (absent) {
      kh_value(mcmap, k) = 1;
    } else {
      kh_value(mcmap, k) += 1;
    }
  }
  mm_count_to_vec(mcmap, cp);
}

void mm_count_to_vec(khash_t(MMC) * mcmap, mm_count_v *cp) {
  mm_count_t val;
  khint_t __i;
  for (__i = kh_begin(mcmap); __i != kh_end(mcmap); ++__i) {
    if (!kh_exist(mcmap, __i)) continue;
    val.mer = kh_key(mcmap, __i);
    val.count = kh_val(mcmap, __i);
    kv_push(mm_count_t, NULL, *cp, val);
  }
}

void aggregate_mm_count(khash_t(MMC) * mcmap, mm_count_v *p) {
  uint32_t idx;
  khiter_t k;
  uint64_t mhash;
  int32_t absent;

  for (idx = 0; idx < p->n; idx++) {
    mhash = p->a[idx].mer;
    k = kh_put(MMC, mcmap, mhash, &absent);
    if (absent) {
      kh_value(mcmap, k) = 0;
    }
    kh_value(mcmap, k) += p->a[idx].count;
  }
}

void write_mm_count(char *fn, mm_count_v *p) {
  FILE *out_data;
  out_data = fopen(fn, "wb");
  if (!out_data) {
    fprintf(stderr, "file '%s' open error: %s\n", fn, strerror(errno));
    exit(1);
  }
  fwrite(&(kv_size(*p)), sizeof(size_t), 1, out_data);
  fwrite(p->a, sizeof(mm_count_t), kv_size(*p), out_data);
  fclose(out_data);
};

mm_count_v read_mm_count(char *fn) {
  mm_count_v p = {0, 0, 0};
  FILE *in_data;
  in_data = fopen(fn, "rb");
  if (!in_data) {
    fprintf(stderr, "file '%s' open error: %s\n", fn, strerror(errno));
    exit(1);
  }
  fread(&(kv_size(p)), sizeof(size_t), 1, in_data);
  kv_resize(mm_count_t, NULL, p, kv_size(p));
  fread(p.a, sizeof(mm_count_t), kv_size(p), in_data);
  fclose(in_data);
  return p;
};

khash_t(RID) *
    build_read_index(char *fpath, seq_data_v *seq_data, khash_t(RLEN) * rlmap) {
  // map the original read ID to read length and offset
  khash_t(RID) *hmap = kh_init(RID);
  int absent;
  khiter_t k;
  FILE *index_file;
  uint32_t i, rid, rlen;
  size_t offset;
  char name_buf[256];
  char *name;
  rl_t rl;

  index_file = fopen(fpath, "r");
  while (fscanf(index_file, "%u %255s %u %lu", &rid, name_buf, &rlen,
                &offset) != EOF) {
    kv_resize(seq_data_t, NULL, *seq_data, seq_data->n + 8192);
    name = kmalloc(NULL, 256 * sizeof(char));
    strncpy(name, name_buf, 256);
    i = seq_data->n;
    seq_data->a[i].name = name;
    seq_data->a[i].rid = rid;
    seq_data->n++;
    k = kh_put(RID, hmap, seq_data->a[i].name, &absent);
    if (absent) kh_value(hmap, k) = rid;
    k = kh_put(RLEN, rlmap, rid, &absent);
    rl.len = rlen;
    rl.offset = offset;
    kh_value(rlmap, k) = rl;
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

khash_t(RLEN) * get_read_length_map(char *fpath) {
  khash_t(RLEN) *rlmap = kh_init(RLEN);
  int absent;
  khiter_t k;
  FILE *index_file;
  char name_buf[256];
  uint32_t rid, rlen;
  size_t offset;
  rl_t rl;

  index_file = fopen(fpath, "r");
  while (fscanf(index_file, "%u %255s %u %lu", &rid, name_buf, &rlen,
                &offset) != EOF) {
    k = kh_put(RLEN, rlmap, rid, &absent);
    rl.len = rlen;
    rl.offset = offset;
    kh_value(rlmap, k) = rl;
  }
  fclose(index_file);
  return rlmap;
}

char *get_read_seq(FILE *seqdb, uint32_t rid, khash_t(RLEN) * rlmap) {
  char *seq;
  rl_t rl;
  khiter_t k;
  k = kh_get(RLEN, rlmap, rid);
  assert(k != kh_end(rlmap));
  rl = kh_val(rlmap, k);
  seq = kmalloc(NULL, sizeof(char) * rl.len + 1);
  fseek(seqdb, rl.offset, SEEK_SET);
  fread(seq, sizeof(char), rl.len, seqdb);
  decode_biseq((uint8_t *)seq, seq, rl.len, ORIGINAL);
  seq[rl.len] = 0;  // terminate the string
  return seq;
}

inline uint8_t *get_read_seq_mmap_ptr(uint8_t *seq_p, uint32_t rid,
                                      khash_t(RLEN) * rlmap) {
  rl_t rl;
  khiter_t k;
  k = kh_get(RLEN, rlmap, rid);
  assert(k != kh_end(rlmap));
  rl = kh_val(rlmap, k);
  return seq_p + rl.offset;
}

void build_map(mm128_v *mmers, khash_t(MMER0) * mmer0_map,
               khash_t(RLEN) * rlmap, khash_t(MMC) * mcmap, uint32_t chunk,
               uint32_t total_chunk, uint32_t lowerbound, uint32_t upperbound) {
  uint64_t mhash;
  mm128_t mmer0, mmer1;
  mp128_v *mpv;
  uint32_t rid;
  uint32_t pos, rpos;
  uint32_t span;
  uint32_t mcount = 0;
  int32_t absent;
  mp128_t rp;
  khiter_t k;
  khash_t(MMER1) * mmer1_map;
  size_t s = 0;

  for (;;) {
    if (s >= mmers->n) break;
    mmer0 = mmers->a[s];
    mhash = mmer0.x >> 8;
    k = kh_get(MMC, mcmap, mhash);
    assert(k != kh_end(mcmap));
    mcount = kh_val(mcmap, k);
    if (mcount >= lowerbound && mcount < upperbound) break;
    s++;
  }
  for (size_t i = s + 1; i < mmers->n; i++) {
    mmer1 = mmers->a[i];
    mhash = mmer1.x >> 8;
    k = kh_get(MMC, mcmap, mhash);
    assert(k != kh_end(mcmap));
    mcount = kh_val(mcmap, k);
    if (mcount < lowerbound || mcount > upperbound) continue;

    if ((mmer0.y >> 32) == (mmer1.y >> 32)) {  // the pairs are in the same read

      // don't use two minimers that are too close to each other
      if (((mmer1.y >> 1) & 0xFFFFFFF) - ((mmer0.y >> 1) & 0xFFFFFFF) < 100) {
        mmer0 = mmer1;
        continue;
      }

      if ((mmer0.x >> 8) % total_chunk == chunk % total_chunk) {
        k = kh_put(MMER0, mmer0_map, mmer0.x, &absent);
        if (absent) kh_value(mmer0_map, k) = kh_init(MMER1);

        mmer1_map = kh_value(mmer0_map, k);
        k = kh_put(MMER1, mmer1_map, mmer1.x, &absent);
        if (absent) {
          mpv = kmalloc(NULL, sizeof(mp128_v));
          mpv->n = 0;
          mpv->m = 0;
          mpv->a = NULL;
          kh_value(mmer1_map, k) = mpv;
        } else {
          mpv = kh_value(mmer1_map, k);
        }
        rp.y0 = mmer0.y;
        rp.y1 = mmer1.y;

        rp.direction = ORIGINAL;
        kv_push(mp128_t, NULL, *mpv, rp);
        // printf("Y %lu %lu %09u\n", mmer0.x >> 8, mmer1.x >> 8, (uint32_t)
        // (rp.y0 >> 32));
      }

      // reverse direction
      if ((mmer1.x >> 8) % total_chunk == chunk % total_chunk) {
        k = kh_put(MMER0, mmer0_map, mmer1.x, &absent);
        if (absent) kh_value(mmer0_map, k) = kh_init(MMER1);

        mmer1_map = kh_value(mmer0_map, k);
        k = kh_put(MMER1, mmer1_map, mmer0.x, &absent);
        if (absent) {
          mpv = kmalloc(NULL, sizeof(mp128_v));
          mpv->n = 0;
          mpv->m = 0;
          mpv->a = NULL;
          kh_value(mmer1_map, k) = mpv;
        } else {
          mpv = kh_value(mmer1_map, k);
        }
        rp.y0 = mmer1.y;
        span = mmer1.x & 0xFF;
        rid = rp.y0 >> 32;
        pos = ((rp.y0 & 0xFFFFFFFF) >> 1) + 1;
        k = kh_get(RLEN, rlmap, rid);
        assert(k != kh_end(rlmap));
        rpos = kh_val(rlmap, k).len - pos + span - 1;
        rp.y0 = ((rp.y0 & 0xFFFFFFFF00000001) | (rpos << 1)) ^
                0x1;  // ^0x1 -> switch strand

        rp.y1 = mmer0.y;
        span = mmer0.x & 0xFF;
        rid = rp.y1 >> 32;
        pos = ((rp.y1 & 0xFFFFFFFF) >> 1) + 1;
        k = kh_get(RLEN, rlmap, rid);
        assert(k != kh_end(rlmap));
        rpos = kh_val(rlmap, k).len - pos + span - 1;
        rp.y1 = ((rp.y1 & 0xFFFFFFFF00000001) | (rpos << 1)) ^
                0x1;  // ^0x1 -> switch strand
        rp.direction = REVERSED;

        kv_push(mp128_t, NULL, *mpv, rp);
        // printf("Y %lu %lu %09u\n", mmer1.x >> 8, mmer0.x >> 8, rid);
      }
    }
    mmer0 = mmer1;
  }
}

void build_map4py(mm128_v *mmers, void *mmer0_map, void *rlmap, void *mcmap,
                  uint32_t chunk, uint32_t total_chunk, uint32_t lowerbound,
                  uint32_t upperbound) {
  build_map(mmers, (khash_t(MMER0) *)mmer0_map, (khash_t(RLEN) *)rlmap,
            (khash_t(MMC) *)mcmap, chunk, total_chunk, lowerbound, upperbound);
}

uint32_t mmer_pos(mm128_t *mmer) { return (mmer->y & 0xFFFFFFFF) >> 1; }

void get_ridmm(khash_t(RIDMM) * ridmm, mm128_v *mmers) {
  mm128_t mmer0;
  uint32_t rid;
  khiter_t k;
  mm128_v *_v;
  int32_t absent;
  size_t s = 0;

  for (;;) {
    if (s >= mmers->n) break;
    mmer0 = mmers->a[s];
    rid = (uint32_t)(mmer0.y >> 32);

    k = kh_get(RIDMM, ridmm, rid);
    if (k == kh_end(ridmm)) {
      _v = calloc(sizeof(mm128_v), 1);
      _v->a = mmers->a + s;
      _v->n++;
      _v->m++;
      k = kh_put(RIDMM, ridmm, rid, &absent);
      kh_val(ridmm, k) = _v;
    } else {
      _v = kh_val(ridmm, k);
      _v->n++;
      _v->m++;
    }
    s++;
  }
}
