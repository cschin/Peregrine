#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define __STDC_LIMIT_MACROS
#include "kvec.h"
#include "shimmer.h"

typedef struct {
  uint8_t size;
  uint8_t head;
  mm128_t *mers;
} small_m_buffer_t;

static inline uint64_t hash64(uint64_t key, uint64_t mask) {
  key = (~key + (key << 21)) & mask;  // key = (key << 21) - key - 1;
  key = key ^ key >> 24;
  key = ((key + (key << 3)) + (key << 8)) & mask;  // key * 265
  key = key ^ key >> 14;
  key = ((key + (key << 2)) + (key << 4)) & mask;  // key * 21
  key = key ^ key >> 28;
  key = (key + (key << 31)) & mask;
  return key;
}

void pop_push(small_m_buffer_t *smb, mm128_t mer) {
  smb->mers[smb->head] = mer;
  smb->head++;
  smb->head %= smb->size;
}

void find_minimizer(small_m_buffer_t *smb, mm128_t *mmer) {
  uint32_t i = 0;
  uint64_t min_val = UINT64_MAX;
  uint64_t h;

  mmer->x = smb->mers[0].x;
  mmer->y = smb->mers[0].y;
  min_val = smb->mers[0].x >> 8;

  for (i = 1; i < smb->size; i++) {
    h = smb->mers[i].x >> 8;
    if (h < min_val) {
      min_val = h;
      mmer->x = smb->mers[i].x;
      mmer->y = smb->mers[i].y;
    }
  }
}

/* rs: reduction size */
void mm_reduce(mm128_v *p, mm128_v *p_out, uint8_t rs) {
  uint32_t idx;
  uint32_t rid;
  uint32_t rid_ = UINT32_MAX;
  uint32_t r_offset = 0;
  mm128_t mmer, mmer_;
  small_m_buffer_t smb;

  kv_resize(mm128_t, NULL, *p_out, p->n);

  smb.size = rs;
  smb.head = 0;
  smb.mers = (mm128_t *)alloca(sizeof(mm128_t) * smb.size);
  memset(smb.mers, UINT8_MAX, rs * 16);

  mmer_.y = UINT64_MAX;

  for (idx = 0; idx < p->n; idx++, r_offset++) {
    rid = p->a[idx].y >> 32;
    if (rid != rid_) {
      r_offset = 0;
      memset(smb.mers, UINT8_MAX, rs * 16);
      smb.head = 0;
      rid_ = rid;
    }
    pop_push(&smb, p->a[idx]);
    if (r_offset < rs - 1) {
      continue;
    }
    find_minimizer(&smb, &mmer);
    if (mmer.y != mmer_.y) {
      // printf("%lu\n", mmer.x >> 8);
      kv_push(mm128_t, NULL, *p_out, mmer);
      mmer_.x = mmer.x;
      mmer_.y = mmer.y;
    }
  }
}
