#include <assert.h>
#include <errno.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <wordexp.h>
#include "kalloc.h"
#include "khash.h"
#include "kvec.h"
#include "shimmer.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>

#define MAX_SMALL_ALNS 4800

shmr_aln_v *shmr_aln(mm128_v *mmers0, mm128_v *mmers1, uint8_t direction,
                     uint32_t max_diff, uint32_t max_dist,
                     uint32_t max_repeat) {
  /* generate a list of co-aligned mimimizer from two
   * minimizer lists: mv1 and mv2
   */
  uint64_t mhash;
  mm128_t mmer0, mmer1;
  khash_t(MMIDX) *mmidx_map = kh_init(MMIDX);
  shmr_aln_v *alns;
  khiter_t k;
  int32_t absent;

  mm_idx_v *idx_tmp;

  alns = calloc(sizeof(shmr_aln_v), 1);

  mm_idx_t s = 0;
  /* build a hasmap from mhash to the index of the minimizer array */
  for (;;) {
    if (s >= mmers0->n) break;
    mmer0 = mmers0->a[s];
    mhash = mmer0.x >> 8;

    k = kh_put(MMIDX, mmidx_map, mhash, &absent);
    if (absent) {
      idx_tmp = calloc(sizeof(mm_idx_v), 1);
      kv_push(mm_idx_t, 0, *idx_tmp, s);
      kh_val(mmidx_map, k) = idx_tmp;
    } else {
      k = kh_get(MMIDX, mmidx_map, mhash);
      assert(k != kh_end(mmidx_map));
      idx_tmp = kh_val(mmidx_map, k);
      kv_push(mm_idx_t, 0, *idx_tmp, s);
    }
    s++;
  }

  /* loop through 2nd shimmer list to build alginements */
  mm_idx_t ss = 0;
  uint32_t small_aln_count = 0;
  for (;;) {
    if (ss >= mmers1->n) break;
    if (direction == 1) {  // reversed
      s = mmers1->n - ss;
    } else {
      s = ss;
    }
    mmer1 = mmers1->a[s];
    mhash = mmer1.x >> 8;
    k = kh_get(MMIDX, mmidx_map, mhash);
    if (k == kh_end(mmidx_map)) {
      ss++;
      continue;
    }
    idx_tmp = kh_val(mmidx_map, k);
    if (idx_tmp->n > max_repeat) {
      ss++;
      continue;
    }

    for (uint32_t i = 0; i < idx_tmp->n; i++) {
      mmer0 = mmers0->a[idx_tmp->a[i]];
      int64_t delta0, delta1;
      int64_t mm_dist;
      if (direction == 0 && (mmer0.y & 0x1) != (mmer1.y & 0x1)) {
        continue;
      }

      if (direction == 1 && (mmer0.y & 0x1) == (mmer1.y & 0x1)) {
        continue;
      }

      if (direction == 1) {
        delta0 = abs(mmer_pos(&mmer0) + mmer_pos(&mmer1));
      } else {
        delta0 = abs(mmer_pos(&mmer0) - mmer_pos(&mmer1));
      }
      uint32_t best_aln_idx = UINT32_MAX;
      double min_diff = max_diff;
      uint8_t best_found = 0;
      small_aln_count = 0;
      for (uint32_t aln_idx = 0; aln_idx < alns->n; aln_idx++) {
        mm128_t m0, m1;
        shmr_aln_t *aln;
        size_t n;
        aln = alns->a + aln_idx;
        n = aln->idx0.n;

        if (n < 3) small_aln_count++;

        if (idx_tmp->a[i] < aln->idx0.a[n - 1]) continue;

        m0 = mmers0->a[aln->idx0.a[n - 1]];
        m1 = mmers1->a[aln->idx1.a[n - 1]];

        mm_dist = abs(mmer_pos(&mmer0) - mmer_pos(&m0));
        if (mm_dist >= max_dist) continue;

        if (direction == 1) {
          delta1 = abs(mmer_pos(&m0) + mmer_pos(&m1));
        } else {
          delta1 = abs(mmer_pos(&m0) - mmer_pos(&m1));
        }
        // double diff = (double) abs(delta0 - delta1) / (double) (mm_dist);
        uint32_t diff = (uint32_t)abs((int32_t)delta0 - (int32_t)delta1);
        if (diff < max_diff && diff < min_diff && mm_dist < max_dist) {
          min_diff = diff;
          best_aln_idx = aln_idx;
          best_found = 1;
        }
      }
      if (best_found == 1) {
        shmr_aln_t *aln;
        aln = alns->a + best_aln_idx;
        kv_push(mm_idx_t, 0, aln->idx0, idx_tmp->a[i]);
        kv_push(mm_idx_t, 0, aln->idx1, s);
        // printf("best %d %d %d\n", best_aln_idx, idx_tmp->a[i], s);
      } else {
        shmr_aln_t *aln;
        aln = calloc(sizeof(shmr_aln_t), 1);
        kv_push(mm_idx_t, 0, aln->idx0, idx_tmp->a[i]);
        kv_push(mm_idx_t, 0, aln->idx1, s);
        kv_push(shmr_aln_t, 0, *alns, *aln);
        // printf("new %d %d %d\n", best_aln_idx, idx_tmp->a[i], s);
      }
    }
    ss++;
    if (small_aln_count > MAX_SMALL_ALNS) break;
  }

  for (khiter_t __i = kh_begin(mmidx_map); __i != kh_end(mmidx_map); ++__i) {
    if (!kh_exist(mmidx_map, __i)) continue;
    idx_tmp = kh_val(mmidx_map, __i);
    kv_destroy(*idx_tmp);
    free(idx_tmp);
  }
  kh_destroy(MMIDX, mmidx_map);
  return alns;
}

void free_shmr_alns(shmr_aln_v *alns) {
  for (uint32_t aln_idx = 0; aln_idx < alns->n; aln_idx++) {
    kv_destroy(alns->a[aln_idx].idx0);
    kv_destroy(alns->a[aln_idx].idx1);
  }
  kv_destroy(*alns);
  free(alns);
}
