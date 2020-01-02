#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define __STDC_LIMIT_MACROS
#include "khash.h"
#include "kvec.h"
#include "shimmer.h"

void mm_end_filter(mm128_v *p, mm128_v *p_out_5, mm128_v *p_out_3,
                   khash_t(RLEN) * rlmap, uint32_t end_length) {
  uint32_t idx;
  uint32_t rid;
  uint32_t rlen;
  uint32_t pos, r_pos, span;
  khiter_t k;
  mm128_t mmer;

  for (idx = 0; idx < p->n; idx++) {
    mmer = p->a[idx];
    rid = mmer.y >> 32;
    span = mmer.x & 0xFF;
    k = kh_get(RLEN, rlmap, rid);
    // is_missing = (k == kh_end(hmap));
    rlen = kh_value(rlmap, k).len;
    pos = ((mmer.y & 0xFFFFFFFF) >> 1) + 1;
    r_pos = rlen - pos + span;
    if (pos < end_length) {
      kv_push(mm128_t, NULL, *p_out_5, mmer);
    };
    if (r_pos < end_length) {
      kv_push(mm128_t, NULL, *p_out_3, mmer);
    };
  }
}
