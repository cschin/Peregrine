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

extern char *optarg;
extern int optind, opterr, optopt;

#define handle_error(msg) \
  do {                    \
    perror(msg);          \
    exit(EXIT_FAILURE);   \
  } while (0)

#define MMER_COUNT_LOWER_BOUND 1
#define MMER_COUNT_UPPER_BOUND 240
#ifndef ORIGINAL
#define ORIGINAL 0
#endif
#ifndef REVERSED
#define REVERSED 1
#endif
#define READ_END_FUZZINESS 48
#define LOCAL_OVERLAP_UPPERBOUND 120
#define ALNBANDSIZE 100

KHASH_MAP_INIT_INT64(RPAIR, uint8_t);

int mp128_comp(const void *a, const void *b) {
  mp128_t *a0 = (mp128_t *)a;
  mp128_t *b0 = (mp128_t *)b;
  return ((a0->y0 & 0xFFFFFFFF) >> 1) < ((b0->y0 & 0xFFFFFFFF) >> 1);
}

void process_map(char *refdb_file_path, char *seqdb_file_path,
                 mm128_v *ref_mmers, khash_t(RLEN) * ref_lmap,
                 khash_t(MMER0) * mmer0_map, khash_t(RLEN) * rlmap,
                 khash_t(MMC) * mcmap, uint32_t lowerbound,
                 uint32_t upperbound) {
  int rfd, sfd;
  struct stat rsb, ssb;
  uint8_t *rseq_p, *seq_p;
  mp128_v *mpv;

  khash_t(MMER1) * mmer1_map;

  rfd = open(refdb_file_path, O_RDONLY);
  if (rfd == -1) handle_error("open");

  if (fstat(rfd, &rsb) == -1) /* To obtain file size */
    handle_error("fstat");

  rseq_p =
      (uint8_t *)mmap((void *)0, rsb.st_size, PROT_READ, MAP_SHARED, rfd, 0);

  sfd = open(seqdb_file_path, O_RDONLY);
  if (sfd == -1) handle_error("open");

  if (fstat(sfd, &ssb) == -1) /* To obtain file size */
    handle_error("fstat");

  seq_p =
      (uint8_t *)mmap((void *)0, ssb.st_size, PROT_READ, MAP_SHARED, sfd, 0);

  // clock_t begin = clock();
  // clock_t end;
  mm128_t mmer0, mmer1;
  khiter_t k;

  size_t s = 0;
  assert(ref_mmers->n > 0);
  for (;;) {
    mmer0 = ref_mmers->a[s];
    if (s >= ref_mmers->n) break;
    k = kh_get(MMER0, mmer0_map, mmer0.x);
    if (k != kh_end(mmer0_map)) break;
    s++;
  }

  for (size_t i = s + 1; i < ref_mmers->n; i++) {
    mmer1 = ref_mmers->a[i];
    uint64_t mhash = mmer1.x >> 8;
    k = kh_get(MMC, mcmap, mhash);
    if (k == kh_end(mcmap)) continue;
    uint32_t mcount = kh_val(mcmap, k);
    if (mcount < lowerbound || mcount > upperbound) continue;

    if ((mmer0.y >> 32) != (mmer1.y >> 32)) {
      mmer0 = mmer1;
      continue;  // the pairs are in the same read
    }

    k = kh_get(MMER0, mmer0_map, mmer0.x);
    if (k == kh_end(mmer0_map)) {
      mmer0 = mmer1;
      continue;
    }

    mmer1_map = kh_val(mmer0_map, k);
    k = kh_get(MMER1, mmer1_map, mmer1.x);
    if (k == kh_end(mmer1_map)) {
      mmer0 = mmer1;
      continue;
    }

    if (((mmer1.y >> 1) & 0xFFFFFFF) - ((mmer0.y >> 1) & 0xFFFFFFF) < 100) {
      mmer0 = mmer1;
      continue;
    }

    mpv = kh_val(mmer1_map, k);

    uint32_t ref_id;
    uint32_t ref_bgn;
    uint32_t ref_end;
    ref_id = (uint32_t)(mmer0.y >> 32);
    ref_bgn = (uint32_t)((mmer0.y & 0xFFFFFFFF) >> 1);
    ref_end = (uint32_t)((mmer1.y & 0xFFFFFFFF) >> 1);

    for (int j = 0; j < mpv->n; j++) {
      uint32_t read_id;
      uint32_t read_bgn;
      uint32_t read_end;
      uint8_t read_direction;

      read_id = mpv->a[j].y0 >> 32;
      read_bgn = (uint32_t)((mpv->a[j].y0 & 0xFFFFFFFF) >> 1);
      read_end = (uint32_t)((mpv->a[j].y1 & 0xFFFFFFFF) >> 1);
      read_direction = mpv->a[j].direction;
      assert(read_bgn < read_end);

      uint64_t mhash = mmer0.x >> 8;
      k = kh_get(MMC, mcmap, mhash);
      assert(k != kh_end(mcmap));
      uint32_t mcount0 = kh_val(mcmap, k);
      mhash = mmer1.x >> 8;
      k = kh_get(MMC, mcmap, mhash);
      assert(k != kh_end(mcmap));
      uint32_t mcount1 = kh_val(mcmap, k);
      printf("%u %u %u %u %u %u %d %u %u\n", ref_id, ref_bgn, ref_end, read_id,
             read_bgn, read_end, read_direction, mcount0, mcount1);
    }
    mmer0 = mmer1;
  }

  munmap(rseq_p, rsb.st_size);
  munmap(seq_p, ssb.st_size);
}

int main(int argc, char *argv[]) {
  char *refdb_prefix = NULL;
  char *seqdb_prefix = NULL;
  char *ref_shimmer_prefix = NULL;
  char *shimmer_prefix = NULL;

  char mmc_file_path[8192];
  char mmer_file_path[8192];
  char ref_idx_file_path[8192];
  char refdb_file_path[8192];
  char seq_idx_file_path[8192];
  char seqdb_file_path[8192];
  int c;
  uint32_t total_chunk = 1, mychunk = 1;

  uint32_t mc_upper = MMER_COUNT_UPPER_BOUND;
  uint32_t mc_lower = MMER_COUNT_LOWER_BOUND;

  wordexp_t p;
  char **mmc_fns;
  char **shimmer_fns;

  mm128_v ref_mmers = {0, 0, 0};
  mm128_v mmers = {0, 0, 0};
  mm128_v mmers_;
  mm_count_v mmc;

  khash_t(RLEN) * ref_lmap;
  khash_t(RLEN) * rlmap;
  khash_t(MMC) *mcmap = kh_init(MMC);

  khash_t(MMER0) * mmer0_map;
  khash_t(MMER1) * mmer1_map;

  mp128_v *mpv;

  opterr = 0;

  while ((c = getopt(argc, argv, "r:m:p:l:M:n:t:c:b:")) != -1) {
    switch (c) {
      case 'r':
        refdb_prefix = optarg;
        break;
      case 'm':
        ref_shimmer_prefix = optarg;
        break;
      case 'p':
        seqdb_prefix = optarg;
        break;
      case 'l':
        shimmer_prefix = optarg;
        break;
      case 'M':
        mc_upper = atoi(optarg);
        break;
      case 'n':
        mc_lower = atoi(optarg);
        break;
      case 't':
        total_chunk = atoi(optarg);
        break;
      case 'c':
        mychunk = atoi(optarg);
        break;
      case '?':
        if (optopt == 'r') {
          fprintf(stderr,
                  "Option -%c not specified, using 'ref' as the ref sequence "
                  "db prefix\n",
                  optopt);
        }
        if (optopt == 'p') {
          fprintf(stderr,
                  "Option -%c not specified, using 'seq_dataset' as the "
                  "sequence db prefix\n",
                  optopt);
        }
        if (optopt == 'l') {
          fprintf(stderr,
                  "Option -%c not specified, using 'shimmer-L2' as the L2 "
                  "index prefix\n",
                  optopt);
        }
        return 1;
      default:
        abort();
    }
  }

  assert(total_chunk > 0);
  assert(mychunk > 0 && mychunk <= total_chunk);

  if (refdb_prefix == NULL) {
    refdb_prefix = (char *)calloc(8192, 1);
    snprintf(refdb_prefix, 8191, "ref");
  }

  if (ref_shimmer_prefix == NULL) {
    ref_shimmer_prefix = (char *)calloc(8192, 1);
    snprintf(ref_shimmer_prefix, 8191, "ref-L2");
  }

  if (seqdb_prefix == NULL) {
    seqdb_prefix = (char *)calloc(8192, 1);
    snprintf(seqdb_prefix, 8191, "seq_dataset");
  }

  if (shimmer_prefix == NULL) {
    shimmer_prefix = (char *)calloc(8192, 1);
    snprintf(shimmer_prefix, 8191, "shimmer-L2");
  }

  int written;
  written = snprintf(ref_idx_file_path, sizeof(ref_idx_file_path), "%s.idx",
                     refdb_prefix);
  assert(written < sizeof(ref_idx_file_path));
  fprintf(stderr, "using ref index file: %s\n", ref_idx_file_path);

  ref_lmap = get_read_length_map(ref_idx_file_path);

  written = snprintf(refdb_file_path, sizeof(seqdb_file_path), "%s.seqdb",
                     refdb_prefix);
  assert(written < sizeof(refdb_file_path));
  fprintf(stderr, "using ref seqdb file: %s\n", refdb_file_path);

  written = snprintf(mmer_file_path, sizeof(mmer_file_path),
                     "%s-[0-9]*-of-[0-9]*.dat", ref_shimmer_prefix);
  assert(written < sizeof(mmer_file_path));
  wordexp(mmer_file_path, &p, 0);
  shimmer_fns = p.we_wordv;
  for (int i = 0; i < p.we_wordc; i++) {
    fprintf(stderr, "using ref shimmer data file: %s\n", shimmer_fns[i]);
    mmers_ = read_mmlist(shimmer_fns[i]);
    fprintf(stderr, "number of shimmers load: %lu\n", mmers_.n);
    append_mmlist(&ref_mmers, &mmers_);
    kv_destroy(mmers_);
  }
  wordfree(&p);

  written = snprintf(seq_idx_file_path, sizeof(seq_idx_file_path), "%s.idx",
                     seqdb_prefix);
  assert(written < sizeof(seq_idx_file_path));
  fprintf(stderr, "using index file: %s\n", seq_idx_file_path);

  rlmap = get_read_length_map(seq_idx_file_path);

  written = snprintf(seqdb_file_path, sizeof(seqdb_file_path), "%s.seqdb",
                     seqdb_prefix);
  assert(written < sizeof(seqdb_file_path));
  fprintf(stderr, "using seqdb file: %s\n", seqdb_file_path);

  written = snprintf(mmer_file_path, sizeof(mmer_file_path),
                     "%s-[0-9]*-of-[0-9]*.dat", shimmer_prefix);

  assert(written < sizeof(mmer_file_path));
  wordexp(mmer_file_path, &p, 0);
  shimmer_fns = p.we_wordv;
  for (int i = 0; i < p.we_wordc; i++) {
    fprintf(stderr, "using shimmer data file: %s\n", shimmer_fns[i]);
    mmers_ = read_mmlist(shimmer_fns[i]);
    fprintf(stderr, "number of shimmers load: %lu\n", mmers_.n);
    append_mmlist(&mmers, &mmers_);
    kv_destroy(mmers_);
  }
  wordfree(&p);

  char buffer[32768];

  setvbuf(stdout, buffer, _IOFBF, sizeof(buffer));

  written = snprintf(mmc_file_path, sizeof(mmc_file_path),
                     "%s-MC-[0-9]*-of-[0-9]*.dat", shimmer_prefix);

  assert(written < sizeof(mmc_file_path));
  wordexp(mmc_file_path, &p, 0);
  mmc_fns = p.we_wordv;
  for (int i = 0; i < p.we_wordc; i++) {
    fprintf(stderr, "using shimmer count file: %s\n", mmc_fns[i]);
    mmc = read_mm_count(mmc_fns[i]);
    aggregate_mm_count(mcmap, &mmc);
    kv_destroy(mmc);
  }

  wordfree(&p);

  mmer0_map = kh_init(MMER0);

  build_map(&mmers, mmer0_map, rlmap, mcmap, mychunk, total_chunk, mc_lower,
            mc_upper);

  process_map(refdb_file_path, seqdb_file_path, &ref_mmers, ref_lmap, mmer0_map,
              rlmap, mcmap, mc_lower, mc_upper);

  for (khiter_t __i = kh_begin(mmer0_map); __i != kh_end(mmer0_map); ++__i) {
    if (!kh_exist(mmer0_map, __i)) continue;
    mmer1_map = kh_val(mmer0_map, __i);
    for (khiter_t __j = kh_begin(mmer1_map); __j != kh_end(mmer1_map); ++__j) {
      if (!kh_exist(mmer1_map, __j)) continue;
      mpv = kh_val(mmer1_map, __j);
      kv_destroy(*mpv);
    }
    kh_destroy(MMER1, mmer1_map);
  }

  kh_destroy(MMER0, mmer0_map);
  kh_destroy(MMC, mcmap);
  kh_destroy(RLEN, rlmap);
  kv_destroy(mmers);
  kv_destroy(ref_mmers);
  fflush(stdout);
}
