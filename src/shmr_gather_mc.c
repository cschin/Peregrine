#include <assert.h>
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>
#include <zlib.h>
#include "khash.h"
#include "kvec.h"
#include "shimmer.h"

extern char *optarg;
extern int optind, opterr, optopt;

int main(int argc, char *argv[]) {
  char *data_path_prefix = NULL;
  char mc_chunk_file_path[8192];
  char mc_file_path[8192];
  int written;
  int total_chunk = 1;
  int chunk = 1;
  int c;

  opterr = 0;

  while ((c = getopt(argc, argv, "p:t:")) != -1) {
    switch (c) {
      case 'p':
        data_path_prefix = optarg;
        break;
      case 't':
        total_chunk = atoi(optarg);
        break;
      case '?':
        if (optopt == 'd') {
          fprintf(stderr,
                  "Option -%c not specified, please specify a prefix of file "
                  "path the data filis\n",
                  optopt);
        }
        if (optopt == 't') {
          fprintf(stderr,
                  "Option -%c not specified, please specify the total number "
                  "of chunks\n",
                  optopt);
        }
        return 1;
      default:
        abort();
    }
  }

  assert(total_chunk > 0);

  if (data_path_prefix == NULL) {
    data_path_prefix = (char *)calloc(8192, 1);
    snprintf(data_path_prefix, 8191, "shimmer");
  }

  khash_t(MMC) *mcmap = kh_init(MMC);

  for (chunk = 1; chunk <= total_chunk; chunk++) {
    mm_count_v mc = {0, 0, 0};
    written = snprintf(mc_chunk_file_path, sizeof mc_chunk_file_path,
                       "%s-MC-%02d-of-%02d.dat", data_path_prefix, chunk,
                       total_chunk);
    assert(written < sizeof(mc_chunk_file_path));
    fprintf(stderr, "input data file: %s\n", mc_chunk_file_path);
    mc = read_mm_count(mc_chunk_file_path);
    aggregate_mm_count(mcmap, &mc);
    kv_destroy(mc);
  }

  mm_count_v mc_all = {0, 0, 0};
  mm_count_to_vec(mcmap, &mc_all);

  written = snprintf(mc_file_path, sizeof mc_file_path, "%s-MC-all.dat",
                     data_path_prefix);
  assert(written < sizeof(mc_file_path));
  fprintf(stderr, "output data file: %s\n", mc_file_path);

  write_mm_count(mc_file_path, &mc_all);

  kv_destroy(mc_all);
  kh_destroy(MMC, mcmap);

  if (!data_path_prefix) free(data_path_prefix);
  return 0;
}
