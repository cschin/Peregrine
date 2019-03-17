from cffi import FFI

ffi = FFI()

ffi.cdef("""
void decode_biseq(uint8_t * src, char * seq,
                  size_t len, uint8_t strand);

typedef int32_t seq_coor_t;

typedef struct {
    seq_coor_t m_size, dist ;
    seq_coor_t q_bgn, q_end;
    seq_coor_t t_bgn, t_end;
} ovlp_match_t;

ovlp_match_t * ovlp_match(uint8_t * query_seq,
                          seq_coor_t q_len,
                          uint8_t q_strand,
                          uint8_t * target_seq,
                          seq_coor_t t_len,
                          uint8_t t_strand,
                          seq_coor_t band_tolerance);

void free_ovlp_match(ovlp_match_t * match);
""")

ffi.set_source("_shimmer4py",
               """
               #include "../src/shimmer.h"
               """,
               sources=['../src/shimmer4py.c',
                        '../src/DWmatch.c',
                        '../src/shmr_utils.c',
                        '../src/kalloc.c'])   # library name, for the linker

if __name__ == "__main__":
    ffi.compile(verbose=True)
