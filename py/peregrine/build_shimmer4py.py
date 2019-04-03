from cffi import FFI
import os

basedir = os.environ["peregrine_base"]

ffibuilder = FFI()

ffibuilder.cdef("""
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

ffibuilder.set_source("peregrine._shimmer4py",
               f"""
               #include "{basedir}/src/shimmer.h"
               """,
               sources=[f'{basedir}/src/shimmer4py.c',
                        f'{basedir}/src/DWmatch.c',
                        f'{basedir}/src/shmr_utils.c',
                        f'{basedir}/src/kalloc.c'])   # library name, for the linker

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
