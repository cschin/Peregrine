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
    seq_coor_t t_m_end, q_m_end;
} ovlp_match_t;

ovlp_match_t * ovlp_match(uint8_t * query_seq,
                          seq_coor_t q_len,
                          uint8_t q_strand,
                          uint8_t * target_seq,
                          seq_coor_t t_len,
                          uint8_t t_strand,
                          seq_coor_t band_tolerance);

void free_ovlp_match(ovlp_match_t * match);

typedef struct { uint64_t x, y; } mm128_t;

typedef struct { size_t n, m; mm128_t *a; } mm128_v;

mm128_v read_mmlist(char *);

void free(void *ptr);

typedef unsigned int khint32_t;

typedef unsigned long khint64_t;

typedef khint32_t khint_t;

typedef struct {
    mm128_v * mmers;
    void * mmer0_map;
    void * rlmap;
    void * mcmap;} py_mmer_t;

void build_shimmer_map4py(py_mmer_t *,
        char *, char *,
        uint32_t, uint32_t, uint32_t, uint32_t);

typedef struct { uint64_t x0, x1, y0, y1; uint8_t direction;} mp256_t;
typedef struct { size_t n, m; mp256_t *a; } mp256_v;

uint32_t get_mmer_count(py_mmer_t * , uint64_t);
void get_shimmer_hits(mp256_v *, py_mmer_t *, uint64_t, uint32_t);


typedef uint32_t mm_idx_t;
typedef struct { size_t n, m; mm_idx_t *a; } mm_idx_v;

typedef struct {
        mm128_v m0;
        mm128_v m1;
        mm_idx_v idx0;
        mm_idx_v idx1;
} shmr_aln_t;

typedef struct { size_t n, m; shmr_aln_t *a; } shmr_aln_v;


shmr_aln_v * shmr_aln( mm128_v *, mm128_v *, uint8_t, double);
void free_shmr_alns(shmr_aln_v *);

""")

ffibuilder.set_source("peregrine._shimmer4py",
               f"""
               #include "{basedir}/src/shimmer.h"
               """,
               sources=[f'{basedir}/src/shimmer4py.c',
                        f'{basedir}/src/DWmatch.c',
                        f'{basedir}/src/shmr_align.c',
                        f'{basedir}/src/shmr_utils.c',
                        f'{basedir}/src/kalloc.c'])   # library name, for the linker

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
