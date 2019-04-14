from cffi import FFI
import os

basedir = os.environ["peregrine_base"]

ffibuilder = FFI()

ffibuilder.cdef("""

typedef int seq_coor_t;

typedef struct {
    seq_coor_t aln_str_size ;
    seq_coor_t dist ;
    seq_coor_t aln_q_s;
    seq_coor_t aln_q_e;
    seq_coor_t aln_t_s;
    seq_coor_t aln_t_e;
    char * q_aln_str;
    char * t_aln_str;

} alignment;

typedef struct {
    seq_coor_t t_pos;
    uint8_t delta;
    char q_base;
    seq_coor_t p_t_pos;   // the tag position of the previous base
    uint8_t p_delta; // the tag delta of the previous base
    char p_q_base;        // the previous base
    unsigned q_id;
} align_tag_t;

typedef struct {
    seq_coor_t len;
    align_tag_t * align_tags;
} align_tags_t;

typedef struct {
    seq_coor_t s1;
    seq_coor_t e1;
    seq_coor_t s2;
    seq_coor_t e2;
    long int score;
} aln_range;

typedef struct {
    char * sequence;
    uint8_t * eqv;
} consensus_data;


align_tags_t * get_align_tags( char * aln_q_seq,
                            char * aln_t_seq,
                            seq_coor_t aln_seq_len,
                            aln_range * range,
                            unsigned q_id,
                            seq_coor_t t_offset);

void free_align_tags( align_tags_t * tags);

consensus_data * get_cns_from_align_tags( align_tags_t ** tag_seqs,
                                          unsigned n_tag_seqs,
                                          unsigned t_len,
                                          unsigned min_cov );

void free_consensus_data( consensus_data * consensus );

alignment * align(char * query_seq, seq_coor_t q_len,
                  char * target_seq, seq_coor_t t_len,
                  seq_coor_t band_tolerance,
                  int get_aln_str);

void free_alignment(alignment *);
void *malloc(size_t size);
void free(void *ptr);
""")

ffibuilder.set_source("peregrine._falcon4py",
               f"""
               #include "{basedir}/falcon/common.h"
               #include "{basedir}/falcon/falcon.h"
               """, sources = [f'{basedir}/falcon/falcon.c',
                               f'{basedir}/falcon/DW_banded.c',
                               f'{basedir}/falcon/kalloc.c'])   # library name, for the linker

if __name__ == "__main__":
    import sys
    ffibuilder.compile(verbose=True)
