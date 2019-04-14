#include "kvec.h"
#include "khash.h"
#include "kalloc.h"

typedef struct {
    seq_coor_t t_pos;
    uint8_t delta;
    char q_base;
    seq_coor_t p_t_pos;   // the tag position of the previous base
    uint8_t p_delta;      // the tag delta of the previous base
    char p_q_base;        // the previous base
    unsigned q_id;
} align_tag_t;

typedef struct {
    seq_coor_t len;
    align_tag_t * align_tags;
} align_tags_t;



typedef struct { size_t n, m; uint64_t *a; } uint64_v;

typedef struct {
    uint64_t ctag_key;
    uint64_t ptag_key;
    uint16_t coverage;
    uint16_t count;
    double score;
} align_edge_t; 

typedef struct { size_t n, m; align_edge_t *a; } align_edge_v;

KHASH_MAP_INIT_INT64(PTAG, uint16_t);
typedef khash_t(PTAG) ptag_to_count_t; 
KHASH_MAP_INIT_INT64(CTAG, khash_t(PTAG) *);
typedef khash_t(CTAG) ctag_to_ptag_t; 

typedef struct {
    uint64_t ctag_key;
    align_edge_t * best_edge;
    double best_score;
} align_node_t;

KHASH_MAP_INIT_INT64(NODE, align_node_t *);
typedef khash_t(NODE) align_node_map_t; 

align_tags_t * get_align_tags( char *, char *, seq_coor_t, aln_range *, unsigned, seq_coor_t);
void free_align_tags( align_tags_t * tags);
consensus_data * get_cns_from_align_tags( align_tags_t **, unsigned, unsigned, unsigned ); 
