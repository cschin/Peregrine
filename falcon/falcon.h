
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
    uint16_t size;
    uint16_t n_link;
    seq_coor_t * p_t_pos;   // the tag position of the previous base
    uint8_t * p_delta; // the tag delta of the previous base
    char * p_q_base;        // the previous base
    uint16_t * link_count;
    uint16_t count;
    seq_coor_t best_p_t_pos;
    uint8_t best_p_delta;
    uint8_t best_p_q_base; // encoded base
    double score;
} align_tag_col_t;

typedef struct {
    align_tag_col_t * base;
} msa_base_group_t;

typedef struct {
    uint8_t size;
    uint8_t max_delta;
    msa_base_group_t * delta;
} msa_delta_group_t;

typedef msa_delta_group_t * msa_pos_t;

align_tags_t * get_align_tags( char *, char *, seq_coor_t, aln_range *, unsigned, seq_coor_t);
void free_align_tags( align_tags_t * tags);
consensus_data * get_cns_from_align_tags( align_tags_t **, unsigned, unsigned, unsigned ); 
