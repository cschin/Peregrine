/*
 * =====================================================================================
 *
 *       Filename:  fastcon.c
 *
 *    Description:
 *
 *        Version:  0.1
 *        Created:  07/20/2013 17:00:00
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jason Chin,
 *        Company:
 *
 * =====================================================================================

 #################################################################################$$
 # Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
 #
 # All rights reserved.
 #
 # Redistribution and use in source and binary forms, with or without
 # modification, are permitted (subject to the limitations in the
 # disclaimer below) provided that the following conditions are met:
 #
 #  * Redistributions of source code must retain the above copyright
 #  notice, this list of conditions and the following disclaimer.
 #
 #  * Redistributions in binary form must reproduce the above
 #  copyright notice, this list of conditions and the following
 #  disclaimer in the documentation and/or other materials provided
 #  with the distribution.
 #
 #  * Neither the name of Pacific Biosciences nor the names of its
 #  contributors may be used to endorse or promote products derived
 #  from this software without specific prior written permission.
 #
 # NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
 # GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
 # BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
 # WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 # OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 # DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
 # CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 # SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 # LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 # USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 # ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 # OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 # OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 # SUCH DAMAGE.
 #################################################################################$$
 */

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include "common.h"
#include "falcon.h"


align_tags_t * get_align_tags( char * aln_q_seq,
                               char * aln_t_seq,
                               seq_coor_t aln_seq_len,
                               aln_range * range,
                               unsigned q_id,
                               seq_coor_t t_offset) {
    char p_q_base;
    align_tags_t * tags;
    seq_coor_t i, j, jj, k, p_j, p_jj;

    tags = calloc( 1, sizeof(align_tags_t) );
    tags->len = aln_seq_len;
    tags->align_tags = calloc( aln_seq_len + 1, sizeof(align_tag_t) );
    i = range->s1 - 1;
    j = range->s2 - 1;
    jj = 0;
    p_j = -1;
    p_jj = 0;
    p_q_base = '.';

    for (k = 0; k < aln_seq_len; k++) {
        if (aln_q_seq[k] != '-') {
            i ++;
			if (jj < 12) {  // we cap the biggest insert gap up to 12 bases
				jj ++;
			}
        }
        if (aln_t_seq[k] != '-') {
            j ++;
            jj = 0;
        }
        //printf("t %d %d %d %c %c\n", q_id, j, jj, aln_t_seq[k], aln_q_seq[k]);

        if ( j + t_offset >= 0 && jj < UINT8_MAX && p_jj < UINT8_MAX) {
			(tags->align_tags[k]).t_pos = j + t_offset;
			(tags->align_tags[k]).delta = jj;
			(tags->align_tags[k]).p_t_pos = p_j + t_offset;
			(tags->align_tags[k]).p_delta = p_jj;
			(tags->align_tags[k]).p_q_base = p_q_base;
			(tags->align_tags[k]).q_base = aln_q_seq[k];
			(tags->align_tags[k]).q_id = q_id;
            p_j = j;
            p_jj = jj;
            p_q_base = aln_q_seq[k];
        } else {
            break; // when there is a big alignment gap > UINT8_MAX, stop to extned the tagging string
        }
    }
    // sentinal at the end
    //k = aln_seq_len;
    tags->len = k;
    (tags->align_tags[k]).t_pos = UINT_MAX;
    (tags->align_tags[k]).delta = UINT8_MAX;
    (tags->align_tags[k]).q_base = '.';
    (tags->align_tags[k]).q_id = UINT_MAX;
    return tags;
}

void free_align_tags( align_tags_t * tags) {
    free( tags->align_tags );
    free( tags );
}


void allocate_aln_col( align_tag_col_t * col) {
    col->p_t_pos = ( seq_coor_t * ) calloc(col->size, sizeof( seq_coor_t ));
    col->p_delta = ( uint8_t * ) calloc(col->size, sizeof( uint8_t ));
    col->p_q_base = ( char * )calloc(col->size, sizeof( char ));
    col->link_count = ( uint16_t * ) calloc(col->size, sizeof( uint16_t ));
}

void realloc_aln_col( align_tag_col_t * col ) {
    col->p_t_pos = (seq_coor_t *) realloc( col->p_t_pos, (col->size) * sizeof( seq_coor_t ));
    col->p_delta = ( uint8_t *)  realloc( col->p_delta, (col->size) * sizeof( uint8_t ));
    col->p_q_base = (char *) realloc( col->p_q_base, (col->size) * sizeof( char ));
    col->link_count = ( uint16_t *) realloc( col->link_count, (col->size) * sizeof( uint16_t ));
}

void free_aln_col( align_tag_col_t * col) {
    free(col->p_t_pos);
    free(col->p_delta);
    free(col->p_q_base);
    free(col->link_count);
}


void allocate_delta_group( msa_delta_group_t * g) {
    int i,j;
    g->max_delta = 0;
    g->delta = (msa_base_group_t *) calloc( g->size, sizeof(msa_base_group_t));
    for (i = 0; i< g->size; i++) {
        g->delta[i].base = ( align_tag_col_t * ) calloc( 5, sizeof(align_tag_col_t ) );
        for (j = 0; j < 5; j++ ) {
             g->delta[i].base[j].size = 8;
             allocate_aln_col(&(g->delta[i].base[j]));
        }
    }
}

void realloc_delta_group( msa_delta_group_t * g, uint16_t new_size ) {
    int i, j, bs, es;
    bs = g->size;
    es = new_size;
    g->delta = (msa_base_group_t *) realloc(g->delta, new_size * sizeof(msa_base_group_t));
    for (i=bs; i < es; i++) {
        g->delta[i].base = ( align_tag_col_t *) calloc( 5, sizeof(align_tag_col_t ) );
        for (j = 0; j < 5; j++ ) {
             g->delta[i].base[j].size = 8;
             allocate_aln_col(&(g->delta[i].base[j]));
        }
    }
    g->size = new_size;
}

void free_delta_group( msa_delta_group_t * g) {
    //manything to do here
    int i, j;
    for (i = 0; i < g->size; i++) {
        for (j = 0; j < 5; j++) {
            free_aln_col( &(g->delta[i].base[j]) );
        }
        free(g->delta[i].base);
    }
    free(g->delta);
}

void update_col( align_tag_col_t * col, seq_coor_t p_t_pos, uint8_t p_delta, char p_q_base) {
    int updated = 0;
    int link;
    col->count += 1;
    for (link = 0; link < col->n_link; link++) {
        if ( p_t_pos == col->p_t_pos[link] &&
             p_delta == col->p_delta[link] &&
             p_q_base == col->p_q_base[link] ) {
            col->link_count[link] ++;
            updated = 1;
            break;
        }
    }
    if (updated == 0) {
        if (col->n_link + 1 > col->size) {
            if (col->size < (UINT16_MAX >> 1)-1) {
                col->size *= 2;
            } else {
                col->size += 256;
            }
            assert( col->size < UINT16_MAX-1 );
            realloc_aln_col(col);
        }
        link = col->n_link;

        col->p_t_pos[link] = p_t_pos;
        col->p_delta[link] = p_delta;
        col->p_q_base[link] = p_q_base;
        col->link_count[link] = 1;
        col->n_link++;
    }
}


msa_pos_t * get_msa_working_sapce(unsigned int max_t_len) {
    msa_pos_t * msa_array;
    unsigned int i;
    msa_array = calloc(max_t_len, sizeof(msa_pos_t *));
    for (i = 0; i < max_t_len; i++) {
        msa_array[i] = calloc(1, sizeof(msa_delta_group_t));
        msa_array[i]->size = 16;
        allocate_delta_group(msa_array[i]);
    }
    return msa_array;
}

void clean_msa_working_space( msa_pos_t * msa_array, unsigned int max_t_len) {
    unsigned int i,k,j;
    align_tag_col_t * col;
    for (i = 0; i < max_t_len; i++) {
        for (j =0; j < msa_array[i]->max_delta + 1; j++) {
            for (k = 0; k < 5; k++ ) {
                col = msa_array[i]->delta[j].base + k;
                /*
                for (c =0; c < col->size; c++) {
                    col->p_t_pos[c] = 0;
                    col->p_delta[c] = 0;
                    col->p_q_base[c] = 0;
                    col->link_count[c] =0;
                }
                */
                col->n_link = 0;
                col->count = 0;
                col->best_p_t_pos = 0;
                col->best_p_delta = 0;
                col->best_p_q_base = 0;
                col->score = 0;
            }
        }
        msa_array[i]->max_delta = 0;
    }
}

#define STATIC_ALLOCATE
//#undef STATIC_ALLOCATE

consensus_data * get_cns_from_align_tags( align_tags_t ** tag_seqs,
                                          unsigned n_tag_seqs,
                                          unsigned t_len,
                                          unsigned min_cov ) {

    seq_coor_t i, j;
    seq_coor_t t_pos = 0;
    unsigned int * coverage;
    unsigned int * local_nbase;

    consensus_data * consensus;
    //char * consensus;
    align_tag_t * c_tag;

    coverage = calloc( t_len, sizeof(unsigned int) );
    local_nbase = calloc( t_len, sizeof(unsigned int) );

#ifndef STATIC_ALLOCATE

    msa_pos_t * msa_array = NULL; // For more efficiency, this should be injected.
    msa_array = calloc(t_len, sizeof(msa_pos_t *));

    for (i = 0; i < t_len; i++) {
        msa_array[i] = calloc(1, sizeof(msa_delta_group_t));
        msa_array[i]->size = 8;
        allocate_delta_group(msa_array[i]);
    }

#else

    static msa_pos_t * msa_array = NULL;
    if ( msa_array == NULL) {
        msa_array = get_msa_working_sapce( 100000 );
    }

    assert(t_len < 100000);

#endif


    // loop through every alignment
    //printf("XX %d\n", n_tag_seqs);
    for (unsigned ii = 0; ii < n_tag_seqs; ii++) {

        // for each alignment position, insert the alignment tag to msa_array
        for (int jj = 0; jj < tag_seqs[ii]->len; jj++) {
            c_tag = tag_seqs[ii]->align_tags + jj;
            unsigned int delta;
            delta = c_tag->delta;
            if (delta == 0) {
                t_pos = c_tag->t_pos;
                coverage[ t_pos ] ++;
            }
            // Assume t_pos was set on earlier iteration.
            // (Otherwise, use its initial value, which might be an error. ~cd)
            if (delta > msa_array[t_pos]->max_delta) {
                msa_array[t_pos]->max_delta = delta;
                if (msa_array[t_pos]->max_delta + 4 > msa_array[t_pos]->size ) {
                    realloc_delta_group(msa_array[t_pos], msa_array[t_pos]->max_delta + 16);
                }
            }

            unsigned int base = -1;
            switch (c_tag->q_base) {
                case 'A': base = 0; break;
                case 'C': base = 1; break;
                case 'G': base = 2; break;
                case 'T': base = 3; break;
                case '-': base = 4; break;
            }
            // Note: On bad input, base may be -1.
            update_col( &(msa_array[t_pos]->delta[delta].base[base]), 
					c_tag->p_t_pos, c_tag->p_delta, c_tag->p_q_base);
            local_nbase[ t_pos ] ++;
        }
    }

    // propogate score throught the alignment links, setup backtracking information
    align_tag_col_t * g_best_aln_col = 0;
    unsigned int g_best_k = 0;
    seq_coor_t g_best_t_pos = 0;
    {
        int k;
		int best_k;
        double score;
        double best_score;
        double g_best_score;

        align_tag_col_t * aln_col;

        g_best_score = -1;

		best_k = 0;
        for (i = 0; (unsigned) i < t_len; i++) {  //loop through every template base
            for (j = 0; j <= msa_array[i]->max_delta; j++) { // loop through every delta position
                for (k = 0; k < 5; k++) {  // loop through diff bases of the same delta posiiton
                    aln_col = msa_array[i]->delta[j].base + k;
                    if (aln_col->count >= 0) {
                        best_score = -1;

                        for (int link = 0; link < aln_col->n_link; link++) { // loop through differnt link to previous column
                            int pi;
                            int pj;
                            int pk;
                            pi = aln_col->p_t_pos[link];
                            pj = aln_col->p_delta[link];
                            switch (aln_col->p_q_base[link]) {
                                case 'A': pk = 0; break;
                                case 'C': pk = 1; break;
                                case 'G': pk = 2; break;
                                case 'T': pk = 3; break;
                                case '-': pk = 4; break;
                                default: pk = 4;
                            }

                            if (aln_col->p_t_pos[link] == -1) {
                                score =  (double) aln_col->link_count[link] - (double) coverage[i] * 0.5;
                            } else {
								//printf("XXX %d %d\n", pj, msa_array[pi]->size);
								assert( pj < msa_array[pi]->size);
                                score = msa_array[pi]->delta[pj].base[pk].score +
                                        (double) aln_col->link_count[link] - (double) coverage[i] * 0.5;
                            }
                            if (score > best_score) {
                                best_score = score;
                                aln_col->best_p_t_pos = pi;
                                aln_col->best_p_delta = pj;
                                aln_col->best_p_q_base = pk;
								best_k = k;
                            }
                        }
                        aln_col->score = best_score;
                        if (best_score > g_best_score) {
                            g_best_score = best_score;
                            g_best_aln_col = aln_col;
                            g_best_t_pos = i;
							g_best_k = best_k;
                        }
                    }
                }
            }
        }
        assert(g_best_score > 0);
    }

    // reconstruct the sequences
    unsigned int index;
    char bb = '$';
    int k;
    char * cns_str;
    int * eqv;
    double score0;

    consensus = calloc( 1, sizeof(consensus_data) );
    consensus->sequence = calloc( t_len * 2 + 1, sizeof(char) );
    consensus->eqv = calloc( t_len * 2 + 1, sizeof(unsigned int) );
    cns_str = consensus->sequence;
    eqv =  consensus->eqv;
    
    index = 0;
    k = g_best_k;
    i = g_best_t_pos;

    while (1) {
        if (coverage[i] > min_cov) {
            switch (k) {
                case 0: bb = 'A'; break;
                case 1: bb = 'C'; break;
                case 2: bb = 'G'; break;
                case 3: bb = 'T'; break;
                case 4: bb = '-'; break;
            }
        } else {
            switch (k) {
                case 0: bb = 'a'; break;
                case 1: bb = 'c'; break;
                case 2: bb = 'g'; break;
                case 3: bb = 't'; break;
                case 4: bb = '-'; break;
            }
        }
        // Note: On bad input, bb will keep previous value, possibly '$'.

        score0 = g_best_aln_col->score;
        i = g_best_aln_col->best_p_t_pos;
        if (i == 0 || index >= t_len * 2) break;
        j = g_best_aln_col->best_p_delta;
        k = g_best_aln_col->best_p_q_base;
        g_best_aln_col = msa_array[i]->delta[j].base + k;

        if (bb != '-') {
            cns_str[index] = bb;
            eqv[index] = (int) score0 - (int) g_best_aln_col->score;
            //printf("C %d %d %c %lf %d %d\n", i, index, bb, g_best_aln_col->score, coverage[i], eqv[index] );
            index ++;
        }
    }

    // reverse the sequence
    for (i = 0; (unsigned) i < index/2; i++) {
        cns_str[i] = cns_str[i] ^ cns_str[index-i-1];
        cns_str[index-i-1] = cns_str[i] ^ cns_str[index-i-1];
        cns_str[i] = cns_str[i] ^ cns_str[index-i-1];
        eqv[i] = eqv[i] ^ eqv[index-i-1];
        eqv[index-i-1] = eqv[i] ^ eqv[index-i-1];
        eqv[i] = eqv[i] ^ eqv[index-i-1];
    }

    cns_str[index] = 0;
    //printf("%s\n", cns_str);
#ifndef STATIC_ALLOCATE
    for (i = 0; i < t_len; i++) {
        free_delta_group(msa_array[i]);
        free(msa_array[i]);
    }

    free(msa_array);
#else
    clean_msa_working_space(msa_array, t_len+1);
#endif

    free(coverage);
    free(local_nbase);
    return consensus;
}


void free_consensus_data( consensus_data * consensus ){
    free(consensus->sequence);
    free(consensus->eqv);
    free(consensus);
}
