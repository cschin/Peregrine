/*
 * =====================================================================================
 *
 *       Filename:  falconX.c
 *
 *    Description:
 *
 *        Version:  0.1
 *        Created:  04/13/2019 07:00:00 
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
#include <ctype.h>
#include "common.h"
#include "falcon.h"


align_tags_t * get_align_tags( 
        char * aln_q_seq,
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
            } else {
                continue;
            }
        }
        if (aln_t_seq[k] != '-') {
            break;
        }
    }
    for (; k < aln_seq_len; k++) {
        if (aln_q_seq[k] != '-') {
            i ++;
            if (jj < 12) {  // we cap the biggest insert gap up to 12 bases
                jj ++;
            } else {
                continue;
            }
        }
        if (aln_t_seq[k] != '-') {
            j ++;
            jj = 0;
        }
        // fprintf(stderr, "t %d %d %d %c %c\n", q_id, j, jj, aln_t_seq[k], aln_q_seq[k]);

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

uint64_t get_tag_key(
        seq_coor_t t_pos, 
        uint8_t delta, 
        char base) {
    assert(t_pos > 0);
    return (((uint64_t) t_pos) << 32) | (delta << 8) | base;   
}


int tag_key_comp(const void * a, const void * b) {
    uint64_t t1, t1p;
    uint64_t t2, t2p;
    align_edge_t * e1 = (align_edge_t *) a;
    align_edge_t * e2 = (align_edge_t *) b;
    t1  = get_tag_key(e1->t_pos, e1->delta, e1->q_base);
    t1p = get_tag_key(e1->p_t_pos, e1->p_delta, e1->p_q_base);
    t2  = get_tag_key(e2->t_pos, e2->delta, e2->q_base);
    t2p = get_tag_key(e2->p_t_pos, e2->p_delta, e2->p_q_base);

    return (t1 == t2) ? (t1p > t2p) : (t1 > t2); 
}


void sort_accumulate_tags(
        align_edge_v * edges,
        align_edge_v * input_edges) {

    align_tag_t c_tag =  {0, 0, '.', 0, 0, '.'};
    align_edge_t align_edge;

    qsort(input_edges->a, input_edges->n, sizeof(align_edge_t), tag_key_comp);

    for (size_t i = 0; i < input_edges->n; i++){
        align_edge_t e;
        e = input_edges->a[i];

        // fprintf(stderr,"E0 %d %d %c %d %d %c %d\n", 
        // e.t_pos, e.delta, e.q_base, e.p_t_pos, 
        // e.p_delta, e.p_q_base, e.count); 

        if (c_tag.t_pos != e.t_pos || 
                c_tag.delta != e.delta || 
                c_tag.q_base != e.q_base || 
                c_tag.p_t_pos != e.p_t_pos || 
                c_tag.p_delta != e.p_delta ||
                c_tag.p_q_base != e.p_q_base) {

            align_edge.t_pos = e.t_pos;
            align_edge.delta = e.delta;
            align_edge.q_base = e.q_base;
            align_edge.p_t_pos = e.p_t_pos;
            align_edge.p_delta = e.p_delta;
            align_edge.p_q_base = e.p_q_base;
            align_edge.count = 1;
            align_edge.score = 0.0;
            kv_push(align_edge_t, NULL, *edges, align_edge); 
        } else {
            edges->a[edges->n-1].count ++;
        }
        e = edges->a[edges->n-1];

        // fprintf(stderr, "E1 %d %d %c %d %d %c %d\n", 
        // e.t_pos, e.delta, e.q_base, 
        // e.p_t_pos, e.p_delta, e.p_q_base, e.count); 

        c_tag.t_pos = e.t_pos;
        c_tag.delta = e.delta;
        c_tag.q_base = e.q_base;
        c_tag.p_t_pos = e.p_t_pos;
        c_tag.p_delta = e.p_delta;
        c_tag.p_q_base = e.p_q_base;
    }
}

uint64_t get_node_score(
        align_node_map_t * node_map,
        align_edge_v * edges, 
        uint16_t * coverage) {

    uint64_t node_key, p_node_key;
    uint64_t best_node_key=0;
    double global_best_score=0;
    int absent;
    khiter_t k;
    align_node_t * node;
    align_node_t * p_node;

    for (size_t i = 0; i < edges->n; i++){
        align_edge_t * e;
        e = edges->a+i;
        e->score = (double) e->count - 0.5 * ((double) coverage[e->t_pos]+1);
        /* 
           fprintf(stderr,"E2 %d %d %c %d %d %c %d %0.3f\n", 
           e->t_pos, e->delta, e->q_base, 
           e->p_t_pos, e->p_delta, e->p_q_base, e->count, e->score); 
           */
        node_key = get_tag_key(e->t_pos, e->delta, e->q_base);

        k = kh_get(NODE, node_map, node_key);
        if (k == kh_end(node_map)) {
            k = kh_put(NODE, node_map, node_key, &absent);
            node = malloc(sizeof(align_node_t));
            node->t_pos = e->t_pos;
            node->delta = e->delta;
            node->q_base = e->q_base;
            node->best_edge = e;
            node->best_score = e->score;
            kh_val(node_map, k) = node;
            // fprintf(stderr, "XX %lu %d %d %c\n", node_key, e.t_pos, e.delta, e.q_base);

        } else {
            node = kh_val(node_map, k);
        }

        if (e->p_q_base == '.') {
            continue;
        }

        p_node_key = get_tag_key(e->p_t_pos, e->p_delta, e->p_q_base);
        k = kh_get(NODE, node_map, p_node_key);
        if (k == kh_end(node_map)) {
            continue;
        };
        p_node = kh_val(node_map, k);
        double new_score;
        new_score = e->score + p_node->best_score;
        if (new_score > node->best_score) {
            node->best_score = new_score;
            node->best_edge = e;
            /*
               fprintf(stderr, "N %d %d %c %d %d %c %0.3f\n", 
               node->t_pos, node->delta, node->q_base, 
               p_node->t_pos, p_node->delta, p_node->q_base, node->best_score); 
               */
            if (new_score > global_best_score) {
                global_best_score = new_score;
                best_node_key = node_key;
            }
        }
    }
    return best_node_key;
}

consensus_data * backtracking(
        align_node_map_t * node_map,
        uint64_t best_node_key,
        align_edge_v *edges, 
        uint16_t * coverage,
        unsigned min_cov) {

    khiter_t k;
    align_node_t * node;
    uint64_t node_key;
    consensus_data * consensus;
    char * cns;
    // uint8_t * eqv;
    uint32_t idx = 0;

    consensus = (consensus_data *) calloc(1, sizeof(consensus_data));
    consensus->sequence = calloc( edges->n, sizeof(char) );
    // consensus->eqv = calloc(1, sizeof(uint8_t)); //not used, reserve for future QV output
    cns = consensus->sequence; // just an alias
    // eqv = consensus->eqv;  // just an alias

    node_key = best_node_key;
    k = kh_get(NODE, node_map, node_key);
    assert( k != kh_end(node_map) );
    node = kh_val(node_map, k);
    idx = 0;
    for (;;) {
        if (node->q_base != '-') {
            if (coverage[node->t_pos] > min_cov) {
                cns[idx] = node->q_base;
            } else {
                cns[idx] = tolower(node->q_base);
            }
            idx++;
        }

        align_edge_t * e;
        e = node->best_edge;
        /*
           fprintf(stderr, "N2 %d %d %c %d %d %c %0.3f\n", 
           e->t_pos, e->delta, e->q_base, 
           e->p_t_pos, e->p_delta, e->p_q_base, node->best_score); 
           */
        if (node->best_edge == NULL || e->p_q_base == '.') break;
        node_key = get_tag_key( e->p_t_pos, e->p_delta, e->p_q_base );
        k = kh_get(NODE, node_map, node_key);
        node = kh_val(node_map, k);
    }

    // reverse the sequence
    for (uint32_t i = 0; i < idx/2; i++) {
        cns[i] = cns[i] ^ cns[idx-i-1];
        cns[idx-i-1] = cns[i] ^ cns[idx-i-1];
        cns[i] = cns[i] ^ cns[idx-i-1];
    }

    cns[idx] = 0; //terminate the string
    return consensus;
}

consensus_data * get_cns_from_align_tags( 
        align_tags_t ** tag_seqs,
        unsigned n_tag_seqs,
        unsigned t_len,
        unsigned min_cov ) {

    seq_coor_t t_pos = 0;
    uint16_t * coverage;

    consensus_data * consensus;
    align_tag_t * c_tag;
    align_edge_t align_edge;
    align_edge_v all_edges = {0, 0, 0};  // initial k_vec with {0, 0, 0}
align_edge_v edges = {0, 0, 0};
align_node_map_t * node_map = kh_init(NODE);

coverage = calloc( t_len, sizeof(unsigned int) );

// loop through every alignment
for (unsigned ii = 0; ii < n_tag_seqs; ii++) {
    // for each alignment position, insert the alignment tag to the table
    for (int jj = 0; jj < tag_seqs[ii]->len; jj++) {
        c_tag = tag_seqs[ii]->align_tags + jj; 
        align_edge.t_pos = c_tag->t_pos;
        align_edge.delta = c_tag->delta;
        align_edge.q_base = c_tag->q_base;
        align_edge.p_t_pos = c_tag->p_t_pos;
        align_edge.p_delta = c_tag->p_delta;
        align_edge.p_q_base = c_tag->p_q_base;
        align_edge.coverage = 1;
        align_edge.count = 1;
        align_edge.score = 0.0;
        kv_push(align_edge_t, NULL, all_edges, align_edge);
        if (c_tag->delta == 0) {
            t_pos = c_tag->t_pos;
            coverage[t_pos] ++;
        }

    }
}

sort_accumulate_tags(&edges, &all_edges);
uint64_t best_node_key;
best_node_key = get_node_score(node_map, &edges, coverage);
consensus = backtracking(node_map, best_node_key, &edges, coverage ,min_cov);

// need to clean the nodes here
for (khiter_t __i = kh_begin(node_map); __i != kh_end(node_map); ++__i) {
    if (!kh_exist(node_map,__i)) continue;
    free(kh_val(node_map, __i));
}
kh_destroy(NODE, node_map);
kv_destroy(edges);
kv_destroy(all_edges);
free(coverage);
return consensus;
}


void free_consensus_data( consensus_data * consensus ){
    free(consensus->sequence);
    // free(consensus->eqv);
    free(consensus);
}

