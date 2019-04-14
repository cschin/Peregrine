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
    
    for (k=0; k < aln_seq_len; k++) {
        if (aln_q_seq[k] != '-') {
            i ++;
            jj ++;
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


int uint64_comp(const void * a, const void * b) {
    return  * ((uint64_t *) a) > * ((uint64_t *) b); 
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
        seq_coor_t t_pos = (seq_coor_t) ((e->ctag_key >> 32) & 0xFFFFFFFF);
        e->score = (double) e->count - 0.5 * ((double) coverage[t_pos]-1);
        /* 
           fprintf(stderr,"E2 %d %d %c %d %d %c %d %0.3f\n", 
           e->t_pos, e->delta, e->q_base, 
           e->p_t_pos, e->p_delta, e->p_q_base, e->count, e->score); 
           */
        node_key = e->ctag_key;

        k = kh_get(NODE, node_map, node_key);
        if (k == kh_end(node_map)) {
            k = kh_put(NODE, node_map, node_key, &absent);
            node = malloc(sizeof(align_node_t));
            node->ctag_key = node_key;
            node->best_edge = e;
            node->best_score = e->score;
            kh_val(node_map, k) = node;
            // fprintf(stderr, "XX %lu %d %d %c\n", node_key, e.t_pos, e.delta, e.q_base);

        } else {
            node = kh_val(node_map, k);
        }

        if ((char) (e->ptag_key & 0xFF) == '.') {
            continue;
        }

        p_node_key = e->ptag_key;
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
               fprintf(stderr, "N0 %d %d %c %d %d %c %0.3f\n", 
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
    char q_base;
    seq_coor_t t_pos;

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
        t_pos = (seq_coor_t) ((node->ctag_key >> 32) & 0xFFFFFFFF);
        q_base = (char)( node->ctag_key & 0xFF);
        if (q_base != '-') {
            if (coverage[t_pos] > min_cov) {
                cns[idx] = q_base;
            } else {
                cns[idx] = tolower(q_base);
            }
            idx++;
        } 

        align_edge_t * e;
        e = node->best_edge;

        char p_q_base = (char) ( e->ptag_key & 0xFF);
        if (node->best_edge == NULL || p_q_base == '.') break;
        /*
           fprintf(stderr, "N2 %d %d %c %d %d %c %0.3f\n", 
           e->t_pos, e->delta, e->q_base, 
           e->p_t_pos, e->p_delta, e->p_q_base, node->best_score); 
           */
        node_key = e->ptag_key; 
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

    align_edge_v edges = {0, 0, 0}; // initial kvec with {0, 0, 0}
    align_node_map_t * node_map = kh_init(NODE);
    ctag_to_ptag_t * ctag_to_ptag = kh_init(CTAG); 
    ptag_to_count_t * ptag_to_count;
    uint64_t ctag_key, ptag_key;

    khiter_t k1, k2;
    int absent;

    coverage = calloc( t_len, sizeof(unsigned int) );

    // loop through every alignment
    for (unsigned ii = 0; ii < n_tag_seqs; ii++) {
        // for each alignment position, insert the alignment tag to the table
        int flag = 0;
        for (int jj = 0; jj < tag_seqs[ii]->len; jj++) {
            c_tag = tag_seqs[ii]->align_tags + jj; 
            if (flag == 0 && c_tag->p_q_base == '-') {
                continue;
            } else {
                flag = 1;
            }

            ctag_key = get_tag_key(c_tag->t_pos,  c_tag->delta,  c_tag->q_base);
            ptag_key = get_tag_key(c_tag->p_t_pos, c_tag->p_delta, c_tag->p_q_base);
            k1 = kh_get(CTAG, ctag_to_ptag, ctag_key);
            if (k1 == kh_end(ctag_to_ptag)) { 
                k1 = kh_put(CTAG, ctag_to_ptag, ctag_key, &absent);
                ptag_to_count = kh_init(PTAG);
                kh_val(ctag_to_ptag, k1) = ptag_to_count;
                k2 = kh_put(PTAG, ptag_to_count, ptag_key, &absent);
                kh_val(ptag_to_count, k2) = 1;
            } else {
                ptag_to_count = kh_val(ctag_to_ptag, k1);
                k2 = kh_get(PTAG, ptag_to_count, ptag_key);
                if (k2 == kh_end(ptag_to_count)) {
                    k2 = kh_put(PTAG, ptag_to_count, ptag_key, &absent);
                    kh_val(ptag_to_count, k2) = 1;
                } else {
                    kh_val(ptag_to_count, k2) += 1;
                }
            }

            if (c_tag->delta == 0) {
                t_pos = c_tag->t_pos;
                coverage[t_pos] ++;
            }

        }
    }

    uint64_v ctag_keys = {0, 0, 0};
    for (k1 = kh_begin(ctag_to_ptag); k1 != kh_end(ctag_to_ptag); k1++) {
        if (!kh_exist(ctag_to_ptag, k1)) continue;
        ctag_key = kh_key(ctag_to_ptag, k1);
        kv_push(uint64_t, NULL, ctag_keys, ctag_key); 
    }

    qsort(ctag_keys.a, ctag_keys.n, sizeof(uint64_t), uint64_comp);
    for (size_t i = 0; i < ctag_keys.n; i++) {
        //fprintf(stderr, "test %lu\n", ctag_keys.a[i]);
        k1 = kh_get(CTAG, ctag_to_ptag, ctag_keys.a[i]);
        ptag_to_count = kh_val(ctag_to_ptag, k1);

        uint64_v ptag_keys = {0, 0, 0};
        for (k2 = kh_begin(ptag_to_count); k2 != kh_end(ptag_to_count); k2++) {
            if (!kh_exist(ptag_to_count, k2)) continue;
            ptag_key = kh_key(ptag_to_count, k2);
            kv_push(uint64_t, NULL, ptag_keys, ptag_key); 
        }
        qsort(ptag_keys.a, ptag_keys.n, sizeof(uint64_t), uint64_comp);
        for (size_t j = 0; j < ptag_keys.n; j++) {
            k2 = kh_get(PTAG, ptag_to_count, ptag_keys.a[j]);
            uint32_t count = kh_val(ptag_to_count, k2);
            // fprintf(stderr, "test %lu %lu %u\n", ctag_keys.a[i], ptag_keys.a[j], count);
            ctag_key = ctag_keys.a[i];
            ptag_key = ptag_keys.a[j];
            align_edge_t align_edge;
            align_edge.ctag_key = ctag_key;
            align_edge.ptag_key = ptag_key;
            align_edge.count = count;
            align_edge.score = 0.0;
            kv_push(align_edge_t, NULL, edges, align_edge); 
        }
        kv_destroy(ptag_keys);
    }

    uint64_t best_node_key;
    best_node_key = get_node_score(node_map, &edges, coverage);
    consensus = backtracking(node_map, best_node_key, &edges, coverage ,min_cov);

    for (size_t i = 0; i < ctag_keys.n; i++) {
        //fprintf(stderr, "test %lu\n", ctag_keys.a[i]);
        k1 = kh_get(CTAG, ctag_to_ptag, ctag_keys.a[i]);
        ptag_to_count = kh_val(ctag_to_ptag, k1);
        kh_destroy(PTAG, ptag_to_count);
    }
    kh_destroy(CTAG, ctag_to_ptag);
    kv_destroy(ctag_keys);
    // need to clean the nodes here
    for (khiter_t __i = kh_begin(node_map); __i != kh_end(node_map); ++__i) {
        if (!kh_exist(node_map,__i)) continue;
        free(kh_val(node_map, __i));
    }
    kh_destroy(NODE, node_map);
    kv_destroy(edges);
    free(coverage);
    return consensus;
}

void free_consensus_data( consensus_data * consensus ){
    free(consensus->sequence);
    // free(consensus->eqv);
    free(consensus);
}

