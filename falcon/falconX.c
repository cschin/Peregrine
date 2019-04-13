/*
 * =====================================================================================
 *
 *       Filename:  falconX.c
 *
 *    Description:
 *
 *        Version:  0.1
 *        Created:  04/12/2019 07:00:00 
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
#include "falconX.h"


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
        // printf("t %d %d %d %c %c\n", q_id, j, jj, aln_t_seq[k], aln_q_seq[k]);

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

uint64_t get_tag_key(seq_coor_t t_pos, 
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



void sort_accumulate_tags(align_edge_v *edges,
		align_edge_map_t * edge_map,
		align_node_v * nodes,
		align_node_map_t * node_map,
		align_edge_v * input_edges) {
	
	align_tag_t c_tag =  {0, 0, '.', 0, 0, '.'};
	align_edge_t align_edge;

	qsort(input_edges->a, input_edges->n, sizeof(align_edge_t), tag_key_comp);
	
	for (size_t i = 0; i < input_edges->n; i++){
		align_edge_t e;
		e = input_edges->a[i];
		// printf("E0 %d %d %c %d %d %c %d\n", e.t_pos, e.delta, e.q_base, e.p_t_pos, e.p_delta, e.p_q_base, e.count); 

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
		//printf("E1 %d %d %c %d %d %c %d\n", e.t_pos, e.delta, e.q_base, e.p_t_pos, e.p_delta, e.p_q_base, e.count); 
	   
		c_tag.t_pos = e.t_pos;
		c_tag.delta = e.delta;
		c_tag.q_base = e.q_base;
		c_tag.p_t_pos = e.p_t_pos;
		c_tag.p_delta = e.p_delta;
		c_tag.p_q_base = e.p_q_base;
	}
	for (size_t i = 0; i < edges->n; i++){
		align_edge_t e;
		e = edges->a[i];
		printf("E2 %d %d %c %d %d %c %d\n", e.t_pos, e.delta, e.q_base, e.p_t_pos, e.p_delta, e.p_q_base, e.count); 
	}
}


consensus_data * get_cns_from_align_tags( align_tags_t ** tag_seqs,
        unsigned n_tag_seqs,
        unsigned t_len,
        unsigned min_cov ) {

    seq_coor_t i, j;
    seq_coor_t t_pos = 0;
	uint32_t total_aln_tags = 0;
    unsigned int * coverage;
    unsigned int * local_nbase;

    consensus_data * consensus;
    align_tag_t * c_tag;
	align_edge_t align_edge;
	align_edge_v all_edges = {0, 0, 0};  // initial k_vec with {0, 0, 0}

	align_edge_v edges = {0, 0, 0};
	align_edge_map_t * edge_map = kh_init(EDGE);
    align_node_v nodes = {0, 0, 0};
	align_node_map_t * node_map = kh_init(NODE);


    coverage = calloc( t_len, sizeof(unsigned int) );
    local_nbase = calloc( t_len, sizeof(unsigned int) );

    // loop through every alignment
	printf("n tag %d\n", n_tag_seqs);
    for (unsigned ii = 0; ii < n_tag_seqs; ii++) {
        // for each alignment position, insert the alignment tag to msa_array
        for (int jj = 0; jj < tag_seqs[ii]->len; jj++) {
			c_tag = tag_seqs[ii]->align_tags + jj; 
			align_edge.t_pos = c_tag->t_pos;
			align_edge.delta = c_tag->delta;
			align_edge.q_base = c_tag->q_base;
			align_edge.p_t_pos = c_tag->p_t_pos;
			align_edge.p_delta = c_tag->p_delta;
			align_edge.p_q_base = c_tag->p_q_base;
			printf("C %d %d %c %d %d %c\n", c_tag->t_pos, c_tag->delta, c_tag->q_base, 
					c_tag->p_t_pos, c_tag->p_delta, c_tag->p_q_base); 
			align_edge.coverage = 1;
			align_edge.count = 1;
			align_edge.score = 0.0;
			kv_push(align_edge_t, NULL, all_edges, align_edge);
            if (c_tag->delta == 0) {
                t_pos = c_tag->t_pos;
                coverage[ t_pos ] ++;
            }

		}
	}

    sort_accumulate_tags(&edges, edge_map,
		&nodes, node_map,
		&all_edges);

    free(coverage);
    return consensus;
}


void free_consensus_data( consensus_data * consensus ){
    free(consensus->sequence);
    free(consensus->eqv);
    free(consensus);
}

