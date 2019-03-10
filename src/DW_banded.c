
/*
 * =====================================================================================
 *
 *       Filename:  DW_banded.c
 *
 *    Description:  A banded version for the O(ND) greedy sequence alignment algorithm
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
#include <stdbool.h>
#include "shimmer.h"

alignment_t * align(char * query_seq, seq_coor_t q_len,
		char * target_seq, seq_coor_t t_len,
		seq_coor_t band_tolerance) {
    seq_coor_t * V;
    seq_coor_t * U;  // array of matched bases for each "k"
    seq_coor_t k_offset;
    seq_coor_t d;
    seq_coor_t k, k2;
    seq_coor_t best_m;  // the best "matches" for each d
    seq_coor_t min_k, new_min_k;
    seq_coor_t max_k, new_max_k;
    seq_coor_t x, y;
    seq_coor_t x1, y1;
    seq_coor_t max_d;
    seq_coor_t band_size;
	bool start = false;

    alignment_t * align_rtn;
    bool aligned = false;

    //printf("debug: %ld %ld\n", q_len, t_len);
    //printf("%s\n", query_seq);

    max_d = (int) (0.3 * (q_len + t_len));

    band_size = band_tolerance * 2;

    V = calloc( max_d * 2 + 1, sizeof(seq_coor_t) );
    U = calloc( max_d * 2 + 1, sizeof(seq_coor_t) );

    k_offset = max_d;

    align_rtn = calloc( 1, sizeof(alignment_t));
    align_rtn->astr_size = 0;
    align_rtn->q_bgn = 0;
    align_rtn->q_end = 0;
    align_rtn->t_bgn = 0;
    align_rtn->t_end = 0;

    //printf("max_d: %lu, band_size: %lu\n", max_d, band_size);
    best_m = -1;
    min_k = 0;
    max_k = 0;
    for (d = 0; d < max_d; d ++ ) {
        if (max_k - min_k > band_size) {
            break;
        }

        for (k = min_k; k <= max_k;  k += 2) {

            if ( (k == min_k) || ((k != max_k) && (V[ k - 1 + k_offset ] < V[ k + 1 + k_offset])) ) {
                x = V[ k + 1 + k_offset];
            } else {
                x = V[ k - 1 + k_offset] + 1;
            }
            y = x - k;
            x1 = x;
			y1 = y;

            while ( x < q_len && y < t_len && query_seq[x] == target_seq[y] ){
                x++;
                y++;
            }

			if ( (x - x1 > 16) && (start == false) ) {
				align_rtn->q_bgn = x1;
				align_rtn->t_bgn = y1;
				start = true;
			}

            V[ k + k_offset ] = x;
            U[ k + k_offset ] = x + y;

            if ( x + y > best_m) {
                best_m = x + y;
            }

            if ( x >= q_len || y >= t_len) {
                aligned = true;
                break;
            }
        }

        // For banding
        new_min_k = max_k;
        new_max_k = min_k;

        for (k2 = min_k; k2 <= max_k;  k2 += 2) {
            if (U[ k2 + k_offset] >= best_m - band_tolerance ) {
                if ( k2 < new_min_k ) {
                    new_min_k = k2;
                }
                if ( k2 > new_max_k ) {
                    new_max_k = k2;
                }
            }
        }

        max_k = new_max_k + 1;
        min_k = new_min_k - 1;

        if (aligned == true) {
            align_rtn->q_end = x;
            align_rtn->t_end = y;
            align_rtn->dist = d;
            align_rtn->astr_size = (align_rtn->q_end - align_rtn->q_bgn + align_rtn->t_end - align_rtn->t_bgn + 2*d) / 2;
            break;
        } 
    }
	if (aligned == false) {
		align_rtn->q_bgn = 0;
		align_rtn->t_bgn = 0;
	}

    free(V);
    free(U);
    return align_rtn;
}

void free_alignment(alignment_t * aln) {
    free(aln);
}
