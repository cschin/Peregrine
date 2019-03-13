#include <errno.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include "shimmer.h"
#include "khash.h"
#include "kvec.h"
#include "kalloc.h"

#define OVERLAP 0
#define CONTAINMENT 1

KHASH_MAP_INIT_INT64(RPAIR, uint8_t);

int main(int argc, char *argv[]) {

//002408115 004118624 -14416 99.6 0 27 14387 15129 1 0 14392 14392 contains	

	uint32_t a_bgn, a_end;
	uint32_t b_bgn, b_end;
	uint32_t rid0;
	uint32_t rid1;
	uint64_t ridp;
	int32_t absent;
	khiter_t k;

	khash_t(RPAIR) * rid_pairs = kh_init(RPAIR);

	while ( !feof(stdin) ) {
	 // while ( fscanf(stdin, "%s %s %s %s %s %s %s %s %s %s %s %s %s\n",
	 //		rid0_s, rid1_s, score, acc,
	 //		a_strand, a_bgn, a_end, a_len,
	 //		b_strand, b_bgn, b_end, b_len,
	 //		ovlp_type) != EOF) {
	    ovlp_t ovlp;
		fread(&ovlp, sizeof(ovlp), 1, stdin);

		rid0 = (uint32_t) (ovlp.y0 >> 32);
		rid1 = (uint32_t) (ovlp.y1 >> 32);

		ridp = rid0 < rid1 ? ( ((uint64_t) rid0) << 32) | ((uint64_t) rid1) : (((uint64_t) rid1) << 32) | ((uint64_t) rid0);
		k = kh_get(RPAIR, rid_pairs, ridp);
		if (k == kh_end(rid_pairs)) {
			uint32_t pos0 = (uint32_t) ((ovlp.y0 & 0xFFFFFFFF) >> 1) + 1;
			uint32_t rlen0 = ovlp.rl0;	
			uint8_t strand0 = ovlp.strand0;

			uint32_t pos1 = (uint32_t) ((ovlp.y1 & 0xFFFFFFFF) >> 1) + 1;
			uint32_t rlen1 = ovlp.rl1;	
			uint8_t strand1 = ovlp.strand1;

			ovlp_match_t match = ovlp.match;
			seq_coor_t q_bgn, q_end, t_bgn, t_end;
			q_bgn = match.q_bgn; q_end = match.q_end; 
			t_bgn = match.t_bgn; t_end = match.t_end;
			q_bgn -= t_bgn;
			t_bgn = 0;
			if (strand0 == ORIGINAL) {
				a_bgn = (seq_coor_t) (pos0-pos1) + q_bgn;
				a_end = (seq_coor_t) (pos0-pos1) + q_end;
				a_bgn = a_bgn < 0 ? 0 : a_bgn;              //this ad-hoc fix, read should be stiched by alignment
				a_end = a_end >= rlen0 ? rlen0 : a_end;
			} else {
				q_bgn -= t_bgn;
				t_bgn = 0;
				a_bgn = (seq_coor_t) rlen0 - (seq_coor_t) (pos0 - pos1) - q_end;
				a_end = (seq_coor_t) rlen0 - (seq_coor_t) (pos0 - pos1) - q_bgn;
				a_bgn = a_bgn < 0 ? 0 : a_bgn;              //this ad-hoc fix
				a_end = a_end >= rlen0 ? rlen0 : a_end;
			}
			if (strand1 == ORIGINAL) {
				b_bgn = t_bgn;
				b_end = t_end;
				b_bgn = b_bgn < 0 ? 0 : b_bgn;              //this ad-hoc fix
				b_end = b_end >= rlen1 ? rlen1 : b_end;
			} else {
				b_bgn = (seq_coor_t) rlen1 - t_end;
				b_end = (seq_coor_t) rlen1 - t_bgn;
				b_bgn = b_bgn < 0 ? 0 : b_bgn;              //this ad-hoc fix
				b_end = b_end >= rlen1 ? rlen1 : b_end;
			}
			double err_est;
			err_est	= 100.0 - 100.0 * (double) (match.dist) / (double) (match.astr_size);
			fprintf(stdout,"%09d %09d %d %0.1f %u %d %d %u %u %d %d %u %s\n",
					rid0, rid1, -(match.astr_size), err_est,
					ORIGINAL, a_bgn, a_end, rlen0,
					(strand0 == ORIGINAL ? strand1 : 1-strand1), b_bgn, b_end, rlen1,
					ovlp.ovlp_type == OVERLAP ? "overlap" : "contains");

			kh_put(RPAIR, rid_pairs, ridp, &absent);
		}
	}
	kh_destroy(RPAIR, rid_pairs);

}
