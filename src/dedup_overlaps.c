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

	char rid0_s[24];
	char rid1_s[24];
	char score[24];
	char acc[24];
	char a_strand[24];
	char a_bgn[24];
	char a_end[24];
	char a_len[24];
	char b_strand[24];
	char b_bgn[24];
	char b_end[24];
	char b_len[24];
	char ovlp_type[24];
	uint32_t rid0;
	uint32_t rid1;
        uint64_t ridp;
	int32_t absent;
	khiter_t k;

	khash_t(RPAIR) * rid_pairs = kh_init(RPAIR);
	while ( fscanf(stdin, "%s %s %s %s %s %s %s %s %s %s %s %s %s\n",
			rid0_s, rid1_s, score, acc,
			a_strand, a_bgn, a_end, a_len,
			b_strand, b_bgn, b_end, b_len,
			ovlp_type) != EOF) {
		rid0 = atoi(rid0_s);
		rid1 = atoi(rid1_s);

		ridp = rid0 < rid1 ? ( ((uint64_t) rid0) << 32) | ((uint64_t) rid1) : (((uint64_t) rid1) << 32) | ((uint64_t) rid0);
		k = kh_get(RPAIR, rid_pairs, ridp);
		if (k == kh_end(rid_pairs)) {
			fprintf(stdout,  "%s %s %s %s %s %s %s %s %s %s %s %s %s\n",
					rid0_s, rid1_s, score, acc,
					a_strand, a_bgn, a_end, a_len,
					b_strand, b_bgn, b_end, b_len,
					ovlp_type);
			kh_put(RPAIR, rid_pairs, ridp, &absent);
		}
	}
	kh_destroy(RPAIR, rid_pairs);

}
