#include <errno.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include "shimmer.h"
#include "khash.h"
#include "kvec.h"


void write_mmlist(char *fn, mm128_v *p) {
	FILE *out_data;
	out_data = fopen(fn, "wb");
	if (!out_data) {
		fprintf(stderr, "file '%s' open error: %s\n", fn, strerror(errno));
		exit(1);
	}
    fwrite(&(kv_size(*p)) , sizeof(size_t), 1, out_data);
    fwrite(p->a, sizeof(mm128_t), kv_size(*p), out_data);
    fclose(out_data);
};


mm128_v read_mmlist(char *fn) {
	mm128_v p = {0, 0, 0};
	FILE *in_data;
	in_data = fopen(fn, "rb");
	if (!in_data) {
		fprintf(stderr, "file '%s' open error: %s\n", fn, strerror(errno));
		exit(1);
	}
    fread(&(kv_size(p)) , sizeof(size_t), 1, in_data);
	kv_resize(mm128_t, NULL, p, kv_size(p));
    fread(p.a, sizeof(mm128_t), kv_size(p), in_data);
    fclose(in_data);
	return p;
};

void mm_count(mm128_v * p, khash_t(MMC) *mcmap, mm_count_v * cp) {
	uint32_t idx;
	khiter_t k;
	mm128_t mmer;
	uint64_t mhash;
	int32_t absent;

	for (idx = 0; idx < p->n; idx++) {
		mmer = p->a[idx];
		mhash = mmer.x >> 8;
		k = kh_put(MMC, mcmap, mhash, &absent);
		if (absent) {
			kh_value(mcmap, k) = 1;
		} else {
			kh_value(mcmap, k) += 1;
		}
	}
	mm_count_to_vec(mcmap, cp);
}

void mm_count_to_vec(khash_t(MMC) *mcmap, mm_count_v * cp) {
	mm_count_t val;
	khint_t __i;
	for (__i = kh_begin(mcmap); __i != kh_end(mcmap); ++__i) {
		if (!kh_exist(mcmap,__i)) continue;
		val.mer = kh_key(mcmap, __i);
		val.count = kh_val(mcmap, __i);
		kv_push(mm_count_t, NULL, *cp, val);
	}
}

void aggregate_mm_count(khash_t(MMC) *mcmap,  mm_count_v * p) {
	uint32_t idx;
	khiter_t k;
	uint64_t mhash;
	int32_t absent;

	for (idx = 0; idx < p->n; idx++) {
		mhash = p->a[idx].mer;
		k = kh_put(MMC, mcmap, mhash, &absent);
		if (absent) {
			kh_value(mcmap, k) = 1;
		} else {
			kh_value(mcmap, k) += p->a[idx].count;
		}
	}
}

void write_mm_count(char *fn, mm_count_v *p) {
	FILE *out_data;
	out_data = fopen(fn, "wb");
	if (!out_data) {
		fprintf(stderr, "file '%s' open error: %s\n", fn, strerror(errno));
		exit(1);
	}
    fwrite(&(kv_size(*p)) , sizeof(size_t), 1, out_data);
    fwrite(p->a, sizeof(mm_count_t), kv_size(*p), out_data);
    fclose(out_data);
};


mm_count_v read_mm_count(char *fn) {
	mm_count_v p = {0, 0, 0};
	FILE *in_data;
	in_data = fopen(fn, "rb");
	if (!in_data) {
		fprintf(stderr, "file '%s' open error: %s\n", fn, strerror(errno));
		exit(1);
	}
    fread(&(kv_size(p)) , sizeof(size_t), 1, in_data);
	kv_resize(mm_count_t, NULL, p, kv_size(p));
    fread(p.a, sizeof(mm_count_t), kv_size(p), in_data);
    fclose(in_data);
	return p;
};
