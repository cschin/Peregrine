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

KHASH_MAP_INIT_INT64(MMC, uint32_t);
void mm_count(mm128_v * p, mm_count_v * cp) {
	uint32_t idx;
	khiter_t k;
	mm128_t mmer;
	uint64_t mhash;
	uint32_t count;
	int32_t absent;
	mm_count_t val;
    khash_t(MMC) *mcmap=kh_init(MMC);

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

	{
		khint_t __i;
		for (__i = kh_begin(mcmap); __i != kh_end(mcmap); ++__i) {
			if (!kh_exist(mcmap,__i)) continue;
			mhash = kh_key(mcmap, __i);
			count = kh_val(mcmap, __i);
			val.mer = mhash;
			val.count = count;
			kv_push(mm_count_t, NULL, *cp, val);
		}
	}
	kh_destroy(MMC, mcmap);
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
