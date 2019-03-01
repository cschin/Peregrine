#include <errno.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include "shimmer.h"
#include "khash.h"
#include "kvec.h"
#include "kalloc.h"


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

void append_mmlist(mm128_v * target, mm128_v * source) {
	kv_resize(mm128_t, NULL, *target, kv_size(*target) + kv_size(*source) );
	memcpy(target->a + target->n, source->a, sizeof(mm128_v) * source->n);
	target->n += source->n;
}

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
			kh_value(mcmap, k) = 0 ;
		} 
		kh_value(mcmap, k) += p->a[idx].count;
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

khash_t(RID) * build_read_index(char *fpath, seq_data_v *seq_data, khash_t(RLEN) *rlmap ) {
	khash_t(RID) *hmap = kh_init(RID);
	int absent;
	khiter_t k;
	FILE *index_file;
	uint32_t i, rid, rlen, offset;
	char name_buf[256];
	char * name;
	rl_t rl;

	index_file = fopen(fpath, "r");
	while (fscanf(index_file, "%u %255s %u %u", &rid, name_buf, &rlen, &offset) != EOF) {
		kv_resize(seq_data_t, NULL, *seq_data, seq_data->n + 8192);
		name = kmalloc(NULL, 256 * sizeof(char));
		strncpy(name, name_buf, 256);
		i = seq_data->n;
		seq_data->a[i].name = name;
		seq_data->a[i].rid = rid;
		seq_data->n ++;
		k=kh_put(RID, hmap, seq_data->a[i].name, &absent);
		if (absent) kh_value(hmap, k) = rid;
		k=kh_put(RLEN, rlmap, rid, &absent);
		rl.len = rlen;
		rl.offset = offset;
		kh_value(rlmap, k) = rl;
		/* 
        { 
			khint_t __i;		
			for (__i = kh_begin(hmap); __i != kh_end(hmap); ++__i) {		
				if (!kh_exist(hmap,__i)) continue;						
				printf("testx:%s %u\n", kh_key(hmap,__i), kh_val(hmap,__i));
			} 
		}
	    */	
	}
	fclose(index_file);
	return hmap;
}

khash_t(RLEN) * get_read_length_map(char *fpath ) {
	khash_t(RLEN) *rlmap = kh_init(RLEN);
	int absent;
	khiter_t k;
	FILE *index_file;
	char name_buf[256];
	uint32_t  rid, rlen;
	size_t offset;
	rl_t rl;

	index_file = fopen(fpath, "r");
	while (fscanf(index_file, "%u %255s %u %lu", &rid, name_buf, &rlen, &offset) != EOF) {
		k=kh_put(RLEN, rlmap, rid, &absent);
		rl.len = rlen;
		rl.offset = offset;
		kh_value(rlmap, k) = rl;
	}
	fclose(index_file);
	return rlmap;
}

char * get_read_seq(FILE * seqdb, uint32_t rid, khash_t(RLEN) *rlmap) {
	char * seq;
	rl_t rl;
	khiter_t k;
	k = kh_get(RLEN, rlmap, rid);
    assert( k != kh_end(rlmap));
	rl = kh_val(rlmap, k);
	seq = kmalloc(NULL, sizeof(char) * rl.len+1);
	fseek(seqdb, rl.offset, SEEK_SET);
	fread(seq, sizeof(char), rl.len, seqdb);
	seq[rl.len] = 0; // terminate the string
    return seq;	
}	

char * get_read_seq_mmap(char * seq_p, uint32_t rid, khash_t(RLEN) *rlmap) {
	char * seq;
	rl_t rl;
	khiter_t k;
	k = kh_get(RLEN, rlmap, rid);
    assert( k != kh_end(rlmap));
	rl = kh_val(rlmap, k);
	seq = kmalloc(NULL, sizeof(char) * rl.len+1);
	strcpy(seq, seq_p + rl.offset);
    return seq;	
}	
