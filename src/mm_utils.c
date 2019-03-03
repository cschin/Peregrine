#include <errno.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include "shimmer.h"
#include "khash.h"
#include "kvec.h"
#include "kalloc.h"

#define ORIGINAL 0
#define REVERSED 1

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
	memcpy(target->a + target->n, source->a, sizeof(mm128_t) * source->n);
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


void build_map(
		mm128_v * mmers,
		khash_t(MMER0) * mmer0_map, 
		khash_t(RLEN) *rlmap, 
		khash_t(MMC) *mcmap,
		uint32_t chunk,
		uint32_t total_chunk,
		uint32_t lowerbound,
		uint32_t upperbound) {

    uint64_t mhash;
    mm128_t mmer0, mmer1;
	rp128_v * rpv;
	uint32_t rid;
    uint32_t pos, rpos;
	uint32_t span;
	uint32_t mcount = 0;
	int32_t absent;
	rp128_t rp;
	khiter_t k;
	khash_t(MMER1) * mmer1_map;
	size_t s=0;

	for(;;) {
		mmer0 = mmers->a[s];
		mhash = mmer0.x >> 8;
		k = kh_get(MMC, mcmap, mhash);
		assert(k != kh_end(mcmap));
        mcount = kh_val(mcmap, k); 
		if (mcount >= lowerbound && mcount < upperbound) break;
		s++;
	}
	for( size_t i=s+1; i < mmers->n; i++ ){
		mmer1 = mmers->a[i];
		mhash = mmer1.x >> 8;
		k = kh_get(MMC, mcmap, mhash);
		assert(k != kh_end(mcmap));
        mcount = kh_val(mcmap, k); 
		if (mcount < lowerbound || mcount > upperbound) continue;

	    if ((mmer0.y >> 32) == (mmer1.y >> 32)) {  // the pairs are in the same read

			// don't use two minimers that are too close to each other
			if (((mmer1.y >> 1) & 0xFFFFFFF) - ((mmer0.y >> 1) & 0xFFFFFFF) < 100) {
				mmer0 = mmer1;
				continue;
			}

			if ((mmer0.x >> 8) % total_chunk == chunk % total_chunk) {
				k = kh_put(MMER0, mmer0_map, mmer0.x, &absent);
				if (absent) kh_value(mmer0_map, k) = kh_init(MMER1);

				mmer1_map = kh_value(mmer0_map, k);
				k = kh_put(MMER1, mmer1_map, mmer1.x, &absent);
				if (absent) {
					rpv = kmalloc(NULL, sizeof(rp128_v));
					rpv->n = 0; rpv->m = 0; rpv->a = NULL;
					kh_value(mmer1_map, k) = rpv;
				} else {
					rpv = kh_value(mmer1_map, k);
				}
				rp.y0 = mmer0.y;
				rp.y1 = mmer1.y;

				rp.direction = ORIGINAL;
				kv_push(rp128_t, NULL, *rpv, rp);
				// printf("Y %lu %lu %09u\n", mmer0.x >> 8, mmer1.x >> 8, (uint32_t) (rp.y0 >> 32));
			}
			
			// reverse direction
			if ((mmer1.x >> 8) % total_chunk == chunk % total_chunk) {
				k = kh_put(MMER0, mmer0_map, mmer1.x, &absent);
				if (absent) kh_value(mmer0_map, k) = kh_init(MMER1);

				mmer1_map = kh_value(mmer0_map, k);
				k = kh_put(MMER1, mmer1_map, mmer0.x, &absent);
				if (absent) {
					rpv = kmalloc(NULL, sizeof(rp128_v));
					rpv->n = 0; rpv->m = 0; rpv->a = NULL;
					kh_value(mmer1_map, k) = rpv;
				} else {
					rpv = kh_value(mmer1_map, k);
				}
				rp.y0 = mmer1.y;
				span = mmer1.x & 0xFF;
				rid = rp.y0 >> 32;
				pos = ((rp.y0 & 0xFFFFFFFF) >> 1) + 1;
				k = kh_get(RLEN, rlmap, rid);
				assert(k != kh_end(rlmap));
				rpos = kh_val(rlmap, k).len - pos + span -1;
				rp.y0 = ((rp.y0 & 0xFFFFFFFF00000001) | (rpos << 1 )) ^ 0x1; // ^0x1 -> switch strand

				rp.y1 = mmer0.y;
				span = mmer0.x & 0xFF;
				rid = rp.y1 >> 32;
				pos = ((rp.y1 & 0xFFFFFFFF) >> 1) + 1;
				k = kh_get(RLEN, rlmap, rid);
				assert(k != kh_end(rlmap));
				rpos = kh_val(rlmap, k).len - pos + span - 1;
				rp.y1 = ((rp.y1 & 0xFFFFFFFF00000001) | (rpos << 1 )) ^ 0x1; // ^0x1 -> switch strand
				rp.direction = REVERSED;

				kv_push(rp128_t, NULL, *rpv, rp);
				// printf("Y %lu %lu %09u\n", mmer1.x >> 8, mmer0.x >> 8, rid);
			}
		}	
		mmer0 = mmer1;
	}

}

void build_map4py(
		mm128_v * mmers,
		void * mmer0_map,
		void * rlmap,
		void * mcmap,
		uint32_t chunk,
		uint32_t total_chunk,
		uint32_t lowerbound,
		uint32_t upperbound) {
	build_map(mmers, 
			(khash_t(MMER0) *) mmer0_map,
			(khash_t(RLEN) *) rlmap,
			(khash_t(MMC) *) mcmap,
			chunk, 
			total_chunk,
			lowerbound,
			upperbound);
}
