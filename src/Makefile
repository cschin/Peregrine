CC=gcc
CFLAGS=-O3 
CFLAGS+=-Wall -Wno-unused-result -Wno-unused-function -Werror
#CFLAGS+=-Wall -Wno-unused-result -Wno-unused-function
LDLIBS=-lz

all: build_read_index build_shimmer_index accumulate_mc shimmer_to_overlap dedup_overlaps shimmer4py.so shimmer_utils.so

kalloc.o: kalloc.c

shimmer_reduce.o: shimmer_reduce.c

mm_sketch.o: mm_sketch.c

shimmer_utils.o: shimmer_utils.c

shimmer_utils.so: shimmer_utils.c
	gcc -O3 -shared -fPIC -Wall kalloc.c shimmer_utils.c -o shimmer_utils.so

shimmer4py.so: shimmer_utils.c shimmer4py.c
	gcc -O3 -shared -fPIC -Wall kalloc.c shimmer_utils.c shimmer4py.c -o shimmer4py.so

shimmer_end_filter.o: shimmer_end_filter.c

accumulate_mc.o: accumulate_mc.c

build_read_index: build_read_index.c kalloc.o shimmer_utils.o

build_shimmer_index: build_shimmer_index.c kalloc.o shimmer_reduce.o mm_sketch.o shimmer_utils.o shimmer_end_filter.o

accumulate_mc: accumulate_mc.o kalloc.o shimmer_utils.o

DWmatch.o: DWmatch.c

shimmer_to_overlap: shimmer_to_overlap.c shimmer_utils.o kalloc.o  DWmatch.o

dedup_overlaps: dedup_overlaps.c kalloc.o

clean:
	rm -f dedup_overlaps accumulate_mc shimmer_to_overlap build_read_index build_shimmer_index accumulate_mc *.o *.so
