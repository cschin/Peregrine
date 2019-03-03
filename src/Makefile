CC=gcc
CFLAGS=-O3 
CFLAGS+=-Wall -Wno-unused-result -Wno-unused-function -Werror
#CFLAGS+=-Wall -Wno-unused-result -Wno-unused-function
LDLIBS=-lz

all: build_read_index build_shimmer_index accumulate_mc mm_utils.so print_overlaps dedup_ovlp

kalloc.o: kalloc.c

mm_reduce.o: mm_reduce.c

mm_sketch.o: mm_sketch.c

mm_utils.o: mm_utils.c

mm_utils.so: mm_utils.c
	gcc -g -shared -fPIC -Wall kalloc.c mm_utils.c -o mm_utils.so

shimmer4py.so: mm_utils.c shimmer4py.c
	gcc -g -shared -fPIC -Wall kalloc.c mm_utils.c shimmer4py.c -o shimmer4py.so

mm_end_filter.o: mm_end_filter.c

accumulate_mc.o: accumulate_mc.c

build_seq_db: build_seqdb.c kalloc.o 

build_shimmer_index: build_shimmer_index.c kalloc.o mm_reduce.o mm_sketch.o mm_utils.o mm_end_filter.o

accumulate_mc: accumulate_mc.o kalloc.o mm_utils.o

DW_banded.o: DW_banded.c

print_overlaps: print_overlaps.c mm_utils.o kalloc.o  DW_banded.o

dedup_ovlp: dedup_ovlp.c kalloc.o

clean:
	rm -f dedup_ovlp accumulate_mc print_overlaps build_read_index build_shimmer_index accumulate_mc *.o *.so
