CC=gcc
CFLAGS=-O3 
CFLAGS+=-Wall -Wno-unused-result -Wno-unused-function -Werror
#CFLAGS+=-Wall -Wno-unused-result -Wno-unused-function
LDLIBS=-lz

all: build_read_index build_hmmer accumulate_mc mm_utils.so

kalloc.o: kalloc.c

mm_reduce.o: mm_reduce.c

mm_sketch.o: mm_sketch.c

mm_utils.o: mm_utils.c

mm_utils.so: mm_utils.c
	gcc -g -shared -fPIC -Wall kalloc.c mm_utils.c -o mm_utils.so

mm_end_filter.o: mm_end_filter.c

accumulate_mc.o: accumulate_mc.c

build_seq_db: build_seqdb.c kalloc.o 

build_hmmer: build_hmmer.c kalloc.o mm_reduce.o mm_sketch.o mm_utils.o mm_end_filter.o

accumulate_mc: accumulate_mc.o kalloc.o mm_utils.o



clean:
	rm -f build_read_index build_hmmer accumulate_mc *.o *.so
