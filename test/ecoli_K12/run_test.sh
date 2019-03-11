#!/bin/bash
set -e
find ./reads/ -name "reads_*.fa" > seq_dataset.lst
SHIMMER=../../..
SHIMMER=$(cd "$(dirname "../../../")"; pwd)/$(basename "$1")
SHIMMERBIN=$SHIMMER/src
WORKDIR=./wd/
INDEX=$WORKDIR/index
OVLOUT=$WORKDIR/ovlp
ASM=$WORKDIR/asm
pushd $SHIMMER
echo SHIMMER revision: $(git rev-parse HEAD)
popd
echo get SHIMMER binaries from $SHIMMER
mkdir -p $INDEX
mkdir -p $OVLOUT
mkdir -p $ASM
echo
echo build read index
time (/usr/bin/time $SHIMMERBIN/build_read_index -p $INDEX/seq_dataset -d seq_dataset.lst 2> build_db.log)
echo
echo build shimmer index
time (for c in `seq 1 12`; do echo "/usr/bin/time $SHIMMERBIN/build_shimmer_index -p $INDEX/seq_dataset -t 12 -c $c -o $INDEX/shimmer 2> build_index.$c.log" ; done | parallel -j 4)
echo
echo build overlaps
time (for c in `seq -f "%02g" 1 8`; do echo "/usr/bin/time $SHIMMERBIN/shimmer_to_overlap -p $INDEX/seq_dataset -l $INDEX/shimmer-L2 -t 8 -c $c > $OVLOUT/out.$c 2> ovlp.$c.log"; done | parallel -j 4)
echo
echo faclon ovlp to graph
cd $ASM
time (cat ../ovlp/out.* | $SHIMMERBIN/dedup_overlaps > preads.ovl; echo "-" >> preads.ovl)
cp $SHIMMER/py/graph_to_contig.py .
cp $SHIMMER/py/ovlp_to_graph.py .
/usr/bin/time pypy ovlp_to_graph.py >& asm.log
ln -sf ../index/seq_dataset.* .
/usr/bin/time pypy graph_to_contig.py >& to_contig.log
