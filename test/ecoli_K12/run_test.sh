#!/bin/bash
set -e
find ./reads/ -name "reads_*.fa" > seq_dataset.lst
SHIMMER=../../..
SHIMMER=$(cd "$(dirname "../../../")"; pwd)/$(basename "$1")
SHIMMERBIN=$SHIMMER/src/bin
WORKDIR=$PWD/wd
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
time (/usr/bin/time $SHIMMERBIN/shmr_mkseqdb -p $INDEX/seq_dataset -d seq_dataset.lst 2> build_db.log)
echo
echo build shimmer index
time (for c in `seq 1 12`; do echo "/usr/bin/time $SHIMMERBIN/shmr_index -p $INDEX/seq_dataset -t 12 -c $c -o $INDEX/shmr 2> build_index.$c.log" ; done | parallel -j 4)
echo
echo build overlaps
time (for c in `seq -f "%02g" 1 8`; do echo "/usr/bin/time $SHIMMERBIN/shmr_overlap -p $INDEX/seq_dataset -l $INDEX/shmr-L2 -t 8 -c $c > $OVLOUT/out.$c 2> ovlp.$c.log"; done | parallel -j 4)
echo
echo faclon ovlp to graph
cd $ASM
time (cat ../ovlp/out.* | $SHIMMERBIN/shmr_dedup > preads.ovl; echo "-" >> preads.ovl)
/usr/bin/time ovlp_to_graph.py >& asm.log
ln -sf ../index/seq_dataset.* .
#/usr/bin/time pypy graph_to_contig.py >& to_contig.log
/usr/bin/time graph_to_path.py >& to_path.log
/usr/bin/time path_to_contig.py $INDEX/seq_dataset p_ctg_tiling_path > p_ctg.fa 2> to_contig.log
echo $PWD/p_ctg.fa > p_ctg.lst
time (/usr/bin/time $SHIMMERBIN/shmr_mkseqdb -p $INDEX/p_ctg -d p_ctg.lst 2> build_p_ctg_db.log)
time (for c in `seq 1 1`; do echo "/usr/bin/time $SHIMMERBIN/shmr_index -p $INDEX/p_ctg -t 1 -c $c -o $INDEX/p_ctg 2> build_p_ctg_index.$c.log" ; done | parallel -j 4)
time (/usr/bin/time $SHIMMERBIN/shmr_map -r $INDEX/p_ctg -m $INDEX/p_ctg-L2 -p $INDEX/seq_dataset -l $INDEX/shmr-L2 -t 1 -c 1 > read_map.txt 2> map.log)
time (/usr/bin/time cns_prototype.py $INDEX/seq_dataset $INDEX/p_ctg read_map.txt 1 1 > p_ctg_cns.fa 2> cns.log)

