#!/bin/bash
set -e
ln -sf ../ecoli_K12/reads/ .
ln -sf ../ecoli_K12/K12MG1655.fa .
find ./reads/ -name "reads_*.fa" > seq_dataset.lst
SHIMMER=../../..
SHIMMER=$(cd "$(dirname "../../../")"; pwd)/$(basename "$1")
SHIMMERBIN=$SHIMMER/src
WORKDIR=./wd/
INDEX=$WORKDIR/index
pushd $SHIMMER
echo SHIMMER revision: $(git rev-parse HEAD)
popd
echo get SHIMMER binaries from $SHIMMER
mkdir -p $INDEX

echo
echo build read index
time (/usr/bin/time $SHIMMERBIN/shmr_mkseqdb -p $INDEX/seq_dataset -d seq_dataset.lst 2> build_db.log)

echo
echo build ref index
echo K12MG1655.fa > ref.lst
time (/usr/bin/time $SHIMMERBIN/shmr_mkseqdb -p $INDEX/ref -d ref.lst 2> build_ref_db.log)

echo build ref shimmer index
time (for c in `seq 1 6`; do echo "/usr/bin/time $SHIMMERBIN/shmr_index -p $INDEX/seq_dataset -t 6 -c $c -o $INDEX/read 2> build_index.$c.log" ; done | parallel -j 4)

echo build ref shimmer index
time (for c in `seq 1 2`; do echo "/usr/bin/time $SHIMMERBIN/shmr_index -p $INDEX/ref -t 2 -c $c -o $INDEX/ref 2> build_ref_index.$c.log" ; done | parallel -j 2)

echo run shimmer_map
$SHIMMERBIN/shmr_map -r $INDEX/ref -m $INDEX/ref-L2 -p $INDEX/seq_dataset  -l $INDEX/read-L2  -t 1 -c 1  >  reads2ref.out

$SHIMMERBIN/shmr_map -r $INDEX/ref -m $INDEX/ref-L2 -p $INDEX/ref  -l $INDEX/ref-L2  -t 1 -c 1 > ref2ref.out
