#!/bin/bash
set -e -o pipefail
. /root/.bashrc
make simreads
make test-pypeflow
mkdir -p /wd
cp -a ./wd-pf/ /wd/ecoli_test_results/
cp K12MG1655.fa /wd/ecoli_test_results/
apt-get install -y mummer
cd /wd/ecoli_test_results/
dnadiff K12MG1655.fa p_ctg_cns.fa 
echo 
echo dnadiff output of the assembled contig to the e. coli genome used for the simulated reads
cat out.report
