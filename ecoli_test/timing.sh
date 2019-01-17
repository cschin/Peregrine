set -x
(time bash assembly.sh  > log)
echo check output contig lengths
fastalength p_ctg.fa
time nucmer K12MG1655.fa p_ctg.fa
time mummerplot -l -t png  out.delta
time gnuplot out.gp
