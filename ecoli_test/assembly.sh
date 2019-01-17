python simread.py > simreads.fa
../samara/build_seq_db
../samara/dump_minimizers
ls mmer*.rs | awk -F "." '{print $2}'  > labels
for label in `cat labels`; do echo "../samara/mm_filter_1 $label | grep L2 > L2.$label"; done | bash
python3 ../samara/check_ovlp2.py L2.0?? > L2.summary.all.2 
python2.7 ovlp_to_graph.py
../samara/dump_seq  > preads4falcon.fasta
fc_graph_to_contig
fastalength p_ctg.fa
