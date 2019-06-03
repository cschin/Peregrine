all: test test-pypeflow
bogus: simreads clean test test-pypeflow

K12MG1655.fa:
	wegt https://s3.amazonaws.com//biologicaldatascience.org/data/ecoli-k12/K12MG1655.fa

simreads: K12MG1655.fa
	mkdir -p ./reads
	python simulate_reads.py

reads_0.fa: simreads

test-pypeflow:
	rm -rf ./wd-pf
	find ${PWD}/reads/ -name "reads_*.fa" > seq_dataset.lst
	/usr/bin/time pg_run.py asm seq_dataset.lst 12 4 8 4 1 1 1 1 1 --with-consensus --output ./wd-pf

clean:
	rm -rf ./wd/ ./logs/ ./reads/ seq_dataset.lst
