all: test test-pypeflow
bogus: simreads clean test test-pypeflow

K12MG1655.fa:
	wget https://www.dropbox.com/s/wqqnzachbdk4d3r/K12MG1655.fa

simreads: K12MG1655.fa
	mkdir -p ./reads
	python simulate_reads.py

reads_0.fa: simreads

test:
	rm -rf ./wd ./logs
	/usr/bin/time ./run_test.sh > all.log 2>&1;  mkdir -p logs;  mv *.log logs 

test-pypeflow:
	rm -rf ./wd-pf
	find ${PWD}/reads/ -name "reads_*.fa" > seq_dataset.lst
	/usr/bin/time pg_run.py asm seq_dataset.lst 12 4 8 4 1 1 1 1 1 --with-consensus --output ./wd-pf

test-pypeflow-with-L0:
	rm -rf ./wd-pf
	find ${PWD}/reads/ -name "reads_*.fa" > seq_dataset.lst
	/usr/bin/time pg_run.py asm seq_dataset.lst 12 4 8 4 1 1 1 1 1 --with-L0-index --with-consensus --output ./wd-pf

test-pypeflow-l1:
	rm -rf ./wd-pf-l1
	find ${PWD}/reads/ -name "reads_*.fa" > seq_dataset.lst
	/usr/bin/time pg_run.py asm seq_dataset.lst 12 4 8 4 1 1 1 1 1 --shimmer-r 24 --with-consensus --shimmer-l 1 --output ./wd-pf-l1

clean:
	rm -rf ./wd/ ./logs/ ./reads/ seq_dataset.lst
