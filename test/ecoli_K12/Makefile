all: test
bogus: simreads clean test

simreads:
	mkdir -p ./reads
	python simulate_reads.py

reads_0.fa: simreads

test:
	/usr/bin/time ./run_test.sh > all.log 2>&1;  mkdir -p logs;  mv *.log logs 

clean:
	rm -rf ./wd/ ./logs/ ./reads/ seq_dataset.lst
