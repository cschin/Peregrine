all: test
bogus: simreads clean test

simreads:
	mkdir -p ./reads
	python simulate_reads.py

reads_0.fa: simreads

test:
	./run_test.sh > all.log;  mkdir -p logs;  mv *.log logs 


clean:
	rm -rf ./wd/ ./logs/ ./reads/ seq_dataset.lst
