all: test
bogus: simreads clean test

K12MG1655.fa:
	wget https://www.dropbox.com/s/wqqnzachbdk4d3r/K12MG1655.fa

simreads: K12MG1655.fa
	mkdir -p ./reads
	python simulate_reads.py

reads_0.fa: simreads

test:
	rm -rf ./wd ./logs
	/usr/bin/time ./run_test.sh > all.log 2>&1;  mkdir -p logs;  mv *.log logs 

clean:
	rm -rf ./wd/ ./logs/ ./reads/ seq_dataset.lst
