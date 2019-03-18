all: test
bogus: clean test

test:
	rm -rf ./wd ./logs
	/usr/bin/time ./run_test.sh > all.log 2>&1;  mkdir -p logs;  mv *.log logs 

clean:
	rm -rf ./wd/ ./logs/ seq_dataset.lst reads2ref.out ref2ref.out
