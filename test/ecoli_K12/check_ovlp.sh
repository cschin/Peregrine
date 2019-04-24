python check_ovlp.py > check_ovlp.out
less check_ovlp.out | grep Z | sort -k 2g -k 4gr | awk '$4 > 10000 && $NF != 1' | wc
less check_ovlp.out | grep Z | sort -k 2g -k 4gr | awk '$4 > 10000' | wc
