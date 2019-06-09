import peregrine
import sys, os
from peregrine._falcon4py import ffi as falcon_ffi
from peregrine._falcon4py import lib as falcon4py
from peregrine._shimmer4py import ffi as shimmer_ffi
from peregrine._shimmer4py import lib as shimmer4py
from peregrine import utils
from FastaReader import FastaReader
f = FastaReader(sys.argv[1])
read_id = None
read_group = []
for r in f:
    crid = "/".join(r.name.split("/")[:2])
    if crid != read_id:
        if len(read_group) >= 5:
            seqs = read_group
            seqs.sort(key = lambda k: len(k))
            m = int(len(seqs)/2)+1
            seqs[0], seqs[m] = seqs[m], seqs[0]
            cns_seq = peregrine.utils.get_cns_from_reads(seqs,
                                                         levels=0,
                                                         k=10,
                                                         w=16 ,
                                                         max_dist=150)
            print(">"+read_id+"/qccs/")
            print(cns_seq.decode("ascii"))
        read_group = []
    read_id = crid
    if len(r.sequence) > 1000:
        read_group.append(r.sequence.encode("ascii"))

