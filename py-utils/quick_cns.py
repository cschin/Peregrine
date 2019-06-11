import peregrine
import sys, os
from peregrine._falcon4py import ffi as falcon_ffi
from peregrine._falcon4py import lib as falcon4py
from peregrine._shimmer4py import ffi as shimmer_ffi
from peregrine._shimmer4py import lib as shimmer4py
from peregrine import utils


read_id = None
read_group = []
total_chunk = int(sys.argv[3])
chunk = int(sys.argv[2])
with open(sys.argv[1], "rb") as f:
    for r in f:
        name, seq = r.strip().split()
        crid = b"/".join(name.split(b"/")[:2])
        if hash(crid) % total_chunk != chunk:
            continue
        if crid != read_id:
            if len(read_group) >= 3:
                seqs = read_group
                if len(seqs) > 20:
                    seqs = seqs[:20]
                cns_seq = peregrine.utils.get_cns_from_reads(seqs,
                                                             levels=0,
                                                             k=10,
                                                             w=16 ,
                                                             max_dist=150)
                print(">"+read_id.decode("ascii")+"/qccs/ {}".format(len(seqs)))
                print(cns_seq.decode("ascii"))
            read_group = []
        read_id = crid
        read_group.append(seq)

