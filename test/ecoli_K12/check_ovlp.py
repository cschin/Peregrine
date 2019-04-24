from intervaltree import Interval, IntervalTree
import glob
tree = IntervalTree()

rname2rid = {}
with open("wd-pf/0-seqdb/seq_dataset.idx") as f:
    for row in f:
        row = row.strip().split()
        rname2rid[row[1]]=row[0]

read_range = {}
for fn in glob.glob("reads/*.bed"):
    with open(fn) as f:
        for row in f:
            row = row.strip().split()
            rname = row[0]
            s = int(row[1])
            e = int(row[2])
            tree.addi(s, e, rname2rid[rname])
            read_range[rname2rid[rname]] = (s, e)
            if s < 40000:
                tree.addi(s+4639694, e+4639694, rname2rid[rname])
readpair = set()
for rid in read_range:
    s, e = read_range[rid]
    for itvl in tree[s:e]:
        if itvl.data == rid:
            continue
        print("X", rid, itvl.data)
        readpair.add( (rid, itvl.data) )
        readpair.add( (itvl.data, rid) )

ovlppair = set()
with open("wd-pf/3-asm/preads.ovl") as f:
    for row in f:
        row = row.strip().split()
        if row[0] == "-":
            continue
        if (row[0], row[1]) in readpair:
            row.append("1")
        else:
            row.append("0")
        print("Y"," ".join(row))
        ovlppair.add( (row[1], row[0]) )
        ovlppair.add( (row[0], row[1]) )

for op in readpair:
    r1 = read_range[op[0]]
    r2 = read_range[op[1]]
    if r1[0] < r2[0]:
        olen = r1[1] - r2[0]
    else:
        olen = r2[1] - r1[0]
    if op in ovlppair:
        op = list(op)

        print("Z {} {} {} 1".format(op[0], op[1], olen))
    else:
        print("Z {} {} {} 0".format(op[0], op[1], olen))




