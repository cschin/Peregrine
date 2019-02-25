from intervaltree import Interval, IntervalTree
tree = IntervalTree()

read_range = {}
with open("reads.bed") as f:
    for row in f:
        row = row.strip().split()
        rname = row[0]
        s = int(row[1])
        e = int(row[2])
        tree.addi(s, e, rname)
        read_range[rname] = (s, e)
        if s < 40000:
            tree.addi(s+4639694, e+4639694, rname)

for rname in read_range:
    s, e = read_range[rname]
    for itvl in tree[s:e]:
        if itvl.data == rname:
            continue
        print(rname, itvl.data)


