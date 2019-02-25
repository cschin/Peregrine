from cffi import FFI
# import redis

ffi = FFI()

ffi.cdef("""
typedef struct { uint64_t x, y; } mm128_t;
typedef struct { size_t n, m; mm128_t *a; } mm128_v;
mm128_v read_mmlist(char *);
void free(void *ptr);
""")

C = ffi.dlopen(None)
mm_utils = ffi.dlopen("../src/mm_utils.so")
# r_conn = redis.Redis(host='127.0.0.1', port=6379, db=0)

rmap = dict(zip(b"ACGT", b"TGCA"))

L0dump = open("L0.txt", "w")

#hmmerL0 = ffi.new("mm128_v *")
#hmmerL2 = ffi.new("mm128_v *")

hmmerL0 = mm_utils.read_mmlist(b"../test/hmmer-L0-01-of-01.dat")

rid2name = {}
rid2len = {}
# rid2seq = {}

with open("../test/seq_dataset.idx") as f:
    for row in f:
        row = row.strip().split()
        rid, rname, rlen = row
        rid = int(rid)
        rlen = int(rlen)
        rid2name[rid] = rname
        rid2len[rid] = rlen

"""
* @param p      minimizers
*               p->a[i].x = kMer<<8 | kmerSpan
*               p->a[i].y = rid<<32 | lastPos<<1 | strand
*               where lastPos is the position of the last base of the i-th minimizer,
*               and strand indicates whether the minimizer comes from the top or the bottom strand.
*               Callers may want to set "p->n = 0"; otherwise results are appended to p
"""

mmer_count = {}
mer_five = {}
mer_three = {}
for i in range(hmmerL0.n):
    span = hmmerL0.a[i].x & 0xFF
    mmer = hmmerL0.a[i].x >> 8
    rid = hmmerL0.a[i].y >> 32
    pos_end = ((hmmerL0.a[i].y & 0xFFFFFFFF) >> 1) + 1
    strand = hmmerL0.a[i].y & 0x1
    mm_str = "{:014X}".format(mmer)
    #
    # mmer_count.setdefault(mm_str, 0)
    # mmer_count[mm_str] += 1
    #
    # kmer = bseq[pos_end-span:pos_end]
    # kmer_r =  bytes([rmap[c] for c in kmer[::-1]])
    r_pos_end = rid2len[rid] - pos_end + span
    name = rid2name[rid]

    if pos_end < 250:
        mer_five.setdefault(mmer, [])
        mer_five[mmer].append(name)
    if r_pos_end < 250:
        mer_three.setdefault(mmer, [])
        mer_three[mmer].append(name)

    print(name, pos_end, r_pos_end,
          strand, mm_str, file=L0dump)

L0dump.close()

dovetail_end = {}
for i in range(hmmerL0.n):
    mmer = hmmerL0.a[i].x >> 8
    rid = hmmerL0.a[i].y >> 32
    rname = rid2name[rid]
    if mmer in mer_five:
        for rname0 in mer_five[mmer]:
            dovetail_end.setdefault(rname0, set())
            dovetail_end[rname0].add( (5, rname) )
    if mmer in mer_three:
        for rname0 in mer_three[mmer]:
            dovetail_end.setdefault(rname0, set())
            dovetail_end[rname0].add( (3, rname) )

dt_file = open("L0_dt.txt","w")
for rname in dovetail_end:
    for e, rname0 in list(dovetail_end[rname]):
        if rname == rname0:
            continue
        intersect = 0
        if e == 5 and (3, rname) in dovetail_end.get(rname0, {}):
            print ( 5, rname0, 3, rname, file=dt_file)
        if e == 5 and (5, rname) in dovetail_end.get(rname0, {}):
            print ( 5, rname0, 5, rname, file=dt_file)
        if e == 3 and (5, rname) in dovetail_end.get(rname0, {}):
            print ( 3, rname0, 5, rname, file=dt_file)
        if e == 3 and (3, rname) in dovetail_end.get(rname0, {}):
            print ( 3, rname0, 3, rname, file=dt_file)

dt_file.close()

C.free(hmmerL0.a)

