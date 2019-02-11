from cffi import FFI
import redis
import FastaReader



ffi = FFI()

ffi.cdef("""
typedef struct { uint64_t x, y; } mm128_t;
typedef struct { size_t n, m; mm128_t *a; } mm128_v;
void mm_sketch(void *km, const char *str, int len, int w, int k,
         uint32_t rid, int is_hpc, mm128_v *p);
void mm_select(mm128_v *p, mm128_v *p_out, uint8_t rs);
void free(void *ptr);
""")

C = ffi.dlopen(None)
sketch = ffi.dlopen("./sketch.so")
mm_select = ffi.dlopen("./mm_select.so")
r_conn = redis.Redis(host='127.0.0.1', port=6379, db=0)

rmap = dict(zip(b"ACGT", b"TGCA"))

f = FastaReader.FastaReader("../ecoli_test/simreads.fa")
#f = FastaReader.FastaReader("./test.fa")
rid = 0
km = ffi.NULL
L0dump = open("L0.txt", "w")
L1dump = open("L1.txt", "w")
L2dump = open("L2.txt", "w")

for r in f:
    bseq=r.sequence.encode("ascii")
    rid += 1

    #seq = ffi.new("char[]", bseq)

    p = ffi.new("mm128_v *")
    p_out = ffi.new("mm128_v *")
    p_out2 = ffi.new("mm128_v *")

    sketch.mm_sketch(km, bseq, len(bseq), 80, 14, rid, 0, p)

    """
    * @param p      minimizers
    *               p->a[i].x = kMer<<8 | kmerSpan
    *               p->a[i].y = rid<<32 | lastPos<<1 | strand
    *               where lastPos is the position of the last base of the i-th minimizer,
    *               and strand indicates whether the minimizer comes from the top or the bottom strand.
    *               Callers may want to set "p->n = 0"; otherwise results are appended to p
    """

    for i in range(p.n):
        span = p.a[i].x & 0xFF
        mmer = p.a[i].x >> 8
        rid = p.a[i].y >> 32
        pos_end =  ((p.a[i].y & 0xFFFFFFFF) >> 1) + 1
        strand = p.a[i].y & 0x1
        #kmer = bseq[pos_end-span:pos_end]
        #kmer_r =  bytes([rmap[c] for c in kmer[::-1]])
        r_pos_end = len(bseq)-pos_end+span
        print(r.name, pos_end, r_pos_end, strand, "{:014X}".format(mmer), file=L0dump)
        #r_conn.rpush(f"rid{rid}:L0", f"{pos_end} {r_pos_end} {strand} {kmer} {kmer_r}, {mmer}")

    mm_select.mm_select(p, p_out, 8)
    for i in range(p_out.n):
        span = p_out.a[i].x & 0xFF
        mmer = p_out.a[i].x >> 8
        rid = p_out.a[i].y >> 32
        pos_end =  ((p_out.a[i].y & 0xFFFFFFFF) >> 1) + 1
        strand = p_out.a[i].y & 0x1
        #kmer = bseq[pos_end-span:pos_end]
        #kmer_r =  bytes([rmap[c] for c in kmer[::-1]])
        r_pos_end = len(bseq)-pos_end+span
        print(r.name, pos_end, r_pos_end, strand, "{:014X}".format(mmer), file=L1dump)
        #print(rid, pos_end, len(bseq)-pos_end+span, strand, kmer, kmer_r, mmer)

    mm_select.mm_select(p_out, p_out2, 8)
    L2list = []
    for i in range(p_out2.n):
        span = p_out2.a[i].x & 0xFF
        mmer = p_out2.a[i].x >> 8
        rid = p_out2.a[i].y >> 32
        pos_end =  ((p_out2.a[i].y & 0xFFFFFFFF) >> 1) + 1
        strand = p_out2.a[i].y & 0x1
        #kmer = bseq[pos_end-span:pos_end]
        #kmer_r =  bytes([rmap[c] for c in kmer[::-1]])
        r_pos_end = len(bseq)-pos_end+span
        #print(rid, pos_end, len(bseq)-pos_end+span, strand, kmer, kmer_r, "{:014X}".format(mmer))
        print(r.name, pos_end, r_pos_end, strand, "{:014X}".format(mmer), file=L2dump)
        L2list.append( (pos_end, r_pos_end, strand, "{:014X}".format(mmer)) )

    if len(L2list) > 2:
        v = L2list[0]
        for w in L2list[1:]:
            print( v[-1], v[2], w[-1], w[2], 0, v[0], w[0], r.name )
            v = w
        v = L2list[-1]
        for w in L2list[-2:0:-1]:
            print( v[-1], 1-v[2], w[-1], 1-w[2], 1, v[1], w[1], r.name )
            v = w

    C.free(p_out2.a)
    C.free(p_out.a)
    C.free(p.a)
L0dump.close()
L1dump.close()
L2dump.close()
