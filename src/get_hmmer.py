from cffi import FFI
# import redis
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
# r_conn = redis.Redis(host='127.0.0.1', port=6379, db=0)

rmap = dict(zip(b"ACGT", b"TGCA"))

f = FastaReader.FastaReader("./test.fa")
rid = 0
km = ffi.NULL
L0dump = open("L0.txt", "w")
L1dump = open("L1.txt", "w")
L2dump = open("L2.txt", "w")

p = ffi.new("mm128_v *")
p_out = ffi.new("mm128_v *")
p_out2 = ffi.new("mm128_v *")

rid2name = {}
rid2len = {}
for r in f:
    bseq=r.sequence.encode("ascii")
    sketch.mm_sketch(km, bseq, len(bseq), 80, 16, rid, 0, p)
    rid2name[rid] = r.name
    rid2len[rid] = len(bseq)
    rid += 1

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
    pos_end = ((p.a[i].y & 0xFFFFFFFF) >> 1) + 1
    strand = p.a[i].y & 0x1
    #
    # kmer = bseq[pos_end-span:pos_end]
    # kmer_r =  bytes([rmap[c] for c in kmer[::-1]])
    r_pos_end = rid2len[rid] - pos_end + span
    print(r.name, pos_end, r_pos_end,
          strand, "{:014X}".format(mmer), file=L0dump)
    # r_conn.rpush(f"rid{rid}:L0", f"{pos_end} {r_pos_end} {strand} {kmer} {kmer_r}, {mmer}")

mm_select.mm_select(p, p_out, 8)
mmer_count = {}
for i in range(p_out.n):
    span = p_out.a[i].x & 0xFF
    mmer = p_out.a[i].x >> 8
    mmer_count.setdefault("{:014X}".format(mmer), 0)
    mmer_count["{:014X}".format(mmer)] += 1
    rid = p_out.a[i].y >> 32
    pos_end = ((p_out.a[i].y & 0xFFFFFFFF) >> 1) + 1
    strand = p_out.a[i].y & 0x1
    #
    # kmer = bseq[pos_end-span:pos_end]
    # kmer_r =  bytes([rmap[c] for c in kmer[::-1]])
    r_pos_end = rid2len[rid] - pos_end + span
    name = rid2name[rid]
    print(name, pos_end, r_pos_end, strand,
          "{:014X}".format(mmer), file=L1dump)
    # print(rid, pos_end, len(bseq)-pos_end+span, strand, kmer, kmer_r, mmer)

mm_select.mm_select(p_out, p_out2, 8)
L2list = {}
for i in range(p_out2.n):
    span = p_out2.a[i].x & 0xFF
    mmer = p_out2.a[i].x >> 8
    rid = p_out2.a[i].y >> 32
    pos_end = ((p_out2.a[i].y & 0xFFFFFFFF) >> 1) + 1
    strand = p_out2.a[i].y & 0x1
    #
    # kmer = bseq[pos_end-span:pos_end]
    # kmer_r =  bytes([rmap[c] for c in kmer[::-1]])
    r_pos_end = rid2len[rid] - pos_end + span
    name = rid2name[rid]
    print(name, pos_end, r_pos_end, strand,
          "{:014X}".format(mmer), file=L2dump)
    L2list.setdefault(rid, [])
    L2list[rid].append((pos_end, r_pos_end, strand,
                      "{:014X}".format(mmer), name))

L2map = {}
rspan = {}
for rid in L2list:
    lst = L2list[rid]
    if len(lst) < 2:
        continue
    rspan[rid] = lst[0][-2], lst[-1][-2]
    v = lst[0]
    for w in lst[1:]:
        v_pos_end, v_r_pos_end, v_strand, v_mmer, v_name = v
        if mmer_count[v_mmer] < 2:
            v = w
            continue
        w_pos_end, w_r_pos_end, w_strand, w_mmer, w_name = w
        if mmer_count[w_mmer] < 2:
            continue
        key = v_mmer, v_strand, w_mmer, w_strand
        L2map.setdefault(key, [])
        L2map[key].append((v_name, rid, 0, v_pos_end, w_pos_end))
        v = w

    v = lst[-1]
    for w in lst[-2::-1]:
        v_pos_end, v_r_pos_end, v_strand, v_mmer, v_name = v
        if mmer_count[v_mmer] < 2:
            v = w
            continue
        w_pos_end, w_r_pos_end, w_strand, w_mmer, w_name = w
        if mmer_count[w_mmer] < 2:
            continue
        key = v_mmer, 1-v_strand, w_mmer, 1-w_strand
        L2map.setdefault(key, [])
        L2map[key].append((v_name, rid, 1, v_r_pos_end, w_r_pos_end))
        v = w

import networkx as nx
G = nx.DiGraph()
for key in L2map.keys():
    v_mmer, v_strand, w_mmer, w_strand = key
    for r in L2map[key]:
        n = len(L2map[key])
        v_name, v_rid, v_strand, v_pos_end, w_pos_end = r
        print(v_mmer, v_strand, w_mmer, w_strand,
              v_name, v_strand, rid2len[v_rid],
              v_pos_end, w_pos_end,
              mmer_count[v_mmer], mmer_count[w_mmer], n)
        if n > 2 and n < 30:
            G.add_edge(v_mmer, w_mmer)

for rid in rspan:
    v, w = rspan[rid]
    if mmer_count[v] >= 30 or mmer_count[v] <= 2:
        continue
    if mmer_count[w] >= 30 or mmer_count[w] <= 2:
        continue
    G.add_edge(v, w)
    G.add_edge(w, v)

nx.write_gexf(G, "test.gexf")
C.free(p_out2.a)
C.free(p_out.a)
C.free(p.a)

L0dump.close()
L1dump.close()
L2dump.close()
