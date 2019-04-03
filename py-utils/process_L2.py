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

L2dump = open("L2.txt", "w")

#hmmerL0 = ffi.new("mm128_v *")
#hmmerL2 = ffi.new("mm128_v *")

hmmerL0 = mm_utils.read_mmlist(b"../test/hmmer-L0-01-of-01.dat")
hmmerL2 = mm_utils.read_mmlist(b"../test/hmmer-L2-01-of-01.dat")

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
for i in range(hmmerL0.n):
    span = hmmerL0.a[i].x & 0xFF
    mmer = hmmerL0.a[i].x >> 8
    rid = hmmerL0.a[i].y >> 32
    pos_end = ((hmmerL0.a[i].y & 0xFFFFFFFF) >> 1) + 1
    strand = hmmerL0.a[i].y & 0x1
    mm_str = "{:014X}".format(mmer)
    mmer_count.setdefault(mm_str, 0)
    mmer_count[mm_str] += 1
    #
    # kmer = bseq[pos_end-span:pos_end]
    # kmer_r =  bytes([rmap[c] for c in kmer[::-1]])
    r_pos_end = rid2len[rid] - pos_end + span
    name = rid2name[rid]

mmer_count_L2 = {}
L2list = {}
for i in range(hmmerL2.n):
    span = hmmerL2.a[i].x & 0xFF
    mmer = hmmerL2.a[i].x >> 8
    rid = hmmerL2.a[i].y >> 32
    pos_end = ((hmmerL2.a[i].y & 0xFFFFFFFF) >> 1) + 1
    strand = hmmerL2.a[i].y & 0x1
    r_pos_end = rid2len[rid] - pos_end + span
    name = rid2name[rid]
    mm_str = "{:014X}".format(mmer)
    mmer_count_L2.setdefault(mm_str, 0)
    mmer_count_L2[mm_str] += 1
    print(name, pos_end, r_pos_end, strand,
          mm_str, mmer_count[mm_str], file=L2dump)
    L2list.setdefault(rid, [])
    L2list[rid].append((pos_end, r_pos_end, strand,
                        mm_str, name))
L2dump.close()

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
        #
        if mmer_count_L2[v_mmer] < 2:
            v = w
            continue
        w_pos_end, w_r_pos_end, w_strand, w_mmer, w_name = w
        #
        if mmer_count_L2[w_mmer] < 2:
            continue
        key = v_mmer, v_strand, w_mmer, w_strand
        L2map.setdefault(key, [])
        L2map[key].append((v_name, rid, 0, v_pos_end, w_pos_end))
        v = w

    v = lst[-1]
    for w in lst[-2::-1]:
        v_pos_end, v_r_pos_end, v_strand, v_mmer, v_name = v
        #
        if mmer_count_L2[v_mmer] < 2:
            v = w
            continue
        w_pos_end, w_r_pos_end, w_strand, w_mmer, w_name = w
        #
        if mmer_count_L2[w_mmer] < 2:
            continue
        key = v_mmer, 1-v_strand, w_mmer, 1-w_strand
        L2map.setdefault(key, [])
        L2map[key].append((v_name, rid, 1, v_r_pos_end, w_r_pos_end))
        v = w

dt_pairs = set()
with open("L0_dt.txt") as f:
    for row in f:
        row = row.strip().split()
        dt_pairs.add( (row[1], row[3]) )
        dt_pairs.add( (row[3], row[1]) )

import networkx as nx
G = nx.DiGraph()

r_dimer_set = {}
read_ovlp = {}
for key in L2map.keys():
    v_mmer, v_strand, w_mmer, w_strand = key
    for r in L2map[key]:
        v_name, v_rid, r_strand, v_pos_end, w_pos_end = r
        if v_name == "ref": continue
        r_dimer_set.setdefault(v_rid, set())
        r_dimer_set[v_rid].add(key)

for key in L2map.keys():
    v_mmer, v_strand, w_mmer, w_strand = key
    rlist = []
    n = len(L2map[key])
    for r in L2map[key]:
        v_name, v_rid, r_strand, v_pos_end, w_pos_end = r
        if v_name == "ref": continue
        left_offset_v = -v_pos_end
        right_offset_v = rid2len[v_rid]-v_pos_end
        dist = w_pos_end - v_pos_end
        print("X", v_mmer, v_strand, w_mmer, w_strand,
              v_name, r_strand, rid2len[v_rid],
              v_pos_end, w_pos_end,
              mmer_count_L2[v_mmer], mmer_count_L2[w_mmer], n,
              left_offset_v, right_offset_v)

        if mmer_count_L2[v_mmer] >= 30 or mmer_count_L2[v_mmer] <= 1:
            continue
        if mmer_count_L2[w_mmer] >= 30 or mmer_count_L2[w_mmer] <= 1:
            continue

        rlist.append((left_offset_v, right_offset_v, v_name,
                      v_rid, r_strand, dist ))

    if len(rlist) == 0: continue
    rlist.sort()

    p_set = set()
    left_offset_v0, right_offset_v0, r_name0, r_id0, r_strand0, dist0 = rlist[0]
    p_set = r_dimer_set[r_id0]
    for left_offset_v, right_offset_v, r_name, r_id, r_strand, dist in rlist[1:]:
        if right_offset_v0 < right_offset_v:
            overlap_count = len(r_dimer_set[r_id] & p_set)
            overlap_len = rid2len[r_id0] - abs(left_offset_v0-left_offset_v)
            dt_match = 1 if (r_name0, r_name) in dt_pairs else 0
            print("Y", v_mmer, v_strand, w_mmer, w_strand,
                  r_name0, r_strand0, r_name, r_strand,
                  overlap_count, overlap_len, left_offset_v0, left_offset_v, dt_match)
            if dt_match == 1:
                read_ovlp.setdefault((r_name0, r_strand0), [])
                read_ovlp.setdefault((r_name, 1-r_strand), [])
                read_ovlp[(r_name0, r_strand0)].append((overlap_len, (r_name, r_strand)))
                read_ovlp[(r_name, 1-r_strand)].append((-overlap_len, (r_name0, 1-r_strand0)))

                # G.add_edge("{}-{}".format(r_name, r_strand), "{}-{}".format(r_name0, r_strand0))
                # G.add_edge("{}-{}".format(r_name0, r_strand0), "{}-{}".format(r_name, r_strand))
        p_set = r_dimer_set[r_id]
        right_offset_v0 = right_offset_v
        left_offset_v0 = left_offset_v
        r_name0 = r_name
        r_strand0 = r_strand

for k in read_ovlp:
    read_ovlp[k].sort()

    for v in read_ovlp[k][:]:
        G.add_edge( "{}-{}".format(*k), "{}-{}".format(*v[-1]))
    #for v in read_ovlp[k][-3:]:
    #    G.add_edge( "{}-{}".format(*k), "{}-{}".format(*v[-1]))

nx.write_gexf(G, "test.gexf")

C.free(hmmerL0.a)
C.free(hmmerL2.a)
