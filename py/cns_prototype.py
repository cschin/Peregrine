import mmap
import sys
from _falcon4py import ffi
from _falcon4py import lib as falcon
from _shimmer4py import lib as shimmer
import numpy as np



f=open("../test/ecoli_K12/wd/index/seq_dataset.seqdb", "rb")
seqdb = mmap.mmap(f.fileno(), 0, flags=mmap.MAP_SHARED, prot=mmap.PROT_READ)

f=open("../test/ecoli_K12/wd/index/p_ctg.seqdb", "rb")
refdb = mmap.mmap(f.fileno(), 0, flags=mmap.MAP_SHARED, prot=mmap.PROT_READ)

read_idx = {}
with open("../test/ecoli_K12/wd/index/seq_dataset.idx") as f:
    for row in f:
        row = row.strip().split()
        rid, rname, rlen, offset = row
        rid = int(rid)
        rlen = int(rlen)
        offset = int(offset)
        read_idx.setdefault(rid, {})
        read_idx[rid]["name"] = rname
        read_idx[rid]["length"] = rlen
        read_idx[rid]["offset"] = offset


ref_idx = {}
with open("../test/ecoli_K12/wd/index/p_ctg.idx") as f:
    for row in f:
        row = row.strip().split()
        rid, rname, rlen, offset = row
        rid = int(rid)
        rlen = int(rlen)
        offset = int(offset)
        ref_idx.setdefault(rid, {})
        ref_idx[rid]["name"] = rname
        ref_idx[rid]["length"] = rlen
        ref_idx[rid]["offset"] = offset


read_map_groups = []
left_anchor = 500
map_group=[]
with open("../test/ecoli_K12/wd/asm/read_map.txt") as f:
    for row in f:
        row = row.strip().split()
        row = tuple(int(c) for c in row)
        ref_p1 = row[1]
        if ref_p1 - left_anchor < 50000:
            map_group.append( row )
        else:
            read_map_groups.append( (left_anchor, ref_p1, map_group) )
            map_group = []
            left_anchor = ref_p1
    read_map_groups.append((left_anchor, ref_p1, map_group))

rng = ffi.new("aln_range[1]")

cns_segments = []
j = 0
for left, right, mapped in read_map_groups:
    print("-", j, left, right, right-left, file=sys.stderr)
    j += 1
    left = left-500
    assert(left>=0)
    rmap = {}
    for d in mapped:
        #print(d)
        read_id = d[3]
        read_offset = d[1] - d[4]
        read_strand = d[6]
        rmap.setdefault((read_id, read_strand),[])
        rmap[(read_id, read_strand)].append(read_offset )


    reads = []
    for (read_id, read_strand), v in rmap.items():
        reads.append((read_id, read_strand, np.min(v)-left, len(v)))

    reads.sort(key = lambda x: x[2])
    s = ref_idx[0]["offset"] + left
    ref_len = right-left

    bseq0 = refdb[s:s+ref_len]

    ref_seq = ffi.new("char[{}]".format(ref_len))

    shimmer.decode_biseq(bseq0, ref_seq, ref_len, 0)

    tags = ffi.new("align_tags_t * [{}]".format(len(reads)))

    aln_count = 0
    for d in reads:
        #print(d)

        read_id = d[0]
        read_strand = d[1]
        read_shift = int(d[2])
        s = read_idx[read_id]["offset"]
        read_len = read_idx[read_id]["length"]
        bseq1 = seqdb[s:s+read_len]
        read_seq = ffi.new("char[{}]".format(read_len))
        shimmer.decode_biseq(bseq1, read_seq , read_len, read_strand)

        aligned = False
        t_offset = 0
        if read_shift < 0:
            aln=falcon.align( read_seq[abs(read_shift):read_len], read_len - abs(read_shift)-5,
                              ref_seq, ref_len-5, 150, 1)

            if abs(abs(aln.aln_q_e-aln.aln_q_s)-(read_len - abs(read_shift))) < 48:
                aligned = True

                rng[0].s1 = aln.aln_q_s
                rng[0].e1 = aln.aln_q_e
                rng[0].s2 = aln.aln_t_s
                rng[0].e2 = aln.aln_t_e
                t_offset = 0


        else:
            aln=falcon.align( read_seq, read_len-5,
                              ref_seq[read_shift:ref_len], ref_len-read_shift-5, 150, 1)

            if abs(abs(aln.aln_q_e-aln.aln_q_s)-read_len) < 48 or \
            abs(ref_len-read_shift-abs(aln.aln_q_e-aln.aln_q_s)) < 48:
                aligned = True
                rng[0].s1 = aln.aln_q_s
                rng[0].e1 = aln.aln_q_e
                rng[0].s2 = aln.aln_t_s
                rng[0].e2 = aln.aln_t_e
                t_offset = read_shift
        if aligned:
            # print(ffi.string(aln.q_aln_str))
            # print(ffi.string(aln.t_aln_str))
            # print(rng[0].s1 , rng[0].e1, rng[0].s2, rng[0].e2 )
            tag=falcon.get_align_tags(aln.q_aln_str, aln.t_aln_str, aln.aln_str_size, rng, 0, t_offset)
            tags[aln_count] = tag
            aln_count+=1
        ffi.release(read_seq)
    #print(aln_count)
    cns = falcon.get_cns_from_align_tags(tags, aln_count, len(ref_seq), 0)
    cns_seq = ffi.string(cns.sequence)
    ##print(cns_seq)

    cns_segments.append(cns_seq)

    falcon.free_consensus_data(cns)
    for i in range(aln_count):
        falcon.free_align_tags(tags[i])
    ffi.release(tags)
    ffi.release(ref_seq)
ffi.release(rng)

s0 = cns_segments[0]
stiched_segments = [s0]
for s1 in cns_segments[1:]:
    aln=falcon.align( s0[-500:], 500,
                       s1[:500], 500, 100, 0)
    #print(aln.aln_q_s, aln.aln_q_e, aln.aln_t_s, aln.aln_t_e, aln.dist)
    stiched_segments.append( s1[aln.aln_t_e:] )
    s0 = s1

contig = b"".join(stiched_segments)
print(">test")
print(contig.decode("ascii"))
