import mmap
import sys
from _shimmer4py import ffi, lib

if __name__ == "__main__":
    seqdb_prefix = sys.argv[1]
    tiling_path_fn = sys.argv[2]

    f = open("{}.seqdb".format(seqdb_prefix), "rb")
    seqdb = mmap.mmap(f.fileno(), 0,
                      flags=mmap.MAP_SHARED, prot=mmap.PROT_READ)

    read_idx = {}
    with open("{}.idx".format(seqdb_prefix)) as f:
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

    tiling_path_data = {}
    with open(tiling_path_fn) as f:
        for row in f:
            row = row.strip().split()
            tiling_path_data.setdefault(row[0], [])
            tiling_path_data[row[0]].append(row)

    sub_seq = []
    for ctg in tiling_path_data:
        segments = []
        # I don't like to have the first read as it breaks the string formulation,
        # but poeple like it for no reason, so I will just do it
        ctg_id, v, w, r, s, e, olen, idt = tiling_path_data[ctg][0]
        v = v.split(":")
        rid0 = int(v[0])
        s0 = read_idx[rid0]["offset"]
        slen0 = read_idx[rid0]["length"]
        e0 = s0 + slen0
        bseq0 = seqdb[s0:e0]
        strand0 = 0 if v[1] == "E" else 1

        seq = ffi.new("char[{}]".format(slen0))
        lib.decode_4bit_bidirection(bseq0, seq, slen0, strand0)
        segments.append(seq)
        for row in tiling_path_data[ctg][1:]:
            ctg_id, v, w, r, s, e, olen, idt = row
            v = v.split(":")
            w = w.split(":")
            s = int(s)
            e = int(e)
            olen = int(olen)
            idt = float(idt)

            rid0 = int(v[0])
            s0 = read_idx[rid0]["offset"]
            slen0 = read_idx[rid0]["length"]
            e0 = s0 + slen0
            bseq0 = seqdb[s0:e0]
            strand0 = 0 if v[1] == "E" else 1

            rid1 = int(w[0])
            s1 = read_idx[rid1]["offset"]
            slen1 = read_idx[rid1]["length"]
            e1 = s1 + slen1
            bseq1 = seqdb[s1:e1]
            strand1 = 0 if w[1] == "E" else 1

            offset1 = slen0 - 500
            offset2 = slen1 - abs(e-s) - 500
            match = lib.ovlp_match(bseq0[offset1:], slen0 - offset1, strand0,
                                   bseq1[offset2:], slen1 - offset2, strand1,
                                   500)
            s -= match.t_end - match.q_end
            lib.free_ovlp_match(match)

            seq = ffi.new("char[{}]".format(abs(e-s)))
            if strand1 == 1:
                s, e = slen1 - s, slen1 - e
            if e > s:
                lib.decode_4bit_bidirection(bseq1[s:e], seq, e-s, strand1)
                segments.append(seq)
            else:
                seq = ffi.new("char[1]", "N")
                segments.append(seq)
                print(row, file=sys.stderr)
        print(">{}".format(ctg_id))
        print(b"".join([ffi.buffer(seq) for seq in segments]).decode("ascii"))
        for seq in segments:
            ffi.release(seq)
