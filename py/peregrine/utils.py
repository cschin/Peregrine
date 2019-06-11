import peregrine
import sys, os
import numpy as np
from peregrine._shimmer4py import ffi as shimmer_ffi
from peregrine._shimmer4py import lib as shimmer4py
from peregrine._falcon4py import ffi as falcon_ffi
from peregrine._falcon4py import lib as falcon4py


rmap = dict(list(zip(b"ACGT", b"TGCA")))


def rc(seq):
    return bytes([rmap[c] for c in seq[::-1]])


def mmer2tuple(mmer):
    x = mmer.x
    y = mmer.y
    span = x & 0xFF
    mmer = x >> 8
    rid = y >> 32
    pos_end = ((y & 0xFFFFFFFF) >> 1) + 1
    direction = y & 0x1
    return (mmer, span, rid, pos_end, direction)


def get_shimmers_from_seq(seq, rid=0,
                          levels=2, reduction_factor=3,
                          k=16, w=80):
    assert levels <= 2
    c_null = shimmer_ffi.NULL
    mmers = shimmer_ffi.new("mm128_v *")
    shimmer4py.mm_sketch(c_null, seq, len(seq), w, k, rid, 0, mmers)
    if levels == 0:
        return mmers
    elif levels == 1:
        mmers_L1 = shimmer_ffi.new("mm128_v *")
        shimmer4py.mm_reduce(mmers, mmers_L1, reduction_factor)
        shimmer_ffi.release(mmers)
        return mmers_L1
    elif levels == 2:
        mmers_L1 = shimmer_ffi.new("mm128_v *")
        mmers_L2 = shimmer_ffi.new("mm128_v *")
        shimmer4py.mm_reduce(mmers, mmers_L1, reduction_factor)
        shimmer4py.mm_reduce(mmers_L1, mmers_L2, reduction_factor)
        shimmer_ffi.release(mmers_L1)
        shimmer_ffi.release(mmers)
        return mmers_L2


def get_shimmer_alns(shimmers0, shimmers1, direction=0,
                     max_diff=100, max_dist=1200,
                     max_repeat=1):
    aln = shimmer4py.shmr_aln(shimmers0, shimmers1, direction,
                              max_diff, max_dist, max_repeat)
    aln_chains = []
    for i in range(aln.n):
        chain = []
        offsets = np.zeros(aln.a[i].idx0.n, dtype=np.float)
        for j in range(aln.a[i].idx0.n):
            idx0, idx1 = aln.a[i].idx0.a[j], aln.a[i].idx1.a[j]
            mmer0 = mmer2tuple(shimmers0.a[idx0])
            mmer1 = mmer2tuple(shimmers1.a[idx1])
            chain.append((mmer0, mmer1))
            if direction == 0:  # same direction
                d = mmer0[3] - mmer1[3]
            else:
                d = mmer0[3] + mmer1[3]
            offsets[j] = d
        aln_chains.append((chain, np.max(d), np.mean(d), np.min(d)))
    shimmer4py.free_shmr_alns(aln)
    return aln_chains

def get_tag_from_seqs(read_seq, ref_seq, read_offset, max_dist=150):
    rng = falcon_ffi.new("aln_range[1]")
    read_len = len(read_seq)
    ref_len = len(ref_seq)
    aligned = False
    if read_offset < 0:
        aln = falcon4py.align(read_seq[abs(read_offset):read_len],
                              read_len - abs(read_offset),
                              ref_seq,
                              len(ref_seq),
                              max_dist, 1)
        if abs(abs(aln.aln_q_e-aln.aln_q_s) -
                (read_len - abs(read_offset))) < 48:
            aligned = True
            rng[0].s1 = aln.aln_q_s
            rng[0].e1 = aln.aln_q_e
            rng[0].s2 = aln.aln_t_s
            rng[0].e2 = aln.aln_t_e
            t_offset = 0
        else:
            falcon4py.free_alignment(aln)
    else:
        aln = falcon4py.align(read_seq,
                              read_len,
                              ref_seq[read_offset:ref_len],
                              ref_len-read_offset,
                              max_dist, 1)
        if abs(abs(aln.aln_q_e-aln.aln_q_s)-read_len) < 48 or \
            abs(ref_len-read_offset-abs(aln.aln_q_e-aln.aln_q_s)) < 48:
            aligned = True
            rng[0].s1 = aln.aln_q_s
            rng[0].e1 = aln.aln_q_e
            rng[0].s2 = aln.aln_t_s
            rng[0].e2 = aln.aln_t_e
            t_offset = read_offset
        else:
            falcon4py.free_alignment(aln)
    tag = None
    if aligned:
        tag = falcon4py.get_align_tags(aln.q_aln_str,
                                       aln.t_aln_str,
                                       aln.aln_str_size,
                                       rng, 0, t_offset)
        falcon4py.free_alignment(aln)
    falcon_ffi.release(rng)

    return tag


def get_cns_from_reads(seqs, levels=2, k=16, w=80, max_dist=150):

    aln_count = 0
    tags = falcon_ffi.new("align_tags_t * [{}]".format(len(seqs)+1))
    seq0 = seqs[0]
    shimmers0 = get_shimmers_from_seq(seq0,
                                      rid=0,
                                      levels=levels,
                                      k=k, w=w)
    alns = get_shimmer_alns(shimmers0, shimmers0, 0)
    aln = alns[0]
    read_offset = aln[0][0][0][3] - aln[0][0][1][3]
    seq = seq0
    tag = get_tag_from_seqs(seq, seq0, read_offset)
    tags[aln_count] = tag
    aln_count += 1

    for i, seq in enumerate(seqs):
        if i == 0:
            continue
        rid = i * 2
        shimmers1 = get_shimmers_from_seq(seq,
                                          rid=rid,
                                          levels=levels,
                                          k=k, w=w)
        alns_0 = get_shimmer_alns(shimmers0, shimmers1, 0)
        alns_0.sort(key=lambda x: -len(x[0]))
        shimmer4py.free(shimmers1.a)
        shimmer_ffi.release(shimmers1)

        rseq = rc(seq)
        rid = i * 2 + 1
        shimmers1 = get_shimmers_from_seq(rseq,
                                          rid=rid,
                                          levels=levels,
                                          k=k, w=w)
        alns_1 = get_shimmer_alns(shimmers0, shimmers1, 0)
        alns_1.sort(key=lambda x: -len(x[0]))
        shimmer4py.free(shimmers1.a)
        shimmer_ffi.release(shimmers1)

        if len(alns_0) > 0 and len(alns_1) > 0:
            if len(alns_0[0][0]) > len(alns_1[0][0]):
                aln = alns_0[0]
            else:
                aln = alns_1[0]
                seq = rseq
        elif len(alns_0) > 0:
            aln = alns_0[0]
        elif len(alns_1) > 0:
            aln = alns_1[0]
            seq = rseq
        else:
            continue

        read_offset = aln[0][0][0][3] - aln[0][0][1][3]
        tag = get_tag_from_seqs(seq, seq0, read_offset,
                                max_dist=max_dist)
        if tag is not None:
            tags[aln_count] = tag
            aln_count += 1

    cns = falcon4py.get_cns_from_align_tags(tags,
                                            aln_count,
                                            len(seq0), 1)
    cns_seq = falcon_ffi.string(cns.sequence)
    falcon4py.free_consensus_data(cns)
    shimmer4py.free(shimmers0.a)
    shimmer_ffi.release(shimmers0)
    for i in range(aln_count):
        falcon4py.free_align_tags(tags[i])
    falcon_ffi.release(tags)

    return cns_seq

