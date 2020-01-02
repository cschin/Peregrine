#!/usr/bin/env pypy3

"""
TODO: (from convo w/ Ivan)
the issue with this script (but would still like to re-read it to refresh my memory). The script loads all edge sequences and tries to do two things at once: create p_ctg and a_ctg sequences, and align the bubbles using those sequences


If we generate:
1. All paths first (as tiling paths) for all p_ctg and all a_ctg without loading sequences - this should not consume much space (take a look at *_tiling_paths files).
2. Load the first read of each tiling path fully, and only edge sequences for every transition, we can generate the output sequences with the same memory/disk consumption.
3. Align bubbles after that.

Our resource consumption should be same

Bubbles?
It aligns them to produce the identity score

After that the dedup_a_tigs.py script is used to deduplicate fake a_ctg.
But that script is simple, and only depends on the alignment info that the previous script stored in the a_ctg header.
"""




#from builtins import zip
#from builtins import range
import argparse
import logging
import sys
import networkx as nx
import contextlib
#from pbcore.io import FastaReader
#from myio import open_progress

LOG = logging.getLogger()
RCMAP = dict(list(zip("ACGTacgtNn-", "TGCAtgcaNn-")))


@contextlib.contextmanager
def open_progress(fn, mode='r', log=LOG.info):
    """
    Usage:
        with open_progress('foo', log=LOG.info) as stream:
            for line in stream:
                use(line)

    That will log progress lines.
    """
    def get_iter(stream, progress):
        for line in stream:
            #progress(len(line))
            yield line

    fp = None
    with open(fn, mode=mode) as stream:
        yield get_iter(stream, fp)
    #fp.finish()

def log(msg):
    sys.stderr.write(msg)
    sys.stderr.write('\n')


def rc(seq):
    return "".join([RCMAP[c] for c in seq[::-1]])


def reverse_end(node_id):
    node_id, end = node_id.split(":")
    new_end = "B" if end == "E" else "E"
    return node_id + ":" + new_end


bit2base = ('N', 'A', 'C', 'N', 'G', 'N', 'N', 'N', 'T')
def decode_seq(encoded_seq, strand):
    ## encoded_seq = bytes(encoded_seq)
    if strand == 0: # ORIGIANL
        seq = "".join([ bit2base[ord(c) & 0x0F] for c in encoded_seq ])
    else: # REVERSED
        seq = "".join([ bit2base[ord(c) >> 4] for c in encoded_seq ])
    return seq


def yield_first_seq(one_path_edges, seqs, read_idx):
    if one_path_edges and one_path_edges[0][0] != one_path_edges[-1][1]:
        # If non-empty, and non-circular,
        # prepend the entire first read.
        (vv, ww) = one_path_edges[0]
        (vv_rid, vv_letter) = vv.split(":")
        s = read_idx[int(vv_rid)]["offset"]
        e = s + read_idx[int(vv_rid)]["length"]
        seq = decode_seq(seqs[s:e], 0)

        if vv_letter == 'E':
            first_seq = seq
        else:
            assert vv_letter == 'B'
            first_seq = "".join([RCMAP[c] for c in seq])
        yield first_seq

def compose_ctg(edge_data, ctg_id, path_edges, proper_ctg):
    total_score = 0
    total_length = 0
    edge_lines = []

    # Splice-in the rest of the path sequence.
    ctg_length = 0
    for vv, ww in path_edges:
        rid, s, t, aln_score, idt = edge_data[(vv, ww)]
        dl = abs(s-t)
        edge_lines.append('%s %s %s %s %d %d %d %0.2f %d %d' % (
            ctg_id, vv, ww, rid, s, t, aln_score, idt, ctg_length, dl))
        ctg_length += dl
        total_score += aln_score
    return edge_lines, total_score, ctg_length


def get_read_index(idx_fn):
    read_idx = {}
    with open(idx_fn) as f:
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
    return read_idx


def run(improper_p_ctg, proper_a_ctg, seqdb_prefix, sg_edges_list_fn, utg_data_fn, ctg_paths_fn):
    """improper==True => Neglect the initial read.
    We used to need that for unzip.
    """
    reads_in_layout = set()
    with open_progress(sg_edges_list_fn) as f:
        for l in f:
            l = l.strip().split()
            """001039799:E 000333411:E 000333411 17524 20167 17524 99.62 G"""
            v, w, rid, s, t, aln_score, idt, type_ = l
            if type_ != "G":
                continue
            r1 = v.split(":")[0]
            reads_in_layout.add(r1)
            r2 = w.split(":")[0]
            reads_in_layout.add(r2)

    edge_data = {}
    with open_progress(sg_edges_list_fn) as f:
        for l in f:
            l = l.strip().split()
            """001039799:E 000333411:E 000333411 17524 20167 17524 99.62 G"""
            v, w, rid, s, t, aln_score, idt, type_ = l

            if type_ != "G":
                continue
            r1, dir1 = v.split(":")
            reads_in_layout.add(r1) # redundant, but harmless
            r2, dir2 = w.split(":")
            reads_in_layout.add(r2) # redundant, but harmless

            s = int(s)
            t = int(t)
            aln_score = int(aln_score)
            idt = float(idt)

            if s < t:
                assert 'E' == dir2
            else:
                assert 'B' == dir2
            edge_data[(v, w)] = (rid, s, t, aln_score, idt)

    utg_data = {}
    with open_progress(utg_data_fn) as f:
        for l in f:
            l = l.strip().split()
            s, v, t, type_, length, score, path_or_edges = l
            if type_ not in ["compound", "simple", "contained"]:
                continue
            length = int(length)
            score = int(score)
            if type_ in ("simple", "contained"):
                path_or_edges = path_or_edges.split("~")
            else:
                path_or_edges = [tuple(e.split("~"))
                                 for e in path_or_edges.split("|")]
            utg_data[(s, v, t)] = type_, length, score, path_or_edges

    p_ctg_t_out = open("p_ctg_tiling_path", "w")
    a_ctg_t_out = open("a_ctg_tiling_path", "w")
    layout_ctg = set()

    with open_progress(ctg_paths_fn) as f:
        for row in f:
            row = row.strip().split()
            ctg_id, c_type_, i_utig, t0, length, score, utgs = row
            ctg_id = ctg_id
            s0 = i_utig.split("~")[0]

            if (reverse_end(t0), reverse_end(s0)) in layout_ctg:
                continue
            else:
                layout_ctg.add((s0, t0))

            length = int(length)
            utgs = utgs.split("|")
            one_path = []
            total_score = 0
            total_length = 0

            a_ctg_group = {}

            for utg in utgs:
                s, v, t = utg.split("~")
                type_, length, score, path_or_edges = utg_data[(s, v, t)]
                total_score += score
                total_length += length
                if type_ == "simple":
                    if len(one_path) != 0:
                        one_path.extend(path_or_edges[1:])
                    else:
                        one_path.extend(path_or_edges)
                if type_ == "compound":

                    c_graph = nx.DiGraph()

                    all_alt_path = []
                    for ss, vv, tt in path_or_edges:
                        type_, length, score, sub_path = utg_data[(ss, vv, tt)]

                        v1 = sub_path[0]
                        for v2 in sub_path[1:]:
                            c_graph.add_edge(
                                v1, v2, e_score=edge_data[(v1, v2)][3])
                            v1 = v2

                    shortest_path = nx.shortest_path(c_graph, s, t, "e_score")
                    score = nx.shortest_path_length(c_graph, s, t, "e_score")
                    all_alt_path.append((score, shortest_path))

                    # a_ctg_data.append( (s, t, shortest_path) ) #first path is the same as the one used in the primary contig
                    while 1:
                        n0 = shortest_path[0]
                        for n1 in shortest_path[1:]:
                            c_graph.remove_edge(n0, n1)
                            n0 = n1
                        try:
                            shortest_path = nx.shortest_path(
                                c_graph, s, t, "e_score")
                            score = nx.shortest_path_length(
                                c_graph, s, t, "e_score")
                            #a_ctg_data.append( (s, t, shortest_path) )
                            all_alt_path.append((score, shortest_path))

                        except nx.exception.NetworkXNoPath:
                            break
                        # if len(shortest_path) < 2:
                        #    break
                    # Is sorting required, if we are appending the shortest paths in order?
                    all_alt_path.sort()
                    all_alt_path.reverse()
                    shortest_path = all_alt_path[0][1]
                    # The longest branch in the compound unitig is added to the primary path.
                    if len(one_path) != 0:
                        one_path.extend(shortest_path[1:])
                    else:
                        one_path.extend(shortest_path)

                    a_ctg_group[(s, t)] = all_alt_path

            if len(one_path) == 0:
                continue

            one_path_edges = list(zip(one_path[:-1], one_path[1:]))

            # Compose the primary contig.
            p_edge_lines, p_total_score, p_total_length = compose_ctg(
                edge_data, ctg_id, one_path_edges, (not improper_p_ctg))

            # Write out the tiling path.
            p_ctg_t_out.write('\n'.join(p_edge_lines))
            p_ctg_t_out.write('\n')


            a_id = 0
            for v, w in a_ctg_group:
                atig_output = []

                # Compose the base sequence.
                for sub_id in range(0, len(a_ctg_group[(v, w)])):
                    if sub_id == 0:  # the alte path is used in the primary
                        continue
                    score, atig_path = a_ctg_group[(v, w)][sub_id]
                    atig_path_edges = list(zip(atig_path[:-1], atig_path[1:]))

                    a_ctg_id = '%s-%03d-%02d' % (ctg_id, a_id + 1, sub_id)
                    a_edge_lines, a_total_score, a_total_length = compose_ctg(
                        edge_data, a_ctg_id, atig_path_edges, proper_a_ctg)

                    # Keep the placeholder for these values for legacy purposes, but mark
                    # them as for deletion.
                    # The base a_ctg will also be output to the same file, for simplicity.
                    delta_len = 0
                    idt = 1.0
                    cov = 1.0
                    atig_output.append((v, w, atig_path, a_total_length,
                                        a_total_score, atig_path_edges,
                                        a_ctg_id, a_edge_lines,
                                        delta_len, idt, cov))

                for data in atig_output:
                    (v, w, tig_path, a_total_length, a_total_score,
                     tig_path_edges, a_ctg_id, a_edge_lines,
                     delta_len, a_idt, cov) = data

                    # Write out the tiling path.
                    a_ctg_t_out.write('\n'.join(a_edge_lines))
                    a_ctg_t_out.write('\n')

                a_id += 1

    a_ctg_t_out.close()
    p_ctg_t_out.close()


def main(argv=sys.argv):
    description = 'Generate the primary and alternate contig fasta files and tiling paths, given the string graph.'
    epilog = """
We write these:

    p_ctg_out = open("p_ctg.fa", "w")
    a_ctg_out = open("a_ctg_all.fa", "w")
    p_ctg_t_out = open("p_ctg_tiling_path", "w")
    a_ctg_t_out = open("a_ctg_tiling_path", "w")
"""
    parser = argparse.ArgumentParser(
            description=description,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog=epilog)
    parser.add_argument('--improper-p-ctg', action='store_true',
            help='Skip the initial read in each p_ctg path.')
    parser.add_argument('--proper-a-ctg', action='store_true',
            help='Skip the initial read in each a_ctg path.')
    parser.add_argument('--seqdb-prefix', type=str,
            default='./seq_dataset',
            help='Input. Preads file, required to construct the contigs.')
    parser.add_argument('--sg-edges-list-fn', type=str,
            default='./sg_edges_list',
            help='Input. File containing string graph edges, produced by ovlp_to_graph.py.')
    parser.add_argument('--utg-data-fn', type=str,
            default='./utg_data',
            help='Input. File containing unitig data, produced by ovlp_to_graph.py.')
    parser.add_argument('--ctg-paths-fn', type=str,
            default='./ctg_paths',
            help='Input. File containing contig paths, produced by ovlp_to_graph.py.')
    args = parser.parse_args(argv[1:])

    run(**vars(args))

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main(sys.argv)
