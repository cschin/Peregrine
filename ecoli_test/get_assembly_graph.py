import networkx as nx


def find_best_edges(edges):
    best_edge = None
    best_weight = 0
    for v, w in edges:
        weight=G.edges[v,w]["weight"]
        if weight > best_weight:
            best_edge = v, w
            best_weight = weight
    return best_edge, best_weight

mdata = {}
with open("L2.summary.all.2") as f:
    for row in f:
        row = row.strip().split()
        m1, m2, rid, pos, _ = (int(c) for c in row)
        k = m1, m2
        mdata.setdefault(k, [])
        mdata[k].append( (rid, pos, _) )

read_length = []
with open("seq_db.idx") as f:
    for row in f:
        row = row.strip().split()
        read_length.append(int(row[2]))

contained = set()
for k in mdata:
    v = mdata[k]
    if len(v) < 2 or len(v) > 45:
        continue

    v.sort(key=lambda x:-x[1])
    v0 = v[0]
    rid = v0[0]//2
    r_pos = read_length[rid] - v0[1]

    for v1 in v[1:]:
        rid = v1[0]//2
        r_pos2 = read_length[rid] - v1[1]
        if  r_pos2 < r_pos:
            contained.add(rid)
        else:
            r_pos = r_pos2

G = nx.DiGraph()

def reverse_end(s):
    s,e = s.split(":")
    e = "B" if e == "E" else "E"
    return s + ":" + e

def get_read_id(s):
    rid0, d0 = s // 2, s % 2
    e0 = "E" if d0 == 0 else "B"
    n0 = "{:09d}:{:s}".format(rid0, e0)
    return n0


for k in mdata:
    v = mdata[k]
    if len(v) < 2 or len(v) > 45:
        continue

    v.sort(key=lambda x:-x[1])
    v = [ x for x in v if x[0]//2 not in contained ]
    if len(v) < 2:
        continue
    v0 = v[0]

    for v1 in v[1:]:
        rid0 = get_read_id(v0[0])
        rid1 = get_read_id(v1[0])
        rrid0 = reverse_end(rid0)
        rrid1 = reverse_end(rid1)
        r_len0 = read_length[v0[0]//2]
        r_len1 = read_length[v1[0]//2]

        if G.has_edge(rid0, rid1):
            G.edges[rid0, rid1]["weight"] += 1
            G.edges[rrid1, rrid0]["weight"] += 1
        else:
            length = r_len0 - abs( v0[1] - v1[1])
            score = length
            G.add_edge(rid0, rid1, weight=1, length=length, score=score)
            G.add_edge(rrid1, rrid0, weight=1, length=length, score=score)
        v0 = v1


subG = nx.DiGraph()
for n in list(G.nodes()):

    best_edge, best_weight = find_best_edges(G.out_edges(n))
    if best_edge is not None:
        length = G.edges[best_edge[0], best_edge[1]]["length"]
        score = G.edges[best_edge[0], best_edge[1]]["score"]
        subG.add_edge(best_edge[0], best_edge[1], weight=best_weight, length=length, score=score)

    best_edge, best_weight = find_best_edges(G.in_edges(n))
    if best_edge is not None:
        length = G.edges[best_edge[0], best_edge[1]]["length"]
        score = G.edges[best_edge[0], best_edge[1]]["score"]
        subG.add_edge(best_edge[0], best_edge[1], weight=best_weight, length=length, score=score)

nodes = set(subG.nodes())
for n in list(nodes):
    n2 = reverse_end(n)
    if n2 in nodes:
        subG.add_edge(n, n2)

nx.write_gexf(subG,"test2.gexf")
