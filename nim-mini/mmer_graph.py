import networkx as nx


G = nx.DiGraph()
m_count = {}

#fn = "preads4falcon_mer"

fn = "H08_mer"
with open(fn) as f:
    for row in f:
        row = row.strip()
        if row[0] == ">":
            continue
        row = row.split()
        m_count.setdefault(row[2], 0)
        m_count[row[2]] += 1


with open(fn) as f:
    for row in f:
        row = row.strip()
        if row[0] == ">":
            v = None
            w = None
        else:
            row = row.split()
            if v is not None:
                w = row[2]
                if m_count[v] > 5 and m_count[v] < 60 and \
                   m_count[w] > 5 and m_count[w] < 60:
                    G.add_edge(v, w)
                    if "count" not in G[v][w]:
                        G[v][w]["count"] = 0
                    G[v][w]["count"] += 1
                v = w
            else:
                v = row[2]

#for v, w in G.edges():
#    print(v, w, G[v][w]["count"], G.out_degree(v), G.in_degree(w))

remove_nodes = set()
for v in G.nodes():
    if G.out_degree(v) > 1 or G.in_degree(v) > 1:
        remove_nodes.add(v)

G2 = G.copy()
for v in list(remove_nodes):
    G2.remove_node(v)

remove_nodes = set()
for subG in nx.weakly_connected_component_subgraphs(G2):
    if len(subG) == 1:
        remove_nodes.update(subG.nodes())

for v in list(remove_nodes):
    G.remove_node(v)

remove_nodes = set()
for v in G.nodes():
    if G.out_degree(v) == 0 or G.in_degree(v) == 0:
        remove_nodes.add(v)

for v in list(remove_nodes):
    G.remove_node(v)

for subG in nx.weakly_connected_component_subgraphs(G):
    subG_size = len(subG.nodes())
    for v in subG.nodes():
        print(subG_size, v, subG.in_degree(v), subG.out_degree(v))


nx.write_gexf(G, "test.gexf")
