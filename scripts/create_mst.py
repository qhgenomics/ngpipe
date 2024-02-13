import sys

import networkx as nx
from matplotlib import pyplot as plt



max_branch = 500

G = nx.Graph()

merged = {}


with open(snakemake.input.summary) as f:
    header = f.readline()
    node_dict = {}
    dist_dict = {}
    profile_dict = {}
    for line in f:
        name, mlst, a1, a2, a3, a4, a5, a6, a7 = line.split("\t")[:9]
        if mlst == '-':
            mlstname_label = (a1, a2, a3, a4, a5, a6, a7)
            mlstname = "Unknown"
        else:
            mlstname_label = mlst
            mlstname = "ST" + mlst
        if not mlstname_label in node_dict:
            node_dict[mlstname_label] = mlstname
            dist_dict[mlstname_label] = {}
            for i in profile_dict:
                dist = 0
                for j, k in zip([a1, a2, a3, a4, a5, a6, a7], profile_dict[i]):
                    if j != k:
                        dist += 1
                dist_dict[mlstname_label][i] = dist
            profile_dict[mlstname_label] = [a1, a2, a3, a4, a5, a6, a7]
        node_dict[mlstname_label] += "\n{}".format(name)


print(dist_dict)
for i in dist_dict:
    for j in dist_dict[i]:
        name = node_dict[i]
        hname = node_dict[j]
        G.add_edge(name, hname, weight=int(dist_dict[i][j]), length=(50 + dist_dict[i][j]))





T = nx.minimum_spanning_tree(G, weight="weight")


sizes = []
for n in list(T.nodes):
    sizes.append(20 * len(n.split("\n")))

pos = nx.kamada_kawai_layout(T, weight="length")  # positions for all nodes - seed for reproducibility


node_label_size=3
edge_label_size=4

elarge = [(u, v) for (u, v, d) in T.edges(data=True) if d["weight"] > 2]
esmall = [(u, v) for (u, v, d) in T.edges(data=True) if d["weight"] <= 2]

# nodes
nx.draw_networkx_nodes(T, pos, node_size=sizes)

nx.draw_networkx_edges(T, pos, edgelist=esmall, width=1)
nx.draw_networkx_edges(
    T, pos, edgelist=elarge, width=1, alpha=0.5, edge_color="b", style="dashed"
)

# node labels
nx.draw_networkx_labels(T, pos, font_size=node_label_size)
# edge weight labels
edge_labels = nx.get_edge_attributes(T, "weight")
nx.draw_networkx_edge_labels(T, pos, edge_labels, font_size=edge_label_size)


ax = plt.gca()
ax.margins(0.08)
plt.axis("off")
plt.rcParams['svg.fonttype'] = 'none'
plt.tight_layout()
plt.savefig(snakemake.output.mst)