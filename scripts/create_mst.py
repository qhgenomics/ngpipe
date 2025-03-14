import sys, os

import networkx as nx
from matplotlib import pyplot as plt
import glob




def colorstr(rgb):
    return "#%02x%02x%02x" % (rgb[0],rgb[1],rgb[2])

def add_profiles(profile_file, add_names=False):
    global cl_index
    with open(profile_file) as f:
        header = f.readline()
        for line in f:
            name, ngstar, a1, a2, a3, a4, a5, a6, a7, ngstarCC, pen_mosaic, mlst, b1, b2, b3, b4, b5, b6, b7, \
                ngmast, c1, c2 = line.split("\t")[:22]
            profile = (a1, a2, a3, a4, a5, a6, a7, b1, b2, b3, b4, b5, b6, b7, c1, c2)
            missing = profile.count('-')
            if missing > snakemake.params.max_missing:
                continue
            if not profile in size_dict:
                size_dict[profile] = 0
            if not profile in color_dict:
                if not mlst in mlst_colors:
                    mlst_colors[mlst] = colorstr(colorlist[cl_index%len(colorlist)])
                    cl_index += 1
                color_dict[profile] = mlst_colors[mlst]
            size_dict[profile] += 1
            if mlst == '-':
                mlst = "unknown"
            if ngstar == '-':
                ngstar = "unknown"
            if ngmast == '-':
                ngmast = "unknown"
            if not profile in label_dict:
                label_dict[profile] = "MLST:{}\nNGSTAR:{}\nNGMAST:{}".format(mlst, ngstar, ngmast)
            if add_names:
                label_dict[profile] += "\n{}".format(name)
            if not profile in dist_dict:
                dist_dict[profile] = {}
                for i in dist_dict:
                    if i != profile:
                        dist = 0
                        for j, k in zip(profile, i):
                            if j != k:
                                dist += 1
                        dist_dict[profile][i] = dist


colorlist = [(240, 163, 255), (0, 117, 220), (153, 63, 0), (76, 0, 92), (25, 25, 25), (0, 92, 49), (43, 206, 72),
              (255, 204, 153), (128, 128, 128), (148, 255, 181), (143, 124, 0), (157, 204, 0), (194, 0, 136), (0, 51, 128),
              (255, 164, 5), (255, 168, 187), (66, 102, 0), (255, 0, 16), (94, 241, 242), (0, 153, 143), (224, 255, 102),
              (116, 10, 255), (153, 0, 0), (255, 255, 128), (255, 255, 0), (255, 80, 5), (0, 0, 0), (50, 50, 50)]

dist_dict = {}
label_dict = {}
size_dict = {}
color_dict = {}
cl_index = 0
mlst_colors = {}



add_profiles(snakemake.input.summary, True)


for summary in glob.glob(os.path.join(snakemake.params.previous_runs, "*.tsv")):
    add_profiles(summary)




max_branch = 500
G = nx.Graph()

for i in dist_dict:
    for j in dist_dict[i]:
        G.add_edge(i, j, weight=int(dist_dict[i][j]), length=(50 + dist_dict[i][j]))

if len (dist_dict) == 1:
     G.add_node(list(dist_dict)[0])



T = nx.minimum_spanning_tree(G, weight="weight")


sizes = []
colors = []
for n in list(T.nodes):
    sizes.append(20 * size_dict[n])
    colors.append(color_dict[n])

pos = nx.kamada_kawai_layout(T, weight="length")  # positions for all nodes - seed for reproducibility


node_label_size=3
edge_label_size=4

elarge = [(u, v) for (u, v, d) in T.edges(data=True) if d["weight"] > 2]
esmall = [(u, v) for (u, v, d) in T.edges(data=True) if d["weight"] <= 2]

# nodes
nx.draw_networkx_nodes(T, pos, node_size=sizes, node_color=colors)

nx.draw_networkx_edges(T, pos, edgelist=esmall, width=1)
nx.draw_networkx_edges(
    T, pos, edgelist=elarge, width=1, alpha=0.5, edge_color="b", style="dashed"
)

# node labels
nx.draw_networkx_labels(T, pos, label_dict, font_size=node_label_size)
# edge weight labels
edge_labels = nx.get_edge_attributes(T, "weight")
nx.draw_networkx_edge_labels(T, pos, edge_labels, font_size=edge_label_size)


ax = plt.gca()
ax.margins(0.08)
plt.axis("off")
plt.rcParams['svg.fonttype'] = 'none'
plt.tight_layout()
plt.savefig(snakemake.output.mst)
