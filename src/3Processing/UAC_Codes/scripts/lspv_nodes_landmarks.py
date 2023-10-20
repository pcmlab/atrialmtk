#!/usr/bin/python3
#
import os
import sys
import argparse
import numpy as np
from sklearn.neighbors import NearestNeighbors
from uac import utils


# Create the parser
my_parser = argparse.ArgumentParser(description="Calculating PV nodes")

# Add the arguments
my_parser.add_argument("target_path", metavar="path", type=str, help="Path for the new target mesh")
my_parser.add_argument("mesh_name", metavar="filename", type=str, help="Mesh name (usually Labelled)")
my_parser.add_argument("seed_landmarks", metavar="filename", type=str, help="Seed file name (e.g. Landmarks.txt)")
my_parser.add_argument("seed_region", metavar="filename", type=str, help="Seed file name (e.g. Region.txt)")
my_parser.add_argument("scale_factor", metavar="scale factor", type=int, help="Scale factor for landmarks")

# Execute parse_args()
args = my_parser.parse_args()

base_dir = args.target_path
mesh_name = args.mesh_name
seed_landmarks = args.seed_landmarks
seed_region = args.seed_region
scale_factor = args.scale_factor

if not os.path.isdir(base_dir):
    print("The target mesh path specified does not exist")
    sys.exit()

sys.path.append(base_dir)


print(base_dir)
pts, elems, *_, surface = utils.read_carp(base_dir, mesh_name)

pts_src = np.loadtxt(base_dir + seed_region, delimiter=",") * scale_factor

# find nodes at PV/LAA
size = 5
nodes_pv = [None] * size
pts_pv = [None] * size
for i in range(size):
    *_, nodes_pv[i] = utils.find_nodes(elems, pts, i, True)

    pts_pv[i] = np.empty((len(nodes_pv[i]), 3), float)
    for ind in range(len(nodes_pv[i])):
        pts_pv[i][ind] = pts[nodes_pv[i][ind]]

pts_src_remap = [pts_src[0], pts_src[1], pts_src[2], pts_src[3], pts_src[5]]
sp = [None] * size
for i in range(size):
    sp[i] = utils.pv_index_remap(pts_src_remap[i], pts_pv, np.mean)

pv_choices = [None] * size
for i in range(size):
    pv_choices[i] = sp.index(i)

nodes_rspv = nodes_pv[pv_choices.index(0)]
nodes_ripv = nodes_pv[pv_choices.index(1)]
nodes_lipv = nodes_pv[pv_choices.index(2)]
nodes_lspv = nodes_pv[pv_choices.index(3)]
nodes_laa = nodes_pv[pv_choices.index(4)]



# find MV nodes
utils.write_vtx_intra(nodes_lspv, base_dir + "LSPV_nodes.vtx")
utils.write_vtx_intra(nodes_lipv, base_dir + "LIPV_nodes.vtx")
utils.write_vtx_intra(nodes_rspv, base_dir + "RSPV_nodes.vtx")
utils.write_vtx_intra(nodes_ripv, base_dir + "RIPV_nodes.vtx")

