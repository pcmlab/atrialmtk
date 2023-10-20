#!/usr/bin/python3
#
import os
import sys
import argparse
import numpy as np
from uac import utils


# Create the parser
my_parser = argparse.ArgumentParser(description="Calculation of labels using laplace solves only")

# Add the arguments
my_parser.add_argument("target_path", metavar="path", type=str, help="Path for the new target mesh")
my_parser.add_argument("fibre_path", metavar="path", type=str, help="Path for the atlas mesh")
my_parser.add_argument("laplace_path", metavar="path", type=str, help="Path for the laplace par files")
my_parser.add_argument("mesh_name", metavar="filename", type=str, help="Mesh name before labelling")
my_parser.add_argument("threshold_value", metavar="float", type=float, help="Threshold for PV e.g. 0.5")
my_parser.add_argument("threshold_value_laa", metavar="float", type=float, help="Threshold for LAA e.g. 0.5")
my_parser.add_argument("region_name", metavar="filename", type=str, help="Region file name (e.g. Regions.txt)")
my_parser.add_argument("scale_factor", metavar="scale factor", type=int, help="Scale factor for landmarks")

# Execute parse_args()
args = my_parser.parse_args()

base_dir = args.target_path
laplace_dir = args.laplace_path
fibre_dir = args.fibre_path
mesh_name = args.mesh_name
region_name = args.region_name
scale_factor = args.scale_factor

if not os.path.isdir(base_dir):
    print("The target mesh path specified does not exist")
    sys.exit()

if not os.path.isdir(laplace_dir):
    print("The atlas path specified does not exist")
    sys.exit()

if not os.path.isdir(fibre_dir):
    print("The laplace files path specified does not exist")
    sys.exit()


print(base_dir)
pts, elems, *_ = utils.read_carp(base_dir, mesh_name, return_surface=False)

pts_src = utils.load_landmarks(base_dir + region_name, scale_factor)
laa_marker = pts_src[4]
nodes_laa = np.argsort(utils.closest_node(laa_marker, pts))[:50]
print(laa_marker)
print(nodes_laa)

size = 3
nodes_pv = [None] * size
pts_pv = [None] * size
for i in range(size):
    *_, nodes_pv[i] = utils.find_nodes(elems, pts, i + 1)

    pts_pv[i] = np.empty((len(nodes_pv[i]), 3), float)
    for ind in range(len(nodes_pv[i])):
        pts_pv[i][ind] = pts[nodes_pv[i][ind]]

pts_src_remap = [pts_src[1], pts_src[3], pts_src[0]]
sp = [None] * size
for i in range(size):
    sp[i] = utils.pv_index_remap(pts_src_remap[i], pts_pv, np.mean)

pv_choices = [None] * size
for i in range(size):
    pv_choices[i] = sp.index(i)

nodes_rspv = nodes_pv[pv_choices.index(0)]
nodes_ripv = nodes_pv[pv_choices.index(1)]
nodes_lipv = nodes_pv[pv_choices.index(2)]

# find MV nodes
*_, nodes_mv = utils.find_nodes(elems, pts)

utils.write_vtx_intra(nodes_mv, base_dir + "MV.vtx")
utils.write_vtx_intra(nodes_lipv, base_dir + "PV2.vtx")
utils.write_vtx_intra(nodes_rspv, base_dir + "PV3.vtx")
utils.write_vtx_intra(nodes_ripv, base_dir + "PV4.vtx")
utils.write_vtx_intra(nodes_laa, base_dir + "LAA.vtx")
