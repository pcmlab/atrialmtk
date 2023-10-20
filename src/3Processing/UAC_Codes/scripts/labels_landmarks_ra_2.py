#!/usr/bin/python3
#
import os
import sys
import argparse
import shutil
import numpy as np
from uac import utils


# Create the parser
my_parser = argparse.ArgumentParser(description="Calculation of labels using laplace solves and landmarks")

# Add the arguments
my_parser.add_argument("target_path", metavar="path", type=str, help="Path for the new target mesh")
my_parser.add_argument("mesh_name", metavar="filename", type=str, help="Mesh name before labelling")
my_parser.add_argument("seed_landmarks", metavar="filename", type=str, help="Seed file name (e.g. Landmarks.txt)")
my_parser.add_argument("seed_region", metavar="filename", type=str, help="Seed file name (e.g. Region.txt)")
my_parser.add_argument("threshold_value", metavar="float", type=float, help="Threshold delta e.g 0.2")
my_parser.add_argument("scale_factor", metavar="scale factor", type=int, help="Scale factor for landmarks")

# Execute parse_args()
args = my_parser.parse_args()

base_dir = args.target_path
mesh_name = args.mesh_name
seed_region = args.seed_region
delta = args.threshold_value
scale_factor = args.scale_factor

if not os.path.isdir(base_dir):
    print("The target mesh path specified does not exist")
    sys.exit()


print(base_dir)
print((base_dir + mesh_name))

pts, elems, *_, surface = utils.read_carp(base_dir, mesh_name)

print(pts)

pts_region = np.loadtxt(base_dir + seed_region, delimiter=",") * scale_factor

size = 3
nodes_pv = [None] * size
pts_pv = [None] * size
for i in range(size):
    *_, nodes_pv[i] = utils.find_nodes(elems, pts, i)

    pts_pv[i] = np.empty((len(nodes_pv[i]), 3), float)
    for ind in range(len(nodes_pv[i])):
        pts_pv[i][ind] = pts[nodes_pv[i][ind]]

pts_src_remap = [pts_region[0], pts_region[3]]
print(pts_src_remap)
size = len(pts_src_remap)
sp = [None] * size
for i in range(size):
    sp[i] = utils.pv_index_remap(pts_src_remap[i], pts_pv, np.mean)

print(sp)

nodes_rspv = nodes_pv[sp[0]]
nodes_lspv = nodes_pv[sp[1]]

print(len(nodes_rspv))
print(len(nodes_lspv))


# set MV as other choice...
considered_choices = [0, 1, 2]
mv_choice = [x for x in considered_choices if x not in sp]
print(mv_choice)


# now labelling for RAA

print(base_dir)

laa_marker = pts_region[4]
nodes_laa = np.argsort(utils.closest_node(laa_marker, pts))[:50]
print(laa_marker)
print(nodes_laa)

laa_marker_base = pts_region[5]
nodes_laa_base = np.argsort(utils.closest_node(laa_marker_base, pts))[0]
print(nodes_laa_base)


# find MV nodes
*_, nodes_mv = utils.find_nodes(elems, pts)

utils.write_vtx_extra(nodes_mv, base_dir + "MV.vtx")
utils.write_vtx_extra(nodes_laa, base_dir + "LAA.vtx")


mv_laa = utils.read_array_igb(base_dir + "MV_LAA/phie.igb")[0]

# combine to allow thresholding of each PV
scaled5 = 1 * mv_laa

utils.write_dat_simple(scaled5, "Scaled5.dat")

# LAA use region landmark
threshold_value_laa_landmark = scaled5[nodes_laa_base]
print(threshold_value_laa_landmark)

node_inds = np.ones(len(pts))
node_inds[scaled5 > delta] = 2
utils.write_dat_simple(node_inds, "Labels.dat")

# move to elements and write a mesh
regs = np.empty(len(elems), int)
for i in range(len(elems)):
    regs[i] = node_inds[min(elems[i].n)]

surface = utils.add_celldata_array(surface, "T")
utils.write_surface_to_carp_mesh(surface, base_dir + "RA_only_RAA", labels=regs)
surface = utils.add_labels(surface, "T", regs)
utils.write_surface_to_vtk(surface, base_dir + "T_labels_2.vtk", write_ascii=True)
shutil.copyfile(base_dir + mesh_name + ".lon", base_dir + "RA_only_RAA.lon")
