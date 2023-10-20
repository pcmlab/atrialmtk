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
threshold_value = args.threshold_value
threshold_value_laa = args.threshold_value_laa

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
pts, elems, *_, surface = utils.read_carp(base_dir, mesh_name)

mv_pv1 = utils.read_array_igb(base_dir + "MV_PV1/phie.igb")[0]
mv_pv2 = utils.read_array_igb(base_dir + "MV_PV2/phie.igb")[0]
mv_pv3 = utils.read_array_igb(base_dir + "MV_PV3/phie.igb")[0]
mv_pv4 = utils.read_array_igb(base_dir + "MV_PV4/phie.igb")[0]
mv_laa = utils.read_array_igb(base_dir + "MV_LAA/phie.igb")[0]

node_inds = np.ones(len(pts), int) * 11
node_inds[mv_pv1 > threshold_value] = 21
node_inds[mv_pv2 > threshold_value] = 23
node_inds[mv_pv3 > threshold_value] = 25
node_inds[mv_pv4 > threshold_value] = 27
node_inds[mv_laa > threshold_value_laa] = 13
utils.write_dat_simple(node_inds, "Labels.dat")

# move to elements and write a mesh
regs = np.empty(len(elems), int)
for i in range(len(elems)):
    regs[i] = node_inds[min(elems[i].n)]

surface = utils.add_celldata_array(surface, "T")
utils.write_surface_to_carp_mesh(surface, base_dir + "Labelled_Labels", labels=regs)
surface = utils.add_labels(surface, "T", regs)
utils.write_surface_to_vtk(surface, base_dir + "T_labels_2.vtk", write_ascii=True)
