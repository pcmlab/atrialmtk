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

if not os.path.isdir(base_dir):
    print("The target mesh path specified does not exist")
    sys.exit()


# read carp, write to vtk file
print(base_dir)
print((base_dir + mesh_name))

pts, elems, *_, surface = utils.read_carp(base_dir, mesh_name)

print(pts)

# move to elements and write a mesh

NodeIndexList = np.ones(len(pts)) * 1
regs = []
for e in elems:
    regs.append(NodeIndexList[min(e.n)])

surface = utils.add_celldata_array(surface, "T")
utils.write_surface_to_carp_mesh(surface, base_dir + "RA_only_RAA", labels=regs)

surface = utils.add_labels(surface, "T", regs)
utils.write_surface_to_vtk(surface, base_dir + "T_labels_2.vtk", write_ascii=True)

shutil.copyfile(base_dir + mesh_name + ".lon", base_dir + "RA_only_RAA.lon")
