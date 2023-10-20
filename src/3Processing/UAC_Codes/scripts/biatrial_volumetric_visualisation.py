#!/usr/bin/python3
#
import os
import sys
import argparse
import shutil
from uac import utils


# Create the parser
my_parser = argparse.ArgumentParser(description="UAC fibre mapping")

# Add the arguments
my_parser.add_argument("target_path", metavar="path", type=str, help="Path for the new target mesh")
my_parser.add_argument("mesh_name", metavar="filename", type=str, help="Mesh name target (usually Labelled)")

# Execute parse_args()
args = my_parser.parse_args()

base_dir_bi = args.target_path
mesh_name = args.mesh_name

if not os.path.isdir(base_dir_bi):
    print("The target mesh path specified does not exist")
    sys.exit()


pts, elems, fiber, _ = utils.read_carp(base_dir_bi, mesh_name, return_surface=False)
M_XYZ_1 = utils.mp_calc_ele_4(pts, elems)
mname = base_dir_bi + "Aux_2"
utils.write_carp(M_XYZ_1, elems, fiber, None, mname)
shutil.copyfile(base_dir_bi + "Aux_2.pts", base_dir_bi + "All_Fibres.vpts")
pts, elems, fiber, _ = utils.read_carp(base_dir_bi, "MergeVol_Threshold", return_surface=False)
utils.write_carp(pts, elems, fiber, None, mname)
shutil.copyfile(base_dir_bi + "Aux_2.lon", base_dir_bi + "All_Fibres.vec")

# repeat for each structure

pts, elems, fiber, _ = utils.read_carp(base_dir_bi, "LA_mesh", return_surface=False)
M_XYZ_1 = utils.mp_calc_ele_4(pts, elems)
mname = base_dir_bi + "Aux_2"
utils.write_carp(M_XYZ_1, elems, fiber, None, mname)
shutil.copyfile(base_dir_bi + "Aux_2.pts", base_dir_bi + "LA_mesh.vpts")
utils.write_carp(pts, elems, fiber, None, mname)
shutil.copyfile(base_dir_bi + "Aux_2.lon", base_dir_bi + "LA_mesh.vec")


pts, elems, fiber, _ = utils.read_carp(base_dir_bi, "RA_mesh", return_surface=False)
M_XYZ_1 = utils.mp_calc_ele_4(pts, elems)
mname = base_dir_bi + "Aux_2"
utils.write_carp(M_XYZ_1, elems, fiber, None, mname)
shutil.copyfile(base_dir_bi + "Aux_2.pts", base_dir_bi + "RA_mesh.vpts")
utils.write_carp(pts, elems, fiber, None, mname)
shutil.copyfile(base_dir_bi + "Aux_2.lon", base_dir_bi + "RA_mesh.vec")


pts, elems, fiber, _ = utils.read_carp(base_dir_bi, "RA_structures", return_surface=False)
M_XYZ_1 = utils.mp_calc_ele_4(pts, elems)
mname = base_dir_bi + "Aux_2"
utils.write_carp(M_XYZ_1, elems, fiber, None, mname)
shutil.copyfile(base_dir_bi + "Aux_2.pts", base_dir_bi + "RA_structures.vpts")
utils.write_carp(pts, elems, fiber, None, mname)
shutil.copyfile(base_dir_bi + "Aux_2.lon", base_dir_bi + "RA_structures.vec")
