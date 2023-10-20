#!/usr/bin/python3
#
import os
import sys
import argparse
import numpy as np
from sklearn.neighbors import NearestNeighbors
from uac import utils


# Create the parser
my_parser = argparse.ArgumentParser(description="UAC scalar mapping")

# Add the arguments
my_parser.add_argument("target_path", metavar="path", type=str, help="Path for the new target mesh")
my_parser.add_argument("fibre_path", metavar="path", type=str, help="Path for the atlas mesh")
my_parser.add_argument("mesh_name", metavar="filename", type=str, help="Mesh name")
my_parser.add_argument("mesh_name_atlas", metavar="filename", type=str, help="Mesh name atlas")
my_parser.add_argument("mesh_name_atlas_2d", metavar="filename", type=str, help="Mesh name atlas 2D")
my_parser.add_argument("mesh_name_2d", metavar="filename", type=str, help="Mesh name 2D")
my_parser.add_argument("scalar_name", metavar="filename", type=str, help="Scalar file name")
my_parser.add_argument("scalar_output_name", metavar="filename", type=str, help="Scalar output file name")

# Execute parse_args()
args = my_parser.parse_args()

base_dir = args.target_path
fibre_dir = args.fibre_path
scalar_file_name = args.scalar_name
scalar_output_file_name = args.scalar_output_name
mesh_name_atlas = args.mesh_name_atlas
mesh_name_atlas_2d = args.mesh_name_atlas_2d
mesh_name_2d = args.mesh_name_2d

if not os.path.isdir(base_dir):
    print("The target mesh path specified does not exist")
    sys.exit()

if not os.path.isdir(fibre_dir):
    print("The fibre files path specified does not exist")
    sys.exit()


def Scalar_Triangles(scalar, ele, triNo):
    MPX = []
    for e in ele:
        Face = e.n
        FX = scalar[Face[triNo]]
        MPX.append(FX)
    return MPX


# read carp files
print(base_dir)
Pts_ABC_1 = utils.read_pts(base_dir + mesh_name_2d)

print(fibre_dir)
Elems_XYZ_0 = utils.read_elem(fibre_dir + mesh_name_atlas)
Pts_ABC_0 = utils.read_pts(fibre_dir + mesh_name_atlas_2d)

print((fibre_dir + scalar_file_name))
Scalar_XYZ_0 = np.loadtxt(fibre_dir + scalar_file_name)
print(Scalar_XYZ_0)


# 1. atlas mesh mid-points from XYZ to ABC
M_ABC_0 = utils.to_ele(Elems_XYZ_0, Pts_ABC_0)

# 2.

Triabc_0_0 = utils.triangle_calc(Pts_ABC_0, Elems_XYZ_0, 0)
Triabc_0_1 = utils.triangle_calc(Pts_ABC_0, Elems_XYZ_0, 1)
Triabc_0_2 = utils.triangle_calc(Pts_ABC_0, Elems_XYZ_0, 2)

# 4. closest points for ABC1 in mid point ABC0

neigh = NearestNeighbors(n_neighbors=1)
neigh.fit(M_ABC_0)
Closest_MP_1_0 = neigh.kneighbors(Pts_ABC_1, return_distance=False)


# 5. barycentric for each

Scalar_0_0 = Scalar_Triangles(Scalar_XYZ_0, Elems_XYZ_0, 0)
Scalar_0_1 = Scalar_Triangles(Scalar_XYZ_0, Elems_XYZ_0, 1)
Scalar_0_2 = Scalar_Triangles(Scalar_XYZ_0, Elems_XYZ_0, 2)

WS = []

for ind in range(len(Pts_ABC_1)):
    p = Pts_ABC_1[ind]
    Index = Closest_MP_1_0[ind]
    a = Triabc_0_0[Index[0]]
    b = Triabc_0_1[Index[0]]
    c = Triabc_0_2[Index[0]]
    u, v, _ = utils.barycentric_coord_v1(p, a, b, c)
    Weighted_Scaler = Scalar_0_2[Index[0]] * u + Scalar_0_1[Index[0]] * v + Scalar_0_0[Index[0]] * (1 - u - v)
    WS.append(Weighted_Scaler)


utils.write_dat_simple(WS, base_dir + scalar_output_file_name)
