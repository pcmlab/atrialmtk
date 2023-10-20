#!/usr/bin/python3
#
import os
import sys
import argparse
import shutil
import numpy as np
from sklearn.neighbors import NearestNeighbors
from uac import utils
import uac.vars as var


# Create the parser
my_parser = argparse.ArgumentParser(description="UAC fibre mapping")

# Add the arguments
my_parser.add_argument("target_path", metavar="path", type=str, help="Path for the new target mesh")
my_parser.add_argument("fibre_path", metavar="path", type=str, help="Path for the atlas mesh")
my_parser.add_argument("laplace_path", metavar="path", type=str, help="Path for the laplace par files")
my_parser.add_argument("mesh3d_input", metavar="filename", type=str, help="Mesh name 3D")
my_parser.add_argument("mesh2d_input", metavar="filename", type=str, help="Mesh name 2D")
my_parser.add_argument("mesh3d_fibre", metavar="filename", type=str, help="Fibre Mesh name 3D")
my_parser.add_argument("mesh2d_fibre", metavar="filename", type=str, help="Fibre Mesh name 2D")
my_parser.add_argument("fibre_input_file", metavar="filename", type=str, help="Fibre atlas file name")
my_parser.add_argument("output_file_name", metavar="filename", type=str, help="Fibre output file name")

# Execute parse_args()
args = my_parser.parse_args()

base_dir = args.target_path
laplace_dir = args.laplace_path
fibre_dir = args.fibre_path
fibre_file = args.fibre_input_file
output_file_name = args.output_file_name
mesh3d_input = args.mesh3d_input
mesh2d_input = args.mesh2d_input
mesh3d_fibre = args.mesh3d_fibre
mesh2d_fibre = args.mesh2d_fibre

if not os.path.isdir(base_dir):
    print("The target mesh path specified does not exist")
    sys.exit()

if not os.path.isdir(laplace_dir):
    print("The atlas path specified does not exist")
    sys.exit()

if not os.path.isdir(fibre_dir):
    print("The laplace files path specified does not exist")
    sys.exit()


def Fibre_Map(Pts_XYZ_0, Pts_ABC_0, Elems_XYZ_0, Fibre_End_Points_XYZ_0, M_ABC_0):
    Tri_0_0 = utils.triangle_calc(Pts_XYZ_0, Elems_XYZ_0, 0)
    Tri_0_1 = utils.triangle_calc(Pts_XYZ_0, Elems_XYZ_0, 1)
    Tri_0_2 = utils.triangle_calc(Pts_XYZ_0, Elems_XYZ_0, 2)
    Tri_0_3 = utils.triangle_calc(Pts_XYZ_0, Elems_XYZ_0, 3)
    Triabc_0_0 = utils.triangle_calc(Pts_ABC_0, Elems_XYZ_0, 0)
    Triabc_0_1 = utils.triangle_calc(Pts_ABC_0, Elems_XYZ_0, 1)
    Triabc_0_2 = utils.triangle_calc(Pts_ABC_0, Elems_XYZ_0, 2)
    Triabc_0_3 = utils.triangle_calc(Pts_ABC_0, Elems_XYZ_0, 3)

    # now loop over points and add them

    Fibres_ABC_X = []
    Fibres_ABC_Y = []
    Fibres_ABC_Z = []

    for ind in range(len(Fibre_End_Points_XYZ_0)):
        p = Fibre_End_Points_XYZ_0[ind]
        a = Tri_0_0[ind]
        b = Tri_0_1[ind]
        c = Tri_0_2[ind]
        d = (Tri_0_3[ind])

        u, v, w, _ = utils.barycentric_coord_v2(p, a, b, c, d)
        a_abc = Triabc_0_0[ind]
        b_abc = Triabc_0_1[ind]
        c_abc = Triabc_0_2[ind]
        d_abc = (Triabc_0_3[ind])
        Location_ABC = a_abc * u + b_abc * v + c_abc * w + d_abc * (1 - u - v - w)

        Vector_UAC = Location_ABC - M_ABC_0[ind]
        Vector_UAC = Vector_UAC / np.linalg.norm(Vector_UAC)
        Fibres_ABC_X.append(Vector_UAC[0])
        Fibres_ABC_Y.append(Vector_UAC[1])
        Fibres_ABC_Z.append(Vector_UAC[2])

    Fibres_ABC = [Fibres_ABC_X, Fibres_ABC_Y, Fibres_ABC_Z]
    Fibres_ABC = list(zip(*Fibres_ABC))
    return Fibres_ABC, Fibres_ABC_X, Fibres_ABC_Y, Fibres_ABC_Z


def Fibre_Map_2(Pts_XYZ_1, Pts_ABC_1, Elems_XYZ_1, Fibres_ABC_Mesh1, M_ABC_1, M_XYZ_1):
    Tri_1_0 = utils.triangle_calc(Pts_XYZ_1, Elems_XYZ_1, 0)
    Tri_1_1 = utils.triangle_calc(Pts_XYZ_1, Elems_XYZ_1, 1)
    Tri_1_2 = utils.triangle_calc(Pts_XYZ_1, Elems_XYZ_1, 2)
    Tri_1_3 = utils.triangle_calc(Pts_XYZ_1, Elems_XYZ_1, 3)
    Triabc_1_0 = utils.triangle_calc(Pts_ABC_1, Elems_XYZ_1, 0)
    Triabc_1_1 = utils.triangle_calc(Pts_ABC_1, Elems_XYZ_1, 1)
    Triabc_1_2 = utils.triangle_calc(Pts_ABC_1, Elems_XYZ_1, 2)
    Triabc_1_3 = utils.triangle_calc(Pts_ABC_1, Elems_XYZ_1, 3)

    Fibres_XYZ_X1 = []
    Fibres_XYZ_Y1 = []
    Fibres_XYZ_Z1 = []

    for ind in range(len(M_ABC_1)):
        Fibres_ABC_Mesh1_0 = Fibres_ABC_Mesh1[ind]
        Fibres_ABC_Mesh1_0 = np.array(Fibres_ABC_Mesh1_0)
        Fibre_End_Points_ABC_1 = M_ABC_1[ind] + var.EPSILON_ABC * Fibres_ABC_Mesh1_0
        p = Fibre_End_Points_ABC_1
        a = Triabc_1_0[ind]
        b = Triabc_1_1[ind]
        c = Triabc_1_2[ind]
        d = (Triabc_1_3[ind])

        u, v, w, _ = utils.barycentric_coord_v2(p, a, b, c, d)
        a_xyz = Tri_1_0[ind]
        b_xyz = Tri_1_1[ind]
        c_xyz = Tri_1_2[ind]
        d_xyz = (Tri_1_3[ind])

        Location_XYZ = a_xyz * u + b_xyz * v + c_xyz * w + d_xyz * (1 - u - v - w)

        Vector_XYZ = Location_XYZ - M_XYZ_1[ind]
        Vector_XYZ = Vector_XYZ / np.linalg.norm(Vector_XYZ)
        Fibres_XYZ_X1.append(Vector_XYZ[0])
        Fibres_XYZ_Y1.append(Vector_XYZ[1])

        Fibres_XYZ_Z1.append(Vector_XYZ[2])

    Fibres_XYZ_1 = [Fibres_XYZ_X1, Fibres_XYZ_Y1, Fibres_XYZ_Z1]
    Fibres_XYZ_1 = list(zip(*Fibres_XYZ_1))
    return Fibres_XYZ_1, Fibres_XYZ_X1, Fibres_XYZ_Y1, Fibres_XYZ_Z1


# read carp files, change to 3D versions
print(base_dir)
Pts_XYZ_1 = utils.read_pts(base_dir + mesh3d_input)
Elems_XYZ_1 = utils.read_elem(base_dir + mesh3d_input)
Pts_ABC_1 = utils.read_pts(base_dir + mesh2d_input)

print(fibre_dir)
Pts_XYZ_0 = utils.read_pts(fibre_dir + mesh3d_fibre)
Elems_XYZ_0 = utils.read_elem(fibre_dir + mesh3d_fibre)
Pts_ABC_0 = utils.read_pts(fibre_dir + mesh2d_fibre)


Fibres_XYZ_0 = np.loadtxt(fibre_dir + fibre_file)
print(Fibres_XYZ_0)


# 1. atlas mesh mid-points from XYZ to ABC
M_ABC_0 = utils.to_ele(Elems_XYZ_0, Pts_ABC_0)


#  Midpoints in XYZ for both meshes
M_XYZ_0 = utils.mp_calc_ele_4(Pts_XYZ_0, Elems_XYZ_0)

# 2. atlas mesh fibre end points from XYZ to ABC
# use barycentric coordinates
Fibre_End_Points_XYZ_0 = M_XYZ_0 + var.EPSILON * Fibres_XYZ_0

del M_XYZ_0

print((1))

_, Fibres_ABC_X, Fibres_ABC_Y, Fibres_ABC_Z = Fibre_Map(
    Pts_XYZ_0, Pts_ABC_0, Elems_XYZ_0, Fibre_End_Points_XYZ_0, M_ABC_0
)

print((123))

# 3. new mesh mid-points from XYZ to ABC
M_ABC_1 = utils.to_ele(Elems_XYZ_1, Pts_ABC_1)


# 4. closest midpoints for ABC1 in ABC0
neigh = NearestNeighbors(n_neighbors=1)
neigh.fit(M_ABC_0)
Closest_MP_1_0 = neigh.kneighbors(M_ABC_1, return_distance=False)
print((Closest_MP_1_0[0]))

# 5.  calculate fibre in a,b,c on mesh 1
Fibres_ABC_A_1 = []
Fibres_ABC_B_1 = []
Fibres_ABC_C_1 = []

for ind in range(len(Closest_MP_1_0)):
    Index = Closest_MP_1_0[ind]
    Index = Index[0]
    Fibres_ABC_A_1.append(Fibres_ABC_X[Index])
    Fibres_ABC_B_1.append(Fibres_ABC_Y[Index])
    Fibres_ABC_C_1.append(Fibres_ABC_Z[Index])


Fibres_ABC_Mesh1 = [Fibres_ABC_A_1, Fibres_ABC_B_1, Fibres_ABC_C_1]
Fibres_ABC_Mesh1 = list(zip(*Fibres_ABC_Mesh1))

del Closest_MP_1_0, Fibres_ABC_A_1, Fibres_ABC_B_1, Fibres_ABC_C_1

print((2))

# 6.  map fibre back to x, y, z on mesh 1
M_XYZ_1 = utils.mp_calc_ele_4(Pts_XYZ_1, Elems_XYZ_1)

Fibres_XYZ_1, *_ = Fibre_Map_2(Pts_XYZ_1, Pts_ABC_1, Elems_XYZ_1, Fibres_ABC_Mesh1, M_ABC_1, M_XYZ_1)

print((3))


print(base_dir)

# write carp
pts, elems, fiber, data = utils.read_carp(base_dir, mesh3d_input, return_surface=False)
utils.write_vtk(pts, elems, fiber, data, base_dir + "Fibres.vtk")
mname = base_dir + "Fibres"
utils.write_carp(pts, elems, Fibres_XYZ_1, None, mname)

# write carp visualisation
mname = base_dir + "Aux_2"
utils.write_carp(M_XYZ_1, elems, Fibres_XYZ_1, None, mname)
shutil.copyfile(mname + ".pts", base_dir + output_file_name + ".vpts")
shutil.copyfile(mname + ".lon", base_dir + output_file_name + ".vec")

mname = base_dir + "UAC_Angles"
utils.write_carp(M_XYZ_1, elems, Fibres_ABC_Mesh1, None, mname)
