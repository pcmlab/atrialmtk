#!/usr/bin/python3
#
import os
import sys
import argparse
import shutil
from copy import deepcopy
import numpy as np
from sklearn.neighbors import NearestNeighbors
from uac import utils


# Create the parser
my_parser = argparse.ArgumentParser(description="UAC fibre mapping")

# Add the arguments
my_parser.add_argument("target_path", metavar="path", type=str, help="Path for the new target mesh")
my_parser.add_argument("mesh3d_input", metavar="filename", type=str, help="Mesh name 3D")
my_parser.add_argument("mesh2d_input", metavar="filename", type=str, help="Mesh name 2D")
my_parser.add_argument("endo_dir", metavar="path", type=str, help="Path for the endo")
my_parser.add_argument("endo_mesh", metavar="filename", type=str, help="Mesh name endo")
my_parser.add_argument("endo_fibre_file", metavar="filename", type=str, help="Mesh name endo fibre")
my_parser.add_argument("epi_dir", metavar="path", type=str, help="Path for the epi")
my_parser.add_argument("epi_mesh", metavar="filename", type=str, help="Mesh name epi")
my_parser.add_argument("epi_fibre_file", metavar="filename", type=str, help="Mesh name epi fibre")
my_parser.add_argument("bb_dir", metavar="path", type=str, help="Path for BB")
my_parser.add_argument("bb_mesh", metavar="filename", type=str, help="Mesh name BB")
my_parser.add_argument("bb_fibre_file", metavar="filename", type=str, help="Mesh name BB fibre")
my_parser.add_argument("output_file_name", metavar="filename", type=str, help="Fibre output file name")

# Execute parse_args()
args = my_parser.parse_args()

base_dir = args.target_path
mesh3d_input = args.mesh3d_input
mesh2d_input = args.mesh2d_input
endo_dir = args.endo_dir
endo_mesh = args.endo_mesh
endo_fibre_file = args.endo_fibre_file
epi_dir = args.epi_dir
epi_mesh = args.epi_mesh
epi_fibre_file = args.epi_fibre_file
bb_dir = args.bb_dir
bb_mesh = args.bb_mesh
bb_fibre_file = args.bb_fibre_file
output_file_name = args.output_file_name

if not os.path.isdir(base_dir):
    print("The target mesh path specified does not exist")
    sys.exit()


# read carp files, change to 3D versions
print(base_dir)
pts_vol, elems_vol, fiber_vol, data_vol = utils.read_carp(base_dir, mesh3d_input, return_surface=False)
Pts_ABC_1 = utils.read_pts(base_dir + mesh2d_input)

print(endo_dir)
Pts_XYZ_0_Endo = utils.read_pts(endo_dir + endo_mesh)
Elems_XYZ_0_Endo = utils.read_elem(endo_dir + endo_mesh)
print(endo_fibre_file)
Fibres_XYZ_0_Endo = np.loadtxt(endo_dir + endo_fibre_file)


Fibres_XYZ_0_Endo_X = []
Fibres_XYZ_0_Endo_Y = []
Fibres_XYZ_0_Endo_Z = []

for loop in Fibres_XYZ_0_Endo:
    Fibres_XYZ_0_Endo_X.append(loop[0])
    Fibres_XYZ_0_Endo_Y.append(loop[1])
    Fibres_XYZ_0_Endo_Z.append(loop[2])

print(epi_dir)
Pts_XYZ_0_Epi = utils.read_pts(epi_dir + epi_mesh)
Elems_XYZ_0_Epi = utils.read_elem(epi_dir + epi_mesh)
Fibres_XYZ_0_Epi = np.loadtxt(epi_dir + epi_fibre_file)

Fibres_XYZ_0_Epi_X = []
Fibres_XYZ_0_Epi_Y = []
Fibres_XYZ_0_Epi_Z = []

for loop in Fibres_XYZ_0_Epi:
    Fibres_XYZ_0_Epi_X.append(loop[0])
    Fibres_XYZ_0_Epi_Y.append(loop[1])
    Fibres_XYZ_0_Epi_Z.append(loop[2])


M_XYZ_0_2D = utils.mp_calc_ele_4(Pts_ABC_1, elems_vol)

# need to import BB too
Pts_XYZ_0_BB = utils.read_pts(bb_dir + bb_mesh)
print(epi_mesh)
Elems_XYZ_0_BB = utils.read_elem(bb_dir + bb_mesh)
Fibres_XYZ_0_BB = np.loadtxt(bb_dir + bb_fibre_file)


M_XYZ_0_BB = utils.mp_calc_ele_3(Pts_XYZ_0_BB, Elems_XYZ_0_BB)


Fibres_XYZ_0_BB_X = []
Fibres_XYZ_0_BB_Y = []
Fibres_XYZ_0_BB_Z = []

for loop in Fibres_XYZ_0_BB:
    Fibres_XYZ_0_BB_X.append(loop[0])
    Fibres_XYZ_0_BB_Y.append(loop[1])
    Fibres_XYZ_0_BB_Z.append(loop[2])


#  Midpoints in XYZ for both meshes
M_XYZ_0_Endo = utils.mp_calc_ele_3(Pts_XYZ_0_Endo, Elems_XYZ_0_Endo)
M_XYZ_0_Epi = utils.mp_calc_ele_3(Pts_XYZ_0_Epi, Elems_XYZ_0_Epi)

# now for 3D
M_XYZ_0_Vol = utils.mp_calc_ele_4(pts_vol, elems_vol)


# find nearest neighbours and assign fibre depending on endo (<0.5) or epi (>0.5)

# update to KD tree as much faster
neigh = NearestNeighbors(n_neighbors=1)
neigh.fit(M_XYZ_0_BB)
Closest_BB = neigh.kneighbors(M_XYZ_0_Vol, return_distance=False)
neigh.fit(M_XYZ_0_Endo)
Closest_Endo = neigh.kneighbors(M_XYZ_0_Vol, return_distance=False)
neigh.fit(M_XYZ_0_Epi)
Closest_Epi = neigh.kneighbors(M_XYZ_0_Vol, return_distance=False)


Scalar_Label_BB = np.loadtxt(epi_dir + "MappedScalar_BB_LA.dat")

Scalar_Label_Elem = []
Scalar_Label_Elem_Labels = []
for ele_P in Elems_XYZ_0_Epi:
    Value_BB = (Scalar_Label_BB[ele_P.n[0]] + Scalar_Label_BB[ele_P.n[1]] + Scalar_Label_BB[ele_P.n[2]]) / 3
    ToA = 0
    ToB = 12
    if 1 < Value_BB:
        ToA = 10
        ToB = 10
    Scalar_Label_Elem.append(ToA)
    Scalar_Label_Elem_Labels.append(ToB)

print(Scalar_Label_Elem_Labels)
print(len(Scalar_Label_Elem_Labels))

Elems_Labelled = []
for ind in range(len(elems_vol)):
    print(Closest_Epi[ind][0])
    Label_Found = Scalar_Label_Elem_Labels[Closest_Epi[ind][0]]
    ele_P = deepcopy(elems_vol[ind])
    ele_P.region(10 if M_XYZ_0_2D[ind][2] > 0.6 and Label_Found == 10 else 11)
    Elems_Labelled.append(ele_P)


# 5.  calculate fibre in a,b,c on mesh 1
Fibres_ABC_A_1 = []
Fibres_ABC_B_1 = []
Fibres_ABC_C_1 = []

print(len(Closest_Endo))

for ind in range(len(Closest_Endo)):
    Label_Found = Scalar_Label_Elem_Labels[Closest_Epi[ind][0]]
    print(Label_Found)
    if M_XYZ_0_2D[ind][2] < 0.5:
        Index = Closest_Endo[ind][0]
        Fibres_ABC_A_1.append(Fibres_XYZ_0_Endo_X[Index])
        Fibres_ABC_B_1.append(Fibres_XYZ_0_Endo_Y[Index])
        Fibres_ABC_C_1.append(Fibres_XYZ_0_Endo_Z[Index])
    elif M_XYZ_0_2D[ind][2] > 0.6 and Label_Found == 10:
        Index = Closest_BB[ind][0]
        Fibres_ABC_A_1.append(Fibres_XYZ_0_BB_X[Index])
        Fibres_ABC_B_1.append(Fibres_XYZ_0_BB_Y[Index])
        Fibres_ABC_C_1.append(Fibres_XYZ_0_BB_Z[Index])
    else:
        Index = Closest_Epi[ind][0]
        Fibres_ABC_A_1.append(Fibres_XYZ_0_Epi_X[Index])
        Fibres_ABC_B_1.append(Fibres_XYZ_0_Epi_Y[Index])
        Fibres_ABC_C_1.append(Fibres_XYZ_0_Epi_Z[Index])
    if ind % 500 == 0:
        print(ind)

Fibres_ABC_Mesh1 = [Fibres_ABC_A_1, Fibres_ABC_B_1, Fibres_ABC_C_1]
Fibres_ABC_Mesh1 = list(zip(*Fibres_ABC_Mesh1))


# write carp
utils.write_vtk(pts_vol, Elems_Labelled, fiber_vol, data_vol, base_dir + "Fibres_Threshold.vtk")

mname = base_dir + "Fibres_Threshold"
utils.write_carp(pts_vol, Elems_Labelled, Fibres_ABC_Mesh1, None, mname)

# write carp visualisation
mname = base_dir + "Aux_2"
utils.write_carp(M_XYZ_0_Vol, elems_vol, Fibres_ABC_Mesh1, None, mname)
shutil.copyfile(mname + ".pts", base_dir + output_file_name + ".vpts")
shutil.copyfile(mname + ".lon", base_dir + output_file_name + ".vec")
