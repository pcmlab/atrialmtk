#!/usr/bin/python3
#
import os
import sys
import argparse
import shutil
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
my_parser.add_argument("final_dir", metavar="path", type=str, help="Path for RA final")
my_parser.add_argument("output_file_name", metavar="filename", type=str, help="Fibre output file name")

# Execute parse_args()
args = my_parser.parse_args()

base_dir = args.target_path
mesh3d_input = args.mesh3d_input
mesh2d_input = args.mesh2d_input
endo_dir = args.endo_dir
endo_fibre_file = args.endo_fibre_file
epi_dir = args.epi_dir
epi_mesh = args.epi_mesh
epi_fibre_file = args.epi_fibre_file
final_dir = args.final_dir
output_file_name = args.output_file_name

if not os.path.isdir(base_dir):
    print("The target mesh path specified does not exist")
    sys.exit()


# read carp files, change to 3D versions
print(base_dir)
Pts_XYZ_1 = utils.read_pts(base_dir + mesh3d_input)
Elems_XYZ_1 = utils.read_elem(base_dir + mesh3d_input)
Pts_ABC_1 = utils.read_pts(base_dir + mesh2d_input)

print(endo_dir)
print(endo_fibre_file)
Fibres_XYZ_0_Endo = np.loadtxt(final_dir + endo_fibre_file + ".vec")


Fibres_XYZ_0_Endo_X = []
Fibres_XYZ_0_Endo_Y = []
Fibres_XYZ_0_Endo_Z = []

for loop in Fibres_XYZ_0_Endo:
    Fibres_XYZ_0_Endo_X.append(loop[0])
    Fibres_XYZ_0_Endo_Y.append(loop[1])
    Fibres_XYZ_0_Endo_Z.append(loop[2])

print(epi_dir)
Pts_XYZ_0_Epi = utils.read_pts(epi_dir + epi_mesh)
print(epi_mesh)
Elems_XYZ_0_Epi = utils.read_elem(epi_dir + epi_mesh)
Fibres_XYZ_0_Epi = np.loadtxt(final_dir + epi_fibre_file + ".vec")

Fibres_XYZ_0_Epi_X = []
Fibres_XYZ_0_Epi_Y = []
Fibres_XYZ_0_Epi_Z = []

for loop in Fibres_XYZ_0_Epi:
    Fibres_XYZ_0_Epi_X.append(loop[0])
    Fibres_XYZ_0_Epi_Y.append(loop[1])
    Fibres_XYZ_0_Epi_Z.append(loop[2])

M_XYZ_0_Epi = utils.mp_calc_ele_3(Pts_XYZ_0_Epi, Elems_XYZ_0_Epi)
M_XYZ_0_2D = utils.mp_calc_ele_4(Pts_ABC_1, Elems_XYZ_1)

# now for 3D
M_XYZ_0_Vol = utils.mp_calc_ele_4(Pts_XYZ_1, Elems_XYZ_1)


# find nearest neighbours and assign fibre depending on endo (<0.5) or epi (>0.5)

# update to KD tree as much faster

print((len(M_XYZ_0_Epi)))
print((len(Fibres_XYZ_0_Endo_X)))

neigh = NearestNeighbors(n_neighbors=1)
neigh.fit(M_XYZ_0_Epi)
_, Closest_Epi = neigh.kneighbors(M_XYZ_0_Vol, return_distance=True)

# 5.  calculate fibre in a,b,c on mesh 1


# need to add in check for if close enough to Endo

# only add if label
Scalar_Label_PM = np.loadtxt(epi_dir + "MappedScalar_PM.dat")
Scalar_Label_SAN = np.loadtxt(epi_dir + "MappedScalar_SAN.dat")
Scalar_Label_CT = np.loadtxt(epi_dir + "MappedScalar_CT.dat")

_, elems, fiber, data = utils.read_carp(epi_dir, epi_mesh, return_surface=False)

Scalar_Label_Elem = []
Scalar_Label_Elem_Labels = []
for ele_P in elems:
    Value_SAN = (Scalar_Label_SAN[ele_P.n[0]] + Scalar_Label_SAN[ele_P.n[1]] + Scalar_Label_SAN[ele_P.n[2]]) / 3
    Value_CT = (Scalar_Label_CT[ele_P.n[0]] + Scalar_Label_CT[ele_P.n[1]] + Scalar_Label_CT[ele_P.n[2]]) / 3
    Value_PM = (Scalar_Label_PM[ele_P.n[0]] + Scalar_Label_PM[ele_P.n[1]] + Scalar_Label_PM[ele_P.n[2]]) / 3
    ToA = 0
    ToB = 4
    if 1 < Value_SAN:  # relabelled SAN as 5 so different to LA body
        ToA = 5
        ToB = 5
    if 6 < Value_CT:
        ToA = 8
        ToB = 8
    if 7 < Value_PM:
        ToA = 9
        ToB = 9
    Scalar_Label_Elem.append(ToA)
    Scalar_Label_Elem_Labels.append(ToB)


# now add region labels

Pts_Vol, Elems_Vol, fiber, data = utils.read_carp(base_dir, mesh3d_input, return_surface=False)

Scalar_Label_Vol = []
for ind in range(len(Elems_Vol)):
    ele_P = Elems_Vol[ind]
    Index_E = Closest_Epi[ind][0]
    Label_Found = (Scalar_Label_Elem_Labels[Index_E] if M_XYZ_0_2D[ind][2] < 0.75 else 4)
    ele_P.region(Label_Found)
    Scalar_Label_Vol.append(Label_Found)


Fibres_ABC_A_1 = []
Fibres_ABC_B_1 = []
Fibres_ABC_C_1 = []
for ind in range(len(Closest_Epi)):
    Label_Found = Scalar_Label_Vol[ind]
    Index = Closest_Epi[ind][0]
    if M_XYZ_0_2D[ind][2] < 0.75 and Label_Found > 4:
        Fibres_ABC_A_1.append(Fibres_XYZ_0_Endo_X[Index])
        Fibres_ABC_B_1.append(Fibres_XYZ_0_Endo_Y[Index])
        Fibres_ABC_C_1.append(Fibres_XYZ_0_Endo_Z[Index])
    else:
        Fibres_ABC_A_1.append(Fibres_XYZ_0_Epi_X[Index])
        Fibres_ABC_B_1.append(Fibres_XYZ_0_Epi_Y[Index])
        Fibres_ABC_C_1.append(Fibres_XYZ_0_Epi_Z[Index])

Fibres_ABC_Mesh1 = [Fibres_ABC_A_1, Fibres_ABC_B_1, Fibres_ABC_C_1]
Fibres_ABC_Mesh1 = list(zip(*Fibres_ABC_Mesh1))


# write carp
utils.write_vtk(Pts_Vol, Elems_Vol, fiber, data, base_dir + "Fibres_Threshold.vtk")
mname = base_dir + "Fibres_Threshold"
utils.write_carp(Pts_Vol, Elems_Vol, Fibres_ABC_Mesh1, None, mname)

# write carp visualisation
mname = base_dir + "Aux_2"
utils.write_carp(M_XYZ_0_Vol, elems, Fibres_ABC_Mesh1, None, mname)
shutil.copyfile(mname + ".pts", base_dir + output_file_name + ".vpts")
shutil.copyfile(mname + ".lon", base_dir + output_file_name + ".vec")

mname = base_dir + "UAC_Angles_Threshold"
utils.write_carp(M_XYZ_0_Vol, elems, Fibres_ABC_Mesh1, None, mname)
