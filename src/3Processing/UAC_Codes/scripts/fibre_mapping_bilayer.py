#!/usr/bin/python3
#
import os
import sys
import argparse
import shutil
from copy import deepcopy
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from sklearn.neighbors import NearestNeighbors
from uac.element import Element
from uac import utils
import uac.vars as var


# Create the parser
my_parser = argparse.ArgumentParser(description="UAC fibre mapping")

# Add the arguments
my_parser.add_argument("target_path", metavar="path", type=str, help="Path for the new target mesh")
my_parser.add_argument("fibre_path", metavar="path", type=str, help="Path for the atlas mesh")
my_parser.add_argument("fibre_epi_path", metavar="path", type=str, help="Path for the LA epi atlas mesh")
my_parser.add_argument("laplace_path", metavar="path", type=str, help="Path for the laplace par files")
my_parser.add_argument("mesh_name", metavar="filename", type=str, help="Mesh name (usually Labelled)")
my_parser.add_argument("fibre_input_file", metavar="filename", type=str, help="Fibre atlas file name")
my_parser.add_argument("fibre_input_file_epi", metavar="filename", type=str, help="Fibre atlas Epi file name")
my_parser.add_argument("output_file_name", metavar="filename", type=str, help="Fibre output file prefix name")
my_parser.add_argument("bilayer_projection_distance", metavar="value", type=int, help="Bilayer projection distance")
my_parser.add_argument(
    "labels_tag",
    metavar="value",
    type=int,
    help="Tag to indicate if endo with regions (11), epi with regions (12), endo no regions (3), epi no regions (4)",
)

# Execute parse_args()
args = my_parser.parse_args()

base_dir = args.target_path
laplace_dir = args.laplace_path
fibre_dir = args.fibre_path
fibre_epi_dir = args.fibre_epi_path
fibre_file = args.fibre_input_file
fibre_file_epi = args.fibre_input_file_epi
mesh_name = args.mesh_name
bilayer_projection_distance = args.bilayer_projection_distance
labels_tag = args.labels_tag
output_file_name = args.output_file_name

if not os.path.isdir(base_dir):
    print("The target mesh path specified does not exist")
    sys.exit()

if not os.path.isdir(laplace_dir):
    print("The atlas path specified does not exist")
    sys.exit()

if not os.path.isdir(fibre_dir):
    print("The laplace files path specified does not exist")
    sys.exit()


# read carp files
print(base_dir)
Pts_XYZ_1 = utils.read_pts(base_dir + mesh_name)
Elems_XYZ_1 = utils.read_elem(base_dir + mesh_name)
Pts_ABC_1 = utils.read_pts(base_dir + "Labelled_Coords_2D_Rescaling_v3_C")

print(fibre_dir)
Pts_XYZ_0 = utils.read_pts(fibre_dir + "Labelled")
Elems_XYZ_0 = utils.read_elem(fibre_dir + "Labelled")
Pts_ABC_0 = utils.read_pts(fibre_dir + "Labelled_Coords_2D_Rescaling_v3_C")


Fibres_XYZ_0 = np.loadtxt(fibre_dir + fibre_file)
print(Fibres_XYZ_0)


# 1. atlas mesh mid-points from XYZ to ABC
M_ABC_0 = utils.to_ele(Elems_XYZ_0, Pts_ABC_0)

#  Midpoints in XYZ for both meshes
M_XYZ_0 = utils.mp_calc_ele_3(Pts_XYZ_0, Elems_XYZ_0)
M_XYZ_1 = utils.mp_calc_ele_3(Pts_XYZ_1, Elems_XYZ_1)

# 2. atlas mesh fibre end points from XYZ to ABC
# use barycentric coordinates
Fibre_End_Points_XYZ_0 = M_XYZ_0 + var.EPSILON * Fibres_XYZ_0
print(Fibre_End_Points_XYZ_0)


Tri_0_0 = utils.triangle_calc(Pts_XYZ_0, Elems_XYZ_0, 0)
Tri_0_1 = utils.triangle_calc(Pts_XYZ_0, Elems_XYZ_0, 1)
Tri_0_2 = utils.triangle_calc(Pts_XYZ_0, Elems_XYZ_0, 2)

Triabc_0_0 = utils.triangle_calc(Pts_ABC_0, Elems_XYZ_0, 0)
Triabc_0_1 = utils.triangle_calc(Pts_ABC_0, Elems_XYZ_0, 1)
Triabc_0_2 = utils.triangle_calc(Pts_ABC_0, Elems_XYZ_0, 2)

# now loop over points and add them

Fibres_ABC_X = []
Fibres_ABC_Y = []
Fibres_ABC_Z = []

for ind in range(len(Fibre_End_Points_XYZ_0)):
    p = Fibre_End_Points_XYZ_0[ind]
    a = Tri_0_0[ind]
    b = Tri_0_1[ind]
    c = Tri_0_2[ind]
    u, v, _ = utils.barycentric_coord_v1(p, a, b, c)
    Location_XYZ = c * u + b * v + a * (1 - u - v)
    a_abc = Triabc_0_0[ind]
    b_abc = Triabc_0_1[ind]
    c_abc = Triabc_0_2[ind]
    Location_ABC = c_abc * u + b_abc * v + a_abc * (1 - u - v)
    Vector_UAC = Location_ABC - M_ABC_0[ind]
    Vector_UAC = Vector_UAC / np.linalg.norm(Vector_UAC)
    Fibres_ABC_X.append(Vector_UAC[0])
    Fibres_ABC_Y.append(Vector_UAC[1])
    Fibres_ABC_Z.append(Vector_UAC[2])


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


# 6.  map fibre back to x, y, z on mesh 1


Tri_1_0 = utils.triangle_calc(Pts_XYZ_1, Elems_XYZ_1, 0)
Tri_1_1 = utils.triangle_calc(Pts_XYZ_1, Elems_XYZ_1, 1)
Tri_1_2 = utils.triangle_calc(Pts_XYZ_1, Elems_XYZ_1, 2)

Triabc_1_0 = utils.triangle_calc(Pts_ABC_1, Elems_XYZ_1, 0)
Triabc_1_1 = utils.triangle_calc(Pts_ABC_1, Elems_XYZ_1, 1)
Triabc_1_2 = utils.triangle_calc(Pts_ABC_1, Elems_XYZ_1, 2)


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
    u, v, _ = utils.barycentric_coord_v1(p, a, b, c)
    Location_XYZ = c * u + b * v + a * (1 - u - v)
    a_xyz = Tri_1_0[ind]
    b_xyz = Tri_1_1[ind]
    c_xyz = Tri_1_2[ind]
    Location_XYZ = c_xyz * u + b_xyz * v + a_xyz * (1 - u - v)
    Vector_XYZ = Location_XYZ - M_XYZ_1[ind]
    Vector_XYZ = Vector_XYZ / np.linalg.norm(Vector_XYZ)
    A = Vector_XYZ[0]
    B = Vector_XYZ[1]
    C = Vector_XYZ[2]
    A = np.array(A)
    B = np.array(B)
    C = np.array(C)
    where_are_NaNs = np.isnan(A)
    A[where_are_NaNs] = 1
    B[where_are_NaNs] = 0
    C[where_are_NaNs] = 0
    Fibres_XYZ_X1.append(A)
    Fibres_XYZ_Y1.append(B)
    Fibres_XYZ_Z1.append(C)

Fibres_XYZ_1 = [Fibres_XYZ_X1, Fibres_XYZ_Y1, Fibres_XYZ_Z1]
Fibres_XYZ_1 = list(zip(*Fibres_XYZ_1))

Fibres_XYZ_1_Endo = [Fibres_XYZ_X1, Fibres_XYZ_Y1, Fibres_XYZ_Z1]
Fibres_XYZ_1_Endo = list(zip(*Fibres_XYZ_1_Endo))


# repeat for Epi


# read carp files
print(fibre_epi_dir)
Pts_XYZ_0 = utils.read_pts(fibre_epi_dir + "Labelled")
Elems_XYZ_0 = utils.read_elem(fibre_epi_dir + "Labelled")
Pts_ABC_0 = utils.read_pts(fibre_epi_dir + "Labelled_Coords_2D_Rescaling_v3_C")


Fibres_XYZ_0 = np.loadtxt(fibre_epi_dir + fibre_file_epi)
print(Fibres_XYZ_0)


# 1. atlas mesh mid-points from XYZ to ABC
M_ABC_0 = utils.to_ele(Elems_XYZ_0, Pts_ABC_0)

#  Midpoints in XYZ for both meshes
M_XYZ_0 = utils.mp_calc_ele_3(Pts_XYZ_0, Elems_XYZ_0)

xlist_1 = []
ylist_1 = []
zlist_1 = []

for loop in Pts_XYZ_1:
    xlist_1.append(loop[0])
    ylist_1.append(loop[1])
    zlist_1.append(loop[2])

# 2. atlas mesh fibre end points from XYZ to ABC
# use barycentric coordinates
Fibre_End_Points_XYZ_0 = M_XYZ_0 + var.EPSILON * Fibres_XYZ_0
print(Fibre_End_Points_XYZ_0)

Tri_0_0 = utils.triangle_calc(Pts_XYZ_0, Elems_XYZ_0, 0)
Tri_0_1 = utils.triangle_calc(Pts_XYZ_0, Elems_XYZ_0, 1)
Tri_0_2 = utils.triangle_calc(Pts_XYZ_0, Elems_XYZ_0, 2)

Triabc_0_0 = utils.triangle_calc(Pts_ABC_0, Elems_XYZ_0, 0)
Triabc_0_1 = utils.triangle_calc(Pts_ABC_0, Elems_XYZ_0, 1)
Triabc_0_2 = utils.triangle_calc(Pts_ABC_0, Elems_XYZ_0, 2)

# now loop over points and add them

Fibres_ABC_X = []
Fibres_ABC_Y = []
Fibres_ABC_Z = []

for ind in range(len(Fibre_End_Points_XYZ_0)):
    p = Fibre_End_Points_XYZ_0[ind]
    a = Tri_0_0[ind]
    b = Tri_0_1[ind]
    c = Tri_0_2[ind]
    u, v, _ = utils.barycentric_coord_v1(p, a, b, c)
    Location_XYZ = c * u + b * v + a * (1 - u - v)
    a_abc = Triabc_0_0[ind]
    b_abc = Triabc_0_1[ind]
    c_abc = Triabc_0_2[ind]
    Location_ABC = c_abc * u + b_abc * v + a_abc * (1 - u - v)
    Vector_UAC = Location_ABC - M_ABC_0[ind]
    Vector_UAC = Vector_UAC / np.linalg.norm(Vector_UAC)
    Fibres_ABC_X.append(Vector_UAC[0])
    Fibres_ABC_Y.append(Vector_UAC[1])
    Fibres_ABC_Z.append(Vector_UAC[2])


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


# 6.  map fibre back to x, y, z on mesh 1

Fibres_XYZ_X1_Epi = []
Fibres_XYZ_Y1_Epi = []
Fibres_XYZ_Z1_Epi = []

for ind in range(len(M_ABC_1)):
    Fibres_ABC_Mesh1_0 = Fibres_ABC_Mesh1[ind]
    Fibres_ABC_Mesh1_0 = np.array(Fibres_ABC_Mesh1_0)
    Fibre_End_Points_ABC_1 = M_ABC_1[ind] + var.EPSILON_ABC * Fibres_ABC_Mesh1_0

    p = Fibre_End_Points_ABC_1
    a = Triabc_1_0[ind]
    b = Triabc_1_1[ind]
    c = Triabc_1_2[ind]
    u, v, _ = utils.barycentric_coord_v1(p, a, b, c)
    Location_XYZ = c * u + b * v + a * (1 - u - v)
    a_xyz = Tri_1_0[ind]
    b_xyz = Tri_1_1[ind]
    c_xyz = Tri_1_2[ind]
    Location_XYZ = c_xyz * u + b_xyz * v + a_xyz * (1 - u - v)
    Vector_XYZ = Location_XYZ - M_XYZ_1[ind]
    Vector_XYZ = Vector_XYZ / np.linalg.norm(Vector_XYZ)
    A = Vector_XYZ[0]
    B = Vector_XYZ[1]
    C = Vector_XYZ[2]
    A = np.array(A)
    B = np.array(B)
    C = np.array(C)
    where_are_NaNs = np.isnan(A)
    A[where_are_NaNs] = 1
    B[where_are_NaNs] = 0
    C[where_are_NaNs] = 0
    Fibres_XYZ_X1.append(A)
    Fibres_XYZ_Y1.append(B)
    Fibres_XYZ_Z1.append(C)
    Fibres_XYZ_X1_Epi.append(A)
    Fibres_XYZ_Y1_Epi.append(B)
    Fibres_XYZ_Z1_Epi.append(C)


Fibres_XYZ_1 = [Fibres_XYZ_X1, Fibres_XYZ_Y1, Fibres_XYZ_Z1]
Fibres_XYZ_1 = list(zip(*Fibres_XYZ_1))

Fibres_XYZ_1_Epi = [Fibres_XYZ_X1_Epi, Fibres_XYZ_Y1_Epi, Fibres_XYZ_Z1_Epi]
Fibres_XYZ_1_Epi = list(zip(*Fibres_XYZ_1_Epi))


# now add line fibers

for loop in Pts_XYZ_1:
    Fibres_XYZ_X1.append(1)
    Fibres_XYZ_Y1.append(0)
    Fibres_XYZ_Z1.append(0)


Fibres_XYZ_1 = [Fibres_XYZ_X1, Fibres_XYZ_Y1, Fibres_XYZ_Z1]
Fibres_XYZ_1 = list(zip(*Fibres_XYZ_1))


# duplicate mesh to make a bilayer
pts, elems, fiber, data = utils.read_carp(base_dir, mesh_name, return_surface=False)
utils.write_vtk(pts, elems, fiber, data, base_dir + "Fibres.vtk")
surface = utils.get_vtk_from_file(base_dir + "Fibres.vtk")
surface = utils.convert_unstructureddata_to_polydata(surface)
_, selems = utils.poly2carp_with_labels(surface)

normalFilt = vtk.vtkPolyDataNormals()
if vtk.VTK_MAJOR_VERSION < 6:
    normalFilt.SetInput(surface)
else:
    normalFilt.SetInputData(surface)

normalFilt.GetOutput().GetPointData().RemoveArray("Normals")
normalFilt.GetOutput().GetCellData().RemoveArray("Normals")
normalFilt.ComputeCellNormalsOff()
normalFilt.ComputePointNormalsOn()
normalFilt.SplittingOff()
# normalFilt.ConsistencyOn()
# normalFilt.AutoOrientNormalsOn()
normalFilt.Update()

pointNormals = vtk_to_numpy(normalFilt.GetOutput().GetPointData().GetNormals())
print((len(pointNormals)))


# now project outwards
x_norm = []
y_norm = []
z_norm = []

for loop in pointNormals:
    x_norm.append(np.multiply(loop[0] / np.sqrt(loop[0] ** 2 + loop[1] ** 2 + loop[2] ** 2), bilayer_projection_distance))
    y_norm.append(np.multiply(loop[1] / np.sqrt(loop[0] ** 2 + loop[1] ** 2 + loop[2] ** 2), bilayer_projection_distance))
    z_norm.append(np.multiply(loop[2] / np.sqrt(loop[0] ** 2 + loop[1] ** 2 + loop[2] ** 2), bilayer_projection_distance))

xlist_Epi = []
ylist_Epi = []
zlist_Epi = []

for ind in range(len(xlist_1)):
    xlist_Epi.append(x_norm[ind] + xlist_1[ind])
    ylist_Epi.append(y_norm[ind] + ylist_1[ind])
    zlist_Epi.append(z_norm[ind] + zlist_1[ind])


print((len(xlist_1)))
print((len(xlist_Epi)))

x_list_All = np.concatenate([xlist_1, xlist_Epi])
y_list_All = np.concatenate([ylist_1, ylist_Epi])
z_list_All = np.concatenate([zlist_1, zlist_Epi])

Pts_All = [x_list_All, y_list_All, z_list_All]
Pts_All = list(zip(*Pts_All))


Elems_Endo = deepcopy(np.copy(elems))

# need to be for new points...

Elems_Epi = []
for ind in range(len(elems)):
    ele_P = deepcopy(elems[ind])
    if labels_tag == 11:
        ele_P.reg += 1
    elif labels_tag == 12:
        ele_P.reg -= 1
    else:
        ele_P.region(11)
    for i in range(len(ele_P.n)):
        ele_P.n[i] += len(pts)
    Elems_Epi.append(ele_P)


# now add line elements


Elems_Ln = []
for ind in range(len(pts)):
    Elems_Ln.append(Element("Ln", [ind, ind + len(pts)], 32))

Elems_2 = np.concatenate([Elems_Endo, Elems_Epi, Elems_Ln])

if labels_tag == 3:
    for i in range(len(Elems_Endo)):
        Elems_2[i].region(12)


# write mesh
mname = base_dir + output_file_name + "Bilayer"
utils.write_carp(Pts_All, Elems_2, Fibres_XYZ_1, None, mname)

mname = base_dir + output_file_name + "Labelled_Endo"
utils.write_carp(pts, elems, Fibres_XYZ_1_Endo, None, mname)

mname = base_dir + output_file_name + "Labelled_Epi"
utils.write_carp(pts, elems, Fibres_XYZ_1_Epi, None, mname)

# write fibres for endo visualisation need to be at midpoints
mname = base_dir + "Aux_2"
utils.write_carp(M_XYZ_1, selems, Fibres_XYZ_1, None, mname)
shutil.copyfile(base_dir + "Aux_2.pts", base_dir + output_file_name + "Labelled_Endo.vpts")
shutil.copyfile(base_dir + mesh_name + ".elem", base_dir + output_file_name + "Labelled_Endo.elem")
shutil.copyfile(base_dir + output_file_name + "Labelled_Endo.lon", base_dir + output_file_name + "Labelled_Endo.vec")

# write fibres for epi visualisation
M_XYZ_1 = utils.mp_calc_ele_3(list(zip(*[xlist_Epi, ylist_Epi, zlist_Epi])), Elems_XYZ_1)
mname = base_dir + "Aux_2"
utils.write_carp(M_XYZ_1, selems, Fibres_XYZ_1, None, mname)
shutil.copyfile(base_dir + "Aux_2.pts", base_dir + output_file_name + "Labelled_Epi.vpts")
shutil.copyfile(base_dir + mesh_name + ".elem", base_dir + output_file_name + "Labelled_Epi.elem")
shutil.copyfile(base_dir + output_file_name + "Labelled_Epi.lon", base_dir + output_file_name + "Labelled_Epi.vec")
