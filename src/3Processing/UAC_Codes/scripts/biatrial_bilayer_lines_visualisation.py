#!/usr/bin/python3
#
import os
import sys
import argparse
import shutil
import numpy as np
from sklearn.neighbors import NearestNeighbors
from uac.element import Element
from uac import utils


# Create the parser
my_parser = argparse.ArgumentParser(description="UAC fibre mapping")

# Add the arguments
my_parser.add_argument("target_path", metavar="path", type=str, help="Path for the new target mesh")
my_parser.add_argument("mesh_name", metavar="filename", type=str, help="Mesh name target (usually Labelled)")
my_parser.add_argument("output_file_name", metavar="filename", type=str, help="Fibre output file name")

# Execute parse_args()
args = my_parser.parse_args()

base_dir_bi = args.target_path
mesh_name = args.mesh_name

if not os.path.isdir(base_dir_bi):
    print("The target mesh path specified does not exist")
    sys.exit()


# Import LA.

pts, elems, fiber, _ = utils.read_carp(base_dir_bi, mesh_name, return_surface=False)

# add lines for two LA surfaces, and RA structures...

# first define scalar label for which surface

F1 = []
F2 = []
F3 = []
Color = []

for ele_P in elems:
    F1.append(ele_P.n[0])
    F2.append(ele_P.n[1])
    F3.append(ele_P.n[2])
    Color.append(ele_P.reg)

F1 = np.array(F1)
F2 = np.array(F2)
F3 = np.array(F3)
Color = np.array(Color)

Scalar_Color = []
for ind in range(len(pts)):
    Found1 = np.argwhere(F1 == ind)
    Found2 = np.argwhere(F2 == ind)
    Found3 = np.argwhere(F3 == ind)
    Fall = np.append(Found1, Found2)
    Fall = np.append(Fall, Found3)
    # print(max(Fall))
    Scalar_Color.append(Color[max(Fall)])
    if ind % 100 == 0:
        print(ind)

Scalar_Color = np.array(Scalar_Color)
print(Scalar_Color)

utils.write_dat_simple(Scalar_Color, base_dir_bi + "Scalar_Labels.dat")

# loop over if node label =12, find closest node label =11.
LA_epi = np.argwhere(Scalar_Color == 12)
LA_endo = np.argwhere(Scalar_Color == 11)

print(LA_epi)
print(LA_endo)

xlist = []
ylist = []
zlist = []

for loop in pts:
    xlist.append(loop[0])
    ylist.append(loop[1])
    zlist.append(loop[2])

P_x_epi = []
P_y_epi = []
P_z_epi = []

for ind1 in range(len(LA_epi)):
    P1 = LA_epi[ind1]
    P_x_epi.append(xlist[P1[0]])
    P_y_epi.append(ylist[P1[0]])
    P_z_epi.append(zlist[P1[0]])

PTS_LA_epi = [P_x_epi, P_y_epi, P_z_epi]
PTS_LA_epi = list(zip(*PTS_LA_epi))


P_x_endo = []
P_y_endo = []
P_z_endo = []

for ind1 in range(len(LA_endo)):
    P1 = LA_endo[ind1]
    P_x_endo.append(xlist[P1[0]])
    P_y_endo.append(ylist[P1[0]])
    P_z_endo.append(zlist[P1[0]])

PTS_LA_endo = [P_x_endo, P_y_endo, P_z_endo]
print(PTS_LA_endo)
PTS_LA_endo = list(zip(*PTS_LA_endo))
print(PTS_LA_endo)


neigh = NearestNeighbors(n_neighbors=1)
neigh.fit(PTS_LA_endo)
_, Closest_LA = neigh.kneighbors(PTS_LA_epi, return_distance=True)


Elems_Ln = []
for ind in range(len(LA_epi)):
    Elems_Ln.append(Element("Ln", [LA_epi[ind][0], LA_endo[Closest_LA[ind]][0][0]], 32))

Extra = len(Elems_Ln)

Elems_2 = np.concatenate([elems, Elems_Ln])


# now add SAN, CT, PM


RA_epi = np.argwhere(Scalar_Color < 3)
SAN = np.argwhere(Scalar_Color == 3)


P_x_epi = []
P_y_epi = []
P_z_epi = []

for ind1 in range(len(SAN)):
    P1 = SAN[ind1]
    P_x_epi.append(xlist[P1[0]])
    P_y_epi.append(ylist[P1[0]])
    P_z_epi.append(zlist[P1[0]])

PTS_LA_epi = [P_x_epi, P_y_epi, P_z_epi]
PTS_LA_epi = list(zip(*PTS_LA_epi))


P_x_endo = []
P_y_endo = []
P_z_endo = []

for ind1 in range(len(RA_epi)):
    P1 = RA_epi[ind1]
    P_x_endo.append(xlist[P1[0]])
    P_y_endo.append(ylist[P1[0]])
    P_z_endo.append(zlist[P1[0]])

PTS_LA_endo = [P_x_endo, P_y_endo, P_z_endo]
print(PTS_LA_endo)
PTS_LA_endo = list(zip(*PTS_LA_endo))
print(PTS_LA_endo)


neigh = NearestNeighbors(n_neighbors=1)
neigh.fit(PTS_LA_endo)
_, Closest_LA = neigh.kneighbors(PTS_LA_epi, return_distance=True)


Elems_Ln = []
for ind in range(len(SAN)):
    Elems_Ln.append(Element("Ln", [SAN[ind][0], RA_epi[Closest_LA[ind]][0][0]], 33))

Extra += len(Elems_Ln)

Elems_2 = np.concatenate([Elems_2, Elems_Ln])

CT = np.argwhere(Scalar_Color == 8)


P_x_epi = []
P_y_epi = []
P_z_epi = []

for ind1 in range(len(CT)):
    P1 = CT[ind1]
    P_x_epi.append(xlist[P1[0]])
    P_y_epi.append(ylist[P1[0]])
    P_z_epi.append(zlist[P1[0]])

PTS_LA_epi = [P_x_epi, P_y_epi, P_z_epi]
PTS_LA_epi = list(zip(*PTS_LA_epi))


neigh = NearestNeighbors(n_neighbors=1)
neigh.fit(PTS_LA_endo)
_, Closest_LA = neigh.kneighbors(PTS_LA_epi, return_distance=True)


Elems_Ln = []
for ind in range(len(CT)):
    Elems_Ln.append(Element("Ln", [CT[ind][0], RA_epi[Closest_LA[ind]][0][0]], 34))

Extra += len(Elems_Ln)

Elems_2 = np.concatenate([Elems_2, Elems_Ln])


PM = np.argwhere(Scalar_Color == 9)


P_x_epi = []
P_y_epi = []
P_z_epi = []

for ind1 in range(len(PM)):
    P1 = PM[ind1]
    P_x_epi.append(xlist[P1[0]])
    P_y_epi.append(ylist[P1[0]])
    P_z_epi.append(zlist[P1[0]])

PTS_LA_epi = [P_x_epi, P_y_epi, P_z_epi]
PTS_LA_epi = list(zip(*PTS_LA_epi))


neigh = NearestNeighbors(n_neighbors=1)
neigh.fit(PTS_LA_endo)
_, Closest_LA = neigh.kneighbors(PTS_LA_epi, return_distance=True)


Elems_Ln = []
for ind in range(len(PM)):
    Elems_Ln.append(Element("Ln", [PM[ind][0], RA_epi[Closest_LA[ind]][0][0]], 35))

Extra += len(Elems_Ln)

Elems_2 = np.concatenate([Elems_2, Elems_Ln])


# add fibres too

# now add line fibers

Fibres_XYZ_X1 = []
Fibres_XYZ_Y1 = []
Fibres_XYZ_Z1 = []

for loop in fiber:
    Fibres_XYZ_X1.append(loop[0])
    Fibres_XYZ_Y1.append(loop[1])
    Fibres_XYZ_Z1.append(loop[2])


Fibres_XYZ_1 = [Fibres_XYZ_X1, Fibres_XYZ_Y1, Fibres_XYZ_Z1]
Fibres_XYZ_1 = list(zip(*Fibres_XYZ_1))

for ind in range(Extra):
    Fibres_XYZ_X1.append(1)
    Fibres_XYZ_Y1.append(0)
    Fibres_XYZ_Z1.append(0)


Fibres_XYZ_1 = [Fibres_XYZ_X1, Fibres_XYZ_Y1, Fibres_XYZ_Z1]
Fibres_XYZ_1 = list(zip(*Fibres_XYZ_1))


mname = base_dir_bi + "Bilayer_Combined_all_Lines"
utils.write_carp(pts, Elems_2, Fibres_XYZ_1, None, mname)

pts, elems, fiber, _ = utils.read_carp(base_dir_bi, "Bilayer_Combined", return_surface=False)
M_XYZ_1 = utils.mp_calc_ele_3(pts, elems)
mname = base_dir_bi + "Aux_2"
utils.write_carp(M_XYZ_1, elems, fiber, None, mname)
shutil.copyfile(base_dir_bi + "Aux_2.pts", base_dir_bi + "All_Fibres.vpts")
utils.write_carp(pts, elems, fiber, None, mname)
shutil.copyfile(base_dir_bi + "Aux_2.lon", base_dir_bi + "All_Fibres.vec")


# repeat for each structure

pts, elems, fiber, _ = utils.read_carp(base_dir_bi, "LA_endo", return_surface=False)
M_XYZ_1 = utils.mp_calc_ele_3(pts, elems)
mname = base_dir_bi + "Aux_2"
utils.write_carp(M_XYZ_1, elems, fiber, None, mname)
shutil.copyfile(base_dir_bi + "Aux_2.pts", base_dir_bi + "LA_endo.vpts")
utils.write_carp(pts, elems, fiber, None, mname)
shutil.copyfile(base_dir_bi + "Aux_2.lon", base_dir_bi + "LA_endo.vec")


pts, elems, fiber, _ = utils.read_carp(base_dir_bi, "LA_epi", return_surface=False)
M_XYZ_1 = utils.mp_calc_ele_3(pts, elems)
mname = base_dir_bi + "Aux_2"
utils.write_carp(M_XYZ_1, elems, fiber, None, mname)
shutil.copyfile(base_dir_bi + "Aux_2.pts", base_dir_bi + "LA_epi.vpts")
utils.write_carp(pts, elems, fiber, None, mname)
shutil.copyfile(base_dir_bi + "Aux_2.lon", base_dir_bi + "LA_epi.vec")


pts, elems, fiber, _ = utils.read_carp(base_dir_bi, "RA_epi_s", return_surface=False)
M_XYZ_1 = utils.mp_calc_ele_3(pts, elems)
mname = base_dir_bi + "Aux_2"
utils.write_carp(M_XYZ_1, elems, fiber, None, mname)
shutil.copyfile(base_dir_bi + "Aux_2.pts", base_dir_bi + "RA_epi_s.vpts")
utils.write_carp(pts, elems, fiber, None, mname)
shutil.copyfile(base_dir_bi + "Aux_2.lon", base_dir_bi + "RA_epi_s.vec")


pts, elems, fiber, _ = utils.read_carp(base_dir_bi, "RA_structures", return_surface=False)
M_XYZ_1 = utils.mp_calc_ele_3(pts, elems)
mname = base_dir_bi + "Aux_2"
utils.write_carp(M_XYZ_1, elems, fiber, None, mname)
shutil.copyfile(base_dir_bi + "Aux_2.pts", base_dir_bi + "RA_structures.vpts")
utils.write_carp(pts, elems, fiber, None, mname)
shutil.copyfile(base_dir_bi + "Aux_2.lon", base_dir_bi + "RA_structures.vec")
