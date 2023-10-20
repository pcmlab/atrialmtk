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
my_parser.add_argument("la_path", metavar="path", type=str, help="Path for the LA mesh")
my_parser.add_argument("ra_path", metavar="path", type=str, help="Path for the RA mesh")
my_parser.add_argument("target_path", metavar="path", type=str, help="Path for the new target mesh")
my_parser.add_argument("mesh_name", metavar="filename", type=str, help="Mesh name target (usually Labelled)")
my_parser.add_argument("mesh_name_ra", metavar="filename", type=str, help="Mesh name RA (usually Labelled)")
my_parser.add_argument("output_file_name", metavar="filename", type=str, help="Fibre output file name")
my_parser.add_argument("seed_landmarks", metavar="filename", type=str, help="Seed file name (e.g. Landmarks.txt)")
my_parser.add_argument("seed_region", metavar="filename", type=str, help="Seed file name (e.g. Regions.txt)")
my_parser.add_argument("scale_factor", metavar="scale factor", type=int, help="Scale factor for landmarks")
my_parser.add_argument("fo_index", metavar="fo_index", type=int, help="Index of FO landmark")

# Execute parse_args()
args = my_parser.parse_args()

base_dir = args.la_path
base_dir_bi = args.target_path
mesh_name = args.mesh_name
seed_landmarks = args.seed_landmarks
scale_factor = args.scale_factor
fo_index = args.fo_index

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
# CS_epi = np.argwhere(Scalar_Color == 5)

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

# add BB
BB = np.argwhere(Scalar_Color == 10)

P_x_BB = []
P_y_BB = []
P_z_BB = []

for ind1 in range(len(BB)):
    P1 = BB[ind1]
    P_x_BB.append(xlist[P1[0]])
    P_y_BB.append(ylist[P1[0]])
    P_z_BB.append(zlist[P1[0]])

PTS_BB = [P_x_BB, P_y_BB, P_z_BB]
PTS_BB = list(zip(*PTS_BB))

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

neigh_LA = NearestNeighbors(n_neighbors=1)
neigh_LA.fit(PTS_LA_epi)
LA_D, Closest_LA = neigh_LA.kneighbors(PTS_BB, return_distance=True)

P_x_epi = []
P_y_epi = []
P_z_epi = []

for ind1 in range(len(RA_epi)):
    P1 = RA_epi[ind1]
    P_x_epi.append(xlist[P1[0]])
    P_y_epi.append(ylist[P1[0]])
    P_z_epi.append(zlist[P1[0]])

PTS_RA_epi = [P_x_epi, P_y_epi, P_z_epi]
PTS_RA_epi = list(zip(*PTS_RA_epi))


neigh_RA = NearestNeighbors(n_neighbors=1)
neigh_RA.fit(PTS_RA_epi)
RA_D, Closest_RA = neigh_RA.kneighbors(PTS_BB, return_distance=True)


Elems_Ln = []
BB_LA = []
BB_RA = []

for ind in range(len(BB)):
    if LA_D[ind] > RA_D[ind]:
        B = RA_epi[Closest_RA[ind]]
        BB_RA.append(BB[ind])
    else:
        B = LA_epi[Closest_LA[ind]]
        BB_LA.append(BB[ind])
    Elems_Ln.append(Element("Ln", [BB[ind][0], B[0][0]], 41))

Extra += len(Elems_Ln)

Elems_2 = np.concatenate([Elems_2, Elems_Ln])

# now add lines for closest 1000 points
P_x_BB_LA = []
P_y_BB_LA = []
P_z_BB_LA = []

for ind1 in range(len(BB_LA)):
    P1 = BB_LA[ind1]
    P_x_BB_LA.append(xlist[P1[0]])
    P_y_BB_LA.append(ylist[P1[0]])
    P_z_BB_LA.append(zlist[P1[0]])

PTS_BB_LA = [P_x_BB_LA, P_y_BB_LA, P_z_BB_LA]
PTS_BB_LA = list(zip(*PTS_BB_LA))

P_x_BB_RA = []
P_y_BB_RA = []
P_z_BB_RA = []

for ind1 in range(len(BB_RA)):
    P1 = BB_RA[ind1]
    P_x_BB_RA.append(xlist[P1[0]])
    P_y_BB_RA.append(ylist[P1[0]])
    P_z_BB_RA.append(zlist[P1[0]])

PTS_BB_RA = [P_x_BB_RA, P_y_BB_RA, P_z_BB_RA]
PTS_BB_RA = list(zip(*PTS_BB_RA))

# find distances, and store closest
# 2. join each node to closest node on LA_epi
neigh_RA = NearestNeighbors(n_neighbors=1)
neigh_RA.fit(PTS_BB_RA)

D_Store = []

for ind in range(len(P_x_BB_LA)):
    LA_D, _ = neigh_RA.kneighbors([PTS_BB_LA[ind]], return_distance=True)
    D_Store.append(LA_D[0][0])

print(D_Store)
List_To_Check = np.argsort(D_Store)
ToKeepLength = 100
print(List_To_Check)

Elems_Ln = []
Set_RA = []
Set_LA = []

for indover in range(ToKeepLength):
    ind = List_To_Check[indover]
    print(ind)
    _, Closest_RA = neigh_RA.kneighbors([PTS_BB_LA[ind]], return_distance=True)
    index = BB_RA[Closest_RA[0][0]][0]
    if index not in Set_RA:
        Elems_Ln.append(Element("Ln", [index, BB_LA[ind][0]], 45))
        Set_RA.append(index)
        Set_LA.append(BB_LA[ind][0])

# reverse it too
neigh_LA = NearestNeighbors(n_neighbors=1)
neigh_LA.fit(PTS_BB_LA)

D_Store = []

for ind in range(len(P_x_BB_RA)):
    LA_D, _ = neigh_LA.kneighbors([PTS_BB_RA[ind]], return_distance=True)
    D_Store.append(LA_D[0][0])

print(D_Store)
List_To_Check = np.argsort(D_Store)
print(List_To_Check)


Set_RA = []
Set_LA = []

for indover in range(ToKeepLength):
    ind = List_To_Check[indover]
    print(ind)
    _, Closest_RA = neigh_LA.kneighbors([PTS_BB_RA[ind]], return_distance=True)
    index = BB_LA[Closest_RA[0][0]][0]
    if index not in Set_RA and BB_RA[ind][0] not in Set_LA:
        Elems_Ln.append(Element("Ln", [index, BB_RA[ind][0]], 45))
        Set_RA.append(index)
        Set_LA.append(BB_RA[ind][0])

Extra += len(Elems_Ln)

Elems_2 = np.concatenate([Elems_2, Elems_Ln])


# add fibres too


# now add FO line connections

# now add FO line connections

sourcePt = np.loadtxt(base_dir + seed_landmarks, delimiter=",") * scale_factor

LSPV_Marker = sourcePt[fo_index]

# find LA points within a distance radii
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

neigh = NearestNeighbors(n_neighbors=100)
neigh.fit(PTS_LA_epi)
_, Closest_LA = neigh.kneighbors([LSPV_Marker], return_distance=True)
Closest_LA = Closest_LA[0]

P_x_epi = []
P_y_epi = []
P_z_epi = []

for ind1 in range(len(RA_epi)):
    P1 = RA_epi[ind1]
    P_x_epi.append(xlist[P1[0]])
    P_y_epi.append(ylist[P1[0]])
    P_z_epi.append(zlist[P1[0]])

PTS_RA_epi = [P_x_epi, P_y_epi, P_z_epi]
PTS_RA_epi = list(zip(*PTS_RA_epi))


neigh = NearestNeighbors(n_neighbors=1)
neigh.fit(PTS_RA_epi)

Elems_Ln = []

for ind in range(100):
    print(Closest_LA[ind])
    _, Closest_LA1 = neigh.kneighbors([PTS_LA_epi[Closest_LA[ind]]], return_distance=True)
    print(LA_epi[Closest_LA[ind]][0])
    Elems_Ln.append(Element("Ln", [RA_epi[Closest_LA1[0][0]][0], LA_epi[Closest_LA[ind]][0]], 39))

Extra += len(Elems_Ln)

Elems_2 = np.concatenate([Elems_2, Elems_Ln])


# now add CS connections...

# 1. find rim of CS

# add fibres too

# now add line fibers


# 3. find half that is closest to LA_epi - only...


# now add line fibers

Fibres_XYZ_X1 = []
Fibres_XYZ_Y1 = []
Fibres_XYZ_Z1 = []

for loop in fiber:
    Fibres_XYZ_X1.append(loop[0])
    Fibres_XYZ_Y1.append(loop[1])
    Fibres_XYZ_Z1.append(loop[2])


for ind in range(Extra):
    Fibres_XYZ_X1.append(1)
    Fibres_XYZ_Y1.append(0)
    Fibres_XYZ_Z1.append(0)


Fibres_XYZ_1 = [Fibres_XYZ_X1, Fibres_XYZ_Y1, Fibres_XYZ_Z1]
Fibres_XYZ_1 = list(zip(*Fibres_XYZ_1))

utils.write_carp(pts, Elems_2, Fibres_XYZ_1, None, base_dir_bi + "Bilayer_Combined_all_Lines_IAC")


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


pts, elems, fiber, _ = utils.read_carp(base_dir_bi, "BB_mesh", return_surface=False)
M_XYZ_1 = utils.mp_calc_ele_3(pts, elems)
mname = base_dir_bi + "Aux_2"
utils.write_carp(M_XYZ_1, elems, fiber, None, mname)
shutil.copyfile(base_dir_bi + "Aux_2.pts", base_dir_bi + "BB_mesh.vpts")
utils.write_carp(pts, elems, fiber, None, mname)
shutil.copyfile(base_dir_bi + "Aux_2.lon", base_dir_bi + "BB_mesh.vec")
