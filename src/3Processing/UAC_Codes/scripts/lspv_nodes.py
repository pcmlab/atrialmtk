#!/usr/bin/python3
#
import os
import sys
import argparse
import numpy as np
from sklearn.neighbors import NearestNeighbors
from uac import utils


# Create the parser
my_parser = argparse.ArgumentParser(description="Calculation of labels using laplace solves only")

# Add the arguments
my_parser.add_argument("target_path", metavar="path", type=str, help="Path for the new target mesh")
my_parser.add_argument("mesh_name", metavar="filename", type=str, help="Mesh name (usually Labelled)")
my_parser.add_argument("la_label", metavar="LA number", type=int, help="LA label e.g. 11")
my_parser.add_argument("laa_label", metavar="LAA number", type=int, help="LAA label e.g. 13")
my_parser.add_argument("lspv_label", metavar="LSPV number", type=int, help="LSPV label e.g. 21")
my_parser.add_argument("lipv_label", metavar="LIPV number", type=int, help="LIPV label e.g. 23")
my_parser.add_argument("rspv_label", metavar="RSPV number", type=int, help="RSPV label e.g. 25")
my_parser.add_argument("ripv_label", metavar="RIPV number", type=int, help="RIPV label e.g. 27")

# Execute parse_args()
args = my_parser.parse_args()

base_dir = args.target_path
lspv_label = args.lspv_label
lipv_label = args.lipv_label
rspv_label = args.rspv_label
ripv_label = args.ripv_label

if not os.path.isdir(base_dir):
    print("The target mesh path specified does not exist")
    sys.exit()


# read carp, write to vtk file
print(base_dir)

pts = utils.read_pts(base_dir + "Labelled")
elems = utils.read_elem(base_dir + "Labelled")

xlist = []
ylist = []
zlist = []

for loop in pts:
    xlist.append(loop[0])
    ylist.append(loop[1])
    zlist.append(loop[2])

# find nodes at PV
*_, Nodes_PV1 = utils.find_nodes(elems, pts, 1)
*_, Nodes_PV2 = utils.find_nodes(elems, pts, 2)
*_, Nodes_PV3 = utils.find_nodes(elems, pts, 3)
*_, Nodes_PV4 = utils.find_nodes(elems, pts, 4)


xlist_PV_1 = []
ylist_PV_1 = []
zlist_PV_1 = []

for ind in range(len(Nodes_PV1)):
    xlist_PV_1.append(xlist[Nodes_PV1[ind]])
    ylist_PV_1.append(ylist[Nodes_PV1[ind]])
    zlist_PV_1.append(zlist[Nodes_PV1[ind]])


xlist_PV_2 = []
ylist_PV_2 = []
zlist_PV_2 = []

for ind in range(len(Nodes_PV2)):
    xlist_PV_2.append(xlist[Nodes_PV2[ind]])
    ylist_PV_2.append(ylist[Nodes_PV2[ind]])
    zlist_PV_2.append(zlist[Nodes_PV2[ind]])


xlist_PV_3 = []
ylist_PV_3 = []
zlist_PV_3 = []

for ind in range(len(Nodes_PV3)):
    xlist_PV_3.append(xlist[Nodes_PV3[ind]])
    ylist_PV_3.append(ylist[Nodes_PV3[ind]])
    zlist_PV_3.append(zlist[Nodes_PV3[ind]])


xlist_PV_4 = []
ylist_PV_4 = []
zlist_PV_4 = []

for ind in range(len(Nodes_PV4)):
    xlist_PV_4.append(xlist[Nodes_PV4[ind]])
    ylist_PV_4.append(ylist[Nodes_PV4[ind]])
    zlist_PV_4.append(zlist[Nodes_PV4[ind]])


# use region labels to relabel each one...


ELE_TAG = []
for e in elems:
    ELE_TAG.append(e.reg)

MidPoints = [e.centre(pts) for e in elems]
print((MidPoints[0]))


PTS_PV1 = [xlist_PV_1, ylist_PV_1, zlist_PV_1]
PTS_PV1 = list(zip(*PTS_PV1))
PTS_PV2 = [xlist_PV_2, ylist_PV_2, zlist_PV_2]
PTS_PV2 = list(zip(*PTS_PV2))
PTS_PV3 = [xlist_PV_3, ylist_PV_3, zlist_PV_3]
PTS_PV3 = list(zip(*PTS_PV3))
PTS_PV4 = [xlist_PV_4, ylist_PV_4, zlist_PV_4]
PTS_PV4 = list(zip(*PTS_PV4))


# for each, find closest element mid-point
neigh = NearestNeighbors(n_neighbors=1)
neigh.fit(MidPoints)
Closest_PV1 = neigh.kneighbors(PTS_PV1, return_distance=False)

# look up tags
ELE_TAG_PV1 = []
for ind in range(len(Closest_PV1)):
    Toc = Closest_PV1[ind]
    Toc = Toc[0]
    ELE_TAG_PV1.append(ELE_TAG[Toc])
print(ELE_TAG_PV1)


# for each, find closest element mid-point
Closest_PV2 = neigh.kneighbors(PTS_PV2, return_distance=False)

# look up tags
ELE_TAG_PV2 = []
for ind in range(len(Closest_PV2)):
    Toc = Closest_PV2[ind]
    Toc = Toc[0]
    ELE_TAG_PV2.append(ELE_TAG[Toc])
print(ELE_TAG_PV2)

# for each, find closest element mid-point
Closest_PV3 = neigh.kneighbors(PTS_PV3, return_distance=False)

# look up tags
ELE_TAG_PV3 = []
for ind in range(len(Closest_PV3)):
    Toc = Closest_PV3[ind]
    Toc = Toc[0]
    ELE_TAG_PV3.append(ELE_TAG[Toc])
print(ELE_TAG_PV3)


# for each, find closest element mid-point
Closest_PV4 = neigh.kneighbors(PTS_PV4, return_distance=False)

# look up tags
ELE_TAG_PV4 = []
for ind in range(len(Closest_PV4)):
    Toc = Closest_PV4[ind]
    Toc = Toc[0]
    ELE_TAG_PV4.append(ELE_TAG[Toc])
print(ELE_TAG_PV4)


# use to order regions


Ind_PVs = [Nodes_PV1, Nodes_PV2, Nodes_PV3, Nodes_PV4]
SP = [ELE_TAG_PV1[0], ELE_TAG_PV2[0], ELE_TAG_PV3[0], ELE_TAG_PV4[0]]
PV_choices = [SP.index(lspv_label), SP.index(lipv_label), SP.index(rspv_label), SP.index(ripv_label)]
print(PV_choices)

Ind_LSPV = Ind_PVs[PV_choices[0]]
Ind_LIPV = Ind_PVs[PV_choices[1]]
Ind_RSPV = Ind_PVs[PV_choices[2]]
Ind_RIPV = Ind_PVs[PV_choices[3]]

Ind_LSPV_S = []
for ind in range(len(Ind_LSPV)):
    Toc = Ind_LSPV[ind] + len(xlist)
    Ind_LSPV_S.append(Toc)

Ind_LIPV_S = []
for ind in range(len(Ind_LIPV)):
    Toc = Ind_LIPV[ind] + len(xlist)
    Ind_LIPV_S.append(Toc)

Ind_RSPV_S = []
for ind in range(len(Ind_RSPV)):
    Toc = Ind_RSPV[ind] + len(xlist)
    Ind_RSPV_S.append(Toc)

Ind_RIPV_S = []
for ind in range(len(Ind_RIPV)):
    Toc = Ind_RIPV[ind] + len(xlist)
    Ind_RIPV_S.append(Toc)


Ind_LSPV_T = np.concatenate([Ind_LSPV, Ind_LSPV_S])
Ind_LIPV_T = np.concatenate([Ind_LIPV, Ind_LIPV_S])
Ind_RSPV_T = np.concatenate([Ind_RSPV, Ind_RSPV_S])
Ind_RIPV_T = np.concatenate([Ind_RIPV, Ind_RIPV_S])

# find MV nodes
utils.write_vtx_intra(Ind_LSPV, base_dir + "LSPV_nodes.vtx")
utils.write_vtx_intra(Ind_LIPV, base_dir + "LIPV_nodes.vtx")
utils.write_vtx_intra(Ind_RSPV, base_dir + "RSPV_nodes.vtx")
utils.write_vtx_intra(Ind_RIPV, base_dir + "RIPV_nodes.vtx")

utils.write_vtx_intra(Ind_LSPV_T, base_dir + "Bi_LSPV_nodes.vtx")
utils.write_vtx_intra(Ind_LIPV_T, base_dir + "Bi_LIPV_nodes.vtx")
utils.write_vtx_intra(Ind_RSPV_T, base_dir + "Bi_RSPV_nodes.vtx")
utils.write_vtx_intra(Ind_RIPV_T, base_dir + "Bi_RIPV_nodes.vtx")
