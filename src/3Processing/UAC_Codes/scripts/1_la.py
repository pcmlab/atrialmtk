#!/usr/bin/python3
#
import os
import sys
import argparse
import numpy as np
from uac import utils


# Create the parser
my_parser = argparse.ArgumentParser(description="Calculation of old atrial coordinates")

# Add the arguments
my_parser.add_argument("target_path", metavar="path", type=str, help="Path for the new target mesh")
my_parser.add_argument("fibre_path", metavar="path", type=str, help="Path for the atlas mesh")
my_parser.add_argument("laplace_path", metavar="path", type=str, help="Path for the laplace par files")
my_parser.add_argument("mesh_name", metavar="filename", type=str, help="Mesh name (usually Labelled)")
my_parser.add_argument("la_label", metavar="LA number", type=int, help="LA label e.g. 11")
my_parser.add_argument("laa_label", metavar="LAA number", type=int, help="LAA label e.g. 13")
my_parser.add_argument("lspv_label", metavar="LSPV number", type=int, help="LSPV label e.g. 21")
my_parser.add_argument("lipv_label", metavar="LIPV number", type=int, help="LIPV label e.g. 23")
my_parser.add_argument("rspv_label", metavar="RSPV number", type=int, help="RSPV label e.g. 25")
my_parser.add_argument("ripv_label", metavar="RIPV number", type=int, help="RIPV label e.g. 27")
my_parser.add_argument(
    "seed_name", metavar="filename", type=str, help="Seed file name (e.g. seedsfileOUT_Landmarks.vtk)"
)
my_parser.add_argument("scale_factor", metavar="scale factor", type=int, help="Scale factor for landmarks")

# Execute parse_args()
args = my_parser.parse_args()

base_dir = args.target_path
laplace_dir = args.laplace_path
fibre_dir = args.fibre_path
mesh_name = args.mesh_name
la_label = args.la_label
lspv_label = args.lspv_label
lipv_label = args.lipv_label
rspv_label = args.rspv_label
ripv_label = args.ripv_label
seed_name = args.seed_name
scale_factor = args.scale_factor

if not os.path.isdir(base_dir):
    print("The target mesh path specified does not exist")
    sys.exit()

if not os.path.isdir(laplace_dir):
    print("The laplace directory specified does not exist")
    sys.exit()

if not os.path.isdir(fibre_dir):
    print("The fibre file directory path specified does not exist")
    sys.exit()


pts_src = utils.load_landmarks(base_dir + seed_name, scale_factor)

print("read carp, write to vtk file")
print(base_dir)
pts, elems, *_, surface = utils.read_carp(base_dir, mesh_name)

print("find nodes on LA/PV junctions")
nodes_lspv = utils.find_junc_nodes(elems, la_label, lspv_label)
nodes_rspv = utils.find_junc_nodes(elems, la_label, rspv_label)
nodes_lipv = utils.find_junc_nodes(elems, la_label, lipv_label)
nodes_ripv = utils.find_junc_nodes(elems, la_label, ripv_label)

print("find MV nodes")
*_, nodes_mv = utils.find_nodes(elems, pts)

print("Read in surfaces")
tot_pts, surface_filter = utils.read_vtk_surface(surface)

print("1. Find closest points of superior veins to inferior veins. Do this using boundary_markers function.")
_, lipv_lspv_2 = utils.boundary_markers(nodes_lipv, nodes_lspv, tot_pts)
_, ripv_rspv_2 = utils.boundary_markers(nodes_ripv, nodes_rspv, tot_pts)

pts_lspv = pts[nodes_lspv]
p1 = pts[nodes_lspv[lipv_lspv_2]]

print("find furthest point as its conjugate")

pointd = ((pts_lspv - p1) ** 2).sum(1)
found_index = utils.find_index(pointd, max)

print("2. LSPV find opposite point in rim")
surface_filter_sub_lspv = utils.get_sub_surface(base_dir, pts, elems, nodes_lspv)
p2 = pts[nodes_lspv[found_index]]
path_pts, _ = utils.geod_lspv_rspv(surface_filter_sub_lspv, p1, p2, tot_pts)
mid_point = path_pts[len(path_pts) // 2]

increment = 1

index = utils.geod_midpoint(surface_filter_sub_lspv, mid_point, nodes_lspv, increment, tot_pts)
mp = pts[index]
tpts1, _ = utils.geod_lspv_rspv(surface_filter_sub_lspv, p1, mp, tot_pts)
tpts2, _ = utils.geod_lspv_rspv(surface_filter_sub_lspv, mp, p2, tot_pts)
at = np.concatenate([tpts1, tpts2])

print("keep path closest to marker")
# code to select correct path
lspv_marker_pts = pts_src[2]
m1 = ((path_pts - lspv_marker_pts) ** 2).sum(1).min()
m2 = ((at - lspv_marker_pts) ** 2).sum(1).min()

print("reassign to rims.")
# for each node on the boundary, find the closest node point and assign.
# first - boundary nodes: Nodes_LA_LIPV, Nodes_LA_RIPV

path_pts_ls = [path_pts, at]
lspv_indices_path = []
lspv_indices_at = []
d = np.empty(len(path_pts_ls), float)

for ind in range(0, len(pts_lspv)):
    for i in range(len(path_pts_ls)):
        index = utils.closest_node(pts_lspv[ind], path_pts_ls[i]).argmin()
        d[i] = ((pts_lspv[ind] - path_pts_ls[i][index]) ** 2).sum()

    if d[0] < d[1]:
        lspv_indices_path.append(nodes_lspv[ind])
    elif d[0] == d[1]:
        lspv_indices_path.append(nodes_lspv[ind])
        lspv_indices_at.append(nodes_lspv[ind])
    else:
        lspv_indices_at.append(nodes_lspv[ind])

if m1 < m2:
    path_lspv = lspv_indices_path
else:
    path_lspv = lspv_indices_at

# also update RSPV - take path closest to LIPV
print("3. RSPV find subdividers in rims")

p3 = pts[nodes_rspv[ripv_rspv_2]]

# find the furthest point as it's conjugate
pts_rspv = pts[nodes_rspv]
pointd = ((pts_rspv - p3) ** 2).sum(1)
found_index = utils.find_index(pointd, max)

p4 = pts[nodes_rspv[found_index]]

# now find two subdividers in rim and join these...
surface_filter_sub_rspv = utils.get_sub_surface(base_dir, pts, elems, nodes_rspv)
path_pts, _ = utils.geod_lspv_rspv(surface_filter_sub_rspv, p3, p4, tot_pts)
mid_point = path_pts[len(path_pts) // 2]
print(mid_point)

found_index = utils.geod_midpoint(surface_filter_sub_rspv, mid_point, nodes_rspv, increment, tot_pts)
mp = tot_pts[found_index]
tpts1, _ = utils.geod_lspv_rspv(surface_filter_sub_rspv, p3, mp, tot_pts)
tpts2, _ = utils.geod_lspv_rspv(surface_filter_sub_rspv, mp, p4, tot_pts)
at = np.concatenate([tpts1, tpts2])

rspv_marker_pts = pts_src[3]
m3 = ((path_pts - rspv_marker_pts) ** 2).sum(1).min()
m4 = ((at - rspv_marker_pts) ** 2).sum(1).min()

path_pts_ls = [path_pts, at]
rspv_indices_path = []
rspv_indices_at = []
d = np.empty(len(path_pts_ls), float)

for ind in range(0, len(pts_rspv)):
    for i in range(len(path_pts_ls)):
        index = utils.closest_node(pts_rspv[ind], path_pts_ls[i]).argmin()
        d[i] = ((pts_rspv[ind] - path_pts_ls[i][index]) ** 2).sum()

    if d[0] < d[1]:
        rspv_indices_path.append(nodes_rspv[ind])
    elif d[0] == d[1]:
        rspv_indices_path.append(nodes_rspv[ind])
        rspv_indices_at.append(nodes_rspv[ind])
    else:
        rspv_indices_at.append(nodes_rspv[ind])

if m3 < m4:
    path_rspv = rspv_indices_path
else:
    path_rspv = rspv_indices_at

print("4.Geodesic on roof")
_, path_roof_inds = utils.geod_lspv_rspv(surface_filter, p2, p4, tot_pts)

print("5. Add LSPV and RSPV paths")
total_roof_inds = np.concatenate([path_lspv, path_roof_inds, path_rspv])

print("6. Save Output Files for PA and LR")
utils.write_vtx_extra(nodes_mv, base_dir + "PAbc1.vtx")
utils.write_vtx_extra(total_roof_inds, base_dir + "PAbc2.vtx")

lspv_marker = pts_src[0]
rspv_marker = pts_src[1]

_, path_lspv_inds = utils.geod_lspv_rspv(surface_filter, lspv_marker, p1, tot_pts)
_, path_rspv_inds = utils.geod_lspv_rspv(surface_filter, rspv_marker, p3, tot_pts)

index_lspv_mv_marker = utils.geod_marker_mv(surface_filter, lspv_marker, increment, tot_pts, nodes_mv)
index_rspv_mv_marker = utils.geod_marker_mv(surface_filter, rspv_marker, increment, tot_pts, nodes_mv)

index_lspv_mv = np.concatenate((index_lspv_mv_marker, path_lspv_inds), axis=None)
index_rspv_mv = np.concatenate((index_rspv_mv_marker, path_rspv_inds), axis=None)

left_nodes = index_lspv_mv
right_nodes = np.flip(index_rspv_mv, 0)
roof_nodes = total_roof_inds

utils.write_vtx_extra(left_nodes, base_dir + "LSbc1.vtx")
utils.write_vtx_extra(right_nodes, base_dir + "LSbc2.vtx")

nodes_border = np.concatenate([left_nodes, roof_nodes, right_nodes])
utils.write_vtx_extra(nodes_border, base_dir + "BorderNodes.vtx")
