#!/usr/bin/python3
#
import os
import sys
import argparse
import numpy as np
from uac import utils


# Create the parser
my_parser = argparse.ArgumentParser(description="Calculation of old atrial coordinates for RA")

# Add the arguments
my_parser.add_argument("target_path", metavar="path", type=str, help="Path for the new target mesh")
my_parser.add_argument("mesh_name", metavar="filename", type=str, help="Mesh name (usually Labelled)")
my_parser.add_argument("seed_landmarks", metavar="filename", type=str, help="Seed file name (e.g. Landmarks.txt)")
my_parser.add_argument("seed_region", metavar="filename", type=str, help="Seed file name (e.g. Region.txt)")
my_parser.add_argument("scale_factor", metavar="scale factor", type=int, help="Scale factor for landmarks")

# Execute parse_args()
args = my_parser.parse_args()

base_dir = args.target_path
mesh_name = args.mesh_name
seed_landmarks = args.seed_landmarks
seed_region = args.seed_region
scale_factor = args.scale_factor

if not os.path.isdir(base_dir):
    print("The target mesh path specified does not exist")
    sys.exit()


print(base_dir)
pts, elems, *_, surface = utils.read_carp(base_dir, mesh_name)

print(pts)

pts_src = np.loadtxt(base_dir + seed_landmarks, delimiter=",") * scale_factor
pts_region = np.loadtxt(base_dir + seed_region, delimiter=",") * scale_factor

size = 3
nodes_pv = [None] * size
pts_pv = [None] * size
for i in range(size):
    *_, nodes_pv[i] = utils.find_nodes(elems, pts, i, True)

    pts_pv[i] = np.empty((len(nodes_pv[i]), 3), float)
    for ind in range(len(nodes_pv[i])):
        pts_pv[i][ind] = pts[nodes_pv[i][ind]]

pts_src_remap = [pts_region[0], pts_region[3]]
print(pts_src_remap)
size = len(pts_src_remap)
sp = [None] * size
for i in range(size):
    sp[i] = utils.pv_index_remap(pts_src_remap[i], pts_pv, np.mean)

print(sp)

nodes_rspv = nodes_pv[sp[0]]
nodes_lspv = nodes_pv[sp[1]]

print(len(nodes_rspv))
print(len(nodes_lspv))


# set MV as other choice...
considered_choices = [0, 1, 2]
mv_choice = [x for x in considered_choices if x not in sp]
print(mv_choice)
nodes_mv = nodes_pv[mv_choice[0]]


# add in LAA...
pts_lspv = pts[nodes_lspv]
index = ((pts_lspv - pts_src[0]) ** 2).sum(1).argmin()
p1 = pts[nodes_lspv[index]]
index = ((pts_lspv - pts_src[4]) ** 2).sum(1).argmin()
p2 = pts[nodes_lspv[index]]


# 2. LSPV find opposite point in rim
tot_pts, surface_filter = utils.read_vtk_surface(surface)
surface_filter_sub_lspv = utils.get_sub_surface(base_dir, pts, elems, nodes_lspv)
path_pts, _ = utils.geod_lspv_rspv(surface_filter_sub_lspv, p2, p1, tot_pts)
mid_point = path_pts[len(path_pts) // 2]

increment = 1

index = utils.geod_midpoint(surface_filter_sub_lspv, mid_point, nodes_lspv, increment, tot_pts)
mp = pts[index]
tpts1, _ = utils.geod_lspv_rspv(surface_filter_sub_lspv, mp, p1, tot_pts)
tpts2, _ = utils.geod_lspv_rspv(surface_filter_sub_lspv, p2, mp, tot_pts)
at = np.concatenate([tpts1, tpts2])


# keep path closest to marker

print(p1)
print(mp)
print(p2)
print(at[0])
print(at[len(at) // 2])
print(at[-1])


# code to select correct path
lspv_marker_pts = pts_region[3]
m1 = ((path_pts - lspv_marker_pts) ** 2).sum(1).min()
m2 = ((at - lspv_marker_pts) ** 2).sum(1).min()


# reassign to rims.

# for each node on the boundary, find closest node point and assign.
# first - boundary nodes: Nodes_LA_LIPV, Nodes_LA_RIPV

if m1 < m2:
    path_lspv = path_pts
else:
    path_lspv = at


# now do for subsets and build up - increase number of paths

path_inds = [None] * 8
for i in range(len(path_inds)):
    pt1 = p2 if i + 1 == len(path_inds) else path_lspv[(i + 1) * len(path_lspv) // 8]
    pt2 = p1 if i == 0 else path_lspv[i * len(path_lspv) // 8]
    _, path_inds[i] = utils.geod_lspv_rspv(surface_filter, pt1, pt2, tot_pts)

path_lspv_inds = np.concatenate(path_inds)


# 3. RSPV find subdividers in rims

pts_rspv = pts[nodes_rspv]

index = ((pts_rspv - pts_src[1]) ** 2).sum(1).argmin()
p3 = pts[nodes_rspv[index]]
print(p3)

index = ((pts_rspv - pts_src[5]) ** 2).sum(1).argmin()
p4 = pts[nodes_rspv[index]]
print(p4)


# now find two subdividers in rim and join these...
surface_filter_sub_rspv = utils.get_sub_surface(base_dir, pts, elems, nodes_rspv)
path_pts, _ = utils.geod_lspv_rspv(surface_filter_sub_rspv, p4, p3, tot_pts)
mid_point = path_pts[len(path_pts) // 2]
print(mid_point)

found_index = utils.geod_midpoint(surface_filter_sub_rspv, mid_point, nodes_rspv, increment, tot_pts)
mp = tot_pts[found_index]
tpts1, _ = utils.geod_lspv_rspv(surface_filter_sub_rspv, mp, p3, tot_pts)
tpts2, _ = utils.geod_lspv_rspv(surface_filter_sub_rspv, p4, mp, tot_pts)
at = np.concatenate([tpts1, tpts2])

rpsv_marker_pts = pts_region[0]
m3 = ((path_pts - rpsv_marker_pts) ** 2).sum(1).min()
m4 = ((at - rpsv_marker_pts) ** 2).sum(1).min()

if m3 < m4:
    path_rspv = path_pts
else:
    path_rspv = at

path_pts = [None] * 8
path_inds = [None] * 8
for i in range(len(path_inds)):
    pt1 = p4 if i + 1 == len(path_inds) else path_rspv[(i + 1) * len(path_rspv) // 8]
    pt2 = p3 if i == 0 else path_rspv[i * len(path_rspv) // 8]
    _, path_inds[i] = utils.geod_lspv_rspv(surface_filter, pt1, pt2, tot_pts)

path_rspv_inds = np.concatenate(path_inds)


# 4. geodesic on roof
_, path_roof_inds = utils.geod_lspv_rspv(surface_filter, p2, p4, tot_pts)


# also save alternate paths to ensure other nodes are in other path
# enforce that lspv, rspv paths are on


# 5. Add LSPV and RSPV paths
total_roof_inds = np.concatenate([path_lspv_inds, path_roof_inds, path_rspv_inds])

# 6. Laplace solve for UD

utils.write_vtx_extra(nodes_mv, base_dir + "PAbc1.vtx")
utils.write_vtx_extra(total_roof_inds, base_dir + "PAbc2.vtx")

lspv_marker_pts = pts_src[2]
rspv_marker_pts = pts_src[3]

_, path_lspv_inds = utils.geod_lspv_rspv(surface_filter, lspv_marker_pts, p1, tot_pts)
_, path_rspv_inds = utils.geod_lspv_rspv(surface_filter, rspv_marker_pts, p3, tot_pts)

lspv_marker_inds = utils.geod_marker_mv(surface_filter, lspv_marker_pts, increment, tot_pts, nodes_mv)
rspv_marker_inds = utils.geod_marker_mv(surface_filter, rspv_marker_pts, increment, tot_pts, nodes_mv)

lspv_inds = np.concatenate((lspv_marker_inds, path_lspv_inds), axis=None)
rspv_inds = np.concatenate((rspv_marker_inds, path_rspv_inds), axis=None)

left_nodes = lspv_inds
right_nodes = np.flip(rspv_inds, 0)
roof_nodes = total_roof_inds

utils.write_vtx_extra(left_nodes, base_dir + "LSbc1.vtx")
utils.write_vtx_extra(right_nodes, base_dir + "LSbc2.vtx")

nodes_border = np.concatenate([left_nodes, roof_nodes, right_nodes])
utils.write_vtx_extra(nodes_border, base_dir + "BorderNodes.vtx")
