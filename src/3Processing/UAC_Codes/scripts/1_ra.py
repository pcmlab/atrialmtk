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
my_parser.add_argument("fibre_path", metavar="path", type=str, help="Path for the atlas mesh")
my_parser.add_argument("laplace_path", metavar="path", type=str, help="Path for the laplace par files")
my_parser.add_argument("mesh_name", metavar="filename", type=str, help="Mesh name (usually Labelled)")
my_parser.add_argument("ra_label", metavar="RA number", type=int, help="RA body label (default 1)")
my_parser.add_argument("svc_label", metavar="SVC number", type=int, help="SVC label (default 6)")
my_parser.add_argument("ivc_label", metavar="IVC number", type=int, help="IVC label (default 7)")
my_parser.add_argument("cs_label", metavar="CS number", type=int, help="CS label (default 5)")
my_parser.add_argument("raa_label", metavar="RAA number", type=int, help="RAA label (default 2)")
my_parser.add_argument(
    "seed_landmarks", metavar="filename", type=str, help="Seed file name (e.g. seedsfileOUT_Landmarks.vtk)"
)
my_parser.add_argument(
    "seed_region", metavar="filename", type=str, help="Seed file name (e.g. seedsfileOUT_Region.vtk)"
)
my_parser.add_argument("scale_factor", metavar="scale factor", type=int, help="Scale factor for landmarks")

# Execute parse_args()
args = my_parser.parse_args()

base_dir = args.target_path
laplace_dir = args.laplace_path
fibre_dir = args.fibre_path
mesh_name = args.mesh_name
ra_label = args.ra_label  # default = 1
svc_label = args.svc_label  # default = 6
ivc_label = args.ivc_label  # default = 7
seed_landmarks = args.seed_landmarks
seed_region = args.seed_region
scale_factor = args.scale_factor

if not os.path.isdir(base_dir):
    print("The target mesh path specified does not exist")
    sys.exit()

if not os.path.isdir(laplace_dir):
    print("The atlas path specified does not exist")
    sys.exit()

if not os.path.isdir(fibre_dir):
    print("The laplace files path specified does not exist")
    sys.exit()


print("Reading seeds file")
pts_src = utils.load_landmarks(base_dir + seed_landmarks, scale_factor)

print("read carp, write to vtk file")
print(base_dir)
pts, elems, *_, surface = utils.read_carp(base_dir, mesh_name)

print(pts)

print("find nodes on RA/PV junctions")
nodes_lspv = utils.find_junc_nodes(elems, ra_label, svc_label)  # (elems,1,6)
nodes_rspv = utils.find_junc_nodes(elems, ra_label, ivc_label)  # (elems,1,7)

# find MV nodes
*_, nodes_mv = utils.find_nodes(elems, pts)

# Read in surfaces
tot_pts, surface_filter = utils.read_vtk_surface(surface)

print(tot_pts[0])
print(pts[0])
print((tot_pts - pts).max())

# instead use markers for this one. pts_src
pts_lspv = pts[nodes_lspv]
index = ((pts_lspv - pts_src[0]) ** 2).sum(1).argmin()
p1 = pts[nodes_lspv[index]]
index = ((pts_lspv - pts_src[4]) ** 2).sum(1).argmin()
p2 = pts[nodes_lspv[index]]


# 2. LSPV find opposite point in rim
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

print("Reading seeds file")
pts_region = utils.load_landmarks(base_dir + seed_region, scale_factor)

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

print(p1)
print(path_lspv[0])
print(path_lspv[len(path_lspv) // 8])
print(path_lspv[7 * len(path_lspv) // 8])
print(path_lspv[-1])
print(p2)

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
tpts1, _ = utils.geod_lspv_rspv(surface_filter_sub_rspv, p3, mp, tot_pts)
tpts2, _ = utils.geod_lspv_rspv(surface_filter_sub_rspv, mp, p4, tot_pts)
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
    path_pts[i], path_inds[i] = utils.geod_lspv_rspv(surface_filter, pt1, pt2, tot_pts)

path_rspv = np.concatenate(path_pts)
path_rspv_inds = np.concatenate(path_inds)

print(path_rspv[0])
print(path_rspv[-1])
print(path_rspv[len(path_rspv) // 2])


# 4. geodesic on roof
_, path_roof_inds = utils.geod_lspv_rspv(surface_filter, p2, p4, tot_pts)


# also save alternate paths to ensure other nodes are in other path
# enforce that lspv, rspv paths are on


# 5. Add LSPV and RSPV paths
total_roof_inds = np.concatenate([path_lspv_inds, path_roof_inds, path_rspv_inds])

# 6. Laplace solve for UD

utils.write_vtx_intra(nodes_mv, base_dir + "PAbc1.vtx")
utils.write_vtx_intra(total_roof_inds, base_dir + "PAbc2.vtx")

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

utils.write_vtx_intra(left_nodes, base_dir + "LSbc1.vtx")
utils.write_vtx_intra(right_nodes, base_dir + "LSbc2.vtx")

nodes_border = np.concatenate([left_nodes, roof_nodes, right_nodes])
utils.write_vtx_intra(nodes_border, base_dir + "BorderNodes.vtx")
