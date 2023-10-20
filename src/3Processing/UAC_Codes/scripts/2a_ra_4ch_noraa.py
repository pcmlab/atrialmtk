#!/usr/bin/python3
#
import os
import sys
import argparse
import numpy as np
from sklearn.neighbors import NearestNeighbors
from uac import utils
import uac.vars as var


# remove inner circle for lspv, rspv. Removed ripv completely

# Create the parser
my_parser = argparse.ArgumentParser(description="Calculation of new atrial coordinates for RA")

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

sys.path.append(base_dir)


# read carp, write to vtk file

# change here onwards - other script for nodes..

print(base_dir)
pts, elems, fiber, data, surface = utils.read_carp(base_dir, mesh_name)

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
mv_choice = [i for i in considered_choices if i not in sp]
print(mv_choice)

"""import igb files"""
pa_ud = utils.read_array_igb(base_dir + "PA_UAC_N2/phie.igb")[0]
ls_lr = utils.read_array_igb(base_dir + "LR_UAC_N2/phie.igb")[0]

# 1. isolate PV mesh
# Edge list for LA/PV junction
tot_pts, surface_filter = utils.read_vtk_surface(surface)

ant_nodes_pa = []
ant_nodes_ls = []
ant_strength_pa = []
ant_strength_ls = []

# rescale y coord of mesh
nodes_border = utils.read_pts(base_dir + "BorderNodes", vtx=True)

elems_left = utils.element2delete(elems, nodes_border, excluded=True)
utils.write_vtk(pts, elems_left, None, None, base_dir + "Test_Split.vtk")
surface = utils.get_vtk_from_file(base_dir + "Test_Split.vtk")


# split into connected components
ripv_marker = pts_region[1]
output = utils.extract_region_by_marker(surface, ripv_marker)
pts_split, _ = utils.read_vtk_surface(output)
surface = utils.convert_unstructureddata_to_polydata(output)
spts_post, selems_post = utils.poly2carp(surface)
utils.write_carp(spts_post, selems_post, None, None, base_dir + "PosteriorMesh")

# for each point in posterior, find in labelled. Take this list and remove to give anterior mesh
# KD tree as much faster
neigh = NearestNeighbors(n_neighbors=1)
neigh.fit(pts)
nodes_ls = neigh.kneighbors(pts_split, return_distance=False)

# remove these
elems_left = utils.element2delete(elems, nodes_ls, excluded=True)
utils.write_vtk(pts, elems_left, None, None, base_dir + "Test_ant.vtk")

# write anterior mesh
surface = utils.get_vtk_from_file(base_dir + "Test_ant.vtk")
surface = utils.convert_unstructureddata_to_polydata(surface)
surface = utils.clean_polydata(surface)
spts_ant, selems_ant = utils.poly2carp(surface)
utils.write_carp(spts_ant, selems_ant, None, None, base_dir + "AnteriorMesh")


# now calculate UAC

pa_ud_r = 1 - pa_ud * 0.5

for i in nodes_ls:
    pa_ud_r[i] = pa_ud[i] * 0.5

pts2d_r = np.empty((len(pts), 3), float)
for i in range(len(pts)):
    pts2d_r[i] = [1 - ls_lr[i], pa_ud_r[i], 0]

mname = base_dir + "Labelled_Coords_2D_Rescaling_N3.vtk"
utils.write_vtk(pts2d_r, elems, fiber, data, mname)

# write 2D carp too
surface = utils.get_vtk_from_file(mname)
surface = utils.convert_unstructureddata_to_polydata(surface)
spts, selems = utils.poly2carp(surface)
mname = mname[:-4]
utils.write_carp(spts, selems, None, None, mname)


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
lspv_marker_pts = pts_region[3]
m1 = ((path_pts - lspv_marker_pts) ** 2).sum(1).min()
m2 = ((at - lspv_marker_pts) ** 2).sum(1).min()

path_lspv_pts = path_pts if m1 < m2 else at
semi_circle_ap_lspv = at[len(at) // 2] if m1 < m2 else path_pts[len(path_pts) // 2]

# now do for subsets and build up - increase number of paths

path_pts = [None] * 8
path_inds = [None] * len(path_pts)
for i in range(len(path_inds)):
    pt1 = p2 if i + 1 == len(path_inds) else path_lspv_pts[(i + 1) * len(path_lspv_pts) // 8]
    pt2 = p1 if i == 0 else path_lspv_pts[i * len(path_lspv_pts) // 8]
    path_pts[i], path_inds[i] = utils.geod_lspv_rspv(surface_filter, pt1, pt2, tot_pts)

path_lspv_pts = np.concatenate(path_pts)
path_lspv_inds = np.concatenate(path_inds)

start_lspv = path_lspv_pts[0]
end_lspv = path_lspv_pts[-1]


nodes_line, strengths_line = utils.lspv2line_ra(var.XSHIFT_21, var.R1_21_25, path_lspv_inds)

post_nodes_ls = nodes_line
post_strength_ls = strengths_line

utils.write_vtx_strength_intra(nodes_line, strengths_line, base_dir + "LSPV_LS_Line")

nodes, strengths_pa, strengths_ls, _ = utils.lspv2circle_ra(
    surface_filter_sub_lspv,
    tot_pts,
    end_lspv,
    semi_circle_ap_lspv,
    start_lspv,
    var.XSHIFT_21,
    var.R1_21_25,
    var.R2_21_25,
)

utils.write_vtx_strength_intra(nodes, strengths_ls, base_dir + "LSPV_LS")
utils.write_vtx_strength_intra(nodes, strengths_pa, base_dir + "LSPV_PA")

nodes_pa = nodes
nodes = np.concatenate([nodes, nodes_line])
strengths_ls = np.concatenate([strengths_ls, strengths_line])
utils.write_vtx_strength_intra(nodes, strengths_ls, base_dir + "LSPV_LS_All")


# for each node on the boundary, find closest node point and assign.
# first - boundary nodes: Nodes_LA_LIPV, Nodes_LA_RIPV

pts_ud2 = set(utils.read_pts(base_dir + "PAbc2", n=1, vtx=True, item_type=int))
unassigneds = [i for i in nodes_lspv if i not in pts_ud2]
pts_lspv_r = pts[unassigneds]
pts_lspv_f = pts[nodes_pa]

pts_lspv_r_len = len(pts_lspv_r)
post_strength_pa = np.empty(pts_lspv_r_len, float)

for i in range(pts_lspv_r_len):
    index = utils.closest_node(pts_lspv_r[i], pts_lspv_f).argmin()
    post_strength_pa[i] = strengths_pa[index]

post_nodes_pa = unassigneds
post_nodes_ls = np.concatenate([post_nodes_ls, nodes, nodes_line])
post_strength_ls = np.concatenate([post_strength_ls, strengths_ls, strengths_line])


# 3. RSPV find subdividers in rims
pts_rspv = pts[nodes_rspv]
index = ((pts_rspv - pts_src[1]) ** 2).sum(1).argmin()
p3 = pts[nodes_rspv[index]]
index = ((pts_rspv - pts_src[5]) ** 2).sum(1).argmin()
p4 = pts[nodes_rspv[index]]

print(p3)
print(p4)

# now find two subdividers in rim and join these...
surface_filter_sub_rspv = utils.get_sub_surface(base_dir, pts, elems, nodes_rspv)
path_pts, _ = utils.geod_lspv_rspv(surface_filter_sub_rspv, p4, p3, tot_pts)
mid_point = path_pts[len(path_pts) // 2]

print(mid_point)

index = utils.geod_midpoint(surface_filter_sub_rspv, mid_point, nodes_rspv, increment, tot_pts)
mp = tot_pts[index]
tpts1, _ = utils.geod_lspv_rspv(surface_filter_sub_rspv, mp, p3, tot_pts)
tpts2, _ = utils.geod_lspv_rspv(surface_filter_sub_rspv, p4, mp, tot_pts)
at = np.concatenate([tpts1, tpts2])

rspv_marker_pts = pts_region[0]
m3 = ((path_pts - rspv_marker_pts) ** 2).sum(1).min()
m4 = ((at - rspv_marker_pts) ** 2).sum(1).min()

path_rspv_pts = path_pts if m3 < m4 else at
semi_circle_ap_rspv = at[len(at) // 2] if m3 < m4 else path_pts[len(path_pts) // 2]

path_pts = [None] * 8
path_inds = [None] * len(path_pts)
for i in range(len(path_inds)):
    pt1 = p4 if i + 1 == len(path_inds) else path_rspv_pts[(i + 1) * len(path_rspv_pts) // 8]
    pt2 = p3 if i == 0 else path_rspv_pts[i * len(path_rspv_pts) // 8]
    path_pts[i], path_inds[i] = utils.geod_lspv_rspv(surface_filter, pt1, pt2, tot_pts)

path_rspv_pts = np.concatenate(path_pts)
path_rspv_inds = np.concatenate(path_inds)

start_rspv = path_rspv_pts[0]
end_rspv = path_rspv_pts[-1]

nodes_line, strengths_line = utils.rspv2line_ra(var.XSHIFT_25, path_rspv_inds, var.R1_21_25)

post_nodes_ls = np.concatenate([post_nodes_ls, nodes_line])
post_strength_ls = np.concatenate([post_strength_ls, strengths_line])

utils.write_vtx_strength_intra(nodes_line, strengths_line, base_dir + "RSPV_Line")

nodes_ls, nodes_pa, strengths_pa, strengths_ls = utils.rspv2circle_ra(
    surface_filter_sub_rspv,
    var.XSHIFT_25,
    tot_pts,
    var.R1_21_25,
    var.R2_21_25,
    semi_circle_ap_rspv,
    end_rspv,
    start_rspv,
)

utils.write_vtx_strength_intra(nodes_ls, strengths_ls, base_dir + "RSPV_LS")
utils.write_vtx_strength_intra(nodes_pa, strengths_pa, base_dir + "RSPV_PA")

nodes_ls = np.concatenate([nodes_ls, nodes_line])
strengths_ls = np.concatenate([strengths_ls, strengths_line])


# for each node on the boundary, find closest node point and assign.
# first - boundary nodes: Nodes_LA_LIPV, Nodes_LA_RIPV

unassigneds = [i for i in nodes_rspv if i not in pts_ud2]
pts_rspv_r = pts[unassigneds]
pts_rspv_f = pts[nodes_ls]
pts_rspv_pf = pts[nodes_pa]

pts_rspv_r_len = len(pts_rspv_r)
lsi = len(post_strength_ls)
pai = len(post_strength_pa)
post_strength_ls.resize(lsi + pts_rspv_r_len)
post_strength_pa.resize(pai + pts_rspv_r_len)

for i in range(pts_rspv_r_len):
    index = utils.closest_node(pts_rspv_r[i], pts_rspv_f).argmin()
    post_strength_ls[lsi + i] = strengths_ls[index]
    index = utils.closest_node(pts_rspv_r[i], pts_rspv_pf).argmin()
    post_strength_pa[pai + i] = strengths_pa[index]

post_nodes_pa = np.concatenate([post_nodes_pa, unassigneds])
post_nodes_ls = np.concatenate([post_nodes_ls, unassigneds])


_, path_top_inds = utils.geod_lspv_rspv(surface_filter, semi_circle_ap_lspv, semi_circle_ap_rspv, tot_pts)

post_nodes_pa = np.concatenate([post_nodes_pa, path_top_inds])
post_strength_pa = np.concatenate([post_strength_pa, np.ones(len(path_top_inds)) * (1 - var.R2_21_25)])

# Split section...
# separate node lists in anteior and posterior for all except old UAC boundaries - must be in both

nodes = [None] * 4
strengths = [None] * len(nodes)
nodes_src = [post_nodes_pa, post_nodes_ls, ant_nodes_pa, ant_nodes_ls]
strengths_src = [post_strength_pa, post_strength_ls, ant_strength_pa, ant_strength_ls]
spts_src = [spts_post, spts_ant]

for i in range(len(nodes)):
    print(len(nodes_src[i]))
    strengths[i] = strengths_src[i][: len(nodes_src[i])]
    nodes[i] = np.empty(len(nodes_src[i]), int)
    for j in range(len(nodes_src[i])):
        pt = tot_pts[nodes_src[i][j]]
        nodes[i][j] = utils.find_index(((spts_src[i // 2] - pt) ** 2).sum(1), min)


filenames = ["LSbc1", "LSbc2", "PAbc1", "PAbc2"]
pts_src = [None] * len(filenames)
for i in range(len(pts_src)):
    pts_src[i] = utils.read_pts(base_dir + filenames[i], n=1, vtx=True, item_type=int)

for i in range(len(pts_src)):
    pts = tot_pts[pts_src[i]]
    strength_values = np.full(len(pts_src[i]), i % 2, int)
    ind = 1 - (i // 2)
    for j in range(2):
        ind += j * 2
        neigh = NearestNeighbors(n_neighbors=1)
        neigh.fit(spts_src[j])
        nodes_res = neigh.kneighbors(pts, return_distance=False)
        nodes[ind] = np.append(nodes[ind], nodes_res)
        strengths[ind] = np.append(strengths[ind], strength_values)


filenames = ["Post_Strength_Test_PA1", "Post_Strength_Test_LS1", "Ant_Strength_Test_PA1", "Ant_Strength_Test_LS1"]
for i in range(len(filenames)):
    utils.write_vtx_strength_intra(nodes[i], strengths[i], base_dir + filenames[i])
