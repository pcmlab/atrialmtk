#!/usr/bin/python3
#
import os
import sys
import argparse
import numpy as np
from sklearn.neighbors import NearestNeighbors
from uac import utils
import uac.vars as var


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
laa_label = args.laa_label
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

sys.path.append(base_dir)


# read carp, write to vtk file
pts, elems, fiber, data, surface = utils.read_carp(base_dir, mesh_name)


# find nodes on LA/PV junctions
nodes_lspv = utils.find_junc_nodes(elems, la_label, lspv_label)
nodes_rspv = utils.find_junc_nodes(elems, la_label, rspv_label)
nodes_lipv = utils.find_junc_nodes(elems, la_label, lipv_label)
nodes_ripv = utils.find_junc_nodes(elems, la_label, ripv_label)
nodes_laa = utils.find_junc_nodes(elems, la_label, laa_label)


# find MV nodes
*_, nodes_mv = utils.find_nodes(elems, pts)

# Read in surfaces
tot_pts, surface_filter = utils.read_vtk_surface(surface)

pa_ud = utils.read_array_igb(base_dir + "PA_UAC_N2/phie.igb")[0]
ls_lr = utils.read_array_igb(base_dir + "LR_UAC_N2/phie.igb")[0]

# 1. isolate PV mesh
# Edge list for LA/PV junction

surface_filter_sub_lipv = utils.get_sub_surface(base_dir, pts, elems, nodes_lipv)
surface_filter_sub_ripv = utils.get_sub_surface(base_dir, pts, elems, nodes_ripv)


# map inferior veins to circles
(
    nodes_lipv_f,
    strengths_pa_li,
    strengths_ls_li,
    max_ud_li_ind,
    min_ud_li_ind,
    max_lr_li_ind,
    min_lr_li_ind,
) = utils.vein2circle(
    nodes_lipv,
    pa_ud,
    ls_lr,
    tot_pts,
    surface_filter_sub_lipv,
    var.XSHIFT_23,
    var.YSHIFT_23_27,
    var.R1_23_27,
    var.R2_23_27,
)
(
    nodes_ripv_f,
    strengths_pa_ri,
    strengths_ls_ri,
    max_ud_ri_ind,
    min_ud_ri_ind,
    max_lr_ri_ind,
    min_lr_ri_ind,
) = utils.vein2circle(
    nodes_ripv,
    pa_ud,
    ls_lr,
    tot_pts,
    surface_filter_sub_ripv,
    var.XSHIFT_27,
    var.YSHIFT_23_27,
    var.R1_23_27,
    var.R2_23_27,
)


# for each node on the boundary, find closest node point and assign.
# first - boundary nodes: nodes_lipv, nodes_ripv

pts_lipv = pts[nodes_lipv]
pts_ripv = pts[nodes_ripv]
pts_lipv_f = pts[nodes_lipv_f]
pts_ripv_f = pts[nodes_ripv_f]

post_strength_pa = np.empty(len(pts_lipv) + len(pts_ripv), float)
post_strength_ls = np.empty_like(post_strength_pa)

for i in range(len(pts_lipv)):
    index = utils.closest_node(pts_lipv[i], pts_lipv_f).argmin()
    post_strength_pa[i] = strengths_pa_li[index]
    post_strength_ls[i] = strengths_ls_li[index]

ind = len(pts_lipv)

for i in range(len(pts_ripv)):
    index = utils.closest_node(pts_ripv[i], pts_ripv_f).argmin()
    post_strength_pa[ind + i] = strengths_pa_ri[index]
    post_strength_ls[ind + i] = strengths_ls_ri[index]

post_nodes_pa = np.append(nodes_lipv, nodes_ripv)
post_nodes_ls = np.append(nodes_lipv, nodes_ripv)


# use region labels to relabel each one...
pts_pv = [None] * 4
nodes_pv = [None] * 4
for i in range(len(pts_pv)):
    pts_pv[i], nodes_pv[i] = utils.pv_boundary(i, elems, pts)


mid_pts = np.empty((len(elems), 3), float)
for i in range(len(elems)):
    mid_pts[i] = elems[i].centre(pts)

print(mid_pts[0])


# for each, find closest element mid-point
neigh = NearestNeighbors(n_neighbors=1)
neigh.fit(mid_pts)

sp = [None] * len(pts_pv)
for i in range(len(pts_pv)):
    closest_pv = neigh.kneighbors(pts_pv[i], return_distance=False)
    ele_tag_pv = np.empty(len(closest_pv), int)
    for j in range(len(closest_pv)):
        ele_tag_pv[j] = elems[closest_pv[j][0]].reg
    sp[i] = ele_tag_pv[0]
    print(ele_tag_pv)


# use to order regions
pv_choices = [sp.index(lspv_label), sp.index(lipv_label), sp.index(rspv_label), sp.index(ripv_label)]
print(pv_choices)

ind_lspv = nodes_pv[pv_choices[0]]
ind_lipv = nodes_pv[pv_choices[1]]
ind_rspv = nodes_pv[pv_choices[2]]
ind_ripv = nodes_pv[pv_choices[3]]

surface_filter_sub_lipv = utils.get_sub_surface(base_dir, pts, elems, ind_lipv)
surface_filter_sub_ripv = utils.get_sub_surface(base_dir, pts, elems, ind_ripv)

nodes_lipv_f, strengths_pa_li, strengths_ls_li, *_ = utils.vein2circle(
    ind_lipv,
    pa_ud,
    ls_lr,
    tot_pts,
    surface_filter_sub_lipv,
    var.XSHIFT_23,
    var.YSHIFT_23_27,
    var.R1_23_27_INNER,
    var.R2_23_27_INNER,
)
nodes_ripv_f, strengths_pa_ri, strengths_ls_ri, *_ = utils.vein2circle(
    ind_ripv,
    pa_ud,
    ls_lr,
    tot_pts,
    surface_filter_sub_ripv,
    var.XSHIFT_27,
    var.YSHIFT_23_27,
    var.R1_23_27_INNER,
    var.R2_23_27_INNER,
)

post_nodes_pa = np.concatenate([post_nodes_pa, nodes_lipv_f, nodes_ripv_f])
post_nodes_ls = np.concatenate([post_nodes_ls, nodes_lipv_f, nodes_ripv_f])

post_strength_pa = np.concatenate([post_strength_pa, strengths_pa_li, strengths_pa_ri])
post_strength_ls = np.concatenate([post_strength_ls, strengths_ls_li, strengths_ls_ri])


surface_filter_sub_lspv = utils.get_sub_surface(base_dir, pts, elems, ind_lspv)
surface_filter_sub_rspv = utils.get_sub_surface(base_dir, pts, elems, ind_rspv)

nodes_lipv_f, strengths_pa_li, strengths_ls_li, *_ = utils.vein2circle(
    ind_lspv,
    pa_ud,
    ls_lr,
    tot_pts,
    surface_filter_sub_lspv,
    var.XSHIFT_21,
    var.YSHIFT_21_25,
    var.R1_21_25_INNER,
    var.R2_21_25_INNER,
)
nodes_ripv_f, strengths_pa_ri, strengths_ls_ri, *_ = utils.vein2circle(
    ind_rspv,
    pa_ud,
    ls_lr,
    tot_pts,
    surface_filter_sub_rspv,
    var.XSHIFT_25,
    var.YSHIFT_21_25,
    var.R1_21_25_INNER,
    var.R2_21_25_INNER,
)

ant_nodes_pa = np.append(nodes_lipv_f, nodes_ripv_f)
ant_nodes_ls = np.append(nodes_lipv_f, nodes_ripv_f)

ant_strength_pa = np.append(strengths_pa_li, strengths_pa_ri)
ant_strength_ls = np.append(strengths_ls_li, strengths_ls_ri)

surface_filter_sub_laa = utils.get_sub_surface(base_dir, pts, elems, nodes_laa)

(
    nodes_laa_f,
    strengths_pa_laa,
    strengths_ls_laa,
    max_ud_laa_ind,
    min_ud_laa_ind,
    max_lr_laa_ind,
    min_lr_laa_ind,
) = utils.vein2circle(
    nodes_laa,
    pa_ud,
    ls_lr,
    tot_pts,
    surface_filter_sub_laa,
    var.XSHIFT_LAA,
    var.YSHIFT_LAA,
    var.R1_LAA,
    var.R2_LAA,
)

# for each node on the boundary, find closest node point and assign.
# first - boundary nodes: Nodes_LA_LIPV, Nodes_LA_RIPV

pts_laa = pts[nodes_laa]
pts_laa_f = pts[nodes_laa_f]

pai = len(ant_strength_pa)
lsi = len(ant_strength_ls)
laa_len = len(pts_laa)
ant_strength_pa.resize(pai + laa_len)
ant_strength_ls.resize(lsi + laa_len)

for i in range(laa_len):
    index = utils.closest_node(pts_laa[i], pts_laa_f).argmin()
    ant_strength_pa[pai + i] = strengths_pa_laa[index]
    ant_strength_ls[lsi + i] = strengths_ls_laa[index]

ant_nodes_pa = np.append(ant_nodes_pa, nodes_laa)
ant_nodes_ls = np.append(ant_nodes_ls, nodes_laa)


# rescale y coord of mesh
nodes_border = utils.read_pts(base_dir + "BorderNodes", vtx=True)

elems_left = utils.element2delete(elems, nodes_border, excluded=True)
utils.write_vtk(pts, elems_left, None, None, base_dir + "Test_Split.vtk")
surface = utils.get_vtk_from_file(base_dir + "Test_Split.vtk")

ripv_marker = pts_ripv[0]
output = utils.extract_region_by_marker(surface, ripv_marker)
pts_split, _ = utils.read_vtk_surface(output)
surface = utils.convert_unstructureddata_to_polydata(output)
spts_post, selems_post = utils.poly2carp_with_labels(surface)
utils.write_carp(spts_post, selems_post, None, None, base_dir + "PosteriorMesh")


# for each point in posterior, find in labelled. Take this list and remove to give anterior mesh
# KD tree as much faster
neigh = NearestNeighbors(n_neighbors=1)
neigh.fit(pts)
nodes = neigh.kneighbors(pts_split, return_distance=False)

elems_left = utils.element2delete(elems, nodes, excluded=True)
utils.write_vtk(pts, elems_left, None, None, base_dir + "Test_ant.vtk")

# write anterior mesh
surface = utils.get_vtk_from_file(base_dir + "Test_ant.vtk")
surface = utils.convert_unstructureddata_to_polydata(surface)
surface = utils.clean_polydata(surface)
spts_ant, selems_ant = utils.poly2carp_with_labels(surface)
utils.write_carp(spts_ant, selems_ant, None, None, base_dir + "AnteriorMesh")


# now calculate UAC

pa_ud_r = 1 - pa_ud * 0.5

for i in nodes:
    pa_ud_r[i] = pa_ud[i] * 0.5

pts2d_r = np.empty((len(pts), 3), float)
for i in range(len(pts)):
    pts2d_r[i] = [1 - ls_lr[i], pa_ud_r[i], 0]

mname = base_dir + "Labelled_Coords_2D_Rescaling_N3.vtk"
utils.write_vtk(pts2d_r, elems, fiber, data, mname)

# write 2D carp too
surface = utils.get_vtk_from_file(mname)
surface = utils.convert_unstructureddata_to_polydata(surface)
spts, selems = utils.poly2carp_with_labels(surface)
mname = mname[:-4]
utils.write_carp(spts, selems, None, None, mname)

ylist2d = utils.read_pts(mname)[:, 1]

# define above and below inferior veins - these only go outwards (with geodesic in between)
# isolines for above and below inferior veins

nl_pa_06_l, diff_pa_06_l = utils.get_iso_points_lt_post_ll(
    pa_ud, pa_ud[min_ud_li_ind], ls_lr, ls_lr[min_ud_li_ind], 0.01, ylist2d
)

nl_pa_06_r, diff_pa_06_r = utils.get_iso_points_lt_post_gl(
    pa_ud, pa_ud[min_ud_ri_ind], ls_lr, ls_lr[min_ud_ri_ind], 0.01, ylist2d
)

nl_pa_09_l, diff_pa_09_l = utils.get_iso_points_gt_post_ll(
    pa_ud, pa_ud[max_ud_li_ind], ls_lr, ls_lr[max_ud_li_ind], 0.01, ylist2d
)

nl_pa_09_r, diff_pa_09_r = utils.get_iso_points_gt_post_gl(
    pa_ud, pa_ud[max_ud_ri_ind], ls_lr, ls_lr[max_ud_ri_ind], 0.01, ylist2d
)

nl_09, diff_09 = utils.get_iso_points_gt_post_ll(
    ls_lr, ls_lr[max_lr_ri_ind], pa_ud, pa_ud[max_lr_ri_ind], 0.01, ylist2d
)

nl_075, diff_075 = utils.get_iso_points_lt_post_ll(
    ls_lr, ls_lr[min_lr_ri_ind], pa_ud, pa_ud[min_lr_ri_ind], 0.01, ylist2d
)

nl_01, diff_01 = utils.get_iso_points_lt_post_ll(
    ls_lr, ls_lr[min_lr_li_ind], pa_ud, pa_ud[min_lr_li_ind], 0.01, ylist2d
)

nl_025, diff_025 = utils.get_iso_points_gt_post_ll(
    ls_lr, ls_lr[max_lr_li_ind], pa_ud, pa_ud[max_lr_li_ind], 0.01, ylist2d
)

post_nodes_ls = np.concatenate([post_nodes_ls, nl_01, nl_025, nl_075, nl_09])
post_strength_ls = np.concatenate(
    [
        post_strength_ls,
        var.LS_INF_PV_1 + diff_01,
        var.LS_INF_PV_2 + diff_025,
        var.LS_INF_PV_3 + diff_075,
        var.LS_INF_PV_4 + diff_09,
    ]
)


# calculate paths between top and bottom inferior veins

top_li = tot_pts[max_ud_li_ind]
bot_li = tot_pts[min_ud_li_ind]
top_ri = tot_pts[max_ud_ri_ind]
bot_ri = tot_pts[min_ud_ri_ind]

_, path_top_inds = utils.geod_lspv_rspv(surface_filter, top_li, top_ri, tot_pts)
_, path_bot_inds = utils.geod_lspv_rspv(surface_filter, bot_li, bot_ri, tot_pts)

path_top_inds = np.flip(path_top_inds, 0)

post_nodes_pa = np.concatenate(
    [post_nodes_pa, path_top_inds, path_bot_inds, nl_pa_06_l, nl_pa_06_r, nl_pa_09_l, nl_pa_09_r]
)
post_strength_pa = np.concatenate(
    [
        post_strength_pa,
        np.ones(len(path_top_inds)) * var.PA_INF_PV_2,
        np.ones(len(path_bot_inds)) * var.PA_INF_PV_1,
        var.PA_INF_PV_1 + diff_pa_06_l,
        var.PA_INF_PV_1 + diff_pa_06_r,
        var.PA_INF_PV_2 + diff_pa_09_l,
        var.PA_INF_PV_2 + diff_pa_09_r,
    ]
)


# now for LAA - add top and bottom paths (isolines for top and bottom)
# isolines and geodesics for left and right

nodes_la = utils.find_region_nodes(elems, la_label)

nl, diff = utils.get_iso_points_gt_ant(pa_ud, pa_ud[max_ud_laa_ind], 0.01, ylist2d)
indices = np.arange(nl.shape[0])[np.in1d(nl, nodes_la)]
nl_s = nl[indices]
diff_s = diff[indices]
strength_pts = var.PA_ANT_LAA + diff_s

ant_nodes_pa = np.append(ant_nodes_pa, nl_s)
ant_strength_pa = np.append(ant_strength_pa, strength_pts)

nl, diff = utils.get_iso_points_lt_ant(pa_ud, pa_ud[min_ud_laa_ind], 0.01, ylist2d)
indices = np.arange(nl.shape[0])[np.in1d(nl, nodes_la)]
nl_s = nl[indices]
diff_s = diff[indices]
strength_pts = var.PA_ANT_LAA_BOT + diff_s

ant_nodes_pa = np.append(ant_nodes_pa, nl_s)
ant_strength_pa = np.append(ant_strength_pa, strength_pts)

nl, diff = utils.get_iso_points_gt_ant_ll(ls_lr, ls_lr[max_lr_laa_ind], pa_ud, pa_ud[max_lr_laa_ind], 0.01, ylist2d)
indices = np.arange(nl.shape[0])[np.in1d(nl, nodes_la)]
nl_s = nl[indices]
diff_s = diff[indices]
strength_pts = var.LS_ANT_LAA_2 + diff_s

ant_nodes_ls = np.append(ant_nodes_ls, nl_s)
ant_strength_ls = np.append(ant_strength_ls, strength_pts)

nl, diff = utils.get_iso_points_lt_ant_ll(ls_lr, ls_lr[min_lr_laa_ind], pa_ud, pa_ud[min_lr_laa_ind], 0.01, ylist2d)
indices = np.arange(nl.shape[0])[np.in1d(nl, nodes_la)]
nl_s = nl[indices]
diff_s = diff[indices]
strength_pts = var.LS_ANT_LAA_1 + diff_s

ant_nodes_ls = np.append(ant_nodes_ls, nl_s)
ant_strength_ls = np.append(ant_strength_ls, strength_pts)

# 1. Find closest points of superior veins to inferior veins. Do this using boundary_markers function.
_, lipv_lspv = utils.boundary_markers(nodes_lipv, nodes_lspv, tot_pts)
_, ripv_rspv = utils.boundary_markers(nodes_ripv, nodes_rspv, tot_pts)

pts_lspv = pts[nodes_lspv]
p1 = pts[nodes_lspv[lipv_lspv]]
index = utils.find_index(((pts_lspv - p1) ** 2).sum(1), max)
p2 = pts[nodes_lspv[index]]

# 2. LSPV find opposite point in rim
surface_filter_sub_lspv = utils.get_sub_surface(base_dir, pts, elems, nodes_lspv)
path_pts, _ = utils.geod_lspv_rspv(surface_filter_sub_lspv, p1, p2, tot_pts)
mid_point = path_pts[len(path_pts) // 2]

increment = 1

index = utils.geod_midpoint(surface_filter_sub_lspv, mid_point, nodes_lspv, increment, tot_pts)
mp = pts[index]
tpts1, _ = utils.geod_lspv_rspv(surface_filter_sub_lspv, p1, mp, tot_pts)
tpts2, _ = utils.geod_lspv_rspv(surface_filter_sub_lspv, mp, p2, tot_pts)
at = np.concatenate([tpts1, tpts2])

pts_src = utils.load_landmarks(base_dir + seed_name, scale_factor)

lspv_marker_pts = pts_src[2]
m1 = ((path_pts - lspv_marker_pts) ** 2).sum(1).min()
m2 = ((at - lspv_marker_pts) ** 2).sum(1).min()

if m1 < m2:
    semi_circle_mp_lspv = mp
    semi_circle_ap_lspv = mid_point
else:
    semi_circle_mp_lspv = mid_point
    semi_circle_ap_lspv = mp

start_seed_lspv = nodes_lspv[lipv_lspv]
semi_circle_end_lspv = p2


# also update RSPV - take path closest to LIPV
# 3. RSPV find subdividers in rims
# find furthest point as it's conjugate
p3 = pts[nodes_rspv[ripv_rspv]]
pts_rspv = pts[nodes_rspv]
index = utils.find_index(((pts_rspv - p3) ** 2).sum(1), max)
p4 = pts[nodes_rspv[index]]

# now find two subdividers in rim and join these...
surface_filter_sub_rspv = utils.get_sub_surface(base_dir, pts, elems, nodes_rspv)
path_pts, _ = utils.geod_lspv_rspv(surface_filter_sub_rspv, p3, p4, tot_pts)
mid_point = path_pts[len(path_pts) // 2]

index = utils.geod_midpoint(surface_filter_sub_rspv, mid_point, nodes_rspv, increment, tot_pts)
mp = tot_pts[index]
tpts1, _ = utils.geod_lspv_rspv(surface_filter_sub_rspv, p3, mp, tot_pts)
tpts2, _ = utils.geod_lspv_rspv(surface_filter_sub_rspv, mp, p4, tot_pts)
at = np.concatenate([tpts1, tpts2])

rspv_marker_pts = pts_src[3]
m3 = ((path_pts - rspv_marker_pts) ** 2).sum(1).min()
m4 = ((at - rspv_marker_pts) ** 2).sum(1).min()

if m3 < m4:
    semi_circle_mp_rpsv = mp
    semi_circle_ap_rspv = mid_point
else:
    semi_circle_mp_rpsv = mid_point
    semi_circle_ap_rspv = mp

start_seed_rspv = nodes_rspv[ripv_rspv]
semi_circle_end_rspv = p4

nodes, strengths = utils.lspv2line(
    surface_filter_sub_lspv,
    tot_pts,
    semi_circle_ap_lspv,
    semi_circle_end_lspv,
    var.XSHIFT_21,
    start_seed_lspv,
    var.R1_21_25,
)

ant_nodes_ls = np.append(ant_nodes_ls, nodes)
ant_strength_ls = np.append(ant_strength_ls, strengths)

nodes, strengths_pa, strengths_ls, _ = utils.lspv2circle(
    surface_filter_sub_lspv,
    tot_pts,
    semi_circle_mp_lspv,
    semi_circle_end_lspv,
    var.XSHIFT_21,
    start_seed_lspv,
    var.R1_21_25,
    var.R2_21_25,
)

utils.write_vtx_strength_extra(nodes, strengths_ls, base_dir + "LSPV_LS")
utils.write_vtx_strength_extra(nodes, strengths_pa, base_dir + "LSPV_PA")


# for each node on the boundary, find closest node point and assign.
# first - boundary nodes: Nodes_LA_LIPV, Nodes_LA_RIPV

pts_ud2 = set(utils.read_pts(base_dir + "PAbc2", n=1, vtx=True, item_type=int))
unassigneds = [i for i in nodes_lspv if i not in pts_ud2]
pts_lspv_r = pts[unassigneds]
pts_lspv_f = pts[nodes]

pts_lspv_r_len = len(pts_lspv_r)
lsi = len(ant_strength_ls)
pai = len(ant_strength_pa)
ant_strength_ls.resize(lsi + pts_lspv_r_len)
ant_strength_pa.resize(pai + pts_lspv_r_len)

for i in range(pts_lspv_r_len):
    index = utils.closest_node(pts_lspv_r[i], pts_lspv_f).argmin()
    ant_strength_ls[lsi + i] = strengths_ls[index]
    ant_strength_pa[pai + i] = strengths_pa[index]

ant_nodes_pa = np.append(ant_nodes_pa, unassigneds)
ant_nodes_ls = np.append(ant_nodes_ls, unassigneds)

nodes, strengths = utils.rspv2line(
    surface_filter_sub_rspv,
    var.XSHIFT_25,
    tot_pts,
    start_seed_rspv,
    var.R1_21_25,
    semi_circle_ap_rspv,
    semi_circle_end_rspv,
)

ant_nodes_ls = np.append(ant_nodes_ls, nodes)
ant_strength_ls = np.append(ant_strength_ls, strengths)

nodes_ls, nodes_pa, strengths_pa, strengths_ls = utils.rspv2circle(
    surface_filter_sub_rspv,
    var.XSHIFT_25,
    tot_pts,
    start_seed_rspv,
    var.R1_21_25,
    var.R2_21_25,
    semi_circle_end_rspv,
    semi_circle_mp_rpsv,
)

utils.write_vtx_strength_extra(nodes_ls, strengths_ls, base_dir + "RSPV_LS")
utils.write_vtx_strength_extra(nodes_pa, strengths_pa, base_dir + "RSPV_PA")

# for each node on the boundary, find closest node point and assign.
# first - boundary nodes: Nodes_LA_LIPV, Nodes_LA_RIPV

unassigneds = [i for i in nodes_rspv if i not in pts_ud2]
pts_rspv_r = pts[unassigneds]
pts_rspv_f = pts[nodes_ls]
pts_rspv_pf = pts[nodes_pa]

pts_rspv_r_len = len(pts_rspv_r)
lsi = len(ant_strength_ls)
pai = len(ant_strength_pa)
ant_strength_ls.resize(lsi + pts_rspv_r_len)
ant_strength_pa.resize(pai + pts_rspv_r_len)

for i in range(pts_rspv_r_len):
    index = utils.closest_node(pts_rspv_r[i], pts_rspv_f).argmin()
    ant_strength_ls[lsi + i] = strengths_ls[index]
    index = utils.closest_node(pts_rspv_r[i], pts_rspv_pf).argmin()
    ant_strength_pa[pai + i] = strengths_pa[index]

ant_nodes_pa = np.append(ant_nodes_pa, unassigneds)
ant_nodes_ls = np.append(ant_nodes_ls, unassigneds)


_, path_top_inds = utils.geod_lspv_rspv(surface_filter, semi_circle_mp_lspv, semi_circle_mp_rpsv, tot_pts)

ant_nodes_pa = np.append(ant_nodes_pa, path_top_inds)
ant_strength_pa = np.append(ant_strength_pa, np.ones(len(path_top_inds)) * (1 - var.R2_21_25))

# geodesic LAA R. Semi_Circle_End_LSPV. IndFound_max_LR_LAA
_, path_top_inds = utils.geod_lspv_rspv(surface_filter, semi_circle_end_lspv, tot_pts[max_lr_laa_ind], tot_pts)

ant_nodes_ls = np.append(ant_nodes_ls, path_top_inds)
ant_strength_ls = np.append(ant_strength_ls, np.ones(len(path_top_inds)) * var.LS_ANT_LAA_2)

# change to geodesic
indices = utils.geod_marker_mv(surface_filter, semi_circle_end_rspv, increment, tot_pts, nodes_mv)

ant_nodes_ls = np.append(ant_nodes_ls, indices)
ant_strength_ls = np.append(ant_strength_ls, np.ones(len(indices)) * var.ANT_RSPV_LS)

# check strengths so far:
utils.write_vtx_strength_extra(post_nodes_pa, post_strength_pa, base_dir + "P_Checker_PA")
utils.write_vtx_strength_extra(post_nodes_ls, post_strength_ls, base_dir + "P_Checker_LS")
utils.write_vtx_strength_extra(ant_nodes_pa, ant_strength_pa, base_dir + "A_Checker_PA")
utils.write_vtx_strength_extra(ant_nodes_ls, ant_strength_ls, base_dir + "A_Checker_LS")

# code to separate nodes for anterior and posterior solves
# separate node lists in anteior and posterior for all except old UAC boundaries - must be in both

nodes = [None] * 4
strengths = [None] * len(nodes)
nodes_src = [post_nodes_pa, post_nodes_ls, ant_nodes_pa, ant_nodes_ls]
strengths_src = [post_strength_pa, post_strength_ls, ant_strength_pa, ant_strength_ls]
spts_src = [spts_post, spts_ant]

for i in range(len(nodes)):
    pts = tot_pts[nodes_src[i]]
    strengths[i] = strengths_src[i][: len(nodes_src[i])]
    neigh = NearestNeighbors(n_neighbors=1)
    neigh.fit(spts_src[i // 2])
    nodes[i] = neigh.kneighbors(pts, return_distance=False)
    print(i)
    print(len(nodes[i]))
    print(len(strengths[i]))


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
    utils.write_vtx_strength_extra(nodes[i], strengths[i], base_dir + filenames[i])
