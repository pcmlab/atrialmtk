#!/usr/bin/python3
#
import os
import sys
import argparse
from sklearn.neighbors import NearestNeighbors
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
my_parser.add_argument("scale_factor", metavar="scale factor", type=int, help="Scale factor for landmarks")

# Execute parse_args()
args = my_parser.parse_args()

base_dir = args.target_path
laplace_dir = args.laplace_path
fibre_dir = args.fibre_path
mesh_name = args.mesh_name

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
pts = utils.read_pts(base_dir + mesh_name)

xlist = []
ylist = []
zlist = []

for loop in pts:
    xlist.append(loop[0])
    ylist.append(loop[1])
    zlist.append(loop[2])


# SPLITTING
# catching up to this point (read in meshes [xlist, ylist, zlist] so the next code works)

# import igb files
PA_UD_Post = utils.read_array_igb(base_dir + "UD_Post_UAC/phie.igb")[0]
LS_LR_Post = utils.read_array_igb(base_dir + "LR_Post_UAC/phie.igb")[0]
PA_UD_Ant = utils.read_array_igb(base_dir + "UD_Ant_UAC/phie.igb")[0]
LS_LR_Ant = utils.read_array_igb(base_dir + "LR_Ant_UAC/phie.igb")[0]


Pts_Labelled = [xlist, ylist, zlist]
Pts_Labelled = list(zip(*Pts_Labelled))


spts_post = utils.read_pts(base_dir + "PosteriorMesh")
Pxlist = []
Pylist = []
Pzlist = []

for loop in spts_post:
    Pxlist.append(loop[0])
    Pylist.append(loop[1])
    Pzlist.append(loop[2])

Axlist = []
Aylist = []
Azlist = []

spts_ant = utils.read_pts(base_dir + "AnteriorMesh")
for loop in spts_ant:
    Axlist.append(loop[0])
    Aylist.append(loop[1])
    Azlist.append(loop[2])

Posterior_Pts = [Pxlist, Pylist, Pzlist]
Anterior_Pts = [Axlist, Aylist, Azlist]
Posterior_Pts = list(zip(*Posterior_Pts))
Anterior_Pts = list(zip(*Anterior_Pts))

# update to KD tree as much faster
neigh = NearestNeighbors(n_neighbors=1)
neigh.fit(Posterior_Pts)
Closest_Posterior = neigh.kneighbors(Pts_Labelled, return_distance=False)
neigh.fit(Anterior_Pts)
Closest_Anterior = neigh.kneighbors(Pts_Labelled, return_distance=False)


Pts2D = []
ylist_2D = utils.read_pts(base_dir + "Labelled_Coords_2D_Rescaling_N3")[:, 1]

counter1 = 0
counter2 = 0
for inds in range(0, len(Pts_Labelled), 1):
    if ylist_2D[inds] < 0.5:
        Index = Closest_Posterior[inds]
        Index = Index[0]
        Pts2D.append([1 - LS_LR_Post[Index], 0.5 * PA_UD_Post[Index], 0])
        counter1 = counter1 + 1
    else:
        Index = Closest_Anterior[inds]
        Index = Index[0]
        Pts2D.append([1 - LS_LR_Ant[Index], 1 - 0.5 * PA_UD_Ant[Index], 0])
        counter2 = counter2 + 1

    if inds % 500 == 0:
        print(inds)


print("Saving output Labelled_Coords_2D_Rescaling_v3_C")
_, elems, fiber, data = utils.read_carp(base_dir, mesh_name, return_surface=False)
utils.write_vtk(Pts2D, elems, fiber, data, base_dir + "Labelled_Coords_2D_Rescaling_v3_C.vtk")

# write 2D carp too
surface = utils.get_vtk_from_file(base_dir + "Labelled_Coords_2D_Rescaling_v3_C.vtk")
surface = utils.convert_unstructureddata_to_polydata(surface)
spts, selems = utils.poly2carp_with_labels(surface)
mname = base_dir + "Labelled_Coords_2D_Rescaling_v3_C"
utils.write_carp(spts, selems, None, None, mname)
