#!/usr/bin/python3
#
import os
import sys
import argparse
import shutil
import numpy as np
from sklearn.neighbors import NearestNeighbors
from uac import utils


# Create the parser
my_parser = argparse.ArgumentParser(description="Calculation of new atrial coordinates")

# Add the arguments
my_parser.add_argument("target_path", metavar="path", type=str, help="Path for the new target mesh")
my_parser.add_argument("folder1", metavar="path", type=str, help="LA_1")
my_parser.add_argument("folder2", metavar="path", type=str, help="LA_2")
my_parser.add_argument("folder3", metavar="path", type=str, help="RA_1")
my_parser.add_argument("folder4", metavar="path", type=str, help="RA_2")

# Execute parse_args()
args = my_parser.parse_args()

base_dir = args.target_path
folder1 = args.folder1
folder2 = args.folder2
folder3 = args.folder3
folder4 = args.folder4

if not os.path.isdir(base_dir):
    print("The target mesh path specified does not exist")
    sys.exit()

sys.path.append(base_dir)


print(base_dir)
pts_LA_1 = utils.read_pts(base_dir + folder1 + "/LA_1")
pts_LA_2 = utils.read_pts(base_dir + folder2 + "/LA_2")
pts_RA_1 = utils.read_pts(base_dir + folder3 + "/RA_1")
pts_RA_2 = utils.read_pts(base_dir + folder4 + "/RA_2")


os.makedirs(base_dir + "/LA_endo",exist_ok=True)
os.makedirs(base_dir + "/LA_epi",exist_ok=True)
os.makedirs(base_dir + "/RA_endo", exist_ok=True)
os.makedirs(base_dir + "/RA_epi", exist_ok=True)


# update to KD tree as much faster
neigh = NearestNeighbors(n_neighbors=1)
neigh.fit(pts_LA_1)
LA_1_1, _ = neigh.kneighbors(pts_RA_1, return_distance=True)
LA_1_2, _ = neigh.kneighbors(pts_RA_2, return_distance=True)

neigh = NearestNeighbors(n_neighbors=1)
neigh.fit(pts_LA_2)
LA_2_1, _ = neigh.kneighbors(pts_RA_1, return_distance=True)
LA_2_2, _ = neigh.kneighbors(pts_RA_2, return_distance=True)

print((np.min(LA_1_1)))

Min_1 = np.min(LA_1_1)
Min_2 = np.min(LA_1_2)
Min_3 = np.min(LA_2_1)
Min_4 = np.min(LA_2_2)

Min_All = np.append(Min_1, Min_2)
Min_All = np.append(Min_All, Min_3)
Min_All = np.append(Min_All, Min_4)
print(Min_All)

index_min = np.argmin(Min_All)
print(index_min)

# save depending on min...
if index_min == 0:
    # LA epi
    src = base_dir + folder1 + "/LA_1.pts"
    dst = base_dir + "LA_epi" + "/LA_epi.pts"
    shutil.copyfile(src, dst)
    src = base_dir + folder1 + "/LA_1.lon"
    dst = base_dir + "LA_epi" + "/LA_epi.lon"
    shutil.copyfile(src, dst)
    src = base_dir + folder1 + "/LA_1.elem"
    dst = base_dir + "LA_epi" + "/LA_epi.elem"
    shutil.copyfile(src, dst)
    src = base_dir + folder1 + "/LA_1.vtk"
    dst = base_dir + "LA_epi" + "/LA_epi.vtk"
    shutil.copyfile(src, dst)
    # RA epi
    src = base_dir + folder3 + "/RA_1.pts"
    dst = base_dir + "RA_epi" + "/RA_epi.pts"
    shutil.copyfile(src, dst)
    src = base_dir + folder3 + "/RA_1.lon"
    dst = base_dir + "RA_epi" + "/RA_epi.lon"
    shutil.copyfile(src, dst)
    src = base_dir + folder3 + "/RA_1.elem"
    dst = base_dir + "RA_epi" + "/RA_epi.elem"
    shutil.copyfile(src, dst)
    src = base_dir + folder3 + "/RA_1.vtk"
    dst = base_dir + "RA_epi" + "/RA_epi.vtk"
    shutil.copyfile(src, dst)
    # LA endo
    src = base_dir + folder2 + "/LA_2.pts"
    dst = base_dir + "LA_endo" + "/LA_endo.pts"
    shutil.copyfile(src, dst)
    src = base_dir + folder2 + "/LA_2.lon"
    dst = base_dir + "LA_endo" + "/LA_endo.lon"
    shutil.copyfile(src, dst)
    src = base_dir + folder2 + "/LA_2.elem"
    dst = base_dir + "LA_endo" + "/LA_endo.elem"
    shutil.copyfile(src, dst)
    src = base_dir + folder2 + "/LA_2.vtk"
    dst = base_dir + "LA_endo" + "/LA_endo.vtk"
    shutil.copyfile(src, dst)
    # RA endo
    src = base_dir + folder4 + "/RA_2.pts"
    dst = base_dir + "RA_endo" + "/RA_endo.pts"
    shutil.copyfile(src, dst)
    src = base_dir + folder4 + "/RA_2.lon"
    dst = base_dir + "RA_endo" + "/RA_endo.lon"
    shutil.copyfile(src, dst)
    src = base_dir + folder4 + "/RA_2.elem"
    dst = base_dir + "RA_endo" + "/RA_endo.elem"
    shutil.copyfile(src, dst)
    src = base_dir + folder4 + "/RA_2.vtk"
    dst = base_dir + "RA_endo" + "/RA_endo.vtk"
    shutil.copyfile(src, dst)


if index_min == 1:
    # LA epi
    src = base_dir + folder1 + "/LA_1.pts"
    dst = base_dir + "LA_epi" + "/LA_epi.pts"
    shutil.copyfile(src, dst)
    src = base_dir + folder1 + "/LA_1.lon"
    dst = base_dir + "LA_epi" + "/LA_epi.lon"
    shutil.copyfile(src, dst)
    src = base_dir + folder1 + "/LA_1.elem"
    dst = base_dir + "LA_epi" + "/LA_epi.elem"
    shutil.copyfile(src, dst)
    src = base_dir + folder1 + "/LA_1.vtk"
    dst = base_dir + "LA_epi" + "/LA_epi.vtk"
    shutil.copyfile(src, dst)
    # RA epi
    src = base_dir + folder4 + "/RA_2.pts"
    dst = base_dir + "RA_epi" + "/RA_epi.pts"
    shutil.copyfile(src, dst)
    src = base_dir + folder4 + "/RA_2.lon"
    dst = base_dir + "RA_epi" + "/RA_epi.lon"
    shutil.copyfile(src, dst)
    src = base_dir + folder4 + "/RA_2.elem"
    dst = base_dir + "RA_epi" + "/RA_epi.elem"
    shutil.copyfile(src, dst)
    src = base_dir + folder4 + "/RA_2.vtk"
    dst = base_dir + "RA_epi" + "/RA_epi.vtk"
    shutil.copyfile(src, dst)
    # LA endo
    src = base_dir + folder2 + "/LA_2.pts"
    dst = base_dir + "LA_endo" + "/LA_endo.pts"
    shutil.copyfile(src, dst)
    src = base_dir + folder2 + "/LA_2.lon"
    dst = base_dir + "LA_endo" + "/LA_endo.lon"
    shutil.copyfile(src, dst)
    src = base_dir + folder2 + "/LA_2.elem"
    dst = base_dir + "LA_endo" + "/LA_endo.elem"
    shutil.copyfile(src, dst)
    src = base_dir + folder2 + "/LA_2.vtk"
    dst = base_dir + "LA_endo" + "/LA_endo.vtk"
    shutil.copyfile(src, dst)
    # RA endo
    src = base_dir + folder3 + "/RA_1.pts"
    dst = base_dir + "RA_endo" + "/RA_endo.pts"
    shutil.copyfile(src, dst)
    src = base_dir + folder3 + "/RA_1.lon"
    dst = base_dir + "RA_endo" + "/RA_endo.lon"
    shutil.copyfile(src, dst)
    src = base_dir + folder3 + "/RA_1.elem"
    dst = base_dir + "RA_endo" + "/RA_endo.elem"
    shutil.copyfile(src, dst)
    src = base_dir + folder3 + "/RA_1.vtk"
    dst = base_dir + "RA_endo" + "/RA_endo.vtk"
    shutil.copyfile(src, dst)


if index_min == 2:
    # LA epi
    src = base_dir + folder2 + "/LA_2.pts"
    dst = base_dir + "LA_epi" + "/LA_epi.pts"
    shutil.copyfile(src, dst)
    src = base_dir + folder2 + "/LA_2.lon"
    dst = base_dir + "LA_epi" + "/LA_epi.lon"
    shutil.copyfile(src, dst)
    src = base_dir + folder2 + "/LA_2.elem"
    dst = base_dir + "LA_epi" + "/LA_epi.elem"
    shutil.copyfile(src, dst)
    src = base_dir + folder2 + "/LA_2.vtk"
    dst = base_dir + "LA_epi" + "/LA_epi.vtk"
    shutil.copyfile(src, dst)
    # RA epi
    src = base_dir + folder3 + "/RA_1.pts"
    dst = base_dir + "RA_epi" + "/RA_epi.pts"
    shutil.copyfile(src, dst)
    src = base_dir + folder3 + "/RA_1.lon"
    dst = base_dir + "RA_epi" + "/RA_epi.lon"
    shutil.copyfile(src, dst)
    src = base_dir + folder3 + "/RA_1.elem"
    dst = base_dir + "RA_epi" + "/RA_epi.elem"
    shutil.copyfile(src, dst)
    src = base_dir + folder3 + "/RA_1.vtk"
    dst = base_dir + "RA_epi" + "/RA_epi.vtk"
    shutil.copyfile(src, dst)
    # LA endo
    src = base_dir + folder1 + "/LA_1.pts"
    dst = base_dir + "LA_endo" + "/LA_endo.pts"
    shutil.copyfile(src, dst)
    src = base_dir + folder1 + "/LA_1.lon"
    dst = base_dir + "LA_endo" + "/LA_endo.lon"
    shutil.copyfile(src, dst)
    src = base_dir + folder1 + "/LA_1.elem"
    dst = base_dir + "LA_endo" + "/LA_endo.elem"
    shutil.copyfile(src, dst)
    src = base_dir + folder1 + "/LA_1.vtk"
    dst = base_dir + "LA_endo" + "/LA_endo.vtk"
    shutil.copyfile(src, dst)
    # RA endo
    src = base_dir + folder4 + "/RA_2.pts"
    dst = base_dir + "RA_endo" + "/RA_endo.pts"
    shutil.copyfile(src, dst)
    src = base_dir + folder4 + "/RA_2.lon"
    dst = base_dir + "RA_endo" + "/RA_endo.lon"
    shutil.copyfile(src, dst)
    src = base_dir + folder4 + "/RA_2.elem"
    dst = base_dir + "RA_endo" + "/RA_endo.elem"
    shutil.copyfile(src, dst)
    src = base_dir + folder4 + "/RA_2.vtk"
    dst = base_dir + "RA_endo" + "/RA_endo.vtk"
    shutil.copyfile(src, dst)


if index_min == 3:
    # LA epi
    src = base_dir + folder2 + "/LA_2.pts"
    dst = base_dir + "LA_epi" + "/LA_epi.pts"
    shutil.copyfile(src, dst)
    src = base_dir + folder2 + "/LA_2.lon"
    dst = base_dir + "LA_epi" + "/LA_epi.lon"
    shutil.copyfile(src, dst)
    src = base_dir + folder2 + "/LA_2.elem"
    dst = base_dir + "LA_epi" + "/LA_epi.elem"
    shutil.copyfile(src, dst)
    src = base_dir + folder2 + "/LA_2.vtk"
    dst = base_dir + "LA_epi" + "/LA_epi.vtk"
    shutil.copyfile(src, dst)
    # RA epi
    src = base_dir + folder4 + "/RA_2.pts"
    dst = base_dir + "RA_epi" + "/RA_epi.pts"
    shutil.copyfile(src, dst)
    src = base_dir + folder4 + "/RA_2.lon"
    dst = base_dir + "RA_epi" + "/RA_epi.lon"
    shutil.copyfile(src, dst)
    src = base_dir + folder4 + "/RA_2.elem"
    dst = base_dir + "RA_epi" + "/RA_epi.elem"
    shutil.copyfile(src, dst)
    src = base_dir + folder4 + "/RA_2.vtk"
    dst = base_dir + "RA_epi" + "/RA_epi.vtk"
    shutil.copyfile(src, dst)
    # LA endo
    src = base_dir + folder1 + "/LA_1.pts"
    dst = base_dir + "LA_endo" + "/LA_endo.pts"
    shutil.copyfile(src, dst)
    src = base_dir + folder1 + "/LA_1.lon"
    dst = base_dir + "LA_endo" + "/LA_endo.lon"
    shutil.copyfile(src, dst)
    src = base_dir + folder1 + "/LA_1.elem"
    dst = base_dir + "LA_endo" + "/LA_endo.elem"
    shutil.copyfile(src, dst)
    src = base_dir + folder1 + "/LA_1.vtk"
    dst = base_dir + "LA_endo" + "/LA_endo.vtk"
    shutil.copyfile(src, dst)
    # RA endo
    src = base_dir + folder3 + "/RA_1.pts"
    dst = base_dir + "RA_endo" + "/RA_endo.pts"
    shutil.copyfile(src, dst)
    src = base_dir + folder3 + "/RA_1.lon"
    dst = base_dir + "RA_endo" + "/RA_endo.lon"
    shutil.copyfile(src, dst)
    src = base_dir + folder3 + "/RA_1.elem"
    dst = base_dir + "RA_endo" + "/RA_endo.elem"
    shutil.copyfile(src, dst)
    src = base_dir + folder3 + "/RA_1.vtk"
    dst = base_dir + "RA_endo" + "/RA_endo.vtk"
    shutil.copyfile(src, dst)
