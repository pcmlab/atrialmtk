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
my_parser.add_argument(
    "lat_file_name", metavar="filename", type=str, help="LAT field file name (e.g. LAT_Spiral4_B.dat)"
)
my_parser.add_argument("mesh_name", metavar="filename", type=str, help="Mesh Name(e.g. Fibres)")

# Execute parse_args()
args = my_parser.parse_args()

base_dir = args.target_path
lat_file_name = args.lat_file_name

if not os.path.isdir(base_dir):
    print("The target mesh path specified does not exist")
    sys.exit()


# read carp, write to vtk file
print(base_dir)
Pts_2D, *_ = utils.read_carp(base_dir, "Labelled_Coords_2D_Rescaling_v3_C", return_surface=False)

print(Pts_2D)

xcen = 0.5
ycen = 0.5

xcen2 = -0.5
ycen2 = 0.5

Pts_2D = Pts_2D * 2
Pts_2D = Pts_2D - 1

X = []
Y = []

for loop in Pts_2D:
    X.append(loop[0])
    Y.append(loop[1])

X = np.array(X)
Y = np.array(Y)
Xshift = X - xcen
Yshift = Y - ycen
Xshift2 = X - xcen2
Yshift2 = Y - ycen2


kappa = 1
omega = -1
D_spiral = 1
t = (
    -2 * np.pi / (D_spiral * omega) * np.sqrt(Xshift**2 + Yshift**2)
    + 1 / omega * np.arctan2(Yshift, Xshift)
    + 2 * np.pi * kappa / omega
)

omega = -1
kappa = 1
D_spiral = -1
time2 = (
    -2 * np.pi / (D_spiral * omega) * np.sqrt(Xshift2**2 + Yshift2**2)
    + 1 / omega * np.arctan2(Yshift2, Xshift2)
    + 2 * np.pi * kappa / omega
)

time2 = -time2 - 14


time = t
TimeAll = time
TimeAll_new = np.where(time > time2, time2, TimeAll)
TimeAll = TimeAll_new

xcen = -0.5
ycen = -0.5

xcen2 = 0.5
ycen2 = -0.5
Xshift = X - xcen
Yshift = Y - ycen
Xshift2 = X - xcen2
Yshift2 = Y - ycen2

kappa = 1
omega = -1
D_spiral = 1
t = (
    -2 * np.pi / (D_spiral * omega) * np.sqrt(Xshift**2 + Yshift**2)
    + 1 / omega * np.arctan2(Yshift, Xshift)
    + 2 * np.pi * kappa / omega
)


TimeAll_new = np.where(TimeAll > t, t, TimeAll)
TimeAll = TimeAll_new

omega = -1
kappa = 1
D_spiral = -1

time2 = (
    -2 * np.pi / (D_spiral * omega) * np.sqrt(Xshift2**2 + Yshift2**2)
    + 1 / omega * np.arctan2(Yshift2, Xshift2)
    + 2 * np.pi * kappa / omega
)
time2 = -time2 - 14


TimeAll_new = np.where(TimeAll > time2, time2, TimeAll)
TimeAll = TimeAll_new


TT = TimeAll - np.min(TimeAll)
TT = TT * 34
print(TT)

utils.write_dat_simple(TT, base_dir + lat_file_name)

TT = np.concatenate([TT, TT])

# bilayer version
utils.write_dat_simple(TT, base_dir + "Bi" + lat_file_name)
