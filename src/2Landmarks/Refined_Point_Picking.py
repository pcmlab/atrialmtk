"""
.. _surface_picking_example:

Picking a Point on the Surface of a Mesh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This example demonstrates how to pick meshes using
:func:`surface_mesh_picking() <pyvista.Plotter.enable_surface_picking>`.

This allows you to pick points on the surface of a mesh.

"""
############ IMPORTANT!!!! ############
# Before starting make sure you select points in the correct order
# The order is as follows:
	# 1. Lateral wall 
	# 2. Septal wall
	# 3. LSPV body junction
	# 4. RSPV body junction

# For no PV version and RA
        # 1. LSPV roof, SVC underneath
        # 2. RSPV roof, IVC underneath
        # 3. LAA lateral wall, RAA lateral wall
        # 4. FO septal wall, FO septal wall
        # 5. LSPV path choice, SVC roof
        # 6. RSPV path choice, IVC roof 
	
# Importing libraries

import pyvista as pv
import numpy as np
import pandas as pd

# Define path
DataPath='/Volumes/Elements_CR/atrialmtk/Examples/Example-Biatrial-MRI/2Landmarks/RA_Mesh1/'
MeshName='Clipped.stl'
OutputName='Landmarks.txt'


###############################################################################
# Reading in personal mesh (add stl file)
# Make sure path is correct one to the .stl file

# LA
mesh = pv.read(DataPath + MeshName)


###############################################################################
# Enable a callback that creates a cube at the clicked point and add a label at
# the point as well it.

Allpoints = []

def callback(point):
    """Create a cube and a label at the click point."""
    mesh = pv.read(DataPath + MeshName)
    pl.add_mesh(mesh, style='wireframe', color='r')
    pl.add_point_labels(point, [f"{point[0]:.2f}, {point[1]:.2f}, {point[2]:.2f}"])
    print(point)
    Allpoints.append(point)
	


pl = pv.Plotter()
pl.add_mesh(mesh, show_edges=True)
pl.enable_surface_picking(callback=callback, left_clicking=False, show_point=False)
pl.show()

print(Allpoints)

# Turning array into dataframe and saving the df as a .txt file
Allpoints_df = pd.DataFrame(Allpoints)
Allpoints_df.to_csv(DataPath + OutputName, header=None, 
                    index=None, sep=',') # check path is correct
