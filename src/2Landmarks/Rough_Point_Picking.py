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
	# 1. RSPV, IVC path choice
	# 2. RIPV, CS
	# 3. LIPV, IVC/SVC/roof
	# 4. LSPV, SVC path choice 
	# 5. Appendage tip, RAA tip
	# 6. Appendage base, RAA base


# Importing libraries

import pyvista as pv
import numpy as np
import pandas as pd

# Define path
DataPath='/Volumes/Elements_CR/atrialmtk/Examples/Example-Biatrial-MRI/2Landmarks/RA_Mesh1/'
MeshName='Clipped.stl'
OutputName='Regions.txt'


###############################################################################
# Reading in personal mesh (add stl file)
# Make sure path is correct one to the .stl file and has correct name
# IMPORTANT! If using Linux make sure you use .stl file, NOT .vtk

# LA 
mesh = pv.read(DataPath + MeshName)

###############################################################################
# Enable a callback that creates a cube at the clicked point and add a label at
# the point as well it.

Allpoints = []

def callback(point):
    """Create a cube and a label at the click point."""
    # LA
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

# Turning array into dataframe and saving the df as a .txt filepythobn
Allpoints_df = pd.DataFrame(Allpoints)
Allpoints_df.to_csv(DataPath + OutputName, 
header=None, 
                    index=None, sep=',') # check path is correct

