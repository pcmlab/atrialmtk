# **Usage** 

Each of the examples has a separate README file detailing all steps, please check instead the folders for these. 

## **Example-LeftAtrium**

This example generates a left atrial surface model from the Atrial Challenge dataset. All steps are explained in the following categories:

0ImagingData

1Clipping

2Landmarks

3Processing

4Simulation



## **Example-Biatrial-MRI**

This example generates a biatrial bilayer model from the Atrial Challenge dataset. This example can be run after the Example-LeftAtrium to add the right atrial component of the model, and interatrial connections as a biatrial bilayer model. The steps are detailed in the same categories as for Example-LeftAtrium. 



## **Example-Biatrial-CT-shape-model**

This example generates a biatrial bilayer model and a biatrial volumetric model from a CT-derived statistical shape model. In this case, we start with the atrial surfaces, and progress through the steps to generate two versions of the model in a bilayer format (triangular mesh) and a volumetric format (tetrahedral mesh). The steps are as follows:

0SurfaceMeshData 

2Landmarks

3Processing

4Simulation

(these meshes do not require any clipping). 
