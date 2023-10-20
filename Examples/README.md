# **Usage** 

Each of the examples has a separate README file detailing all steps, please check instead the folders for these. 

## **[Example-LeftAtrium](/Examples/Example-LeftAtrium/README.md)**
This example generates a left atrial surface model from the Atrial Challenge dataset.  All steps are explained in the following categories:

0ImagingData

1Clipping

2Landmarks

3Processing

4Simulation

![LA model](https://github.com/pcmlab/atrialmtk/blob/main/images/laP.png?raw=true)


## **[Example-Biatrial-MRI](/Examples/Example-Biatrial-MRI/README.md)**

This example generates a biatrial bilayer model from the Atrial Challenge dataset. This example can be run after the Example-LeftAtrium to add the right atrial component of the model, and interatrial connections as a biatrial bilayer model. The steps are detailed in the same categories as for Example-LeftAtrium. 

![Biatrial model](https://github.com/pcmlab/atrialmtk/blob/main/images/biatrial3.png?raw=true)

## **[Example-Biatrial-CT-shape-model](/Examples/Example-Biatrial-CT-shape-model/README.md)**

This example generates a biatrial bilayer model and a biatrial volumetric model from a CT-derived statistical shape model. In this case, we start with the atrial surfaces, and progress through the steps to generate two versions of the model in a bilayer format (triangular mesh) and a volumetric format (tetrahedral mesh). The steps are as follows:

0SurfaceMeshData 

2Landmarks

3Processing

4Simulation

(these meshes do not require any clipping). 

![Biatrial CT model](https://github.com/pcmlab/atrialmtk/blob/main/images/ct-biatrial.png?raw=true)
