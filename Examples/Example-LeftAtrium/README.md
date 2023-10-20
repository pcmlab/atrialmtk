## Steps for LA mesh generation:

**Imaging Data:**

1. Start with an Atrial Challenge Imaging dataset (see ImagingData/Dataset.txt). You should have the raw MRI file lgemri.nrrd and corresponding segmentation mask laendo.nrrd. See Examples/Example-LeftAtrium/0ImagingData
    
    These files can be visualised using different software including 3DSlicer, CemrgApp, itksnap (e.g. in itksnap: do “File/load main image” - lgemri.nrrd, “Segmentation/load segmentation” - laendo.nrrd, “update”)
    
    They can also be exported as surface meshes (.vtk format) using any of these software. 
    
    E.g. itksnap:

   ![itk-snap](https://github.com/pcmlab/atrialmtk/blob/main/images/itk-snap-LA.png?raw=true)

    

**Clipping**

1. Set up the folder Examples/Example-LeftAtrium/1Clipping/LA_mesh1, and copy across the data from Examples/Example-LeftAtrium/0ImagingData
2. Load segmented.vtk in paraview and open the mesh by clipping it at the pulmonary veins and mitral valve, extracting the surface after each clip. If you open the state file Clips.pvsm with the provided segmented.vtk file (use the 'Search names under specified directory' option); you should be able to see the sphere clippers that we used to clip the veins/valve in the example (and in the figures below). 
3. Triangulate and extract the final surface, and save it as Clipped.stl, which you will also be able to open and view in paraview. This will be the input for the landmarking step.
        
    Clip1 (Clips.pvsm): right supervior pulmonary vein
   ![Clip1](https://github.com/pcmlab/atrialmtk/blob/main/images/rspv.png?raw=true)
        
    Clip2 (Clips.pvsm): right inferior pulmonary vein
    ![Clip2](https://github.com/pcmlab/atrialmtk/blob/main/images/ripv.png?raw=true)
       
    Clip3 (Clips.pvsm): left superior pulmonary vein
    ![Clip3](https://github.com/pcmlab/atrialmtk/blob/main/images/lspv.png?raw=true)
      
    Clip4 (Clips.pvsm): left inferior pulmonary vein
    ![Clip4](https://github.com/pcmlab/atrialmtk/blob/main/images/lipv.png?raw=true)
       
    Clip5 (Clips.pvsm): mitral valve
    ![Clip5](https://github.com/pcmlab/atrialmtk/blob/main/images/mv.png?raw=true)

   
    Clipped.stl:
    ![Clipped](https://github.com/pcmlab/atrialmtk/blob/main/images/clipped.png?raw=true)

    
**Landmarking**

1. Make a new folder Examples/Example-LeftAtrium/2Landmarking/LA_Mesh1 and copy Clipped.stl across.
2. You will need to select several points on the surface in order to inform the model, so that it can correctly separate the surface into its different atrial regions in the next step. First you will need to roughly select some points on each key anatomical landmark.
3. If you are following these steps for the first time you will first need to create a conda environment for point selection (follow the steps for “1a. Point picking environment” below). If you have already created this environment, skip this step.
4. Activate the point picking environment by doing:

    ```
    conda activate pointpicking
    ```
    
5. Update the file path in the script **Rough_Point_Picking.py** in the src/2Landmarks folder. DataPath should be the path of the folder containing your Clipped.stl.
6. cd to src/2Landmarking and run:

    ```
    python Rough_Point_Picking.py
    ```
![pyvista-regions](https://github.com/pcmlab/atrialmtk/blob/main/images/la-regions-pyvista.png?raw=true)    
    
7. Follow the steps for “general point picking” from the “Instructions for landmark selection” below. Select the 6 landmarks by right clicking or pressing P with the cursor at the desired location,  making sure to select the points in the same order as listed in the python script, and close the window. The coordinates will be saved automatically in a text file Regions.txt.
8. Next, you will need to provide the model with some more specific information. Update DataPath in Refined_Point_Picking.py to match the folder containing Clipped.stl.
9. With the point picking environment still active, run:
    
    ```
    python Refined_Point_Picking.py
    ```
    and follow the steps for “specific point selection” in the “Instructions for landmark selection”below to select these 4 additional points. Again make sure to do this in the same order as in the script. These points will need to be selected more carefully, so make sure to read through the steps fully, and make use of the reference images below. When you are finished close the window. You will find these coordinates saved in the text file Landmarks.txt.
   
   ![pyvista-landmarks](https://github.com/pcmlab/atrialmtk/blob/main/images/la-landmarks-pyvista.png?raw=true)       
  
11. Once you have finished with the point selection, deactivate the conda environment using:

    ```
    conda deactivate
    ```
    
    
    

**Processing**

1. Copy the surface Clipped.stl and the coordinate files Regions.txt and Landmarks.txt from Examples/Example-LeftAtrium/2Landmarking/LA_Mesh1 to the processing folder Examples/Example-LeftAtrium/3Processing/LA_Mesh1.
2. Update the paths in the file **mri-la.sh**, found in the folder src/3Processing, so that PROJECT is the path to the folder containing the universal atrial coordinate (UAC) codes (src/3Processing/UAC_Codes), DATA contains the path to the processing step for the example (Examples/Example-LeftAtrium/3Processing), and ConvertFilesLoc is the path to the python script trans_mod.py (src/3Processing/UAC_Codes).
3. Set index=1 to correspond to LA_Mesh1 of the Example case, or change it to 2 if you have set up a second version of the example, LA_Mesh2 etc.
4. Follow the steps on the general README for the docker installation of openCARP, in order to run the script mri-la.sh. (You only need to do this step once.)
5. In order to run mri-la.sh, which automatically separates the surface into its different atrial regions, adds atrial fibres and generates a bilayer model, you will need to create and activate the second conda environment for UAC (follow the steps “1b. Universal Atrial Coordinates (UAC) environment” in the general README). If you have already created this environment, you will only need to activate it:

    ```
    conda activate uac
    ```
    
    
    
7. With the uac environment active, cd to the src folder and run mri-la.sh to generate a mesh with regions, fibres, atrial coordinates, and initial conditions. Each stage will produce a number of output files. These can be checked by running the commands from mri-la.sh step-by-step and comparing with the example outputs. (Note you may need to change the permissions for mri-la.sh to allow it to run as an executable e.g. in Linux).
8. When this code has finished, type: conda deactivate
9. You can open the meshes in Paraview or meshalyzer and compare to the examples below. 

 ![lauac](https://github.com/pcmlab/atrialmtk/blob/main/images/lauac.png?raw=true) 

Please also see notes here if you would like more information on the different stages of the code: 
[UAC_README](/src/3Processing/README.md)


**Simulation**

1. From the outputs of the Processing step, you will need the mesh Fibre_l.pts, Fibre_l.elem, the fibre file, Fibre_l.lon, and the LAT field, LAT_Spiral4_B.dat. Copy them from Examples/Example-LeftAtrium/3Processing/LA_Mesh1 to the simulation folder Examples/Example-LeftAtrium/4Simulation/LA_Mesh1
2. Also copy AF_Simulation.par from src to the Simulation folder.
3. Use the following command, updated to have the path to your simulation folder there, to run the simulation step of the model:
docker run --rm --volume=/Volumes/Elements_CR/atrialmtk/Examples/Example-LeftAtrium/4Simulation/LA_Mesh1:/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F AF_Simulation.par -simID AF

 
---

**********************************************************************Instructions for landmark selection**********************************************************************

1. **********************************Landmark selection for the left atrium**********************************
    
    **General point selection**
    
    First, roughly select a point somewhere on each anatomical landmark that is listed in Rough_Point_Picking.py. Make sure to do this in the same order as in the script:
    
    i) Right superior pulmonary vein (RSPV)
    
    ii) Right inferior pulmonary vein (RIPV)
    
    iii) Left inferior pulmonary vein (LIPV)
    
    iv) Left superior pulmonary vein (LSPV)
    
    v) Left atrial appendage tip (LAA tip)
    
    vi) Left atrial appendage base (LAA base)
    
    View 1 (general points):
    ![LA_Regions_roof](https://github.com/pcmlab/atrialmtk/blob/main/images/LA_Regions_roof.png?raw=true) 
    
    View 2 (general points):
    ![LA_Regions_LAA](https://github.com/pcmlab/atrialmtk/blob/main/images/LA_Regions_LAA.png?raw=true) 
    
    **Specific point selection**
    
    Now select points at specific locations in the following order (as listed in the python script Refined_Point_Picking.py) 
    
    i) On the lateral wall, in line with the LSPV, posterior of the LAA; 
  View 1 (specific points):
   ![LA_Landmarks_lateralwall](https://github.com/pcmlab/atrialmtk/blob/main/images/LA_Landmarks_lateralwall.png?raw=true) 

   
    ii) On the septal wall, in line with the RSPV;
   View 2 (specific points):
   ![LA_Landmarks_septalwall](https://github.com/pcmlab/atrialmtk/blob/main/images/LA_Landmarks_septalwall.png?raw=true) 
   
   iii) Approximately at the intersection of the LSPV and the body of the left atrium (LA), level with the roof
   View 3 (specific points):
   ![LA_Landmarks_roof](https://github.com/pcmlab/atrialmtk/blob/main/images/LA_Landmarks_roof.png?raw=true) 
    
    iv) Approximately at the intersection between the RSPV and the LA, at the level of the roof

    
    The atrial regions will then be identified automatically using Laplace solvers in CARPentry

