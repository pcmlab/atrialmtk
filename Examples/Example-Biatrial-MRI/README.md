# Steps for Biatrial bilayer model generation:

Please follow the steps for LA model generation first, and generate all of the LA_Mesh1 subfolders.

Next follow these steps for RA mesh generation, and then for biatrial bilayer and volumetric model generation. 

## Steps for RA mesh generation:

**Imaging Data:**

1. Start with an Atrial Challenge Imaging dataset (see ImagingData/Dataset.txt). You should have the raw MRI file lgemri.nrrd and corresponding segmentation mask raendo.nrrd (note that we have segmented the RA ourselves as the label is not provided in the challenge dataset). See Examples/Example-Biatrial-MRI/0ImagingData
    
    These files can be visualised using different software including 3DSlicer, CemrgApp, itksnap (e.g. in itksnap: do “File/load main image” - lgemri.nrrd, “Segmentation/load segmentation” - raendo.nrrd, “update”)
    
    They can also be exported as surface meshes (.vtk format) using any of these software. 
    
    E.g. itksnap:
    ![itk-snap](https://github.com/pcmlab/atrialmtk/blob/main/images/itk-snap-RA.png?raw=true)

**Clipping**

1. Set up the folder Examples/Example-Biatrial-MRI/1Clipping/RA_mesh1, and copy across the RA data from Examples/Example-Biatrial-MRI0ImagingData
2. Load segmented.vtk in paraview and open the mesh by clipping it at the vena cava, coronary sinus and tricuspid valve, extracting the surface after each clip. If you open the state file Clips.pvsm with the provided segmented.vtk file (use the 'Search names under specified directory' option); you should be able to see the sphere clippers that we used to clip the vena cava/valve in the example (also see the figures below). 
3. Triangulate and extract the final surface, and save it as Clipped.stl, which you will also be able to open and view in paraview. This will be the input for the landmarking step.
    
    **segmented.vtk:**
    
    
    **Clip1 (Clips.pvsm): superior vena cava** 
    ![Clip1](https://github.com/pcmlab/atrialmtk/blob/main/images/svc.png?raw=true)
    
    **Clip2 (Clips.pvsm): inferior vena cava** 
    ![Clip2](https://github.com/pcmlab/atrialmtk/blob/main/images/ivc.png?raw=true)
    
    **Clip3 (Clips.pvsm): coronary sinus** 
    ![Clip3](https://github.com/pcmlab/atrialmtk/blob/main/images/cs.png?raw=true)
    
    **Clip4 (Clips.pvsm): tricuspid valve** 
    ![Clip4](https://github.com/pcmlab/atrialmtk/blob/main/images/tv.png?raw=true)
    
    **Clipped.stl:**
    ![Clipped](https://github.com/pcmlab/atrialmtk/blob/main/images/clippedra.png?raw=true)
    

**Landmarking**

1. Make a new folder Examples/Example-Biatrial-MRI/2Landmarking/RA_Mesh1 and copy Clipped.stl across.
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
    
 ![pyvista-regions](https://github.com/pcmlab/atrialmtk/blob/main/images/ra-regions-pyvista.png?raw=true) 
 
7. Follow the steps for “general point picking” from the “Instructions for landmark selection” below. Select the 6 landmarks by right clicking or pressing P with the cursor at the desired location,  making sure to select the points in the same order as listed in the python script, and close the window. The coordinates will be saved automatically in a text file Regions.txt.
8. Next, you will need to provide the model with some more specific information. Update DataPath in Refined_Point_Picking.py to match the folder containing Clipped.stl.
9. With the point picking environment still active, run:
    
    ```
    python Refined_Point_Picking.py
    ```
    
    and follow the steps for “specific point selection” in the “Instructions for landmark selection”below to select these 6 additional points. Again make sure to do this in the same order as in the script. These points will need to be selected more carefully, so make sure to read through the steps fully, and make use of the reference images below. When you are finished close the window. You will find these coordinates saved in the text file Landmarks.txt.
    
10. Once you have finished with the point selection, deactivate the conda environment using:
    
    ```
    conda deactivate
    ```
    
    

**Processing**

1. Copy the surface Clipped.stl and the coordinate files Regions.txt and Landmarks.txt from Examples/Example-Biatrial-MRI/2Landmarking/RA_Mesh1 to the processing folder Examples/Example-Biatrial-MRI/3Processing/RA_Mesh1.
2. Update the paths in the file **mri-ra.sh**, found in the folder src/3Processing, so that PROJECT is the path to the folder containing the universal atrial coordinate (UAC) codes (src/3Processing/UAC_Codes/UAC_Codes-refactorisation), DATA contains the path to the processing step for the example (Examples/Example--Biatrial-MRI/3Processing), and ConvertFilesLoc is the path to the python script trans_mod.py (src/3Processing/UAC_Codes).
3. Set index=1 to correspond to RA_Mesh1 of the Example case, or change it to 2 if you have set up a second version of the example, RA_Mesh2 etc.
4. Follow the steps on the general README for the docker installation of openCARP, in order to run the script mri-ra.sh. (You only need to do this step once.)
5. In order to run mri-ra.sh, which automatically separates the surface into its different atrial regions, adds atrial fibres and generates a bilayer model, you will need to create and activate the second conda environment for UAC (follow the steps “1b. Universal Atrial Coordinates (UAC) environment” in the general README). If you have already created this environment, you will only need to activate it:
    
    ```
    conda activate uac
    ```
    
6. With the uac environment active, cd to the src folder and run ./mri-ra.sh to generate an RA mesh with regions, fibres, atrial coordinates, and initial conditions. Each stage will produce a number of output files. These can be checked by running the commands from mri-ra.sh step-by-step and comparing with the example outputs. (Note you may need to change the permissions for mri-ra.sh to allow it to run as an executable e.g. in Linux).
7. You can open the meshes in Paraview or meshalyzer and compare to the examples.

 
**Finally, run ./mri-biatrial.**
This code will combine the LA and RA meshes to make a bilayer model with atrial regions and fibres. 

When this code has finished, type: conda deactivate

 ![biatrial](https://github.com/pcmlab/atrialmtk/blob/main/images/biatrial.png?raw=true) 
  ![biatrial2](https://github.com/pcmlab/atrialmtk/blob/main/images/biatrial2.png?raw=true) 
   ![biatrial3](https://github.com/pcmlab/atrialmtk/blob/main/images/biatrial3.png?raw=true) 
    ![biatrial4](https://github.com/pcmlab/atrialmtk/blob/main/images/biatrial4.png?raw=true) 

**Simulation**

1. From the outputs of the Processing step, you will need the mesh Fibre_l.pts, Fibre_l.elem, the fibre file, Fibre_l.lon, and the LAT field, LAT_Spiral4_B.dat. Copy them from Examples/Example-LeftAtrium/3Processing/LA_Mesh1 to the simulation folder Examples/Example-LeftAtrium/4Simulation/LA_Mesh1
2. Also copy AF_Simulation.par from src to the Simulation folder.
3. Use the following command, updated to have the path to your simulation folder there, to run the simulation step of the model:
docker run ...
Initial conditions used for simulation: 
  ![biatriallat](https://github.com/pcmlab/atrialmtk/blob/main/images/biatrial_LAT.png?raw=true) 
    

---

**********************************************************************Instructions for landmark selection**********************************************************************

   
