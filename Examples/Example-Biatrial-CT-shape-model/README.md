##  Steps for Biatrial bilayer and volumetric model generation. 

This example generates a biatrial bilayer model and a biatrial volumetric model from a CT-derived statistical shape model. 
In this case, we start with the atrial surfaces, and progress through the steps to generate two versions of the model in a bilayer format (triangular mesh) and a volumetric format (tetrahedral mesh). 
The steps are as follows:

0SurfaceMeshData

2Landmarks

3Processing

4Simulation

(these meshes do not require any clipping).

In this example, there is no tissue in the model for the pulmonary veins, vena cava, coronary sinus, or LAA. 

## **Surface and volumetric mesh data**


**Landmarking**

1. Make a new folder Examples/Example-Biatrial-CT-shape-model/2Landmarking/LA_Mesh1 and Examples/Example-Biatrial-CT-shape-model/2Landmarking/RA_Mesh1.
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


**Simulation**

1. From the outputs of the Processing step, you will need the mesh Fibre_l.pts, Fibre_l.elem, the fibre file, Fibre_l.lon, and the LAT field, LAT_Spiral4_B.dat for an RA model. For a biatrial simulatiom, use Bilayer_Combined_all_Lines_IAC.pts, Bilayer_Combined_all_Lines_IAC.elem, Bilayer_Combined_all_Lines_IAC.lon and BiatrialcombinedLAT_Spiral4_B.dat. 
Copy them from Examples/Example-Biatrial-MRI/3Processing/RA_Mesh1 to the simulation folder Examples/Example-Biatrial-MRI/4Simulation/RA_Mesh1
2. Also copy AF_Simulation.par from src to the Simulation folder.
3. Use the following command, updated to have the path to your simulation folder there, to run the simulation step of the model:
docker run ...

    

---

**********************************************************************Instructions for landmark selection**********************************************************************

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                                                                                                                          
1. **********************************Landmark selection for the left atrium**********************************
    
    **General point selection**
    
    First, select a point somewhere on each anatomical landmark that is listed in Rough_Point_Picking.py. The RAA tip point is used as a boundary condition location of the Laplace-Dirichlet solve. Make sure to do this in the same order as in the script:
    
    i) **Inferior vena cava path choice** 
    
    ii) **Coronary sinus**
    
    iii) **Roof between SVC and IVC**
    
    iv) **Superior vena cava path choice** 
    
    v) **Right atrial appendage tip (RAA tip)**
    
    vi) **Right atrial appendage base (RAA base)**

Here points i) and iv) are used to say which path is on the lateral wall. The ring of nodes at the vena cava and RA body junction is split into two paths between the lowest and highest points (using specific points i and v or ii and vi described below). 

View 1 (general points):

View 2 (general points):

**Specific point selection**
Now select points at specific locations in the following order (as listed in the python script Refined_Point_Picking.py)    

i) **Lowest point on the SVC/RA ostia.** Two paths are calculated from the SVC/RA ostia point closest to the TV to point (v). I.e. the ring of nodes at the ostia is split into two paths. 
     
  View 1 (specific points):

   ii) **Lowest point on the IVC/RA ostia.** Two paths are calculated from the IVC/RA ostia point closest to the TV to point (vi). 
   
   View 2 (specific points):
   
   iii) **On the septal wall, in line with the SVC, septal of the RAA.*** A geodesic path is calculated from the bottom of the SVC ostia to the TV through this point. You want to ensure this path is septal of the RAA, so that the entire RAA is assigned to the lateral component of the UAC. 
   
   View 3 (specific points):
   
   iv) **On the septal wall, in line with the IVC.** Similar to point iii, a geodesic path is calculated from the bottom of the IVC ostia to the TV through this point. This point is also used to assign line connections in the biatrial bilayer connections at the fossa ovalis, so should be chosen accordingly.  

   v) **Highest point on the SVC/RA ostia (at the roof).**  Two paths are calculated from the SVC/RA ostia point closest to the TV (i) to this point. 

   vi) **Highest point on the IVC/RA ostia (at the roof).** Two paths are calculated from the IVC/RA ostia point closest to the TV (i) to this point.  

                                                                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                                                                          
1. **********************************Landmark selection for the right atrium**********************************
    
    **General point selection**
    
    First, select a point somewhere on each anatomical landmark that is listed in Rough_Point_Picking.py. The RAA tip point is used as a boundary condition location of the Laplace-Dirichlet solve. Make sure to do this in the same order as in the script:
    
    i) **Inferior vena cava path choice** 
    
    ii) **Coronary sinus**
    
    iii) **Roof between SVC and IVC**
    
    iv) **Superior vena cava path choice** 
    
    v) **Right atrial appendage tip (RAA tip)**
    
    vi) **Right atrial appendage base (RAA base)**

Here points i) and iv) are used to say which path is on the lateral wall. The ring of nodes at the vena cava and RA body junction is split into two paths between the lowest and highest points (using specific points i and v or ii and vi described below). 

View 1 (general points):

View 2 (general points):

**Specific point selection**
Now select points at specific locations in the following order (as listed in the python script Refined_Point_Picking.py)    

i) **Lowest point on the SVC/RA ostia.** Two paths are calculated from the SVC/RA ostia point closest to the TV to point (v). I.e. the ring of nodes at the ostia is split into two paths. 
     
  View 1 (specific points):

   ii) **Lowest point on the IVC/RA ostia.** Two paths are calculated from the IVC/RA ostia point closest to the TV to point (vi). 
   
   View 2 (specific points):
   
   iii) **On the septal wall, in line with the SVC, septal of the RAA.*** A geodesic path is calculated from the bottom of the SVC ostia to the TV through this point. You want to ensure this path is septal of the RAA, so that the entire RAA is assigned to the lateral component of the UAC. 
   
   View 3 (specific points):
   
   iv) **On the septal wall, in line with the IVC.** Similar to point iii, a geodesic path is calculated from the bottom of the IVC ostia to the TV through this point. This point is also used to assign line connections in the biatrial bilayer connections at the fossa ovalis, so should be chosen accordingly.  

   v) **Highest point on the SVC/RA ostia (at the roof).**  Two paths are calculated from the SVC/RA ostia point closest to the TV (i) to this point. 

   vi) **Highest point on the IVC/RA ostia (at the roof).** Two paths are calculated from the IVC/RA ostia point closest to the TV (i) to this point.  

    
The atrial regions will then be identified automatically using Laplace solvers in CARPentry
   
