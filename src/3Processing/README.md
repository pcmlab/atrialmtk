# Atrial Model Toolkit Pipelines using Universal Atrial Coordinates

These files represent different pipelines for creating the following: 


**Workflow 1:** 
Assumes there is atrial tissue in the model for the pulmonary veins, vena cava, appendages and coronary sinus. 
* mri-la.sh - left atrial surface model with fibres including four pulmonary veins and the left atrial appendage
* mri-ra.sh - right atrial surface model with fibres including two vena cava, coronary sinus and right atrial appendage
* mri-biatrial.sh - biatrial bilayer surface model with fibres and interatrial connections (mri-la.sh and mri-ra.sh must be run first)

**Workflow 2:** 
Assumes each of the following structures have been clipped and removed from the model: the pulmonary veins, vena cava, left atrial appendage and coronary sinus. 
These files should be executed sequentially to generate bilayer and volumetric models. 
* ct-1-extractsurfaces.sh - file to separate surfaces from CT-derived statistical shape atlas into endocardial and epicardial
* ct-2-bilayer.sh - file to generate a biatrial bilayer model with regions, fibres and interatrial connections 
* ct-3-volumetric.sh - file to generate a biatrial volumetric model with regions, fibres and interatrial connections


Each of the files follow a similar format, calculating Universal Atrial Coordinates, and using these for model construction. We detail here the steps in mri-la.sh.

### Example-LeftAtrium, mri-la.sh

**stage1** (trans_mod.py, from Ed Vigmond):

* Convert surface mesh from clipping step from stl to carp format to remesh it to simulation grade

* inputs: Clipped.stl, Regions.txt, Landmarks.txt

* outputs: Labelled.pts, Labelled.elem, Labelled.fcon, Lablled.lon (carp txt)

**stage2** (meshtool convert - rescale, https://bitbucket.org/aneic/meshtool/src/master/):

* Rescale the mesh into the required units for the simulation

* inputs: Labelled.pts, Labelled.elem, Labelled.fcon, Lablled.lon (carp txt)

* outputs: Labelled.pts, Labelled.elem, Labelled.fcon, Lablled.lon (rescaled 1000 - to micrometres, carp txt) (Note: overwrites previous outputs if they are not copied elsewhere before running this line)

**stage3** (meshtool resample):

* Remesh to the resolution required for simulation

* inputs: Labelled.pts, Labelled.elem, Labelled.fcon, Lablled.lon (rescaled to micrometres, carp txt)

* outputs: Labelled.pts, Labelled.elem, Labelled.fcon, Lablled.lon (rescaled to micrometres, remeshed to average edge length 300 micrometres, carp txt) (Note: overwrites previous outputs if they are not copied elsewhere before running this line)

**stage4** (meshtool clean topology):

* Clean the mesh from bad topology definitions

* inputs: Labelled.pts, Labelled.elem, Labelled.fcon, Lablled.lon (rescaled to micrometres, remeshed to average edge length 300 micrometres, carp txt)

* outputs: Labelled.pts, Labelled.elem, Labelled.fcon, Lablled.lon (rescaled to micrometres, remeshed to average edge length 300 micrometres, topology cleaned, carp txt) (Note: overwrites previous outputs if they are not copied elsewhere before running this line)

**stage5** (meshtool clean quality):

* Deform mesh elements to reach the required quality threshold value for simulations

* inputs: Labelled.pts, Labelled.elem, Labelled.fcon, Lablled.lon (rescaled to micrometres, remeshed to average edge length 300 micrometres, topology cleaned, carp txt)

* outputs: Labelled.pts, Labelled.elem, Labelled.fcon, Lablled.lon (rescaled to micrometres, remeshed to average edge length 300 micrometres, topology cleaned, quality cleaned, carp txt) (Note: overwrites previous outputs if they are not copied elsewhere before running this line) 

So overall…

**stages1-5** (trans_mod,py, meshtool convert, meshtool resample, meshtool clean topology, meshtool clean quality):

* Convert the surface mesh provided from the Clipping step into carp format, and remesh it to provide a surface mesh that meets the standards required for the simulations

* inputs: Clipped.stl

* outputs: Labelled.pts, Labelled.elem, Labelled.fcon, Lablled.lon (mesh)

**stage6** (labels_la_1.py):

* Calculate boundaries required to run Laplace-Dirichlet solves for region definition

* inputs: Labelled.pts, Labelled.elem, Labelled.fcon, Lablled.lon (mesh)

* outputs: MV.txt, PV1.vtx, PV2.vtx, PV3.vtx, PV4.vtx, LAA.vtx

**stage7** (laplace):

* Run Laplace-Dirichlet solves for region definition

* inputs: Labelled.pts, Labelled.elem, Labelled.fcon, Lablled.lon (mesh), MV.txt, PV1.vtx, PV2.vtx, PV3.vtx, PV4.vtx, LAA.vtx

* carpf_laplace_PV1.par, carpf_laplace_PV2.par, carpf_laplace_PV3.par, carpf_laplace_PV4.par, carpf_laplace_LAA.par (copied from UAC folder to run Laplace-Dirichlet simulations)

* outputs:  MV_PV1, MV_PV2, MV_PV3, MV_PV4, MV_LAA (folders containing outputs of the simulations)

**stage8** (labels_la_2.py):

* Threshold Laplace-Dirichlet solves to define regions 

* inputs: Labelled.pts, Labelled.elem, Labelled.fcon, Lablled.lon (mesh), MV_PV1/phie.igb, MV_PV2/phie.igb, MV_PV3/phie.igb, MV_PV4/phie.igb, MV_LAA/phie.igb

* outputs: Labels.dat, Labelled_Labels.elem, Labelled_Labels.pts, T_labels_2.vtk, Labelled.elem (Labelled_Labels.elem copied and renamed to Labelled.elem - overwrites original Labelled.elem if not copied elsewhere beforehand. The mesh now contains region labels)

**stage9** (1_la.py):

* Calculate boundary conditions for Universal Atrial Coordinates approximation

* inputs: Labelled. elem (Labelled_Labels.elem), Labelled.pts, Labelled.lon, Labelled.fcon, Landmarks.txt

* outputs: PAbc1.vtx, PAbc2.vtx, LSbc1.vtx, LSbc2.vtx, BorderNodes.vtx

**stage10** (laplace old UAC):

* Run Laplace-Dirichlet solves for Universal Atrial Coordinates approximation

* inputs: Labelled. elem (Labelled_Labels.elem), Labelled.pts, Labelled.lon, Labelled.fcon, PAbc1.vtx, PAbc2.vtx, LSbc1.vtx, LSbc2.vtx

* carpf_laplace_LS.par, carpf_laplace_PA.par (copied from UAC folder to run Laplace-Dirichlet simulations)

* outputs: PA_UAC_N2, LR_UAC_N2 (folders containing simulation outputs)
  
**stage11** (2a_la.py - new UAC):

* Calculate boundary conditions for Universal Atrial Coordinates calculation

* inputs: Labelled. elem (Labelled_Labels.elem), Labelled.pts, Labelled.lon, Labelled.fcon, PA_UAC_N2/phie.igb, LR_UAC_N2/phie.igb, BorderNodes.vtx, PAbc2.vtx, LSbc1.vtx, LSbc2.vtx, PAbc1.vtx, Landmarks.txt 

* outputs: Test_Split.vtk, PosteriorMesh , Test_Ant.vtk, AnteriorMesh, Lablled_Coords_2D_Rescaling_N3.vtk, LSPV_LS.vtx, LSPV_PA.vtx, P_Checker_PA.vtx, P_Checker_LS.vtx, A_Checker_PA.vtx, A_Checker_LS.vtx,  Post_Strength_Test_PA1.vtx, Post_Strength_Test_LS1.vtx, Ant_Strength_Test_PA1.vtx, Ant_Strength_Test_LS1.vtx

**stage12** (laplace new UAC):

* Run Laplace-Dirichlet solves for Universal Atrial Coordinates calculation

* inputs: PosteriorMesh.elem, PosteriorMesh.pts, PosteriorMesh.lon, AnteriorMesh.elem, AnteriorMesh.pts, AnteriorMesh.lon, Post_Strength_Test_PA1.vtx, Post_Strength_Test_LS1.vtx, Ant_Strength_Test.PA1.vtx, Ant_Strength_Test_LS1.vtx

* carpf_laplace_single_LR_A.par, carpf_laplace_single_LR_P.par, carpf_laplace_single_UD_A.par,  carpf_laplace_single_UD_P.par (copied from UAC folder to run Laplace-Dirichlet simulations)

* outputs: LR_Ant_UAC, LR_Post_UAC, UD_Ant_UAC, UD_Post_UAC (folders containing outputs from simulations)

**stage13** (2b_la.py - new UAC pt.2):

* Calculate Universal Atrial Coordinates 

* inputs: Labelled. elem (Labelled_Labels.elem), Labelled.pts, Labelled.lon, Labelled.fcon, UD_Post_UAC/phie.igb, UD_Post_UAC/phie.igb, LR_Post_UAC/phie.igb, PosteriorMesh.pts, AnteriorMesh.pts, AnteriorMesh.lon,UD_Ant_UAC/phie.igb, LR_Ant_UAC/phie.igb, LabelledCoords_2D_Rescaling_N3.pts, 

* outputs: LabelledCoords_2D_Rescaling_v3_C.pts, LabelledCoords_2D_Rescaling_v3_C.elem, LabelledCoords_2D_Rescaling_v3_C.fcon, LabelledCoords_2D_Rescaling_v3_C.lon, LabelledCoords_2D_Rescaling_v3_C.vtk (UAC representation of mesh)

**stage14** (fibre_mapping.py):

* Use Universal Atrial Coordinates to map fibres from the desired atlas 

* inputs: Labelled. elem (Labelled_Labels.elem), Labelled.pts, Labelled.lon, Labelled.fcon, LabelledCoords_2D_Rescaling_v3_C.pts

* outputs: Fibre_l.vtk, Fibre_l.pts, Fibre_l.elem, Fibre_l.lon, Fibre_l.vec, Fibre_l.vpts (mesh with fibres)

**stage15** (lat_field.py):

* Caculate initial conditions for spiral wave re-entry

* inputs: Fibre_l, LabelledCoords_2D_Rescaling_v3_C.pts, LabelledCoords_2D_Rescaling_v3_C.elem, LabelledCoords_2D_Rescaling_v3_C.lon, LabelledCoords_2D_Rescaling_v3_C.fcon

* outputs: LAT_Spiral4_B.dat, BiLAT_Spiral4.dat (initial conditions used to start spiral wave re-entries)
