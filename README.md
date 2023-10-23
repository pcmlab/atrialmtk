# atrialmtk
**Atrial Modelling Toolkit for Constructing Bilayer and Volumetric Atrial Models at Scale**

To enable large in silico trials and personalised model predictions on clinical timescales, it is imperative that models can be constructed quickly and reproducibly. First, we aimed to overcome the challenges of constructing cardiac models at scale through developing a robust, open-source pipeline for bilayer and volumetric atrial models. Second, we aimed to investigate the effects of fibres, fibrosis and model representation (surface or volumetric) on fibrillatory dynamics. To construct bilayer and volumetric models, we extended our previously developed coordinate system to incorporate transmurality, atrial regions, and fibres (rule-based or data driven diffusion tensor MRI). We demonstrate the generalisation of these methods by creating a cohort of 1000 biatrial bilayer and volumetric models derived from CT data, as well as models from MRI, and electroanatomic mapping. Fibrillatory dynamics diverged between bilayer and volumetric simulations across the CT cohort (correlation co-efficient for phase singularity maps: LA 0.27±0.19, RA 0.41±0.14). Adding fibrotic remodelling stabilised re-entries and reduced the impact of model type in the LA (LA: 0.52±0.20, RA: 0.36±0.18). The choice of fibre field has a small effect on paced activation data (<12ms), but a larger effect on fibrillatory dynamics. Overall, we developed an open-source user-friendly pipeline for generating atrial models from imaging or electroanatomic mapping enabling in silico clinical trials at scale.

![Model pipeline](https://github.com/pcmlab/atrialmtk/blob/main/images/Figure1_Schematicv2.jpg?raw=true)

# Installation instructions

## ************Docker************

You will need to make sure docker is installed and running on your machine to install and run openCARP and meshtool. You will need to install docker with sudo access:

Install docker: https://docs.docker.com/engine/install/

To install with sudo access, please follow these steps post-installation:

1. Create the docker group:

   ```
   sudo groupadd docker
    ```
2. Add the user to the docker group:

    ```
   sudo usermod -aG docker $USER
    ```
    
3. Activate the changes to the groups:

   ```
   newgrp docker
    ```
       
4. You will need to restart your computer to see if docker works without sudo access. 
5. Verify that you can run docker commands without sudo:
    
    ```
   docker run —rm hello-world
    ```
    
    This command downloads a test image and runs it in a container. When the container runs, it prints a message and exits.
    

## ****************openCARP****************

Follow these steps to install openCARP as a docker container on the users local machine (an alternative method can be found here: https://opencarp.org/download/installation for the direct installation of openCARP onto a HPC).

Docker installation: https://opencarp.org/download/installation#installation-of-opencarp-docker-containers

1. Once docker is installed on your machine with sudo access, do: 
    
    ```
   docker pull [docker.opencarp.org/opencarp/opencarp:latest](http://docker.opencarp.org/opencarp/opencarp:latest)
    ```
    
3. Then check you have some output from: 
    
   ```
   docker run -it [docker.opencarp.org/opencarp/opencarp:latest](http://docker.opencarp.org/opencarp/opencarp:latest)
   ```
    
    The command line should change to something like: root@……..:/openCARP#
    
5. Type exit and press enter to close the openCARP interpreter. You will now be able to use openCARP. 


## **meshtool**

meshtool may be included with openCARP depending on the installation method

If you have installed openCARP as a docker container you will already have meshtool. Otherwise you can download it here: https://bitbucket.org/aneic/meshtool/src/master/ and add the location to your bashrc file, or refer back to the openCARP installation pages: https://opencarp.org/download/installation.


## **Instructions for conda environments (environments only need to be created once)**

**To generate conda environments:**
    
1a. **Point picking environment:** 
    
    conda create --name pointpicking python=3.10 pandas numpy
    
    conda activate pointpicking
    
    python -m pip install pyvista
    
(you will then be able to run the code using: python Rough_Point_Picking.py)

1b. **Universal Atrial Coordinates (UAC) environment:** 
    
cd to src/3Processing/UAC_Codes/UAC_Codes-refactorisation
    
    conda env create -f environment.yml
    
(Note: the “Collecting package metadata” and “Solving environment” steps take a few minutes each, so the environment may take up to 10 minutes to create)
    
    conda activate uac
(if you have already installed openCARP, you will now be able to run mri-la.sh from the src folder using ./mri-la.sh)
    

# **Usage** 

We have included the following examples: 

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
