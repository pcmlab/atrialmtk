# UAC_Codes
Atrial coordinates codes for adding fibres to meshes and registering data.

These codes are integrated into Cemrgapp software in docker form, so that they can be used in the atrial model generation pipeline. https://cemrg.com/software/cemrgapp.html

## Methodology
We developed an atrial coordinate system: https://pubmed.ncbi.nlm.nih.gov/31026761/
We updated it to fix the locations of the atrial structures to standard locations: https://pubmed.ncbi.nlm.nih.gov/32458222/


## Pre-processing
The codes assume you have one of the following set-ups:
- Left atrial surface with 4 pulmonary veins labelled and left atrial appendage labelled
- Left atrial surface with mesh clipped at the 4 pulmonary vein openings and at the left atrial appendage opening
- Right atrial surface with superior vena cava, inferior vena cava, coronary sinus and right atrial appendage labelled
- Right atrial surface with mesh clipped at the superior vena cava and inferior vena cava openings

These options are shown in the following Figure:



## Dependencies
To run these atrial coordinates codes, we recommend using a conda virtual environment with the packages installed by pip. From this folder:
```
conda create -n uac -c conda-forge python=3.8 -y
conda activate uac
pip install -r requirements.txt
```
On Apple Silicon Macs, create the environment with `CONDA_SUBDIR=osx-64 conda create -n uac -c conda-forge python=3.8 -y` and run `conda config --env --set subdir osx-64` after activating it. See the main [README](../../../README.md#conda-environments) for details.

To calculate the Laplace solves required for the method, we use opencarp. We suggest downloading a docker and calling as shown in the example usage code. Installation instructions: https://opencarp.org/download/installation#installation-of-opencarp-docker-containers
Tutorial here: https://opencarp.org/documentation/examples/02_ep_tissue/13_laplace

## Running Tests
In order to run tests, call *runtests* script:
```
./runtests.sh
```

To generate coverage report:
```
./runtests.sh --coverage --html
```

## Overview of the codes


## Example Usage
We have included an example in the Example_DTMRI_1. To run it, edit Example_Usage.sh to include your paths.
