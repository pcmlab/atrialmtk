#!/bin/bash

PROJECT=/Volumes/Elements_CR/atrialmtk/src/3Processing/UAC_Codes/UAC_Codes-refactorisation

DATA="/Volumes/Elements_CR/atrialmtk/Examples/Example-Biatrial-CT-shape-model/3Processing"

export PYTHONPATH=$PROJECT

echo $DATA

python $PROJECT/scripts/label_endo_epi_surfaces.py ${DATA}/ LA_1 LA_2 RA_1 RA_2;
