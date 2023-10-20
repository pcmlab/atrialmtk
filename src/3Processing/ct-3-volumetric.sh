#!/bin/bash

PROJECT=/Volumes/Elements_CR/atrialmtk/src/3Processing/UAC_Codes/UAC_Codes-refactorisation
DATA="/Volumes/Elements_CR/atrialmtk/Examples/Example-Biatrial-CT-shape-model/3Processing"
export PYTHONPATH=$PROJECT

UnitsRescale=1000


cp "$PROJECT/laplace_files/carpf_laplace_alpha.par" "$DATA/LA_vol/carpf_laplace_alpha.par"
cp "$PROJECT/laplace_files/carpf_laplace_beta.par" "$DATA/LA_vol/carpf_laplace_beta.par"
cp "$PROJECT/laplace_files/carpf_laplace_EE.par" "$DATA/LA_vol/carpf_laplace_EE.par"


docker run --rm --volume="$DATA/LA_vol":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool convert -imsh=LA_vol -omsh=LA_only -ofmt=carp_txt -scale=1000 -ifmt=carp_txt

docker run --rm --volume="$DATA/LA_vol":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract surface -msh=LA_only -surf=LA_only


python $PROJECT/scripts/volumetric_vtx.py "$DATA/LA_vol/" "$DATA/LA_endo/" "$DATA/LA_epi/" LA_only LA_only LA_only

docker run --rm --volume="$DATA/LA_vol":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_alpha.par -simID Alpha
docker run --rm --volume="$DATA/LA_vol":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_beta.par -simID Beta
docker run --rm --volume="$DATA/LA_vol":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_EE.par -simID EndoEpi

python $PROJECT/scripts/volumetric_write_mesh.py "$DATA/LA_vol/" "$DATA/LA_endo/" "$DATA/LA_epi/" LA_only LA_only LA_only

cp "$DATA/LA_vol/LA_only.elem" "$DATA/LA_vol/Mesh_UAC_3D.elem"
cp "$DATA/LA_vol/LA_only.surf" "$DATA/LA_vol/Mesh_UAC_3D.surf"

cp "$PROJECT/laplace_files/carpf_laplace_alpha_RA.par" "$DATA/RA_vol/carpf_laplace_alpha.par"
cp "$PROJECT/laplace_files/carpf_laplace_beta_RA.par" "$DATA/RA_vol/carpf_laplace_beta.par"
cp "$PROJECT/laplace_files/carpf_laplace_EE_RA.par" "$DATA/RA_vol/carpf_laplace_EE.par"

docker run --rm --volume="$DATA/RA_vol":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool convert -imsh=RA_vol -omsh=RA_only -ofmt=carp_txt -scale=1000 -ifmt=carp_txt

docker run --rm --volume="$DATA/RA_vol":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract surface -msh=RA_only -surf=RA_only

python $PROJECT/scripts/volumetric_vtx.py "$DATA/RA_vol/" "$DATA/RA_endo/" "$DATA/RA_epi/" RA_only RA_only RA_only

docker run --rm --volume="$DATA/RA_vol":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP openCARP +F carpf_laplace_alpha.par -simID Alpha
docker run --rm --volume="$DATA/RA_vol":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP openCARP +F carpf_laplace_beta.par -simID Beta
docker run --rm --volume="$DATA/RA_vol":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP openCARP +F carpf_laplace_EE.par -simID EndoEpi

python $PROJECT/scripts/volumetric_write_mesh.py "$DATA/RA_vol/" "$DATA/RA_endo/" "$DATA/RA_epi/" RA_only RA_only RA_only


cp "$DATA/RA_vol/RA_only.elem" "$DATA/RA_vol/Mesh_UAC_3D.elem"
cp "$DATA/RA_vol/RA_only.surf" "$DATA/RA_vol/Mesh_UAC_3D.surf"


python $PROJECT/scripts/fibre_mapping.py "$DATA/LA_endo/" "$PROJECT/fibre_files/la/endo/l/" "$PROJECT/laplace_files/" LA_only Labelled.lon Fibre_Labarthe

python $PROJECT/scripts/fibre_mapping.py "$DATA/LA_epi/" "$PROJECT/fibre_files/la/epi/l/" "$PROJECT/laplace_files/" Labelled Labelled.lon Fibre_Labarthe

python $PROJECT/scripts/fibre_mapping.py "$DATA/LA_epi/" "$PROJECT/fibre_files/bb/la/" "$PROJECT/laplace_files/" LA_only Labelled.lon Fibre_Labarthe_bb

#now threshold fibres LA

python $PROJECT/scripts/fibre_mapping_volumetric_threshold_bb.py "$DATA/LA_vol/" LA_only Mesh_UAC_3D "$DATA/LA_endo/" LA_only Fibre_Labarthe.vec "$DATA/LA_epi/" LA_only Fibre_Labarthe.vec "$DATA/LA_epi/" LA_only Fibre_Labarthe_bb.vec Fibre_T.vec

#Fibre mapping bilayer for RA

python $PROJECT/scripts/fibre_mapping.py "$DATA/RA_epi/" "$PROJECT/fibre_files/ra/epi/l/" "$PROJECT/laplace_files/" Labelled Labelled.lon Fibre_Labarthe

python $PROJECT/scripts/fibre_mapping.py "$DATA/RA_epi/" "$PROJECT/fibre_files/bb/ra/" "$PROJECT/laplace_files/" Labelled Labelled.lon Fibre_Labarthe_bb


python $PROJECT/scripts/fibre_mapping_volumetric_threshold_ra_bb.py "$DATA/RA_vol/" RA_only Mesh_UAC_3D "$DATA/RA_epi/" RA_only Fibre_Labarthe_Labelled_Endo "$DATA/RA_epi/" RA_only Fibre_Labarthe_Labelled_Epi "$DATA/RA_epi/" RA_only Fibre_Labarthe_Labelled_BB "$DATA/RA_epi/" Fibre_T.vec



cp "$DATA/LA_vol/Fibres_Threshold.pts" "$DATA/RA_vol/Fibres_Threshold_LA.pts"
cp "$DATA/LA_vol/Fibres_Threshold.elem" "$DATA/RA_vol/Fibres_Threshold_LA.elem"
cp "$DATA/LA_vol/Fibres_Threshold.lon" "$DATA/RA_vol/Fibres_Threshold_LA.lon"


docker run --rm --volume="$DATA/RA_vol":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool merge meshes -msh1=Fibres_Threshold_LA -msh2=Fibres_Threshold -ofmt=carp_txt -outmsh=MergeVol_Threshold

docker run --rm --volume="$DATA/RA_vol":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract mesh -msh=MergeVol_Threshold -tags=11 -submsh=LA_mesh -ofmt=carp_txt





docker run --rm --volume="$DATA/RA_vol":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract surface -msh=LA_mesh -surf=LA_mesh

docker run --rm --volume="$DATA/RA_vol":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract mesh -msh=MergeVol_Threshold -tags=4 -submsh=RA_mesh -ofmt=carp_txt

docker run --rm --volume="$DATA/RA_vol":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract surface -msh=RA_mesh -surf=RA_mesh

docker run --rm --volume="$DATA/RA_vol":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract mesh -msh=MergeVol_Threshold -tags=5,8,9 -submsh=RA_structures -ofmt=carp_txt

docker run --rm --volume="$DATA/RA_vol":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract mesh -msh=MergeVol_Threshold -tags=10 -submsh=BB_structures -ofmt=carp_txt


docker run --rm --volume="$DATA/RA_vol":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract surface -msh=RA_structures -surf=RA_structures

##code to visualise fibres
python $PROJECT/scripts/biatrial_volumetric_visualisation_iac.py "$DATA/RA_vol/" MergeVol_Threshold
