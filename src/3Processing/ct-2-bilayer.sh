#!/bin/bash

PROJECT=/Volumes/Elements_CR/atrialmtk/src/3Processing/UAC_Codes
DATA="/Volumes/Elements_CR/atrialmtk/Examples/Example-Biatrial-CT-shape-model/3Processing"
export PYTHONPATH=$PROJECT

UnitsRescale=1000



cp "/Volumes/Elements_CR/atrialmtk/Examples/Example-Biatrial-CT-shape-model/2Landmarks/LA_Mesh1/outputs/Landmarks.txt" "$DATA/LA_endo/Landmarks.txt"
cp "/Volumes/Elements_CR/atrialmtk/Examples/Example-Biatrial-CT-shape-model/2Landmarks/LA_Mesh1/outputs/Regions.txt" "$DATA/LA_endo/Regions.txt"

cp "/Volumes/Elements_CR/atrialmtk/Examples/Example-Biatrial-CT-shape-model/2Landmarks/LA_Mesh1/outputs/Landmarks.txt" "$DATA/LA_epi/Landmarks.txt"
cp "/Volumes/Elements_CR/atrialmtk/Examples/Example-Biatrial-CT-shape-model/2Landmarks/LA_Mesh1/outputs/Regions.txt" "$DATA/LA_epi/Regions.txt"

cp "/Volumes/Elements_CR/atrialmtk/Examples/Example-Biatrial-CT-shape-model/2Landmarks/RA_Mesh1/outputs/Landmarks.txt" "$DATA/RA_epi/Landmarks.txt"
cp "/Volumes/Elements_CR/atrialmtk/Examples/Example-Biatrial-CT-shape-model/2Landmarks/RA_Mesh1/outputs/Regions.txt" "$DATA/RA_epi/Regions.txt"


cp "/Volumes/Elements_CR/atrialmtk/Examples/Example-Biatrial-CT-shape-model/2Landmarks/RA_Mesh1/outputs/Landmarks.txt" "$DATA/RA_endo/Landmarks.txt"
cp "/Volumes/Elements_CR/atrialmtk/Examples/Example-Biatrial-CT-shape-model/2Landmarks/RA_Mesh1/outputs/Regions.txt" "$DATA/RA_endo/Regions.txt"



cp "$PROJECT/laplace_files/carpf_laplace_single_LR_P.par" "$DATA/LA_endo/carpf_laplace_single_LR_P.par"
cp "$PROJECT/laplace_files/carpf_laplace_single_UD_P.par" "$DATA/LA_endo/carpf_laplace_single_UD_P.par"
cp "$PROJECT/laplace_files/carpf_laplace_single_LR_A.par" "$DATA/LA_endo/carpf_laplace_single_LR_A.par"
cp "$PROJECT/laplace_files/carpf_laplace_single_UD_A.par" "$DATA/LA_endo/carpf_laplace_single_UD_A.par"

cp "$PROJECT/laplace_files/carpf_laplace_LS_4Ch.par" "$DATA/LA_endo/carpf_laplace_LS.par"
cp "$PROJECT/laplace_files/carpf_laplace_PA_4Ch.par" "$DATA/LA_endo/carpf_laplace_PA.par"

cp "$PROJECT/laplace_files/carpf_laplace_LAA_4Ch.par" "$DATA/LA_endo/carpf_laplace_LAA.par"


docker run --rm --volume="$DATA/LA_endo":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool convert -imsh=LA_endo.vtk -omsh=LA_only -ofmt=carp_txt -scale=$UnitsRescale

## Old UAC approximation
python $PROJECT/scripts/1_la_4ch.py "$DATA/LA_endo/" LA_only Landmarks.txt Regions.txt $UnitsRescale

docker run --rm --volume="$DATA/LA_endo":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_PA.par -simID PA_UAC_N2
docker run --rm --volume="$DATA/LA_endo":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_LS.par -simID LR_UAC_N2

echo "=====New UAC====="
python $PROJECT/scripts/2a_la_4ch.py "$DATA/LA_endo/" LA_only Landmarks.txt Regions.txt $UnitsRescale

docker run --rm --volume="$DATA/LA_endo":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_single_LR_A.par -simID LR_Ant_UAC
docker run --rm --volume="$DATA/LA_endo":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_single_LR_P.par -simID LR_Post_UAC
docker run --rm --volume="$DATA/LA_endo":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_single_UD_A.par -simID UD_Ant_UAC
docker run --rm --volume="$DATA/LA_endo":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_single_UD_P.par -simID UD_Post_UAC


echo "=====New UAC part 2====="
python $PROJECT/scripts/2b_la.py "$DATA/LA_endo/" "$PROJECT/fibre_files/la/endo/" "$PROJECT/laplace_files/" LA_only 11 13 21 23 25 27 1

echo "=====Add fibres2 ====="
python $PROJECT/scripts/fibre_mapping.py "$DATA/LA_endo/" "$PROJECT/fibre_files/la/endo/l/" "$PROJECT/laplace_files/" LA_only Labelled.lon Fibre_Labarthe




#repeat this for epi!


cp "$PROJECT/laplace_files/carpf_laplace_single_LR_P.par" "$DATA/LA_epi/carpf_laplace_single_LR_P.par"
cp "$PROJECT/laplace_files/carpf_laplace_single_UD_P.par" "$DATA/LA_epi/carpf_laplace_single_UD_P.par"
cp "$PROJECT/laplace_files/carpf_laplace_single_LR_A.par" "$DATA/LA_epi/carpf_laplace_single_LR_A.par"
cp "$PROJECT/laplace_files/carpf_laplace_single_UD_A.par" "$DATA/LA_epi/carpf_laplace_single_UD_A.par"

cp "$PROJECT/laplace_files/carpf_laplace_LS_4Ch.par" "$DATA/LA_epi/carpf_laplace_LS.par"
cp "$PROJECT/laplace_files/carpf_laplace_PA_4Ch.par" "$DATA/LA_epi/carpf_laplace_PA.par"

cp "$PROJECT/laplace_files/carpf_laplace_LAA_4Ch.par" "$DATA/LA_epi/carpf_laplace_LAA.par"

docker run --rm --volume="$DATA/LA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool convert -imsh=LA_epi.vtk -omsh=LA_only -ofmt=carp_txt -scale=$UnitsRescale

## Old UAC approximation
python $PROJECT/scripts/1_la_4ch.py "$DATA/LA_epi/" LA_only Landmarks.txt Regions.txt $UnitsRescale

docker run --rm --volume="$DATA/LA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_PA.par -simID PA_UAC_N2
docker run --rm --volume="$DATA/LA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_LS.par -simID LR_UAC_N2


echo "=====New UAC====="
python $PROJECT/scripts/2a_la_4ch.py "$DATA/LA_epi/" LA_only Landmarks.txt Regions.txt $UnitsRescale

docker run --rm --volume="$DATA/LA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_single_LR_A.par -simID LR_Ant_UAC
docker run --rm --volume="$DATA/LA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_single_LR_P.par -simID LR_Post_UAC
docker run --rm --volume="$DATA/LA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_single_UD_A.par -simID UD_Ant_UAC
docker run --rm --volume="$DATA/LA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_single_UD_P.par -simID UD_Post_UAC


echo "=====New UAC part 2====="
python $PROJECT/scripts/2b_la.py "$DATA/LA_epi/" "$PROJECT/fibre_files/la/endo/" "$PROJECT/laplace_files/" LA_only 11 13 21 23 25 27 1


echo "=====Add fibres2 ====="
python $PROJECT/scripts/fibre_mapping.py "$DATA/LA_epi/" "$PROJECT/fibre_files/la/epi/l/" "$PROJECT/laplace_files/" LA_only Labelled.lon Fibre_Labarthe


#repeat this for RA epi!


cp "$PROJECT/laplace_files/carpf_laplace_single_LR_P.par" "$DATA/RA_epi/carpf_laplace_single_LR_P.par"
cp "$PROJECT/laplace_files/carpf_laplace_single_UD_P.par" "$DATA/RA_epi/carpf_laplace_single_UD_P.par"
cp "$PROJECT/laplace_files/carpf_laplace_single_LR_A.par" "$DATA/RA_epi/carpf_laplace_single_LR_A.par"
cp "$PROJECT/laplace_files/carpf_laplace_single_UD_A.par" "$DATA/RA_epi/carpf_laplace_single_UD_A.par"

cp "$PROJECT/laplace_files/carpf_laplace_LS_4Ch_RA.par" "$DATA/RA_epi/carpf_laplace_LS.par"
cp "$PROJECT/laplace_files/carpf_laplace_PA_4Ch_RA.par" "$DATA/RA_epi/carpf_laplace_PA.par"

cp "$PROJECT/laplace_files/carpf_laplace_LAA_4Ch.par" "$DATA/RA_epi/carpf_laplace_LAA.par"
cp "$PROJECT/laplace_files/carpf_laplace_RAA.par" "$DATA/RA_epi/carpf_laplace_RAA.par"

docker run --rm --volume="$DATA/RA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool convert -imsh=RA_epi.vtk -omsh=RA_only -ofmt=carp_txt -scale=$UnitsRescale

python $PROJECT/scripts/labels_ra_noraa.py "$DATA/RA_epi/" RA_only Landmarks.txt Regions.txt 0.6 $UnitsRescale

## Old UAC approximation
python $PROJECT/scripts/1_ra_4ch.py "$DATA/RA_epi/" RA_only_RAA Landmarks.txt Regions.txt $UnitsRescale

docker run --rm --volume="$DATA/RA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_PA.par -simID PA_UAC_N2
docker run --rm --volume="$DATA/RA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_LS.par -simID LR_UAC_N2

python $PROJECT/scripts/2a_ra_4ch_noraa.py "$DATA/RA_epi/" RA_only_RAA Landmarks.txt Regions.txt $UnitsRescale

docker run --rm --volume="$DATA/RA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_single_LR_P.par -simID LR_Post_UAC
docker run --rm --volume="$DATA/RA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_single_UD_P.par -simID UD_Post_UAC

python $PROJECT/scripts/2b_ra_noraa.py "$DATA/RA_epi/" RA_only_RAA 1 6 7 5 2 1




cp "$PROJECT/laplace_files/carpf_laplace_single_LR_P.par" "$DATA/RA_endo/carpf_laplace_single_LR_P.par"
cp "$PROJECT/laplace_files/carpf_laplace_single_UD_P.par" "$DATA/RA_endo/carpf_laplace_single_UD_P.par"
cp "$PROJECT/laplace_files/carpf_laplace_single_LR_A.par" "$DATA/RA_endo/carpf_laplace_single_LR_A.par"
cp "$PROJECT/laplace_files/carpf_laplace_single_UD_A.par" "$DATA/RA_endo/carpf_laplace_single_UD_A.par"

cp "$PROJECT/laplace_files/carpf_laplace_LS_4Ch_RA.par" "$DATA/RA_endo/carpf_laplace_LS.par"
cp "$PROJECT/laplace_files/carpf_laplace_PA_4Ch_RA.par" "$DATA/RA_endo/carpf_laplace_PA.par"

cp "$PROJECT/laplace_files/carpf_laplace_LAA_4Ch.par" "$DATA/RA_endo/carpf_laplace_LAA.par"
cp "$PROJECT/laplace_files/carpf_laplace_RAA.par" "$DATA/RA_endo/carpf_laplace_RAA.par"

docker run --rm --volume="$DATA/RA_endo":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool convert -imsh=RA_endo.vtk -omsh=RA_only -ofmt=carp_txt -scale=$UnitsRescale

python $PROJECT/scripts/labels_ra_noraa.py "$DATA/RA_endo/" RA_only Landmarks.txt Regions.txt 0.6 $UnitsRescale

## Old UAC approximation
python $PROJECT/scripts/1_ra_4ch.py "$DATA/RA_endo/" RA_only_RAA Landmarks.txt Regions.txt $UnitsRescale

docker run --rm --volume="$DATA/RA_endo":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_PA.par -simID PA_UAC_N2
docker run --rm --volume="$DATA/RA_endo":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_LS.par -simID LR_UAC_N2

echo "=====New UAC====="
python $PROJECT/scripts/2a_ra_4ch_noraa.py "$DATA/RA_endo/" RA_only_RAA Landmarks.txt Regions.txt $UnitsRescale

docker run --rm --volume="$DATA/RA_endo":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_single_LR_P.par -simID LR_Post_UAC
docker run --rm --volume="$DATA/RA_endo":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_single_UD_P.par -simID UD_Post_UAC

python $PROJECT/scripts/2b_ra_noraa.py "$DATA/RA_endo/" RA_only_RAA 1 6 7 5 2 1


#Fibres next

echo "=====Add fibres2 ====="
python $PROJECT/scripts/fibre_mapping.py "$DATA/RA_epi/" "$PROJECT/fibre_files/ra/epi/l/" "$PROJECT/laplace_files/" RA_only_RAA Labelled.lon Fibre_Labarthe


python $PROJECT/scripts/scalar_mapping_bilayer.py  "$DATA/RA_epi/" "$PROJECT/extra_structures/" RA_only_RAA Labelled_Extra_RA_BB Labelled_Coords_2D_Rescaling_v3_C Labelled_Coords_2D_Rescaling_v3_C Extra_SAN.dat MappedScalar_SAN.dat
python $PROJECT/scripts/scalar_mapping_bilayer.py  "$DATA/RA_epi/" "$PROJECT/extra_structures/" RA_only_RAA Labelled_Extra_RA_BB Labelled_Coords_2D_Rescaling_v3_C Labelled_Coords_2D_Rescaling_v3_C Extra_CT.dat MappedScalar_CT.dat
python $PROJECT/scripts/scalar_mapping_bilayer.py  "$DATA/RA_epi/" "$PROJECT/extra_structures/" RA_only_RAA Labelled_Extra_RA_BB Labelled_Coords_2D_Rescaling_v3_C Labelled_Coords_2D_Rescaling_v3_C Extra_PM.dat MappedScalar_PM.dat
python $PROJECT/scripts/scalar_mapping_bilayer.py  "$DATA/RA_epi/" "$PROJECT/extra_structures/" RA_only_RAA Labelled_Extra_RA_BB Labelled_Coords_2D_Rescaling_v3_C Labelled_Coords_2D_Rescaling_v3_C Extra_BB.dat MappedScalar_BB.dat

#LA too
python $PROJECT/scripts/scalar_mapping_bilayer.py "$DATA/LA_epi/" "$PROJECT/extra_structures/" LA_only Labelled_Extra_LA_BB LA_Labelled_Coords_2D_Rescaling_v3_C Labelled_Coords_2D_Rescaling_v3_C Extra_BB_LA.dat MappedScalar_BB_LA.dat


python $PROJECT/scripts/fibre_mapping_bilayer_bb.py "$DATA/LA_epi/" "$PROJECT/fibre_files/la/epi/l/" "$PROJECT/fibre_files/la/endo/l/" "$PROJECT/fibre_files/bb/la/" LA_only Labelled.lon Labelled.lon Labelled.lon Fibre_Labarthe_ -100 100 3


##Fibre mapping bilayer for RA
python $PROJECT/scripts/fibre_mapping_bilayer_ra_bb.py "$DATA/RA_epi/" "$PROJECT/fibre_files/ra/endo/l/" "$PROJECT/fibre_files/ra/epi/l/" "$PROJECT/fibre_files/bb/ra/" Labelled Labelled Labelled RA_only_RAA Labelled.lon Labelled.lon Labelled.lon Fibre_Labarthe_


docker run --rm --volume="$DATA/RA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract mesh -msh=Fibre_Labarthe_Bilayer -tags=1,2,3,8,9 -submsh=Bilayer2 -ofmt=carp_txt


docker run --rm --volume="$DATA/RA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract mesh -msh=Fibre_Labarthe_Bilayer -tags=10 -submsh=Bilayer2_BB -ofmt=carp_txt


docker run --rm --volume="$DATA/RA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract unreachable -msh=Bilayer2_BB -submsh=Bilayer2_BB_c -ofmt=carp_txt -ifmt=carp_txt



if test -f "$DATA/RA_epi/Bilayer2_BB_c.part0.pts"; then
    docker run --rm --volume="$DATA/RA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool convert -imsh=Bilayer2_BB_c.part0 -omsh=Bilayer2_BB_c -ofmt=carp_txt -ifmt=carp_txt
else
    docker run --rm --volume="$DATA/RA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool convert -imsh=Bilayer2_BB -omsh=Bilayer2_BB_c -ofmt=carp_txt -ifmt=carp_txt
fi




## combine LA & RA & add RA elements
docker run --rm --volume="$DATA/LA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract mesh -msh=Fibre_Labarthe_Bilayer -tags=11,12 -submsh=Bilayer2 -ofmt=carp_txt


docker run --rm --volume="$DATA/LA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract mesh -msh=Fibre_Labarthe_Bilayer -tags=10 -submsh=Bilayer2_BB -ofmt=carp_txt


docker run --rm --volume="$DATA/LA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest  meshtool extract unreachable -msh=Bilayer2_BB -submsh=Bilayer2_BB_c -ofmt=vtk -ifmt=carp_txt



if test -f "$DATA/LA_epi/Bilayer2_BB_c.part0.pts"; then
    docker run --rm --volume="$DATA/LA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool convert -imsh=Bilayer2_BB_c.part0 -omsh=Bilayer2_BB_c -ofmt=carp_txt -ifmt=carp_txt
else
    docker run --rm --volume="$DATA/LA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool convert -imsh=Bilayer2_BB -omsh=Bilayer2_BB_c -ofmt=carp_txt -ifmt=carp_txt
fi




cp "$DATA/LA_epi/Bilayer2.pts" "$DATA/RA_epi/Bilayer2_LA.pts"
cp "$DATA/LA_epi/Bilayer2.elem" "$DATA/RA_epi/Bilayer2_LA.elem"
cp "$DATA/LA_epi/Bilayer2.lon" "$DATA/RA_epi/Bilayer2_LA.lon"


## combine LA & RA & add RA elements
docker run --rm --volume="$DATA/RA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool merge meshes -msh1=Bilayer2_LA -msh2=Bilayer2 -ofmt=carp_txt -outmsh=Bilayer_Combined


cp "$DATA/LA_epi/Bilayer2_BB_c.pts" "$DATA/RA_epi/Bilayer2_BB_c_LA.pts"
cp "$DATA/LA_epi/Bilayer2_BB_c.elem" "$DATA/RA_epi/Bilayer2_BB_c_LA.elem"
cp "$DATA/LA_epi/Bilayer2_BB_c.lon" "$DATA/RA_epi/Bilayer2_BB_c_LA.lon"


docker run --rm --volume="$DATA/RA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool merge meshes -msh1=Bilayer_Combined -msh2=Bilayer2_BB_c_LA -ofmt=carp_txt -outmsh=Bilayer_Combined


docker run --rm --volume="$DATA/RA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool merge meshes -msh1=Bilayer_Combined -msh2=Bilayer2_BB_c -ofmt=carp_txt -outmsh=Bilayer_Combined





docker run --rm --volume="$DATA/RA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract mesh -msh=Bilayer_Combined -tags=11 -submsh=LA_endo -ofmt=vtk
docker run --rm --volume="$DATA/RA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool convert -imsh=LA_endo.vtk -omsh=LA_endo -ofmt=carp_txt

docker run --rm --volume="$DATA/RA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract mesh -msh=Bilayer_Combined -tags=12 -submsh=LA_epi -ofmt=vtk
docker run --rm --volume="$DATA/RA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool convert -imsh=LA_epi.vtk -omsh=LA_epi -ofmt=carp_txt

docker run --rm --volume="$DATA/RA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract mesh -msh=Bilayer_Combined -tags=1,2 -submsh=RA_epi_s -ofmt=vtk
docker run --rm --volume="$DATA/RA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool convert -imsh=RA_epi_s.vtk -omsh=RA_epi_s -ofmt=carp_txt






docker run --rm --volume="$DATA/RA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract mesh -msh=Bilayer_Combined -tags=3,8,9 -submsh=RA_structures -ofmt=carp_txt
docker run --rm --volume="$DATA/RA_epi":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract mesh -msh=Bilayer_Combined -tags=10 -submsh=BB_mesh -ofmt=carp_txt



cp "$DATA/LA_epi/LA_only.pts" "$DATA/LA_epi/Labelled.pts"
cp "$DATA/LA_epi/LA_only.elem" "$DATA/LA_epi/Labelled.elem"
cp "$DATA/LA_epi/LA_only.lon" "$DATA/LA_epi/Labelled.lon"

cp "$DATA/RA_epi/RA_only.pts" "$DATA/RA_epi/Labelled.pts"
cp "$DATA/RA_epi/RA_only.elem" "$DATA/RA_epi/Labelled.elem"
cp "$DATA/RA_epi/RA_only.lon" "$DATA/RA_epi/Labelled.lon"

python $PROJECT/scripts/biatrial_bilayer_iac_lines_visualisation_4ch.py "$DATA/LA_epi/" "$DATA/RA_epi/" "$DATA/RA_epi/" Bilayer_Combined Labelled Fibre_Labarthe_Bilayer Landmarks.txt Regions.txt $UnitsRescale 3
