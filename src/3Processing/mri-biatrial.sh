#!/bin/bash


PROJECT=/Volumes/Elements_CR/atrialmtk/src/3Processing/UAC_Codes/UAC_Codes-refactorisation
DATA="/Volumes/Elements_CR/atrialmtk/Examples/Example-Biatrial-MRI/3Processing"
ConvertFilesLoc=/Volumes/Elements_CR/atrialmtk/src/3Processing/UAC_Codes

PVLabelT=0.8
LAALabelT=0.35
UnitsRescale=1000

index=1

export PYTHONPATH=$PROJECT
ES_dir="$PROJECT/extra_structures"
RAendofib_dir="$PROJECT/fibre_files/ra/endo/l/"
RAepifib_dir="$PROJECT/fibre_files/ra/epi/l/"
RAbbfib_dir="$PROJECT/fibre_files/bb/ra/"
LAendofib_dir="$PROJECT/fibre_files/la/endo/l/"
LAepifib_dir="$PROJECT/fibre_files/la/epi/l/"
LAbbfib_dir="$PROJECT/fibre_files/bb/la/"

MeshName="Labelled"
LandmarksRA="Landmarks.txt"
RegionRA="Regions.txt"
LandmarksLA="Landmarks.txt"
RegionLA="Regions.txt"


RAPath="$DATA/RA_Mesh"$index/
LAPath="$DATA/LA_Mesh"$index/

echo $RAPath

'''
python $PROJECT/scripts/scalar_mapping_bilayer.py  "$RAPath" "$ES_dir/" Labelled Labelled_Extra_RA_BB Labelled_Coords_2D_Rescaling_v3_C Labelled_Coords_2D_Rescaling_v3_C Extra_SAN.dat MappedScalar_SAN.dat
python $PROJECT/scripts/scalar_mapping_bilayer.py  "$RAPath" "$ES_dir/" Labelled Labelled_Extra_RA_BB Labelled_Coords_2D_Rescaling_v3_C Labelled_Coords_2D_Rescaling_v3_C Extra_CT.dat MappedScalar_CT.dat
python $PROJECT/scripts/scalar_mapping_bilayer.py  "$RAPath" "$ES_dir/" Labelled Labelled_Extra_RA_BB Labelled_Coords_2D_Rescaling_v3_C Labelled_Coords_2D_Rescaling_v3_C Extra_PM.dat MappedScalar_PM.dat
python $PROJECT/scripts/scalar_mapping_bilayer.py  "$RAPath" "$ES_dir/" Labelled Labelled_Extra_RA_BB Labelled_Coords_2D_Rescaling_v3_C Labelled_Coords_2D_Rescaling_v3_C Extra_BB.dat MappedScalar_BB.dat
python $PROJECT/scripts/scalar_mapping_bilayer.py  "$LAPath" "$ES_dir/" Labelled Labelled_Extra_LA_BB LA_Labelled_Coords_2D_Rescaling_v3_C Labelled_Coords_2D_Rescaling_v3_C Extra_BB_LA.dat MappedScalar_BB_LA.dat


###Fibre mapping bilayer for LA
python $PROJECT/scripts/fibre_mapping_bilayer_bb.py $LAPath "$LAendofib_dir" "$LAepifib_dir" "$LAbbfib_dir" Labelled Labelled.lon Labelled.lon Labelled.lon Fibre_Labarthe_ 100 200 11


##Fibre mapping bilayer for RA
python $PROJECT/scripts/fibre_mapping_bilayer_ra_bb.py $RAPath "$RAendofib_dir" "$RAepifib_dir" "$RAbbfib_dir" Labelled Labelled Labelled Labelled Labelled.lon Labelled.lon Labelled.lon Fibre_Labarthe_


## remove RA endocardial shell

docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract mesh -msh=Fibre_Labarthe_Bilayer -tags=1,2,3,5,6,7,8,9 -submsh=Bilayer2 -ofmt=carp_txt
docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract mesh -msh=Fibre_Labarthe_Bilayer -tags=10 -submsh=Bilayer2_BB -ofmt=carp_txt
docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract unreachable -msh=Bilayer2_BB -submsh=Bilayer2_BB_c -ofmt=carp_txt -ifmt=carp_txt


if test -f "$RAPath/Bilayer2_BB_c.part0.pts"; then
    docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool convert -imsh=Bilayer2_BB_c.part0 -omsh=Bilayer2_BB_c -ofmt=carp_txt -ifmt=carp_txt
else
    docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool convert -imsh=Bilayer2_BB -omsh=Bilayer2_BB_c -ofmt=carp_txt -ifmt=carp_txt
fi



## remove LA endocardial shell
docker run --rm --volume="$LAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract mesh -msh=Fibre_Labarthe_Bilayer -tags=11,12,13,14,21,22,23,24,25,26,27,28 -submsh=Bilayer2 -ofmt=carp_txt

docker run --rm --volume="$LAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract mesh -msh=Fibre_Labarthe_Bilayer -tags=10 -submsh=Bilayer2_BB -ofmt=carp_txt

docker run --rm --volume="$LAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract unreachable -msh=Bilayer2_BB -submsh=Bilayer2_BB_c -ofmt=carp_txt -ifmt=carp_txt



if test -f "$LAPath/Bilayer2_BB_c.part0.pts"; then
    docker run --rm --volume="$LAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool convert -imsh=Bilayer2_BB_c.part0 -omsh=Bilayer2_BB_c -ofmt=carp_txt -ifmt=carp_txt
else
    docker run --rm --volume="$LAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool convert -imsh=Bilayer2_BB -omsh=Bilayer2_BB_c -ofmt=carp_txt -ifmt=carp_txt
fi




## combine LA & RA & add RA elements
cp "$LAPath/Bilayer2.pts" "$RAPath/Bilayer2_LA.pts"
cp "$LAPath/Bilayer2.elem" "$RAPath/Bilayer2_LA.elem"
cp "$LAPath/Bilayer2.lon" "$RAPath/Bilayer2_LA.lon"


docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool merge meshes -msh1=Bilayer2_LA -msh2=Bilayer2 -ofmt=carp_txt -outmsh=Bilayer_Combined


docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool merge meshes -msh1=Bilayer_Combined -msh2=Bilayer2_BB_c -ofmt=carp_txt -outmsh=Bilayer_Combined


cp "$LAPath/Bilayer2_BB_c.pts" "$RAPath/Bilayer2_BB_c_LA.pts"
cp "$LAPath/Bilayer2_BB_c.elem" "$RAPath/Bilayer2_BB_c_LA.elem"
cp "$LAPath/Bilayer2_BB_c.lon" "$RAPath/Bilayer2_BB_c_LA.lon"

docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool merge meshes -msh1=Bilayer_Combined -msh2=Bilayer2_BB_c_LA -ofmt=carp_txt -outmsh=Bilayer_Combined




## add IAC
python $PROJECT/scripts/lat_field_biatrial_bilayer_meshtool.py $LAPath $RAPath LAT_Spiral4_B.dat Labelled Bilayer_Combined



docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract mesh -msh=Bilayer_Combined -tags=11,13,21,23,25,27 -submsh=LA_endo -ofmt=carp_txt


docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract mesh -msh=Bilayer_Combined -tags=12,14,22,24,26,28 -submsh=LA_epi -ofmt=carp_txt


docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract mesh -msh=Bilayer_Combined -tags=1,2,5,6,7 -submsh=RA_epi_s -ofmt=carp_txt

docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract mesh -msh=Bilayer_Combined -tags=3,8,9 -submsh=RA_structures -ofmt=carp_txt


docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool extract mesh -msh=Bilayer_Combined -tags=10 -submsh=BB_mesh -ofmt=carp_txt

'''

#make bilayer mesh with lines & visualisation
python $PROJECT/scripts/biatrial_bilayer_iac_lines_visualisation.py "$LAPath" "$RAPath" "$RAPath" Bilayer_Combined Labelled Fibre_Labarthe_Bilayer $LandmarksLA $RegionRA 1000 1


