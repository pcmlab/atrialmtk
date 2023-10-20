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


RAPath="$DATA/RA_Mesh"$index
LAPath="$DATA/LA_Mesh"$index

echo $RAPath


#stages 1-5 
python $ConvertFilesLoc/trans_mod.py --from-format=stl --to-format=carp $RAPath/Clipped.stl $RAPath/Labelled

docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool convert -imsh=Labelled -omsh=Labelled -ofmt=carp_txt -ifmt=carp_txt -scale=$UnitsRescale

docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool resample surfmesh -msh=Labelled -avrg=300 -outmsh=Labelled -surf_corr=0.95

docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool clean topology -msh=Labelled -outmsh=Labelled

docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool clean quality -msh=Labelled -thr=0.5 -outmsh=Labelled



#stage 6 
python $PROJECT/scripts/labels_ra_1.py "$RAPath/" "$PROJECT/fibre_files/ra/endo/" "$PROJECT/laplace_files/" Labelled $PVLabelT $LAALabelT $RegionRA $UnitsRescale



#stage 7 
cp "$PROJECT/laplace_files/carpf_laplace_PV2.par" "$DATA/RA_Mesh$index/carpf_laplace_PV2.par"
cp "$PROJECT/laplace_files/carpf_laplace_PV3.par" "$DATA/RA_Mesh$index/carpf_laplace_PV3.par"
cp "$PROJECT/laplace_files/carpf_laplace_PV4.par" "$DATA/RA_Mesh$index/carpf_laplace_PV4.par"
cp "$PROJECT/laplace_files/carpf_laplace_LAA.par" "$DATA/RA_Mesh$index/carpf_laplace_LAA.par"

cd "$DATA/RA_Mesh$index/"
docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_PV2.par -simID MV_PV2
docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_PV3.par -simID MV_PV3
docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_PV4.par -simID MV_PV4
docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_LAA.par -simID MV_LAA




#stage 8 
python $PROJECT/scripts/labels_ra_2.py "$RAPath/" "$PROJECT/fibre_files/ra/endo/" "$PROJECT/laplace_files/" Labelled $PVLabelT $LAALabelT $RegionRA $UnitsRescale



cp "$RAPath/Labelled_Labels.elem" "$RAPath/Labelled.elem"

cp "$PROJECT/laplace_files/carpf_laplace_LS.par" "$RAPath/carpf_laplace_LS.par"
cp "$PROJECT/laplace_files/carpf_laplace_PA.par" "$RAPath/carpf_laplace_PA.par"



#stage 9
python $PROJECT/scripts/1_ra.py "$RAPath/" "$PROJECT/fibre_files/ra/endo/" "$PROJECT/laplace_files/" $MeshName 1 6 7 5 2 $LandmarksRA $RegionRA $UnitsRescale


#stage 10
echo "=====Old UAC approximation - openCARP command line ($i)====="
cd "$DATA/RA_Mesh$index/"
docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_PA.par -simID PA_UAC_N2
docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_LS.par -simID LR_UAC_N2



#stage 11
echo "=====New UAC====="
cp "$PROJECT/laplace_files/carpf_laplace_single_LR_P.par" "$RAPath/carpf_laplace_single_LR_P.par"
cp "$PROJECT/laplace_files/carpf_laplace_single_UD_P.par" "$RAPath/carpf_laplace_single_UD_P.par"
cp "$PROJECT/laplace_files/carpf_laplace_single_LR_A.par" "$RAPath/carpf_laplace_single_LR_A.par"
cp "$PROJECT/laplace_files/carpf_laplace_single_UD_A.par" "$RAPath/carpf_laplace_single_UD_A.par"

python $PROJECT/scripts/2a_ra.py "$RAPath/" "$PROJECT/fibre_files/ra/endo/" "$PROJECT/laplace_files/" $MeshName 1 6 7 5 2 $LandmarksRA $RegionRA $UnitsRescale



#stage 12
cd "$DATA/RA_Mesh$index/"
docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_single_LR_A.par -simID LR_Ant_UAC
docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_single_LR_P.par -simID LR_Post_UAC
docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_single_UD_A.par -simID UD_Ant_UAC
docker run --rm --volume="$RAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_single_UD_P.par -simID UD_Post_UAC



#stage 13
echo "=====New UAC part 2====="
python $PROJECT/scripts/2b_ra.py "$RAPath/" $MeshName 1 6 7 5 2 1




#stage 14
echo "=====Add fibres2 ====="
python $PROJECT/scripts/fibre_mapping.py "$RAPath/" "$PROJECT/fibre_files/ra/epi/l/" "$PROJECT/laplace_files/" Labelled Labelled.lon Fibre_l


#stage 15
echo "=====Add LAT field ($i)====="
python $PROJECT/scripts/lat_field.py "$RAPath/" LAT_Spiral4_B.dat Fibre_l

