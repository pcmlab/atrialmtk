#!/bin/bash


PROJECT=/Volumes/Elements_CR/atrialmtk/src/3Processing/UAC_Codes/UAC_Codes-refactorisation
DATA="/Volumes/Elements_CR/atrialmtk/Examples/Example-LeftAtrium/3Processing"
ConvertFilesLoc=/Volumes/Elements_CR/atrialmtk/src/3Processing/UAC_Codes

PVLabelT=0.7
LAALabelT=0.25
UnitsRescale=1000

index=1

export PYTHONPATH=$PROJECT
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

echo $LAPath


#stages 1-5 
python $ConvertFilesLoc/trans_mod.py --from-format=stl --to-format=carp $LAPath/Clipped.stl $LAPath/Labelled

docker run --rm --volume="$LAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool convert -imsh=Labelled -omsh=Labelled -ofmt=carp_txt -ifmt=carp_txt -scale=$UnitsRescale

docker run --rm --volume="$LAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool resample surfmesh -msh=Labelled -avrg=300 -outmsh=Labelled -surf_corr=0.95

docker run --rm --volume="$LAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool clean topology -msh=Labelled -outmsh=Labelled

docker run --rm --volume="$LAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest meshtool clean quality -msh=Labelled -thr=0.5 -outmsh=Labelled

#stage 6 
python $PROJECT/scripts/labels_la_1.py "$LAPath/" "$PROJECT/fibre_files/ra/endo/" "$PROJECT/laplace_files/" Labelled $PVLabelT $LAALabelT $RegionLA $UnitsRescale


#stage 7 
cp "$PROJECT/laplace_files/carpf_laplace_PV1.par" "$DATA/LA_Mesh$index/carpf_laplace_PV1.par"
cp "$PROJECT/laplace_files/carpf_laplace_PV2.par" "$DATA/LA_Mesh$index/carpf_laplace_PV2.par"
cp "$PROJECT/laplace_files/carpf_laplace_PV3.par" "$DATA/LA_Mesh$index/carpf_laplace_PV3.par"
cp "$PROJECT/laplace_files/carpf_laplace_PV4.par" "$DATA/LA_Mesh$index/carpf_laplace_PV4.par"
cp "$PROJECT/laplace_files/carpf_laplace_LAA.par" "$DATA/LA_Mesh$index/carpf_laplace_LAA.par"

cd "$DATA/LA_Mesh$index/"
docker run --rm --volume="$LAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_PV1.par -simID MV_PV1
docker run --rm --volume="$LAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_PV2.par -simID MV_PV2
docker run --rm --volume="$LAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_PV3.par -simID MV_PV3
docker run --rm --volume="$LAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_PV4.par -simID MV_PV4
docker run --rm --volume="$LAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_LAA.par -simID MV_LAA


#stage 8 
python $PROJECT/scripts/labels_la_2.py "$LAPath/" "$PROJECT/fibre_files/ra/endo/" "$PROJECT/laplace_files/" Labelled $PVLabelT $LAALabelT $RegionLA $UnitsRescale


cp "$LAPath/Labelled_Labels.elem" "$LAPath/Labelled.elem"

cp "$PROJECT/laplace_files/carpf_laplace_LS.par" "$LAPath/carpf_laplace_LS.par"
cp "$PROJECT/laplace_files/carpf_laplace_PA.par" "$LAPath/carpf_laplace_PA.par"


#stage 9
python $PROJECT/scripts/1_la.py "$LAPath/" "$PROJECT/fibre_files/la/endo/" "$PROJECT/laplace_files/" Labelled 11 13 21 23 25 27 $LandmarksLA $UnitsRescale


#stage 10
echo "=====Old UAC approximation - openCARP command line ($i)====="
cd "$DATA/LA_Mesh$index/"
docker run --rm --volume="$LAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_PA.par -simID PA_UAC_N2
docker run --rm --volume="$LAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_LS.par -simID LR_UAC_N2


#stage 11
echo "=====New UAC====="
cp "$PROJECT/laplace_files/carpf_laplace_single_LR_P.par" "$LAPath/carpf_laplace_single_LR_P.par"
cp "$PROJECT/laplace_files/carpf_laplace_single_UD_P.par" "$LAPath/carpf_laplace_single_UD_P.par"
cp "$PROJECT/laplace_files/carpf_laplace_single_LR_A.par" "$LAPath/carpf_laplace_single_LR_A.par"
cp "$PROJECT/laplace_files/carpf_laplace_single_UD_A.par" "$LAPath/carpf_laplace_single_UD_A.par"
python $PROJECT/scripts/2a_la.py "$LAPath/" "$PROJECT/fibre_files/la/endo/" "$PROJECT/laplace_files/" Labelled 11 13 21 23 25 27 $LandmarksLA $UnitsRescale


#stage 12
cd "$DATA/LA_Mesh$index/"
docker run --rm --volume="$LAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_single_LR_A.par -simID LR_Ant_UAC
docker run --rm --volume="$LAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_single_LR_P.par -simID LR_Post_UAC
docker run --rm --volume="$LAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_single_UD_A.par -simID UD_Ant_UAC
docker run --rm --volume="$LAPath":/shared:z --workdir=/shared docker.opencarp.org/opencarp/opencarp:latest openCARP +F carpf_laplace_single_UD_P.par -simID UD_Post_UAC


#stage 13
echo "=====New UAC part 2====="
python $PROJECT/scripts/2b_la.py "$LAPath/" "$PROJECT/fibre_files/la/endo/" "$PROJECT/laplace_files/" Labelled 11 13 21 23 25 27 $UnitsRescale



#stage 14
echo "=====Add fibres2 ====="
python $PROJECT/scripts/fibre_mapping.py "$LAPath/" "$PROJECT/fibre_files/la/endo/l/" "$PROJECT/laplace_files/" Labelled Labelled.lon Fibre_l


#stage 15
echo "=====Add LAT field ($i)====="
python $PROJECT/scripts/lat_field.py "$LAPath/" LAT_Spiral4_B.dat Fibre_l


