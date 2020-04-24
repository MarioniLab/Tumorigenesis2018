#!/bin/bash
# 1 Argument home directory of singularity instance: abosulte path to repository (must contain the singularity image
# 2 Argument working directory of singularity instance: abosulte path to /src/Tumorigenesis
hd=$1
wd=$2
singularity exec -c -H $hd ../Container/Tumorigenesis.sif Rscript -e "rmarkdown::render('$wd/01_MakeCountMatrix.Rmd')" &&\
singularity exec -c -H $hd ../Container/Tumorigenesis.sif Rscript -e "rmarkdown::render('$wd/03_QCReport.Rmd')" &&\
singularity exec -c -H $hd ../Container/Tumorigenesis.sif $hd/repo/src/runScript.sh $wd 04_Normalization.R &&\
singularity exec -c -H $hd ../Container/Tumorigenesis.sif $hd/repo/src/runScript.sh $wd 05_computeDoubletScores.R &&\
singularity exec -c -H $hd ../Container/Tumorigenesis.sif $hd/repo/src/runScript.sh $wd 06_fastMNNforDoublets.R &&\
singularity exec -c -H $hd ../Container/Tumorigenesis.sif Rscript -e "rmarkdown::render('$wd/07_DoubletID.Rmd')" &&\
singularity exec -c -H $hd ../Container/Tumorigenesis.sif $hd/repo/src/runScript.sh $wd 08_BatchCorrection.R &&\
singularity exec -c -H $hd ../Container/Tumorigenesis.sif $hd/repo/src/runScript.sh $wd 09_computeUMAP.R &&\
singularity exec -c -H $hd ../Container/Tumorigenesis.sif $hd/repo/src/runScript.sh $wd 10_Clustering.R &&\
singularity exec -c -H $hd ../Container/Tumorigenesis.sif Rscript -e "rmarkdown::render('$wd/11_CellTypeInference.Rmd')" &&\
singularity exec -c -H $hd ../Container/Tumorigenesis.sif Rscript -e "rmarkdown::render('$wd/12_TumorTime.Rmd')" 
