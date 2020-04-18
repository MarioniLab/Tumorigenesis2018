#!/bin/bash
hd=$1
wd=$2
singularity exec -c -H $hd ../Container/Tumorigenesis.sif Rscript -e "rmarkdown::render('$wd/01_MakeCountMatrix.Rmd')" &&\
singularity exec -c -H $hd ../Container/Tumorigenesis.sif Rscript -e "rmarkdown::render('$wd/03_QCReport.Rmd')" &&\
singularity exec -c -H $hd ../Container/Tumorigenesis.sif runScript.sh $wd 04_Normalization.R &&\
singularity exec -c -H $hd ../Container/Tumorigenesis.sif runScript.sh $wd 05_computeDoubletScores.R &&\
singularity exec -c -H $hd ../Container/Tumorigenesis.sif runScript.sh $wd 06_fastMNNforDoublets.R &&\
singularity exec -c -H $hd ../Container/Tumorigenesis.sif Rscript -e "rmarkdown::render('07_DoubletID.Rmd')" &&\
singularity exec -c -H $hd ../Container/Tumorigenesis.sif runScript.sh $wd 08_BatchCorrection.R &&\
singularity exec -c -H $hd ../Container/Tumorigenesis.sif runScript.sh $wd 09_computeUMAP.R &&\
singularity exec -c -H $hd ../Container/Tumorigenesis.sif runScript.sh $wd 10_Clustering.R &&\
singularity exec -c -H $hd ../Container/Tumorigenesis.sif Rscript -e "rmarkdown::render('$wd/11_CellTypeInference.Rmd')" &&\
singularity exec -c -H $hd ../Container/Tumorigenesis.sif Rscript -e "rmarkdown::render('$wd/12_TumorTime.Rmd')"
