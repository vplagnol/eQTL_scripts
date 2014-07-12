#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_vmem=8G
#$ -l tmem=8G
#$ -l h_rt=12:00:00
#$ -o /cluster/project8/vyp/eQTL_integration/cluster/out/
#$ -e /cluster/project8/vyp/eQTL_integration/cluster/error/
#$ -j y
#$ -cwd 


dataset1="monocytes_Knight"
dataset2="Smith_macrophages"
cond1="LPS24logFC"
cond2="logFC"
matchBy="gene"
dataFolder="/cluster/project8/vyp/eQTL_integration/data"

scriptsFolder=
resultsFolder=

Rbin=/share/apps/R-3.0.2/bin/R
script=/cluster/project8/vyp/eQTL_integration/scripts/coloc/main_eqtl_eqtl.R

$Rbin CMD BATCH --no-save --no-restore --scripts.folder=/cluster/project8/vyp/eQTL_integration/scripts  --dataset1=$dataset1 --condition1=$cond1 --dataset2=$dataset2 --condition2=$cond2 --match.by=$matchBy $script /cluster/project8/vyp/kitty/cluster/R/eqtl_${dataset1}_${dataset2}_${cond1}_${cond2}.Rout




