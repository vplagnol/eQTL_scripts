#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_vmem=8G
#$ -l tmem=8G
#$ -l h_rt=12:00:00
#$ -o /cluster/project8/vyp/kitty/cluster/out/ 
#$ -e /cluster/project8/vyp/kitty/cluster/error/
#$ -j y

PATH=$PATH:\/share\/apps\/R-3.0.1\/bin\/
Rscript main_eqtl_eqtl.R 

