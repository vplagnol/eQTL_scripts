

Rbin=/share/apps/R-3.0.2/bin/R

dataset="liver_Schadt"
conditions="Liver"

for condition in $conditions; do
    for chromosome in `seq 1 22`; do
	
	script=cluster/submission/${dataset}_${condition}_${chromosome}_eQTLs.sh
	Rscript=cluster/submission/${dataset}_${condition}_${chromosome}_eQTLs.R
	
	echo "
source('scripts/eQTLs_scripts/find_all_eQTLs.R')
source('scripts/transeQTL_scripts/find_all_trans_eQTLs.R')
source('scripts/eQTLs_scripts/create_eQTL_summary.R')


### define the dataset(s) of interest
dataset <- '$dataset'
condition <- '$condition'

###### cis eQTLs first
#my.tab <- run.eQTL ( dataset, condition, chromosome = as.character($chromosome), start = 1, end = 300*10^6, pvOutputThreshold = 5, force = TRUE, min.MAF = 0.05)

my.sum <- create.eQTL.summary (dataset, condition, min.MAF = 0.03, level = 'probe', pval.threshold = 5,
                                  base.folder = '/cluster/project8/vyp/eQTL_integration',
                                  chromosome = as.character($chromosome), plot = FALSE)



#### and now the trans eQTL modules
#res <- find.all.trans.eQTLs ( choice.sets, min.MAF = 0.05, chromosome = as.character($chromosome), min.gene.module = 6, pval.threshold = 5, with.pca = FALSE, plot = FALSE, run.stepwise = TRUE)

" > $Rscript
	
	echo "
#!/bin/bash
#$ -S /bin/bash
#$ -o cluster/out
#$ -e cluster/error
#$ -cwd
#$ -l tmem=23G,h_vmem=23G
#$ -V
#$ -R y
#$ -l h_rt=36:00:00

$Rbin CMD BATCH --no-save --no-restore  $Rscript cluster/R/${dataset}_${condition}_${chromosome}_eQTLs.out
" > $script
	
	ls -ltrh $script
	qsub $script
    done
done