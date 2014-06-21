

Rbin=/share/apps/R-3.0.2/bin/R


dataset="monocytes_Knight"
conditions="LPS24logFC LPS24 LPS2 normal IFN"
#conditions="LPS24logFC"

for condition in $conditions; do

    script=cluster/submission/${dataset}_${condition}_eQTLs.sh
    Rscript=cluster/submission/${dataset}_${condition}_eQTLs.R
    
    echo "
source('scripts/eQTLs_scripts/find_all_eQTLs.R')
### define the dataset(s) of interest
dataset <- '$dataset'
condition <- '$condition'

### define the slice we are looking at

for (chromosome in as.character(seq(22, 1))) {
  my.tab <- run.eQTL ( dataset, condition, chromosome, start = 1, end = 300*10^6, pvOutputThreshold = 1e-5, force = TRUE)
}

source('scripts/Pickrell/prepare_for_Pickrell.R')
choice.sets <- list()
choice.sets[[ dataset ]] <- condition
test <- prepare.Pickrell.set(choice.sets)

source('scripts/transeQTL_scripts/find_all_trans_eQTLs.R')

res <- find.all.trans.eQTLs ( choice.sets, min.MAF = 0.03, min.gene.module = 6, pval.threshold = 5, with.pca = FALSE, plot = FALSE, run.stepwise = TRUE)

" > $Rscript
    
    echo "
#!/bin/bash
#$ -S /bin/bash
#$ -o cluster/out
#$ -e cluster/error
#$ -cwd
#$ -l tmem=10G,h_vmem=10G
#$ -V
#$ -R y
#$ -l h_rt=36:00:00

$Rbin CMD BATCH --no-save --no-restore  $Rscript cluster/R/${dataset}_${condition}_eQTLs.out



" > $script
    
    ls -ltrh $script
    qsub $script
done