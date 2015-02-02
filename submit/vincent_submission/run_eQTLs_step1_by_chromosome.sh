

Rbin=/share/apps/R-3.0.2/bin/R

#dataset="liver_Schadt"
#conditions="Liver"
#pvOutputThreshold=5

#dataset=monocytes_Knight
#conditions="normal IFN LPS2 LPS24 LPS2logFC LPS24logFC IFNlogFC"
#pvOutputThreshold=5

#dataset="WB_Franke"
#conditions="WB"
#pvOutputThreshold=5

#dataset="brain_UKBEC"
#conditions="probesetCRBL"
#conditions="probesetCRBL probesetFCTX probesetHIPP probesetMEDU probesetOCTX probesetPUTM probesetSNIG probesetTCTX probesetTHAL probesetWHMT"
#conditions="core_AVERAGE"
#pvOutputThreshold=5

dataset="GTex"
conditions="BrainCerebellum_isoform_ratios"
pvOutputThreshold=5

#dataset="GTex"
#conditions="Muscle"
#conditions="AdiposeSubcutaneous  ArteryTibial  Brain  HeartLeftVentricle  Lung  Muscle  Nerve  SkinSunExposedLowerleg  Thyroid  WholeBlood"
#conditions="AdiposeSubcutaneous ArteryTibial BrainCaudatebasalganglia BrainCerebellarHemisphere BrainCerebellum BrainCortex BrainFrontalCortexBA9 BrainHippocampus BrainHypothalamus BrainNucleusaccumbensbasalganglia HeartLeftVentricle Lung Muscle Nerve SkinSunExposedLowerleg Thyroid transcript_AdiposeSubcutaneous WholeBlood"
#pvOutputThreshold=5



step1=TRUE
step2=TRUE
step3=FALSE
memory=1.9



if [[ "$step2" == "TRUE" ]]; then memory=3.9; fi
if [[ "$step1" == "TRUE" ]]; then memory=5.9; fi


#memory=8.9  ##step2 for WB and liver needs that sort of memory


for condition in $conditions; do
    
    script=cluster/submission/${dataset}_${condition}_eQTLs.sh
    Rscript=cluster/submission/${dataset}_${condition}_eQTLs.R
    
    echo "
source('scripts/eQTLs_scripts/find_all_eQTLs.R')
source('scripts/transeQTL_scripts/find_all_trans_eQTLs.R')
source('scripts/eQTLs_scripts/create_eQTL_summary.R')


getArgs <- function() {
  myargs.list <- strsplit(grep(\"=\",gsub(\"--\",\"\",commandArgs()),value=TRUE),\"=\")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}

### define the dataset(s) of interest
dataset <- '$dataset'
condition <- '$condition'


myArgs <- getArgs()
chromosome <- myArgs[[ 'chromosome' ]]

step1 <- $step1
step2 <- $step2
step3 <- $step3

###### cis eQTLs first

if (step1) {
my.tab <- run.eQTL (dataset, 
                     condition, 
                     chromosome = as.character(chromosome), 
                     start = 1, end = 300*10^6, 
                     pvOutputThreshold = $pvOutputThreshold, 
                     force = TRUE, min.MAF = 0.05)
}

if (step2) {
my.sum <- create.eQTL.summary (dataset, condition, min.MAF = 0.03, level = 'probe', pval.threshold = $pvOutputThreshold,
                                  base.folder = '/cluster/project8/vyp/eQTL_integration',
                                  chromosome = as.character(chromosome), plot = TRUE)
}


#### and now the trans eQTL modules
if (step3) {
  res <- find.all.trans.eQTLs (choice.sets, min.MAF = 0.05, 
                               chromosome = as.character(chromosome), min.gene.module = 6, 
                               pval.threshold = $pvOutputThreshold, 
                               with.pca = FALSE, plot = FALSE, run.stepwise = TRUE)
}


" > $Rscript
    
    echo "
#!/bin/bash
#$ -S /bin/bash
#$ -o cluster/out
#$ -e cluster/error
#$ -cwd
#$ -l tmem=${memory}G,h_vmem=${memory}G
#$ -V
#$ -R y
#$ -l h_rt=60:00:0
#$ -pe smp 1
#$ -t 1-22
#$ -tc 22

chromosome=\${SGE_TASK_ID}

$Rbin CMD BATCH --no-save --no-restore --chromosome=\${chromosome}  $Rscript cluster/R/${dataset}_${condition}_\${chromosome}_eQTLs.out

" > $script
    
    ls -ltrh $script
    qsub $script
done
