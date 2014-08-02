##############################################
## Variables to set 
##############################################
biom.dataset <- 'UCLEB' # 'SCZ_Ripke'
eqtl.dataset <- 'liver_Schadt'  # 'brain_UKBEC' # 'GTex'
biom.names="ucleb"
# If UCLEB, leave biom.names as "ucleb" and load later
## condition <- 'AdiposeSubcutaneous'
#cond.all = c("HeartLeftVentricle","AdiposeSubcutaneous","Thyroid","SkinSunExposedLowerleg","Nerve","Muscle","ArteryTibial","WholeBlood")
cond.all="Liver"
type="quant"

########### SCZ and UKBEC
## biom.dataset <- 'SCZ_Ripke'
## eqtl.dataset <- 'brain_UKBEC' 
## biom.names="scz"
## cond.all = c("CRBL","FCTX","HIPP","MEDU","OCTX","PUTM","SNIG","TCTX","THAL","WHMT")
## type = "cc"

#################
## Paths for submission of scripts:
## cluster.folder = '/scratch2/vyp-scratch2/claudia/GTEx/cluster/'
#For now:
cluster.folder = '/cluster/project8/vyp/claudia/cluster/'
path_to_Rscript = '/share/apps/R-3.0.1/bin/Rscript'
#Rbin=/share/apps/R-3.0.2/bin/R
# Create cluster folders: 
cluster.out = paste(cluster.folder, 'out/', sep='')
   if (!file.exists (cluster.out)) dir.create(cluster.out)
cluster.error = paste(cluster.folder, 'error/', sep='')
   if (!file.exists (cluster.error)) dir.create(cluster.error)
cluster.submission = paste(cluster.folder, 'submission/', sep='')
   if (!file.exists (cluster.submission)) dir.create(cluster.submission)

#queue= 'queue14' # Submission queue to use: check which queu to use with: qstat -u "*"

#################
# Base directory of the github repository 
scripts.folder = "/cluster/project8/vyp/claudia/scripts/eQTL_scripts/submit/claudia_submission/"
#scripts.folder = "/cluster/project8/vyp/eQTL_integration/scripts/coloc/"
main_script = paste(scripts.folder, "main_eqtl_biom.R", sep="")

for (condition in cond.all) {
for (chr.name in 1:22) {


script_to_submit = paste(cluster.submission, 'Submit_main_script_', biom.dataset, '_', eqtl.dataset, '_', condition, '_', chr.name, '.sh', sep='') # Root name of scripts that will be submitted for each chr, each biomarker, eqtl

    write(file=script_to_submit, '#!/bin/bash -l', append=F)
    write(file=script_to_submit, '#$ -S /bin/bash', append=T)
    write(file=script_to_submit, paste(path_to_Rscript, main_script, biom.dataset, eqtl.dataset, biom.names, condition, chr.name, type, sep=' '), append=T)

#$Rbin CMD BATCH --no-save --no-restore --scripts.folder=/cluster/project8/vyp/eQTL_integration/scripts  --biom.dataset=$biom.dataset --biom.names=$biom.names --eqtl.dataset=$eqtl.dataset --condition=$condition --chr.name=$chr.name --type=$type $script

    system(paste('qsub -cwd -l  h_vmem=8G -l tmem=8G -l h_rt=20:0:0 -o ', cluster.out, ' -e ', cluster.error, ' ',  script_to_submit, sep=''))


} # END chr
} # END condition
