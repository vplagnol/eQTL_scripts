#################
## Paths for submission of scripts:
cluster.folder = '/scratch2/vyp-scratch2/claudia/GTEx/cluster/'
path_to_Rscript = '/share/apps/R-3.0.1/bin/Rscript'
# Create cluster folders: 
cluster.out = paste(cluster.folder, 'out/', sep='')
   if (!file.exists (cluster.out)) dir.create(cluster.out)
cluster.error = paste(cluster.folder, 'error/', sep='')
   if (!file.exists (cluster.error)) dir.create(cluster.error)
cluster.submission = paste(cluster.folder, 'submission/', sep='')
   if (!file.exists (cluster.submission)) dir.create(cluster.submission)

#queue= 'queue14' # Submission queue to use: check which queu to use with: qstat -u "*"

#################
# Variables to set 
biom.dataset <- 'SCZ_Ripke'
eqtl.dataset <- 'brain_UKBEC' 
# 13,833 schizophrenia cases and 18,310 controls
## condition <- 'AdiposeSubcutaneous'
cond.all = c("CRBL","FCTX","HIPP","MEDU","OCTX","PUTM","SNIG","TCTX","THAL","WHMT")
type = "cc"



main_script = "/cluster/project8/vyp/eQTL_integration/scripts/coloc/main_eqtl_biom_scz.R"


for (condition in cond.all) {
for (chr.name in 1:22) {


script_to_submit = paste(cluster.submission, 'Submit_main_script_', biom.dataset, '_', eqtl.dataset, '_', condition, '_', chr.name, '.sh', sep='') # Root name of scripts that will be submitted for each chr, each biomarker, eqtl

    write(file=script_to_submit, '#!/bin/bash -l', append=F)
    write(file=script_to_submit, '#$ -S /bin/bash', append=T)
    write(file=script_to_submit, paste(path_to_Rscript, main_script, biom.dataset, eqtl.dataset, condition, chr.name, type, sep=' '), append=T)

    system(paste('qsub -cwd -l  h_vmem=8G -l tmem=8G -l h_rt=20:0:0 -o ', cluster.out, ' -e ', cluster.error, ' ',  script_to_submit, sep=''))


} # END chr
} # END condition                                                                                                                                                                 

