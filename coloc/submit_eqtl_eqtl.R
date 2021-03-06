

################
coloc.eqtl.eqtl <- function(dataset1,
                            dataset2,
                            basedir,
                            outdir = '/cluster/project8/vyp/eQTL_integration/coloc',
                            cluster.folder) {

  path_to_Rscript = '/share/apps/R-3.0.1/bin/Rscript'

  eqtl.dataset1 <- names(dataset1)
  eqtl.dataset2 <- names(dataset2)
  cond.1 <- as.character(dataset1[[1]])
  cond.2 <- as.character(dataset2[[1]])

  message('Dataset 1 ', eqtl.dataset1, ' and condition 1 ', cond.1)
  message('Dataset 2 ', eqtl.dataset2, ' and condition 2 ', cond.2)
  
  
  # Create cluster folders: 
  cluster.out = paste(cluster.folder, '/out', sep='')
  if (!file.exists (cluster.out)) dir.create(cluster.out)
  cluster.error = paste(cluster.folder, '/error', sep='')
  if (!file.exists (cluster.error)) dir.create(cluster.error)
  cluster.submission = paste(cluster.folder, '/submission', sep='')
  if (!file.exists (cluster.submission)) dir.create(cluster.submission)
  
  main_script = paste(basedir, "/coloc/main_eqtl_eqtl.R", sep ='')
  outfld.base <-  paste(outdir, "/", eqtl.dataset1, "_", eqtl.dataset2, sep = '')
  outfld = paste(outdir, "/", eqtl.dataset1, "_", eqtl.dataset2, "/", cond.1, "_", cond.2, sep = '') 

  if (!file.exists(outfld.base)) dir.create( outfld.base )
  if (!file.exists(outfld)) dir.create( outfld )
  
  script_to_submit = paste(cluster.submission, '/Submit_main_script_', eqtl.dataset1, '_', eqtl.dataset2, '_', cond.1, '_', cond.2, '.sh', sep='') 
  
  write(file=script_to_submit, '#!/bin/bash -l', append=F)
  write(file=script_to_submit, '#$ -S /bin/bash', append=T)
  
  write(file=script_to_submit, paste(path_to_Rscript, main_script, eqtl.dataset1, eqtl.dataset2, 
          cond.1, cond.2, outfld, basedir, sep=' '), append=TRUE)

  message('R script file written in ', main_script)
  message('Submission script: ', script_to_submit)

  if (!file.exists(main_script)) {stop('The main R submission script does not exist')}
  if (!file.exists(script_to_submit)) {stop('The main bash submission script does not exist')}
  
  system(paste('qsub -cwd -l  h_vmem=8G -l tmem=8G -l h_rt=20:0:0 -o ', cluster.out, ' -e ',  cluster.error, ' ',  script_to_submit, sep=''))  
}




#################
# Variables to set

##kitty's version
#eqtl.dataset1 <- 'monocytes_Knight' 
#eqtl.dataset2 <- 'Smith_macrophages'  
#cond.1="LPS24logFC"
#cond.2="logFC"
#basedir = "/cluster/project8/vyp/kitty/eQTL_scripts/" 
#outdir = "/scratch2/vyp-scratch2/kitty/eQTLs_integration/coloc" 
#cluster.folder = '/scratch2/vyp-scratch2/kitty/eQTLs_integration/coloc/'

##vincent's version
eqtl.dataset1 <- 'LCL_dexamethasone_DiRienzo' 
eqtl.dataset2 <- 'WB_dexamethasone_DiRienzo'  
cond.1="logFC"
cond.2="logFC"
basedir = "/cluster/project8/vyp/eQTL_integration/scripts"
outdir = "/cluster/project8/vyp/eQTL_integration/coloc"
cluster.folder = '/cluster/project8/vyp/eQTL_integration/cluster'



dataset1 <- list ( LCL_dexamethasone_DiRienzo = 'logFC')
dataset2 <- list ( WB_dexamethasone_DiRienzo = 'logFC')
res <- coloc.eqtl.eqtl (dataset1, dataset2, basedir, outdir = '/cluster/project8/vyp/eQTL_integration/coloc', cluster.folder)
