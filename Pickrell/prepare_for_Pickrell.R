



prepare.Pickrell.set <- function ( dataset, condition, min.MAF = 0.03, level = 'probe', base.folder = '/cluster/project8/vyp/eQTL_integration') {

  options(stringsAsFactors = FALSE)


  oFolder <- paste(base.folder, '/data/', dataset, '/eQTLs', sep = '')
  oFolder.1.summary <- paste(oFolder, '/fgwas/fgwas_', condition, '_summary_by_chromosome', sep = '')
  summary.file <- paste(oFolder, '/fgwas/fgwas_', condition, '_summary_eQTLs.csv', sep = '')
  input.fgwas <- paste(oFolder, '/fgwas/fgwas_', condition, '_fine_mapping_input.tab', sep = '')

  #############
  single.chromosome.files <- paste(oFolder.1.summary, '/summary_ciseQTLs_chr', 1:22, '.csv', sep = '')
  if (sum(!file.exists ( single.chromosome.files )) > 0) {
    stop('Some single chromosome files are missing')
  }

  ######## some clean up
  for (file in c(summary.file, input.fgwas)) {
    if (file.exists(file)) file.remove(file)
  }
  

  
  final.table <- data.frame()
  fgwas.table <- data.frame()

  segnumber <- 0
  for (i in 1:22) {
    message('Now reading ', single.chromosome.files[i])
    chrom.table <- read.csv( single.chromosome.files[i],stringsAsFactors = FALSE)
    final.table <- rbind.data.frame (final.table, chrom.table)
                   
    for (j in 1:nrow(chrom.table)) {
      segnumber <- segnumber + 1
      input.file <- chrom.table$output.file[ j ]
      message('Input file : ', input.file)

      #final[, 'SEGNUMBER'] <- as.numeric(factor(final$ProbeID, levels = unique(final$ProbeID)))
    }
    
  }

  message('Writing the summary in ', summary.file)
  write.csv(x = final.table,  file = summary.file, row.names = FALSE)
  stop()
  
###### first input file without the TF annotation


  write.table(x = final, file = input.fgwas, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')
  print(input.fgwas)
  
####### Now annotate the input file with the knowledge of the TF binding sites
  input.fgwas.with.anno <- paste(oFolder, '/fgwas/fgwas_', condition, '_fine_mapping_input_with_annotation.tab', sep = '')
  list.bin.annot <- paste(oFolder, '/fgwas/fgwas_', condition, '_list_bin_annot.tab', sep = '')
  
  my.cmd <- paste('R CMD BATCH --no-save --no-restore  --inputFrame=', input.fgwas, ' --outputFrame=',input.fgwas.with.anno,
                  ' --listBinAnnot=', list.bin.annot,
                  ' --ENCODE=FALSE --TF=TRUE /cluster/project8/vyp/vincent/Software/pipeline/fgwas/fgwas_step1.R cluster/R/fgwas_step1.out', sep = '')
  
  submission.file <- paste('cluster/submission/fgwas_', my.sets[1], '.sh', sep = '')
  cat('#!/bin/bash
#$ -S /bin/bash
#$ -o cluster/out
#$ -e cluster/error
#$ -cwd
#$ -l tmem=15G,h_vmem=15G
#$ -l h_rt=60:0:0
#$ -V
#$ -R y
', my.cmd, file = submission.file)
  
  #system(paste("qsub ", submission.file, sep = ''))
  system(paste("sh ", submission.file, sep = ''))
  return (summary.table)
}



#choice.sets <- list( LCL_dexamethasone_DiRienzo = c('logFC'), WB_dexamethasone_DiRienzo = c('logFC'))  ###key argument
#choice.sets <- list(WB_dexamethasone_DiRienzo = c('logFC'))
#choice.sets <- list(LCL_dexamethasone_DiRienzo = c('logFC'))
#choice.sets <- list(brain_UKBEC = 'lncRNA_cerebellar_cortex')
#choice.sets <- list(liver_Schadt = 'Liver')

test <- prepare.Pickrell.set(dataset = 'liver_Schadt', condition = 'Liver', choice.sets)
