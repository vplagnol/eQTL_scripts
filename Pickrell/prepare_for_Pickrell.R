



prepare.Pickrell.set <- function ( choice.sets, min.MAF = 0.03, level = 'probe' ) {
  base.folder <- '/cluster/project8/vyp/eQTL_integration'
  
  source(paste(base.folder, '/scripts/Pickrell/create_Pickrell_input.R', sep = ''))
  source(paste(base.folder, '/scripts/plotting_functions/gviz_eqtls_v2.R', sep = ''))
  options(stringsAsFactors = FALSE)

  dataset <- names(choice.sets)[1]
  condition <- as.character(unlist(choice.sets))
  final <- data.frame()  ### that will be the final table as an input for fgwas
  summary.table <- data.frame()  ## a nice summary statistics table
  
########### select output folder
  my.sets <- sort(names(choice.sets))
  set.key <-  paste(my.sets, collapse = '_')
  if (length(choice.sets) == 1) {
    oFolder <- paste(base.folder, '/data/', set.key, '/eQTLs', sep = '')
    if (!file.exists(oFolder)) dir.create(oFolder)
  } else {
    oFolder <- paste(base.folder, '/combined_data/', set.key, '/eQTLs', sep = '')
  }

  oFolder.1.figs <- paste(oFolder, '/fgwas/fgwas_', condition, '_individual_figs', sep = '')
  oFolder.1.files <- paste(oFolder, '/fgwas/fgwas_', condition, '_individual_files', sep = '')
  summary.file <- paste(oFolder, '/fgwas/fgwas_', condition, '_summary_eQTLs.csv', sep = '')
  
  for (folder in c(oFolder, oFolder.1.figs, oFolder.1.files)) {
    if (!file.exists(folder)) dir.create(folder)
  }
  
  all.pdfs <- list.files( oFolder.1.figs, full.names = TRUE, pattern = condition)
  all.tabs <- list.files( oFolder.1.files, full.names = TRUE, pattern = condition)

##########################
  expression.file <- paste(base.folder, '/data/', dataset, '/expression_data/expression_', condition, '.RData', sep = '')
  message('Loading ', expression.file)
  load(expression.file)
  expression <- get(condition)
  
################
  message('Remove old files')
  file.remove(c(all.pdfs, all.tabs))
  
  log.file <- paste(oFolder, '/temp.log', sep = '')

  nsignals <- 0
  for (chr in 22:1) {

############### load the genotype data
    genotype.file <- paste(base.folder, '/data/', names(choice.sets)[[1]], '/genotypes/chr', chr, sep = '')
    load(genotype.file)

############ input eQTL 
    input.eQTL.file <- paste(oFolder, '/matrixEQTL/', condition, '/', condition, '_chr', chr, '_pval5.tab', sep = '')
    data <- read.table(file = input.eQTL.file, header = TRUE, stringsAsFactors = FALSE)
    data <- subset(data, cis.eQTL & MAF >= min.MAF)
    
    if (level == 'gene') {list.genes <- unique(data$ensemblID)}
    if (level == 'probe') {list.genes <- unique(data$ProbeID)}
    

    for (loc.gene in list.genes) {  ### find each gene with a cis-eQTL
      nsignals <- nsignals + 1
      
      if (level == 'gene') loc.eQTL <- subset(data, ensemblID == loc.gene)  ### take the subset of the file with the right probe
      if (level == 'probe') loc.eQTL <- subset(data, ProbeID == loc.gene)
      
      best.row <- loc.eQTL[ which.min(loc.eQTL$p.value), ]  ##find the best SNP/probe P-value for that gene
      loc.symbol <- best.row$Gene.name
      ProbeID <- as.character(best.row$ProbeID) ##the character bit is important for probe names that are numbers

      if (is.na(loc.symbol)) loc.symbol <- loc.gene ## to avoid something missing
      Pickrell.table <- create.Pickrell.input.file (choice.sets = choice.sets,
                                                    genotypes = genotypes, expression = expression,
                                                    min.MAF = min.MAF,
                                                    ProbeID = ProbeID, chromosome = best.row$chromosome, snp.name = best.row$SNP)
      
      Pickrell.table$Gene.name <- loc.symbol
      Pickrell.table$ensemblID <- best.row$ensemblID
      
      min.p <- min(Pickrell.table$PVAL, na.rm = TRUE)
      if (min.p < 10^-4) {
        final <- rbind.data.frame(final, Pickrell.table)

        my.list <- as.list(Pickrell.table [ which.min(Pickrell.table$PVAL),])
        my.list[[ 'distance.to.gene' ]] <- min( abs(best.row$position - best.row$gene.position.start), abs(best.row$position - best.row$gene.position.end) )
        
        
        
        output.pdf <- paste(oFolder.1.figs, '/fgwas_', condition, '_', best.row$SNP, '_Probe', best.row$ProbeID, '_', best.row$Gene.name, '.pdf', sep = '')
        annotated <- plot.eQTL (chromosome = chr, positions = Pickrell.table$POS, pvalues = Pickrell.table$PVAL, output.pdf = output.pdf, gene.list = loc.symbol,
                                gene.chromosome = best.row$gene.chromosome,
                                gene.position.start = best.row$gene.position.start, gene.position.end = best.row$gene.position.end,
                                gene.name = loc.symbol, gene.context = TRUE) ##plot a fancy graph
        
        output.file <- paste(oFolder.1.files, '/fgwas_', condition, '_', best.row$SNP, '_Probe', best.row$ProbeID, '_', best.row$Gene.name, '.tab', sep = '')


        my.list[[ 'output.pdf' ]] <- output.pdf
        my.list[[ 'output.file' ]] <- output.file

        
        write.table(x = Pickrell.table, file = output.file, row.names = FALSE, sep = '\t', quote = FALSE)
        print(output.file)
        summary.table <- rbind.data.frame( summary.table, my.list)
        message('Writing the summary in ', summary.file)
        write.csv(x = summary.table,  file = summary.file, row.names = FALSE)
        
      } else {
        cat('Problem with ', best.row$ensemblID, ' ', best.row$SNP, ' ', best.row$Gene.name, file = log.file, append = TRUE)
        
      }
    }

    message('Writing the summary in ', summary.file)
    write.csv(x = summary.table,  file = summary.file, row.names = FALSE)
  }
  
###### first input file without the TF annotation
  input.fgwas <- paste(oFolder, '/fgwas/fgwas_', condition, '_fine_mapping_input.tab', sep = '')
  final[, 'SEGNUMBER'] <- as.numeric(factor(final$ProbeID, levels = unique(final$ProbeID)))
  write.table(x = final, file = input.fgwas, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')
  print(input.fgwas)
  
####### Now annotate the input file with the knowledge of the TF binding sites
  input.fgwas.with.anno <- paste(oFolder, '/fgwas/fgwas_', condition, '_fine_mapping_input_with_annotation.tab', sep = '')
  list.bin.annot <- paste(oFolder, '/fgwas/fgwas_', condition, '_list_bin_annot.tab', sep = '')
  
  my.cmd <- paste('R CMD BATCH --no-save --no-restore  --inputFrame=', input.fgwas, ' --outputFrame=',input.fgwas.with.anno,
                  ' --listBinAnnot=', list.bin.annot,
                  ' --ENCODE=FALSE --TF=TRUE /cluster/project8/vyp/vincent/toolsVarious/fgwas/fgwas_step1.R cluster/R/fgwas_step1.out', sep = '')
  
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
', my.cmd, '\ngzip -f ',  input.fgwas.with.anno, file = submission.file)
  
  #system(paste("qsub ", submission.file, sep = ''))
  system(paste("sh ", submission.file, sep = ''))
  return (summary.table)
}



#choice.sets <- list( LCL_dexamethasone_DiRienzo = c('logFC'), WB_dexamethasone_DiRienzo = c('logFC'))  ###key argument
#choice.sets <- list(WB_dexamethasone_DiRienzo = c('logFC'))
#choice.sets <- list(LCL_dexamethasone_DiRienzo = c('logFC'))
#choice.sets <- list(brain_UKBEC = 'lncRNA_cerebellar_cortex')
#choice.sets <- list(monocytes_Knight= 'LPS2')

#test <- prepare.Pickrell.set(choice.sets)
