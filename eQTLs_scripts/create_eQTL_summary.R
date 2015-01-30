



create.eQTL.summary<- function (dataset, condition, min.MAF = 0.03, level = 'probe', pval.threshold = 5,
                                base.folder = '/cluster/project8/vyp/eQTL_integration',
                                chromosome, plot = TRUE) {
  
  source(paste(base.folder, '/scripts/Pickrell/create_Pickrell_input.R', sep = ''))
  source(paste(base.folder, '/scripts/plotting_functions/gviz_eqtls_v2.R', sep = ''))
  options(stringsAsFactors = FALSE)

  final <- data.frame()  ### that will be the final table as an input for fgwas
  summary.table <- data.frame()  ## a nice summary statistics table
  
########### select output folder
  oFolder <- paste(base.folder, '/data/', dataset, '/eQTLs', sep = '')
  oFolder.1.figs <- paste(oFolder, '/fgwas/fgwas_', condition, '_individual_figs', sep = '')
  oFolder.1.files <- paste(oFolder, '/fgwas/fgwas_', condition, '_individual_files', sep = '')
  oFolder.1.summary <- paste(oFolder, '/fgwas/fgwas_', condition, '_summary_by_chromosome', sep = '')
  
  for (folder in c(oFolder, oFolder.1.figs, oFolder.1.files, oFolder.1.summary)) {
    if (!file.exists(folder)) dir.create(folder)
  }



  log.file <- paste(oFolder, '/temp_chr', chromosome, '.log', sep = '')
  output.file.chrom <- paste(oFolder.1.summary, '/summary_ciseQTLs_chr', chromosome, '.csv', sep = '')

  
##########################
  expression.file <- paste(base.folder, '/data/', dataset, '/expression_data/expression_', condition, '.RData', sep = '')
  message('Loading ', expression.file)
  load(expression.file)
  expression <- get(condition)
  
############### load the genotype data
  genotype.file <- paste(base.folder, '/data/', dataset, '/genotypes/chr', chromosome, sep = '')
  load(genotype.file)

############ input eQTL 
  input.eQTL.file <- paste(oFolder, '/matrixEQTL/', condition, '/', condition, '_chr', chromosome, '_pval', pval.threshold, '.tab', sep = '')
  data <- read.table(file = input.eQTL.file, header = TRUE, stringsAsFactors = FALSE)
  data <- subset(data, cis.eQTL & MAF >= min.MAF)
  
  ###### clean up old files?
  all.pdfs <- grep(list.files( oFolder.1.figs, full.names = TRUE, pattern = condition), pattern = paste('chr', chromosome, '_', sep = ''), value = TRUE)
  all.tabs <- grep(list.files( oFolder.1.files, full.names = TRUE, pattern = condition), pattern = paste('chr', chromosome, '_', sep = ''), value = TRUE)
  message('Remove old files')
  removed <- file.remove(c(all.pdfs, all.tabs))
  message('Done removing old files')
  
  ############################### now start the loop

  if (level == 'gene') {list.genes <- unique(data$ensemblID)}
  if (level == 'probe') {list.genes <- unique(data$ProbeID)}

  nsignals <- 0
  for (loc.gene in list.genes) {  ### find each gene with a cis-eQTL
    nsignals <- nsignals + 1
    
    if (level == 'gene') loc.eQTL <- subset(data, ensemblID == loc.gene)  ### take the subset of the file with the right probe
    if (level == 'probe') loc.eQTL <- subset(data, ProbeID == loc.gene)
    
    best.row <- loc.eQTL[ which.min(loc.eQTL$p.value), ]  ##find the best SNP/probe P-value for that gene based on matrixEQTL
    loc.symbol <- best.row$Gene.name
    ProbeID <- as.character(best.row$ProbeID) ##the character bit is important for probe names that are numbers
    
    if (is.na(loc.symbol)) loc.symbol <- loc.gene ## to avoid something missing
    Pickrell.table <- create.Pickrell.input.file (dataset = dataset,
                                                  condition = condition,
                                                  genotypes = genotypes,
                                                  expression = expression,
                                                  base.folder = base.folder,
                                                  min.MAF = min.MAF,
                                                  ProbeID = ProbeID, chromosome = best.row$chromosome, snp.name = best.row$SNP)
    
    Pickrell.table$Gene.name <- loc.symbol
    Pickrell.table$ensemblID <- best.row$ensemblID
    
    min.p <- min(Pickrell.table$PVAL, na.rm = TRUE)
    #if (min.p < 10^-(pval.threshold -1)) {
    if (min.p < 10^(-2)) {
      final <- rbind.data.frame(final, Pickrell.table)
      
      my.list <- as.list(Pickrell.table [ which.min(Pickrell.table$PVAL),])  ##this is where the best P based on snpStats P-value is selected
      my.list[[ 'distance.to.gene' ]] <- min( abs(best.row$position - best.row$gene.position.start), abs(best.row$position - best.row$gene.position.end) )

      output.pdf <- paste(oFolder.1.figs, '/chr_', chromosome, '_', condition, '_', best.row$SNP, '_Probe', best.row$ProbeID, '_', best.row$Gene.name, '.pdf', sep = '')
      output.file <- paste(oFolder.1.files, '/chr', chromosome, '_', condition, '_', best.row$SNP, '_Probe', best.row$ProbeID, '_', best.row$Gene.name, '.tab', sep = '')
      my.list[[ 'output.pdf' ]] <- output.pdf
      my.list[[ 'output.file' ]] <- output.file
      
      if (plot) {
        annotated <- plot.eQTL (chromosome = chromosome, positions = Pickrell.table$POS, pvalues = Pickrell.table$PVAL, output.pdf = output.pdf, gene.list = loc.symbol,
                                gene.chromosome = best.row$gene.chromosome,
                                gene.position.start = best.row$gene.position.start, gene.position.end = best.row$gene.position.end,
                                gene.name = loc.symbol, gene.context = TRUE) ##plot a fancy graph
      }
      
      
      write.table(x = Pickrell.table, file = output.file, row.names = FALSE, sep = '\t', quote = FALSE)
      print(output.file)
      summary.table <- rbind.data.frame( summary.table, my.list)
    } else {
      cat('Problem with ', best.row$ensemblID, ' ', best.row$SNP, ' ', best.row$Gene.name, file = log.file, append = TRUE)
    }
  }
  
  message('Writing the summary in ', output.file.chrom)
  write.csv(x = summary.table,  file = output.file.chrom, row.names = FALSE)
  return(summary.table)
}
