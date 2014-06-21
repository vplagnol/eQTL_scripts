library(snpStats)
source('/cluster/project8/vyp/eQTL_integration/scripts/plotting_functions/plot_eQTLs_per_SNP.R')
options(stringsAsFactors = FALSE)

plot <- FALSE




find.all.trans.eQTLs <- function( choice.sets, chromosome.list, min.MAF = 0.03, min.gene.module = 6, pval.threshold = 5, with.pca = FALSE, plot = FALSE, run.stepwise = FALSE) {
  
  base.folder <- '/cluster/project8/vyp/eQTL_integration'

  
  if (run.stepwise) {
    script.file <- paste(base.folder, '/scripts/transeQTL_scripts/tools.R', sep= '')   ## this has the stepwise analysis function, which I think is helpful
    source(script.file)
  }
  
  dataset <- names(choice.sets)[[ 1 ]]
  condition <- as.character(choice.sets[1])

  message('Dataset: ', dataset, ' and condition: ', condition)
  
########## prepare the PCA analysis
  expression.file <- paste(base.folder, '/data/', dataset, '/expression_data/expression_', condition, '.RData', sep = '')
  load(expression.file)
  expression <- get(condition)
  
  my.res.trans <- list()
  my.res.cis <- list()
  
###load the expression file
  main.folder <- paste(base.folder, '/data/', dataset, '/eQTLs')
  oFolder <- paste(base.folder, '/data/', dataset, '/modules/', condition, sep = '')
  fig.folder <- paste(oFolder, '/fig', sep = '')
  
  for (folder in c(oFolder, fig.folder)) {
    message('Output folder ', folder)
    if (!file.exists(folder)) dir.create(folder)
  }
  
### do some clean up of old files
  if (length(chromosome.list) == 22) {
    finalFile <- paste(base.folder, '/data/', dataset, '/modules/modules_', condition, '_all_modules.RData', sep = '')
    tableFile <-  paste(base.folder, '/data/', dataset, '/modules/modules_', condition, '_main_table.tab', sep = '')
    for (file in c(finalFile, tableFile)) file.remove(file)
  }

  table.final <- data.frame()
  modules.final <- list()
  
  for (chromosome in chromosome.list) {
    message('Chromosome ', chromosome)

    chromFile <- paste(base.folder, '/data/', dataset, '/modules/modules_', condition, '_chromosome', chromosome, '_all_modules.RData', sep = '')
    tableChrom <-  paste(base.folder, '/data/', dataset, '/modules/modules_', condition, '_chromosome', chromosome, '_main_table.tab', sep = '')
    for (file in c(chromFile, tableChrom)) file.remove(file)

########
    list.modules <- list()
    table.chrom <- list()
    n.modules <- 0
    
    first <- TRUE
   
###### load the genotype file
    geno.file <- paste(base.folder, '/data/', dataset, '/genotypes/chr', chromosome, sep = '')
    load(geno.file)
    
####### take samples in common only between expression and genotypes
    if (first) {
      shared.samples <- intersect(dimnames(genotypes$genotypes)[[1]], dimnames(expression)[[2]])
      expression <- expression[, shared.samples]
      
      if (with.pca)  my.pca <- prcomp(t(expression), center = TRUE)
      first <- FALSE
    }
    genotypes$genotypes <- genotypes$genotypes[ shared.samples, ]
    
#######
    input.file <- paste(base.folder, '/data/', dataset, '/eQTLs/matrixEQTL/', condition, '/', condition, '_chr', chromosome, '_pval', pval.threshold, '.tab', sep = '')
    data <- read.table(input.file, sep = '\t', header = TRUE)
    data$dataset <- dataset
    data <- subset(data, !is.na(ensemblID))
    
############# strategy 1: use the idea of auxiliary genes
    strong.candidates <- subset( data, p.value < 10^-12 & MAF > min.MAF & !cis.eQTL, 'SNP', drop = TRUE)
    
#### strategy 2: look for SNPs with trans effects on many genes
    my.tab <- table(subset(data$SNP, data$MAF > min.MAF & ! data$cis.eQTL))
    good.trans <- names(which(my.tab  >= min.gene.module))
  
    strong.candidates.combined <- subset(data, SNP %in% c(good.trans, strong.candidates))
    
    my.res.trans[[ chromosome ]] <- strong.candidates.combined
    my.res.cis[[ chromosome ]] <- subset(data, cis.eQTL & SNP %in% strong.candidates.combined$SNP)
  
###### Now define the modules below
    loc.table <- my.res.trans[[ chromosome ]]
    my.tab <- table(c('dummy', loc.table$SNP))
    
    while ((nrow(loc.table) > 0) && (max(my.tab) >= min.gene.module)) {
      n.modules <- n.modules + 1
      message('Module number ', n.modules)
      my.sort <- sort(my.tab, decreasing = TRUE)
      best.SNP <- names(my.sort)[[1]]
      best.module <- subset( loc.table, SNP == best.SNP)
      my.position <- best.module$position[ 1 ]
      list.modules[[ n.modules ]] <- best.module
      names(list.modules)[ n.modules ] <- best.SNP
      
      my.res <-  list (chromosome = chromosome,
                       position = my.position,
                       best.SNP = best.SNP,
                       n.genes = length(unique(best.module$ensemblID)),
                       n.cis.eQTLs =  sum(best.module$cis.eQTL, na.rm = TRUE),
                       cis.eQTL.genes = ifelse ( sum(best.module$cis.eQTL, na.rm = TRUE) > 0, paste( subset(best.module, cis.eQTL, 'Gene.name', drop = TRUE), collapse = ','), NA),
                       trans.eQTL.genes = ifelse ( sum(! best.module$cis.eQTL, na.rm = TRUE) > 0, paste( subset(best.module, !cis.eQTL, 'Gene.name', drop = TRUE), collapse = ','), NA),
                       module.number = n.modules)
      
############ Now the stepwise regression analysis
      if (run.stepwise) {
        message('Running the stepwise analysis')
        print(head(best.module))
                                        #browser()
        my.stepwise <- stepwise.analysis (expression.data = expression[ best.module$ProbeID, ] ,  ### need to add covariates
                                          genotype.data = genotypes$genotypes,
                                          snp.name = best.SNP,
                                          gene.names = best.module$Gene.name)
        my.res[['n.independent.signals']] <- my.stepwise$nsteps
      }
      
    
########## Now the PCA analysis
      if (with.pca) {
        geno.num <- as.numeric(genotypes$genotypes[, best.SNP])
        for (i in 1:5) {
          my.res[[  paste('PCA', i, sep = '') ]] <- cor.test( geno.num, my.pca$x[, i])$p.value
        }
      }
      
###########
      tableChrom <- rbind.data.frame(tableChrom, my.res)
                                        #print(head(table.final))
    
################
      loc.table <- subset(loc.table, ! ensemblID %in% best.module$ensemblID) #### remove the genes that have been accounted for already
    
      my.tab <- table(c('dummy', loc.table$SNP))
    
      if (plot) {
        output.pdf <- paste(fig.folder, '/', condition, '_', best.module$SNP[1], '_chr', best.module$chromosome[1], '_',  best.module$position[1], '.pdf', sep = '')
        test <- plot.eQTL (choice.sets, snp.name = best.module$SNP[1], chromosome = best.module$chromosome[1], gene.names = best.module$Gene.name, output.pdf = output.pdf)
        message('Output pdf for this module in ', output.pdf)
      }
      
    }

    if (length(chromosome.list) == 1) save(list = c('modules.final', 'tableChrom'), file = chromFile)
    write.table(x = table.chrom, file = tableChrom, row.names = FALSE, quote = FALSE, sep = '\t')
  }


  ####### Now I save the output below after each chromosome
  if (length(chromosome.list) >= 22) { 
    message('Nb of modules ', n.modules)
    save(list = c('list.modules', 'table.final'), file = finalFile)
    write.table(x = table.final, file = tableFile, row.names = FALSE, quote = FALSE, sep = '\t')
  }
  
  return ( list( modules = list.modules, table = table.final ) )
}
  

#print( subset(my.res.trans[[ 20 ]], SNP == 'rs6027303' ) )
#print ( subset(my.res.trans[[ '22' ]], SNP == 'rs7364123' ) )

#dataset <- 'WB_dexamethasone_DiRienzo'
#condition <- 'logFC'


#dataset <- 'monocytes_Knight'
#condition <- 'LPS24'
#choice.sets <- list()
#choice.sets[[ dataset ]] <- condition


#modules <- find.all.trans.eQTLs ( choice.sets, run.stepwise = TRUE)
