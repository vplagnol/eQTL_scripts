library(snpStats)
source('scripts/plotting_functions/plot_eQTLs_per_SNP.R')
options(stringsAsFactors = FALSE)

plot <- FALSE

pval <- 5

#dataset <- 'WB_dexamethasone_DiRienzo'
#condition <- 'logFC'

dataset <- 'monocytes_Knight'
condition <- 'LPS24'



min.gene.module <- 6

choice.sets <- list()
choice.sets[[ dataset ]] <- condition

extra.datasets <- c()
if (dataset == 'WB_dexamethasone_DiRienzo') {
  extra.datasets <- 'LCL_dexamethasone_DiRienzo'
}

if (dataset == 'LCL_dexamethasone_DiRienzo') {
  extra.datasets <- 'WB_dexamethasone_DiRienzo'
}


with.extra <- length(extra.datasets) > 0
########## prepare the PCA analysis
expression.file <- paste('data/', dataset, '/expression_data/expression_', condition, '.RData', sep = '')
load(expression.file)
expression <- get(condition)


my.res.trans <- list()
my.res.cis <- list()

oFolder <- paste('data/', dataset, '/modules/', condition, sep = '')
if (!file.exists(oFolder)) dir.create(oFolder)
old.files <- list.files(oFolder, pattern = condition, full.names = TRUE)
file.remove(old.files)


table.final <- data.frame()

list.modules <- list()
n.modules <- 0

first <- TRUE
for (chr in as.character(22:1)) {
#for (chr in as.character(1)) {
  
  message('Chromosome ', chr)
  geno.file <- paste('data/', dataset, '/genotypes/chr', chr, sep = '')
  load(geno.file)


  ####### take samples in common only
  if (first) {
    shared.samples <- intersect(dimnames(genotypes$genotypes)[[1]], dimnames(expression)[[2]])
    expression <- expression[, shared.samples]
    my.pca <- prcomp(t(expression), center = TRUE)
    first <- FALSE
  }
  genotypes$genotypes <- genotypes$genotypes[ shared.samples, ]

  
  
  input.file <- paste('data/', dataset, '/eQTLs/matrixEQTL/', condition, '/', condition, '_chr', chr, '_pval', pval, '.tab', sep = '')
  data <- read.table(input.file, sep = '\t', header = TRUE)
  data$dataset <- dataset
  data <- subset(data, !is.na(ensemblID))


  ######## The additional data where we can look for consistent signals
  if (with.extra) {
    for (extra.dataset in extra.datasets) {
      extra.data.file <- paste('data/', extra.dataset, '/eQTLs/matrixEQTL/', condition, '/', condition, '_chr', chr, '_pval', pval, '.tab', sep = '')
      extra.data <- read.table(extra.data.file, sep = '\t', header = TRUE)
      extra.data$dataset <- extra.dataset
    }
  }


  
  ###use the idea of auxiliiary gene
  strong.candidates <- subset( data, p.value < 10^-12 & MAF > 0.03 & !cis.eQTL, 'SNP', drop = TRUE)

  ####now look for SNPs with trans effects on at least 10 genes
  my.tab <- table(subset(data$SNP, data$MAF > 0.03 & ! data$cis.eQTL))
  good.trans <- names(which(my.tab  >= min.gene.module))
  
  strong.candidates.bis <- subset(data, SNP %in% c(good.trans, strong.candidates))

  my.res.trans[[ chr ]] <- strong.candidates.bis
  my.res.cis[[ chr ]] <- subset(data, cis.eQTL & SNP %in% strong.candidates.bis$SNP)

  ###### Now define the modules below
  loc.table <- my.res.trans[[ chr ]]
  my.tab <- table(c('dummy', loc.table$SNP))

  while ((nrow(loc.table) > 0) && (max(my.tab) >= min.gene.module)) {
    n.modules <- n.modules + 1
    my.sort <- sort(my.tab, decreasing = TRUE)
    best.SNP <- names(my.sort)[[1]]
    best.module <- subset( loc.table, SNP == best.SNP)
    my.position <- best.module$position[ 1 ]
    list.modules[[ n.modules ]] <- best.module
    names(list.modules)[ n.modules ] <- best.SNP

    if (with.extra) {
      extra.data.matching <- subset(extra.data, SNP == best.SNP  | position == my.position)
      extra.data.matching$shared.gene <- extra.data.matching$ensemblID %in% best.module$ensemblID
    }
    
    my.res <-  list (chromosome = chr,
                     position = my.position,
                     best.SNP = best.SNP,
                     n.genes = length(unique(best.module$ensemblID)),
                     n.cis.eQTLs =  sum(best.module$cis.eQTL))

    if (with.extra) {
      my.res$n.matching.signals <- nrow(extra.data.matching)
      my.res$n.matching.signals.and.genes <- sum( extra.data.matching$shared.gene )
    }

      
    
########## Now the PCA analysis
    geno.num <- as.numeric(genotypes$genotypes[, best.SNP])
    for (i in 1:5) {
      my.res[[  paste('PCA', i, sep = '') ]] <- cor.test( geno.num, my.pca$x[, i])$p.value
    }
    
###########
    table.final <- rbind.data.frame(table.final, my.res)
    #print(table.final)

################
    loc.table <- subset(loc.table, ! ensemblID %in% best.module$ensemblID)
    
    my.tab <- table(c('dummy', loc.table$SNP))

    if (plot) {
      output.pdf <- paste(oFolder, '/', condition, '_', best.module$SNP[1], '_chr', best.module$chromosome[1], '_',  best.module$position[1], '.pdf', sep = '')
      test <- plot.eQTL (choice.sets, snp.name = best.module$SNP[1], chromosome = best.module$chromosome[1], gene.names = best.module$Gene.name, output.pdf = output.pdf)
      message('Output pdf for this module in ', output.pdf)
    }

  }
  
  message('Nb of modules ', n.modules)
  finalFile <- paste('data/', dataset, '/modules/', condition, '_all_modules.RData', sep = '')
  save(list = c('list.modules', 'table.final'), file = finalFile)
}

write.table(x = table.final, file = paste('data/', dataset, '/modules/', condition, '_main_table.tab', sep = ''), row.names = FALSE, quote = FALSE, sep = '\t')




#print( subset(my.res.trans[[ 20 ]], SNP == 'rs6027303' ) )
print ( subset(my.res.trans[[ '22' ]], SNP == 'rs7364123' ) )
