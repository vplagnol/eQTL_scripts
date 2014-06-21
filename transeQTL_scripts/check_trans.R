source('/cluster/project8/vyp/eQTL_integration/scripts/eQTLs_scripts/tools.R')
source('/cluster/project8/vyp/eQTL_integration/scripts/eQTLs_scripts/find_all_eQTLs.R')

prepare <- TRUE

if (prepare) {
  load('data/WB_dexamethasone_DiRienzo/modules/logFC_all_modules.RData')
  choice.snps <- names(list.modules)
  
  first <- TRUE
  for (chr in 1:22) {
    message('Chromosome ', chr)
    input.file <- paste('data/WB_dexamethasone_DiRienzo/genotypes/chr', chr, sep = '')
    load(input.file)

    inter.snps <- subset( choice.snps, choice.snps %in% dimnames(genotypes$genotypes)[[2]])
    
    if (length(inter.snps) > 0) {
      
      if (first) {
        final.geno <- genotypes$genotypes[, inter.snps ]
        final.map <- genotypes$map[ inter.snps, ]
        first <- FALSE
      } else {
        final.geno <- cbind( final.geno, genotypes$genotypes[, inter.snps ] )
        final.map <- rbind.data.frame(final.map, genotypes$map[ inter.snps, ])
      }
    }
  }
}


load('data/WB_dexamethasone_DiRienzo/expression_data/expression_logFC.RData')
genotypes <- final.geno
min.MAF <- 0.03
                                        #expression.matrix <- logFC
                                        #covariates <- NULL
                                        #pvOutputThreshold <- 10^(-2)
                                        #expression.support <- support.logFC



eQTL.data <- run.eQTL.selected(genotypes = final.geno,
                               expression.matrix = logFC,
                               expression.support = support.logFC,
                               pvOutputThreshold = 10^(-3))


pca <- prcomp(t(logFC))
geno.num <- as(final.geno, 'numeric')
print(table(dimnames(geno.num)[[1]] == dimnames(pca$x)[[1]]))

cor.matrix <- matrix(data = NA,
                     ncol = nrow(geno.num), nrow = ncol(geno.num),
                     dimnames = list( dimnames(geno.num)[[2]], paste('PC', 1:nrow(geno.num), sep = '')))


for (i in 1:ncol(geno.num)) {
  message(dimnames(geno.num)[[2]][i])
  for (j in 1:nrow(geno.num)) {
    cor.matrix[ i, j ] <- cor.test( geno.num[,i], pca$x[,j])$p.value
  }
}


to.check <- subset(eQTL.data, Gene.name == 'OAS2')
print(final.map[ to.check$SNP, ])


too.many <- subset(eQTL.data, SNP == 'rs11192101')$Gene.name
