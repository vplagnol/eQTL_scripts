extra.cov <- read.table('data/WB_dexamethasone_DiRienzo/genotypes/extra_covariates.tab', header = TRUE, stringsAsFactors = FALSE)
condition <- 'logFC'

load('data/WB_dexamethasone_DiRienzo/expression_data/expression_logFC.RData')
load('data/WB_dexamethasone_DiRienzo/modules/logFC_all_modules.RData')
test.module <- my.modules[['15']][[1]]

expression <- get(condition)


my.pca <- prcomp(t(expression), center = TRUE)

chrom <- test.module$chromosome[1]
SNP <- test.module$SNP[1]


library(snpStats)
geno.file <- paste('data/WB_dexamethasone_DiRienzo/genotypes/chr', chrom, sep = '')
load(geno.file)


geno.num <- as.numeric(genotypes$genotype[,SNP])

cov.names <- names(extra.cov)[- (1:2) ]
my.pvals.covar <- rep(NA, length(cov.names))




for (i in 1:length(cov.names)) {
  my.pvals[ i ] <- cor.test( geno.num, extra.cov[, cov.names[ i ]])$p.value
}
