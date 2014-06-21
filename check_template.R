library(snpStats)

#input <- 'data/LCL_dexamethasone_DiRienzo'
#input <- 'data/WB_dexamethasone_DiRienzo'
input <- "data/DC_MTB_Barreiro"

code <- basename(input)

output.file <- paste(input, '/check_automated.txt', sep = '')
if (file.exists(output.file)) file.remove(output.file)
message('Output will be placed in ', output.file)
sink(output.file)




##############
genotypes <- paste(input, '/genotypes/chr', 1:22, sep = '')

if(sum(!file.exists(genotypes)) > 0) {
  stop('Missing genotype files')
}


nSNPs <- 0
for (gfile in genotypes) {
  print(gfile)
  load(gfile)
  nSNPs <- nSNPs + ncol(genotypes$genotypes)

  must.have.names <- c('snp.name', 'chromosome', 'position', 'allele.1', 'allele.2')
  if (sum(!  must.have.names %in% names(genotypes$map)) > 0) {
    print(subset(must.have.names, ! must.have.names %in% names(genotypes$map)))
    stop()
  }
}

nsamples.geno <- nrow(genotypes$fam)
cat('Nb of imputed SNPs ', nSNPs, '\n')


######### is there a covariate file?
covariates.file <- paste(input, '/genotypes/covariates.tab', sep = '')
if (file.exists(covariates.file)) {
  covar <- read.table(covariates.file, header = TRUE)
  if (! 'id' %in% names(covar)) {stop('ID is not in the header of the covariates table')}
  if (sum (covar$id != row.names(genotypes$genotypes)) > 0) {stop('Covariate file ID column is not matching the genotype tables')}
  matrixeQTL.covar <- gsub(pattern = 'tab$', replacement = 'matrixeQTL', covariates.file)
  if (!file.exists(matrixeQTL.covar)) {stop('The matrix eQTL covariates file is missing')}
}



################
expr.files <- list.files(paste(input, '/expression_data', sep = ''), full.names = TRUE)
conditions <- gsub(replacement = '', pattern = 'expression_|.RData', x = basename(expr.files))


cat('Nb of conditions ', length(conditions), '\n')

for (exp in expr.files) {

  cat('Looking at ', exp, '\n')
  message('Looking at ', exp)
  
  load(exp)
  exp.data <- get(conditions[1])
  support <- get(paste('support', conditions[1], sep = '.'))


  ####### check the columns with the gene name information
  required.columns <- c('ensemblID', 'Gene.name')
  for (col in required.columns) {
    if (! col %in% names(support)) {stop('Support file is missing ', col, 'column')}
    message('Non missing column ', col, ' for ', signif(100*sum(!is.na(support[, col]))/nrow(support), 3), '% of the support file')
  }

  
  
  nsamples <- dim(exp.data)[2]

  cat('Number of samples in expression data ', nsamples, '\n')
  cat('Number of samples in genotype data ', nsamples.geno, '\n')

  if (nsamples != nsamples.geno) {
    stop('Non matching number of samples')
  }

  if (sum(row.names(genotypes$fam) != dimnames(exp.data)[[2]])) {
    print(row.names(genotypes$fam))
    print(dimnames(exp.data)[[2]])
    stop('Non matching names between genotypes and expression data')
  }
  
}

sink()

message('Result in ', output.file)
