
create.Pickrell.input.file <- function (dataset,
                                        condition,
                                        genotypes, expression,
                                        ProbeID, chromosome, snp.name, min.MAF = 0.03, min.certain.calls = 0.3,
                                        gene.chromosome, gene.position.start, gene.position.end,
                                        base.folder = '/cluster/project8/vyp/eQTL_integration') {
  library(snpStats)
  
  
####
  shared.samples <- intersect(dimnames(expression)[[2]], dimnames(genotypes$genotypes)[[1]])
  n.samples <- length(shared.samples)
  expression <- expression[, shared.samples]
  genotypes$genotypes <- genotypes$genotypes[ shared.samples, ]
  genotypes$fam <- genotypes$genotypes[ shared.samples, ]
  message('Nb of samples shared by expression and genotypes: ', n.samples)
  
  
########
  expression.frame <- data.frame(samples = shared.samples, exp = NA)
  row.names(expression.frame) <- shared.samples
  
      
############# Now add covariates if there is such a file
  covariates.tab <- paste(base.folder, '/data/', dataset, '/covariates/covariates_', condition, '.tab', sep = '')  ##look for condition specific covariates first
  if (!file.exists(covariates.tab)) {covariates.tab <- paste(base.folder, '/data/', dataset, '/covariates/covariates.tab', sep = '')} ##otherwise the generic one
  
  if (file.exists(covariates.tab)) {
    covariates <- read.table(covariates.tab, header = TRUE, sep = '\t')
    row.names(covariates) <- covariates$id
    covariates <- covariates[ shared.samples, ]
    covar.labels <- subset ( names(covariates), names(covariates) != 'id')
    my.formula <- paste('exp ~ ', paste(covar.labels, collapse = ' + '))
    for (covar.loc in covar.labels) {expression.frame[, covar.loc] <- covariates[, covar.loc]}
  } else {
    my.formula <- 'exp ~ 1'
  }
  

  SNP.position <- genotypes$map[ snp.name, ]$position
  good.SNPs <- which ( genotypes$map$position > SNP.position - 500000 & genotypes$map$position < SNP.position + 500000)
  loc.geno <- genotypes$genotypes[, good.SNPs ]
  my.sum <- col.summary(loc.geno)
  genotypes$map <- genotypes$map[ good.SNPs, ]

  ######### Now subset the good SNPs
  v.good.SNPs <- (my.sum$MAF > min.MAF) & my.sum$Certain.calls > min.certain.calls
  loc.geno <- loc.geno[, v.good.SNPs ]
  my.sum <- my.sum[ v.good.SNPs, ]
  genotypes$map<- genotypes$map[ v.good.SNPs, ]
  
#######
  loc.expression <- expression[ ProbeID, ] ### take the row of the expression data with the matching probe ID
  #expression.frame <- data.frame(exp = as.numeric(loc.expression)/sd(as.numeric(loc.expression)))
  #row.names(expression.frame) <- names(loc.expression)
  expression.frame$exp <- as.numeric(loc.expression)/sd(as.numeric(loc.expression))
  

  ####### Now we can properly compute the P-values
  message('Formula being used: ', my.formula)

  my.tests <- snp.rhs.estimates( snp.data = loc.geno, data = expression.frame, formula = as.formula (my.formula), family = 'gaussian')
  
  effect <- sapply(as(my.tests, 'list'), FUN = function(my.tests) {if (is.null(my.tests)) {return (NA);} else {return((my.tests)$beta)} } )
  SE <- sapply(as(my.tests, 'list'), FUN = function(my.tests) { if (is.null(my.tests)) {return (NA);} else {return(sqrt((my.tests)$Var.beta))} } )

  res <- data.frame (SNPID = dimnames(loc.geno)[[2]],
                     Z = signif(effect/SE, 3),
                     F = signif(my.sum$RAF, 3),
                     PVAL = signif(2*pnorm(-abs(effect/SE)), 3),
                     POS = genotypes$map$position,
                     CHR = chromosome,
                     N = sum(!is.na(expression.frame$exp)),
                     beta = signif(effect, 4),
                     se.beta = signif(SE, 4),
                     ProbeID = ProbeID)

  #### add the information about the imputation now
  info.columns <- grep( names(genotypes$map), pattern = '^info\\.', value = TRUE)
  if (length(info.columns) > 0) {
    for (col in info.columns) {res[, col] <- genotypes$map[, col]}
  }
  
  return(res)

}
