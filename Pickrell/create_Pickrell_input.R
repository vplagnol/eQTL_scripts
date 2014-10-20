
create.Pickrell.input.file <- function (dataset,
                                        condition,
                                        genotypes, expression,
                                        ProbeID, chromosome, snp.name, min.MAF = 0.03, min.certain.calls = 0.,
                                        #gene.chromosome, gene.position.start, gene.position.end,
                                        base.folder = '/cluster/project8/vyp/eQTL_integration') {
  library(snpStats)
  
  
####

  sample.names.in.genotypes <- dimnames(genotypes$genotypes)[[1]]
  
  shared.samples <- intersect(dimnames(expression)[[2]], dimnames(genotypes$genotypes)[[1]])
  n.samples <- length(shared.samples)

                                          
  loc.expression <- as.numeric(expression[ ProbeID, ]) ### take the row of the expression data with the matching probe ID
  names(loc.expression) <- colnames(expression)
  loc.expression <- loc.expression [ sample.names.in.genotypes ]
  
  message('Nb of samples shared by expression and genotypes: ', n.samples)
  
########
  expression.frame <- data.frame(samples = sample.names.in.genotypes, exp = as.numeric(loc.expression))
  row.names(expression.frame) <- sample.names.in.genotypes
  expression.frame$exp <- expression.frame$exp / sd(expression.frame$exp, na.rm = TRUE) 
      
############# Now add covariates if there is such a file
  covariates.tab <- paste(base.folder, '/data/', dataset, '/covariates/covariates_', condition, '.tab', sep = '')  ##look for condition specific covariates first
  if (!file.exists(covariates.tab)) {covariates.tab <- paste(base.folder, '/data/', dataset, '/covariates/covariates.tab', sep = '')} ##otherwise the generic one
  
  if (file.exists(covariates.tab)) {
    covariates <- read.table(covariates.tab, header = TRUE, sep = '\t')
    row.names(covariates) <- covariates$id
    covariates <- covariates[ sample.names.in.genotypes, ]
    covar.labels <- subset ( names(covariates), names(covariates) != 'id')
    my.formula <- paste('exp ~ ', paste(covar.labels, collapse = ' + '))
    for (covar.loc in covar.labels) {expression.frame[, covar.loc] <- covariates[, covar.loc]}
  } else {
    my.formula <- 'exp ~ 1'
  }
  

  SNP.position <- genotypes$map[ snp.name, ]$position
  good.SNPs <- which ( genotypes$map$position > SNP.position - 500000 & genotypes$map$position < SNP.position + 500000)
  loc.geno <- genotypes$genotypes[, good.SNPs ]
  loc.map <- genotypes$map[ good.SNPs, ]
  my.sum <- col.summary(loc.geno)

  ######### Now subset the good SNPs
  v.good.SNPs <- (my.sum$MAF > min.MAF) & my.sum$Certain.calls > min.certain.calls
  loc.geno <- loc.geno[, v.good.SNPs ]
  my.sum <- my.sum[ v.good.SNPs, ]
  loc.map <- loc.map[ v.good.SNPs, ]
  
  ####### Now we can properly compute the P-values
  message('Formula being used: ', my.formula)
  my.tests <- snp.rhs.tests(snp.data = loc.geno, data = expression.frame, formula = as.formula (my.formula), score = TRUE, robust=TRUE, uncertain=TRUE, family = 'gaussian')
  
  score= sapply(my.tests@score, "[[", 1)
  score.var= sapply(my.tests@score, "[[", 2)
  effect = score/score.var
  score.chisq = score^2/score.var
  pval <-  pchisq(score.chisq, df =1, lower.tail = FALSE)
  var.beta <- 1 / score.var
  SE <- sqrt(var.beta)

  res <- data.frame (SNPID = dimnames(loc.geno)[[2]],
                     Z = signif(effect/SE, 3),
                     F = signif(my.sum$RAF, 3),
                     PVAL = p.value ( my.tests),
                     POS = loc.map$position,
                     CHR = chromosome,
                     N = sum(!is.na(expression.frame$exp)),
                     beta = signif(effect, 4),
                     se.beta = signif(SE, 4),
                     ProbeID = ProbeID)

  #### add the information about the imputation now
  info.columns <- grep( names(loc.map), pattern = '^info\\.', value = TRUE)
  if (length(info.columns) > 0) {
    for (col in info.columns) {res[, col] <- loc.map[, col]}
  }
  
  return(res)

}
