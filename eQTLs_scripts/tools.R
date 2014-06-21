pretty.names <- function(dataset, condition) {
  text <- NA
  
  if (dataset == 'LCL_dexamethasone_DiRienzo') {
    if (condition == 'treated') text <- 'LCL, Dex 8h'
    if (condition == 'untreated') text <- 'LCL, Resting'
    if (condition == 'logFC') text <- 'LCL, Dex logFC 8h'
  }

  if (dataset == 'WB_dexamethasone_DiRienzo') {
    if (condition == 'treated') text <- 'WB, Dex 8h'
    if (condition == 'untreated') text <- 'WB, Resting'
    if (condition == 'logFC') text <- 'WB, Dex logFC 8h'
  }

  if (dataset == 'DC_MTB_Barreiro') {
    if (condition == 'infected') text <- 'DC, mTB 16h'
    if (condition == 'normal') text <- 'DC, Resting'
    if (condition == 'logFC') text <- 'DC, mTB infected, logFC at 16h'
  }

  if (dataset == 'monocytes_Knight') {
    if (condition == 'normal') text <- 'Monocytes, rested'
    if (condition == 'LPS2') text <- 'Monocytes, LPS 2h'
    if (condition == 'LPS24') text <- 'Monocytes, LPS 24h'
  }

  if (dataset == 'GTex') {
    text <- paste('GTex, ', condition)
  }
  
  return(text)
}

make.matEQTL.geno.basic <- function ( genotypes, output.file) {
  geno.num <- t(as(genotypes, 'numeric'))
  message('Number of samples: ', ncol(geno.num))
  
  cat('snpid\t', file = output.file)
  write.table(x = geno.num, file = output.file, sep = '\t', row.names = TRUE, quote = FALSE, append = TRUE, col.names = TRUE)  

}

make.matEQTL.geno <- function (genotypes, chromosome, start, end, output.file, min.freq = 0.01) {
  library(snpStats)
  message('Preparing a matrixEQTL genotype dataset')

  good.SNPs <- which ( genotypes$map$position >= start & genotypes$map$position <= end)
  snpStats.format <- genotypes$genotypes[, good.SNPs]
  support.data <- col.summary(snpStats.format)
  support <- genotypes$map[ good.SNPs, ]

######## Now exclude low MAF SNPs
  exclude.low.MAF <- support.data$MAF < min.freq
  snpStats.format <- snpStats.format[, ! exclude.low.MAF]
  support.data <- support.data[ ! exclude.low.MAF, ]
  support <- support[ ! exclude.low.MAF, ]

  support$MAF <- support.data$MAF
  support$Call.rate <- support.data$Call.rate
  
###########
  geno.num <- t(as(snpStats.format, 'numeric'))
  message('Number of samples: ', ncol(geno.num))
  
  cat('snpid\t', file = output.file)
  write.table(x = geno.num, file = output.file, sep = '\t', row.names = TRUE, quote = FALSE, append = TRUE, col.names = TRUE)
 

  return( support )
}


make.matEQTL.expression <- function( expression, output.file) {
  message('Preparing a matrixEQTL expression dataset')
  message('Number of samples: ', ncol(expression))
  
  cat('geneid\t', file = output.file)
  write.table(x = expression, file = output.file, quote = FALSE, row.names = TRUE, col.names = TRUE, append = TRUE, sep ='\t')
  return( 'done' )
}

