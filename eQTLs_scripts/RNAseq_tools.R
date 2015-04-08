# Set of functions to process RNAseq expression data  
# Written by Kitty Lo 
# On 29th Jan, 2015 

##############################################################
# This function log and quantile normalises RPKM similar to the 
# GTEx analysis procedure 
# Inputs: 
# rpkm is a matrix with rows as genes and columsn as samples 
# Outputs: 
# matrix of rpkm values that has been normalised and
# genes with low expression values removed  

normalise.RPKM <- function( rpkm, n.thresh = 10, rpkm.thresh = 0.1 ) { 

  #rpkm <- read.table(rpkm.fname, stringsAsFactor = FALSE, header = TRUE, sep = ",")
  # rownames(rpkm) <- rpkm$ensemblID
  # rpkm <- data.matrix(rpkm[,-c(1,2)]) 
   
   # Filter out low expression genes 
   n.well.expressed <- apply(MAR = 1, rpkm, function(x) { sum(x > rpkm.thresh) } ) 
   rpkm <- rpkm[n.well.expressed > n.thresh, ]  
   
   genes <- rownames(rpkm) 
   sampleIDs <- colnames(rpkm) 
 
   # Log and quantile normalise 
   library(preprocessCore)
   shift    <- 0.01 
   rpkm     <- log2(rpkm + shift) 
   rpkm     <- normalize.quantiles(rpkm) 

   # Outliers correction 
   # For each gene, rank values across samples then map to standard normal 
   rpkm.rank <- t(apply(rpkm, 1, rank, ties.method = "average")  ) 
   rpkm <- qnorm(rpkm.rank/ (ncol(rpkm.rank) + 1) ) 

   rownames(rpkm) <- genes 
   colnames(rpkm) <- sampleIDs

   return(rpkm) 

} 


