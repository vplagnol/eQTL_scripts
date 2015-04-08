# Set of functions to process RNAseq expression data  
# Written by Kitty Lo 
# On 29th Jan, 2015 

##############################################################
# This function writes the PEER factors 
# Taken from Claudia's code: 
# /cluster/project8/vyp/eQTL_integration/data/GTex/scripts/GTEx_process.R

do.peer <- function(x_trans, outFile_pca, samples.id, n.factors = 15) { 
     require(peer)

     # Build the model
     model=PEER()

     # Set the maximum number of unobserved factors to model.
     PEER_setNk(model,n.factors)  #n_unobserved_factors=25?

     # Set expression data
     PEER_setPhenoMean(model, as.matrix(x_trans))   # (NULL response means no error here)
     PEER_update(model)

     factors = PEER_getX(model)
     print(dim(factors)) 

     factors <- as.data.frame(factors)

     names(factors)<- paste("PC", 1:n.factors, sep="")

     factors$id <- samples.id

     write.table(factors,outFile_pca,quote=F,row.names=F,sep="\t")

} 

##############################################################
# This function log and quantile normalises RPKM similar to the 
# GTEx analysis procedure 
# Taken from Claudia's code: 
# /cluster/project8/vyp/eQTL_integration/data/GTex/scripts/GTEx_process.R
# Inputs: 
# rpkm is a matrix with rows as genes and columsn as samples 
# Outputs: 
# matrix of rpkm values that has been normalised and
# genes with low expression values removed  

normalise.RPKM <- function( rpkm, n.thresh = 10, rpkm.thresh = 0.1, log.norm = TRUE ) { 


   if (log.norm) {
     ## First filter out low expression genes 
     n.well.expressed <- apply(MAR = 1, rpkm, function(x) { sum(x > rpkm.thresh) } ) 
     rpkm <- rpkm[n.well.expressed > n.thresh, ]  
     
     shift    <- 0.01 
     rpkm     <- log2(rpkm + shift) 
   }

   genes <- rownames(rpkm) 
   sampleIDs <- colnames(rpkm) 
   
   # Log and quantile normalise 
   library(preprocessCore)
   rpkm     <- normalize.quantiles(rpkm) 

   # Outliers correction 
   # For each gene, rank values across samples then map to standard normal 
   rpkm.rank <- t(apply(rpkm, 1, rank, ties.method = "average")  ) 
   rpkm <- qnorm(rpkm.rank/ (ncol(rpkm.rank) + 1) ) 

   rownames(rpkm) <- genes 
   colnames(rpkm) <- sampleIDs

   return(rpkm) 

} 


