# R Script for running colocalisation analysis of two eQTL datasets 
# Written by Kitty Lo 
# On 9th June 2014 
library(coloc) 
library(VPlib)

##########################
# Variables to set 

setwd("/cluster/project8/vyp/eQTL_integration/") 
f1.name <- "data/monocytes_Knight/eQTLs/fgwas/fgwas_LPS24_summary_eQTLs.csv" 
f2.name <- "data/Smith_macrophages/eQTLs/fgwas/fgwas_logFC_summary_eQTLs.csv" 
p12 <- 1e-6
pdf("/cluster/project8/vyp/kitty/tools/knight_smith.pdf", width = 9, height = 6)
outfname <- "/cluster/project8/vyp/kitty/tools/knight_smith.csv" 
##########################

coloc.eqtl.eqtl <- function(eqtl.f1, eqtl.f2, p12, outfname, pdf.fname) { 

   # pdf of most significant regions with PP4 values 
   pdf(pdf.fname, width = 9, heigth = 6) 

   data1 <- read.csv(eqtl.f1, sep =',', stringsAsFactors = FALSE) 
   data2 <- read.csv(eqtl.f2, sep =',', stringsAsFactors = FALSE)

   # Work out what is the shared list of probes in the two datasets 
   d1.probes <- unique(data1$ProbeID) 
   d2.probes <- unique(data2$ProbeID) 
   common.probes <- d1.probes[which(d1.probes %in% d2.probes)] 

   # Assume that the number of samples is the same for every line in each dataset 
   N1 <- data1$N[1] 
   N2 <- data2$N[1] 

   res.all <- c() 

   for (i in 1:length(common.probes))  { 
      match.1 <- which(data1$ProbeID == common.probes[i]) 
      match.2 <- which(data2$ProbeID == common.probes[i])
    
      region.1 <- read.csv(as.character(data1[match.1[1], "output.file"]) ,
                         sep='\t', stringsAsFactors = FALSE) 
      region.2 <- read.csv(as.character(data2[match.2[1], "output.file"]) ,
                         sep='\t', stringsAsFactors = FALSE) 

      merged.data <- merge(region.1, region.2, by = "SNPID") 

      # Construct the MAF data.frame, since we only want the shared snps, just do it for one dataset  
      maf1   <- data.frame(snp = region.1$SNPID, maf = region.1$F) 
    
      dataset1 = list(snp = merged.data$SNPID, beta = merged.data$beta.x, 
               varbeta = merged.data$se.beta.x^2, N = N1, type = "quant")
      dataset1$MAF <-  maf1[match(merged.data$SNPID, maf1$snp ) ,"maf"]

      dataset2 = list(snp = merged.data$SNPID, beta = merged.data$beta.y, 
               varbeta = merged.data$se.beta.y^2, N = N2, type = "quant")
      dataset2$MAF <-  maf1[match(merged.data$SNPID, maf1$snp ) ,"maf"]
    
      coloc.res <- coloc.abf(dataset1, dataset2, p12 = p12) 
      pp4       <- coloc.res$summary[6]

      gene <- region.1[which(region.1$ProbeID == common.probes[i]), "Gene.name"]

      # Make a plot if the pp4 value is high enough  
      # Do it per probe for each of the two dataset 

      if (pp4 > 0.9) { 
       plot(region.1$POS, -log10(region.1$PVAL), col = "red", cex = 0.7,
            ylab = "-log10(Pval)", xlab = paste("chr", data2[match.2,"CHR"], sep = ""), 
            main = paste(gene[1], ":", common.probes[i], ", PP4 = ", pp4, sep = "")) 
       points(region.2$POS, -log10(region.2$PVAL), col = "blue", cex = 0.7) 
      }
      res.all <- rbind(res.all, c(ProbeID = common.probes[i], 
                       Gene.name = gene[1], pp4) ) 
   } 
   dev.off()  
   res.all <- data.frame(res.all) 
   res.all <- res.all[with(res.all, order(PP.H4.abf)),] 
   write.table(res.all, quote = FALSE, sep = ",", row.names = FALSE, file = outfname) 

} 
