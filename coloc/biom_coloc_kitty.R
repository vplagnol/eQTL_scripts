# R Script for running colocalisation analysis of two eQTL datasets 
# Written by Kitty Lo 
# On 9th June 2014 
library(coloc) 
library(VPlib)

#############################
# For plotting with locuszoom
#refFlat_path = "/cluster/project8/vyp/vincent/toolsVarious/locuszoom/refFlat.RData"
#load(refFlat_path)
#refFlatRaw <- refFlatRaw.VP

##########################
# Variables to set 

setwd("/cluster/project8/vyp/eQTL_integration/") 
chr.name <- "8"
f1.name <- "data/UCLEB/summaryStats/ucleb_biom8.RData" 
f2.name <- "data/brain_UKBEC/eQTLs/fgwas/fgwas_core_CRBL_summary_eQTLs.csv" 
p12 <- 1e-6
#pdf("/cluster/project8/vyp/eQTL_integration/scripts/coloc/ucleb_CRB1_2.pdf", width = 9, height = 6)
outfname <- "/cluster/project8/vyp/eQTL_integration/scripts/coloc/ucleb_CRB.csv" 
biom.data <- "/cluster/project8/vyp/eQTL_integration/scripts/coloc/biom.names.RData" 
load(biom.data) 
# For the eqtl file, this is the +/- region around which we construct the region around the top snp 
##########################

#########################################################
# This function takes in:
#  eqtl.fname -  file name of eQTL file with the summary info 
#  biom.fname -  .RData file name with the biom info 
#  biom.names -  data.frame with the names of the biomarkers  
#  p12        -  probability of trait 1 and trait 2 
#  outfname   -  name of csv file output 
#  chr.name   -  name of chromosome (just the number part) 
 
coloc.eqtl.biom <- function(eqtl.fname, biom.rdata.fname, biom.names, p12, outfname, chr.name) {   
gap <- 500000 
load(biom.rdata.fname)
biom.data <- ucleb 
eqtl.data <- read.csv(eqtl.fname, sep =',', stringsAsFactors = FALSE)

# Subset by chromosome 
eqtl.data <- subset(eqtl.data, CHR == chr.name) 
message("There are ", nrow(eqtl.data), " regions in chr ", chr.name) 

# Work out what is the shared list of probes in the two datasets 
d2.probes <- unique(data2$ProbeID) 
common.probes <- d1.probes[which(d1.probes %in% d2.probes)] 

# Assume that the number of samples is the same for every line in each dataset 
N.eqtl <- eqtl.data$N[1] 

res.all <- c() 

for (i in 1:nrow(eqtl.data))  { 
    pos.start <- eqtl.data[i, "POS"] - gap 
    pos.end   <- eqtl.data[i, "POS"] + gap 
    matches <- which(biom.data$POS > pos.start & biom.data$POS < pos.end ) 
    
    # There can be more than one probe per gene, merge all probes together 
    region.eqtl <- read.csv(as.character(eqtl.data[i, "output.file"]) ,
                         sep='\t', stringsAsFactors = FALSE) 
    region.biom <- biom.data[matches, ] 

    # Loop over each biomarker 

    gene <- eqtl.data[i,"Gene.name"]
    probeID <- eqtl.data[i,"ProbeID"]

    message(gene, ": ", length(matches), " snps in biomarkers. From: ", pos.start, " To: ", pos.end) 
    for (j in 1:length(biom.names)) { 
       colname.pval <-  paste("PVAL.", biom.names[j], sep = "")   
       colname.N    <-  paste("N.", biom.names[j], sep = "")   
       colname.beta <-  paste("BETA.", biom.names[j], sep = "")   
       colname.se   <-  paste("SE.", biom.names[j], sep = "")   
       
       # Construct the MAF data.frame, since we only want the shared snps, just do it for one dataset  
       maf.biom   <- data.frame(snp = region.biom$SNPID, maf = region.biom$F) 
   
       merged.data <- merge(region.biom, region.eqtl, by = "SNPID")  

       # Remove the NA, assuming they are in the beta columns 
       merged.data <- merged.data[!is.na(merged.data$beta),] 
       merged.data <- merged.data[!is.na(merged.data[,colname.beta]),] 

       dataset.biom = list(snp = merged.data$SNPID, beta = merged.data[, colname.beta], 
               varbeta = merged.data[,colname.se]^2, N = merged.data[1,colname.N], type = "quant")
       dataset.biom$MAF <-  maf.biom[match(merged.data$SNPID, maf.biom$snp ) ,"maf"]

       dataset.eqtl = list(snp = merged.data$SNPID, beta = merged.data$beta, 
               varbeta = merged.data$se.beta^2, N = N.eqtl, type = "quant")
       dataset.eqtl$MAF <-  maf.biom[match(merged.data$SNPID, maf.biom$snp ) ,"maf"]

       suppressMessages(capture.output(coloc.res <- coloc.abf(dataset.biom, dataset.eqtl, p12 = p12))) 
       pp3       <- coloc.res$summary[5]
       pp4       <- coloc.res$summary[6]

       res.all <- rbind(res.all, c(ProbeID = probeID, Chr = chr.name, Gene.name = gene, 
                                   biom = biom.names[j],pp3, pp4) ) 
    } 
} 
res.all <- data.frame(res.all)
res.all$PP.H4.abf <- as.numeric(as.character(res.all$PP.H4.abf)) 
res.all$PP.H3.abf <- as.numeric(as.character(res.all$PP.H3.abf)) 
res.all <- res.all[with(res.all, order(PP.H4.abf)),] 
write.table(res.all, quote = FALSE, sep = ",", row.names = FALSE, file = outfname) 

} 
###############################
# Plotting code 
plot.locuszoom <- function(data, snp.chrom, snp.pos, probeID, gene.name) { 
    lz.frame <- data.frame(SNP = data$SNPID,
                           chr_name= data$CHR,
                           chrom_start = data$POS,
                           P.value = data$PVAL)
    lz.frame$locus_snp <- paste('chr', lz.frame$chr_name, ':', lz.frame$chrom_start, sep = '')
    lz.frame.file <- 'PD_lz.tab'
    write.table( x = lz.frame, file = lz.frame.file, row.names = FALSE, quote = FALSE, sep = '\t')
    source('/cluster/project8/vyp/vincent/toolsVarious/locuszoom/call_locuszoom3_temp.R') 
    locuszoom.ugi(metal = lz.frame.file,
              refSnp = paste("chr",snp.chrom, ":", snp.pos, sep = "") , #rs10877835
              title = gene.name,
              pvalCol='P.value',
              legend = 'left',
              markerCol='locus_snp',
              ylab = paste('-log10p for ', gene.name, ' with probe ', probeID, sep=""),
              chrCol= 'chr_name',
              posCol = 'chrom_start',
              chr = snp.chrom,
              showGenes = TRUE,
              show_xlab=FALSE,
              start= min(data$POS) ,
              end = max(data$POS)
              )

}
