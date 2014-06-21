# R Script for running colocalisation analysis of two eQTL datasets 
# Written by Kitty Lo 
# On 9th June 2014 
library(coloc) 
library(VPlib)

 
#############################
# For plotting with locuszoom
refFlat_path = "/cluster/project8/vyp/vincent/toolsVarious/locuszoom/refFlat.RData"
load(refFlat_path)
refFlatRaw <- refFlatRaw.VP

##########################
# Variables to set 

setwd("/cluster/project8/vyp/eQTL_integration/") 
f1.name <- "data/monocytes_Knight/eQTLs/fgwas/fgwas_LPS24_summary_eQTLs.csv" 
f2.name <- "data/Smith_macrophages/eQTLs/fgwas/fgwas_logFC_summary_eQTLs.csv" 
p12 <- 1e-6
pdf("/cluster/project8/vyp/kitty/tools/knight_smith.pdf", width = 9, height = 6)
outfname <- "/cluster/project8/vyp/kitty/tools/knight_smith.csv" 
##########################

data1 <- read.csv(f1.name, sep =',', stringsAsFactors = FALSE) 
data2 <- read.csv(f2.name, sep =',', stringsAsFactors = FALSE)

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
    
    # There can be more than one probe per gene, merge all probes together 
    region.1 <- read.csv(as.character(data1[match.1[1], "output.file"]) ,
                         sep='\t', stringsAsFactors = FALSE) 
    region.2 <- read.csv(as.character(data2[match.2[1], "output.file"]) ,
                         sep='\t', stringsAsFactors = FALSE) 

    merged.data <- merge(region.1, region.2, by = "SNPID") 
    # Construct the MAF data.frame, since we only want the shared snps, just do it for one dataset  
    maf1   <- data.frame(snp = region.1$SNPID, maf = region.1$F) 
    
    dataset1 = list(snp = merged.data$SNPID, beta = merged.data$beta.x, 
               varbeta = merged.data$se.beta.x^2, N = N1, type = "quant")
    #dataset1 = list(snp = merged.data$SNPID, pvalues = merged.data$PVAL.x, 
    #                N = N1, type = "quant")
    dataset1$MAF <-  maf1[match(merged.data$SNPID, maf1$snp ) ,"maf"]
    #dataset1$beta <- ifelse(dataset1$MAF <= 0.5, dataset1$beta , -dataset1$beta )
    #dataset1$MAF  <- ifelse(dataset1$MAF <= 0.5, dataset1$MAF , 1-dataset1$MAF )

    dataset2 = list(snp = merged.data$SNPID, beta = merged.data$beta.y, 
               varbeta = merged.data$se.beta.y^2, N = N2, type = "quant")
    #dataset2 = list(snp = merged.data$SNPID, pvalues = merged.data$PVAL.y, 
    #                N = N2, type = "quant")
    dataset2$MAF <-  maf1[match(merged.data$SNPID, maf1$snp ) ,"maf"]
    #dataset2$beta <- ifelse(dataset2$MAF <= 0.5, dataset2$beta , -dataset2$beta )
    #dataset2$MAF  <- ifelse(dataset2$MAF <= 0.5, dataset2$MAF , 1-dataset2$MAF )
    
    coloc.res <- coloc.abf(dataset1, dataset2, p12 = p12) 
    pp4       <- coloc.res$summary[6]

    gene <- region.1[which(region.1$ProbeID == common.probes[i]), "Gene.name"]

    # Make a plot if the pp4 value is high enough  
    # Do it per probe for each of the two dataset 

    #if (pp4 > 0.8) { 
       #output.pdf <- paste('/cluster/project8/vyp/kitty/tools/', common.probes[i], '.pdf', sep = '') 
       #pdf(output.pdf, width = 9, height = 6)
       #plot.locuszoom(region.1, data1[match.1,"CHR"], data1[match.1, "POS"], 
       #               data1[match.1, "ProbeID"], data1[match.1, "Gene.name"])   
       #system("rm region_ld.*")  
       #plot.locuszoom(region.2, data2[match.2,"CHR"], data2[match.1, "POS"], 
       #               data2[match.2, "ProbeID"], data1[match.2, "Gene.name"])   
       plot(region.1$POS, -log10(region.1$PVAL), col = "red", cex = 0.7,
            ylab = "-log10(Pval)", xlab = paste("chr", data2[match.2,"CHR"], sep = ""), 
            main = paste(gene[1], ":", common.probes[i], ", PP4 = ", pp4, sep = "")) 
       points(region.2$POS, -log10(region.2$PVAL), col = "blue", cex = 0.7) 
       #dev.off()
    #}
    res.all <- rbind(res.all, c(ProbeID = common.probes[i], Gene.name = gene[1], pp4) ) 
} 
dev.off()  
res.all <- data.frame(res.all) 
res.all <- res.all[with(res.all, order(PP.H4.abf)),] 
write.table(res.all, quote = FALSE, sep = ",", row.names = FALSE, file = outfname) 
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
