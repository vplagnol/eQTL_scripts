# Written by Claudia Giambartolomei on 06/10/2014 
# For the eqtl file, this is the +/- region around which we construct the region around the top snp 
##########################

#########################################################
# This function takes in:
#  eqtl.fname -  file name of eQTL file with the summary info 
#  biom.fname -  .RData file name with the biom info 
#  biom.names -  data.frame with the names of the biomarkers  
#  p12        -  probability of trait 1 and trait 2 
#  outfolder   -  name of output folder with csv file output 
#  chr.name   -  name of chromosome (just the number part) 
#  plot.fldr  -  name of folder with locuszoom plots (just if pp4>0.5) 
#  type= "quant" or "cc", if "cc" must have a column with N.cases
#  prefix   - for the output file, usually condition and chromsome


coloc.eqtl.biom <- function(eqtl.fname, biom.rdata.fname, biom.names, p12, chr.name, type="quant", plot=TRUE, outfolder, prefix= "pref") {

suppressPackageStartupMessages({
  require(VPlib);
  require(coloc);
});

# For plotting with locuszoom # Already specified in main script
#refFlat_path = "/cluster/project8/vyp/vincent/toolsVarious/locuszoom/refFlat.RData"
#source('/cluster/project8/vyp/vincent/toolsVarious/locuszoom/call_locuszoom3_temp.R')
#load(refFlat_path)
#refFlatRaw <- refFlatRaw.VP

if (plot) {
   plot.fld = paste(outfolder, "plot/", sep="")
   pval.fld= paste(outfolder, "pval/", sep="")
   dir.create(file.path(plot.fld), showWarnings = FALSE)
   dir.create(file.path(pval.fld), showWarnings = FALSE)
}

   gap <- 500000
   biom.data <- get(load(biom.rdata.fname))
   eqtl.data <- read.csv(eqtl.fname, sep =',', stringsAsFactors = FALSE)

   # Subset by chromosome 
   eqtl.data <- subset(eqtl.data, CHR == chr.name)
   # Filter by imputation quality if column exists
   info.columns <- grep( names(eqtl.data), pattern = '^info\\.', value = TRUE)
   if (length(info.columns) > 0)	{
       eqtl.data = subset(eqtl.data, eqtl.data[,info.columns] > 0.4)
   }
   # Filter by MAF? 
   eqtl.data$MAF  <- ifelse(eqtl.data$F <= 0.5, eqtl.data$F , 1-eqtl.data$F )
   eqtl.data = subset(eqtl.data, eqtl.data$MAF > 0.01)

   message("There are ", nrow(eqtl.data), " regions in chr ", chr.name)
  
   # If Gene.name is missing, use ensemblID instead, then try to retrieve name from biomaRt. 
   if (nrow(eqtl.data[eqtl.data$Gene.name=="",]) > 0) {
     eqtl.data$Gene.name = ifelse(eqtl.data$Gene.name=="", eqtl.data$ensemblID, eqtl.data$Gene.name)
     if (length(eqtl.data$Gene.name[grep("ENSG", eqtl.data$Gene.name)]) >0 ) {
      library(biomaRt)
      mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
      res.gn <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = eqtl.data$Gene.name[grep("ENSG", eqtl.data$Gene.name)], mart = mart)
      res.gn = res.gn[res.gn$hgnc_symbol!="",]
      eqtl.data$Gene.name[which(eqtl.data$Gene.name %in% res.gn$ensembl_gene_id)]= res.gn[match(eqtl.data$Gene.name[which(eqtl.data$Gene.name %in% res.gn$ensembl_gene_id)], res.gn$ensembl_gene_id), "hgnc_symbol"]
     }
   }
   # Assume that the number of samples is the same for every line in each dataset 
   N.eqtl <- eqtl.data$N[1]

   res.all <- data.frame()

   for (i in 1:nrow(eqtl.data))  {
      pos.start <- eqtl.data[i, "POS"] - gap
      pos.end   <- eqtl.data[i, "POS"] + gap
      matches <- which(biom.data$POS > pos.start & biom.data$POS < pos.end )

      # There can be more than one probe per gene 
      region.eqtl <- read.csv(as.character(eqtl.data[i, "output.file"]) ,
                         sep='\t', stringsAsFactors = FALSE)
      region.biom <- biom.data[matches, ]
      # Construct the MAF data.frame, since we only want the shared snps, just do it for one dataset  
      maf.eqtl   <- data.frame(snp = region.eqtl$SNPID, maf = region.eqtl$F)

      # Loop over each biomarker 

      gene <- eqtl.data[i,"Gene.name"]
      probeID <- eqtl.data[i,"ProbeID"]

      message(i, " ",  gene, ": ", length(matches), " snps in biomarkers. From: ", pos.start, " To: ", pos.end)

      for (j in 1:length(biom.names)) {
         colname.pval <-  paste("PVAL.", biom.names[j], sep = "")
         colname.N    <-  paste("N.", biom.names[j], sep = "")
         colname.beta <-  paste("BETA.", biom.names[j], sep = "")
         colname.se   <-  paste("SE.", biom.names[j], sep = "")

         # If have OR instead of BETA (but this is for the whole dataset...correct later!)
         #if (type=="cc") (colname.beta <-  paste("OR.", biom.names[j], sep = ""))
         if (type=="cc") { 
                 colname.Ncases <-  paste("Ncases.", biom.names[j], sep = "")
                 s1 = region.biom[1,colname.Ncases]/region.biom[1,colname.N]
                 } else {
                 s1 = 0.5
                 }
         # Subset biomarker file by biom (this could be made faster with use of list)
         merged.data <- merge(region.biom[, c("SNPID", colname.pval, colname.N, colname.beta, colname.se)], region.eqtl, by = "SNPID")
         nsnps = nrow(merged.data)

         if (nsnps <= 2 ) ("There are not enough common snps in the region")
         if (nsnps > 2 ) {
         # Remove the NA, assuming they are in the beta columns 
         merged.data <- merged.data[!is.na(merged.data$beta),]
         merged.data <- merged.data[!is.na(merged.data[,colname.beta]),]
 
         # For now run with p-values (better for cc data)
         # dataset.biom = list(snp = merged.data$SNPID, beta = merged.data[, colname.beta],varbeta = merged.data[,colname.se]^2,
         dataset.biom = list(snp = merged.data$SNPID, pvalues = merged.data[, colname.pval],
                         N = merged.data[1,colname.N], s=s1, type = "quant")
         dataset.biom$MAF <-  maf.eqtl[match(merged.data$SNPID, maf.eqtl$snp ) ,"maf"]

         # dataset.eqtl = list(snp = merged.data$SNPID, beta = merged.data$beta, varbeta = merged.data$se.beta^2,
         dataset.eqtl = list(snp = merged.data$SNPID, pvalues = merged.data$PVAL,
                          N = N.eqtl, type = "quant")
         dataset.eqtl$MAF <-  maf.eqtl[match(merged.data$SNPID, maf.eqtl$snp ) ,"maf"]

         suppressMessages(capture.output(coloc.res <- coloc.abf(dataset.biom, dataset.eqtl, p12 = p12)))
         pp3       <- as.numeric(coloc.res$summary[5])
         pp4       <- as.numeric(coloc.res$summary[6])         
         snp.biom <- merged.data[which.min(merged.data[,colname.pval]), "SNPID"]
         snp.eqtl <- merged.data[which.min(merged.data$PVAL), "SNPID"]
         min.pval.biom <- min(merged.data[,colname.pval])
         min.pval.eqtl <- min(merged.data$PVAL)
         best.causal = as.character(coloc.res$results$snp[which.max(coloc.res$results$SNP.PP.H4)])

         res.temp = data.frame(ProbeID = probeID, Chr = chr.name, Gene.name = gene, biom = biom.names[j], pos.start=pos.start, pos.end=pos.end, snp.biom=snp.biom, snp.eqtl=snp.eqtl, min.pval.biom=min.pval.biom, min.pval.eqtl=min.pval.eqtl, best.causal=best.causal, pp3, pp4, files=NA)
         

         ############# PLOT
         if (plot & pp4 > 0.5 & nsnps > 2) {

                pvalue_BF_df = as.data.frame(coloc.res[2])            
                region_name <- paste(gene, '.', probeID,'.', biom.names[j], sep= '')
                pvalue_BF_file <- paste(pval.fld, 'pval_', region_name, '.txt', sep="")

                ### LocusZoom arguments:
	        pvalue_BF_df$chr = chr.name
	        pvalue_BF_df$pos = merged.data[match(pvalue_BF_df$results.snp, merged.data$SNPID), "POS"]
	        pvalue_BF_df$locus_snp <- paste(pvalue_BF_df$chr, pvalue_BF_df$pos, sep=":")
	        pvalue_BF_df$locus_snp <- paste("chr", pvalue_BF_df$locus_snp, sep="")
                # INDELS FORMAT FOR LOCUSZOOM: chr1:117930794:AAG_A (not rsid)
                pvalue_BF_df$locus_snp <- ifelse(grepl("*[:][:]*", pvalue_BF_df$results.snp), paste("chr", as.character(pvalue_BF_df$results.snp), sep=""), as.character(pvalue_BF_df$locus_snp)) 
                pvalue_BF_df$results.pvalues.df1 = merged.data[match(pvalue_BF_df$results.snp, merged.data$SNPID), colname.pval]
                pvalue_BF_df$results.pvalues.df2 = merged.data[match(pvalue_BF_df$results.snp, merged.data$SNPID), "PVAL"]
                image.biom = paste(plot.fld, "/", region_name, '_df1.pdf', sep='')
                image.eqtl = paste(plot.fld, "/", region_name, '_df2.pdf', sep='')

                write.table(x =  pvalue_BF_df , file = pvalue_BF_file, row.names = FALSE, quote = FALSE, sep = '\t')
             

                message('Output pdf for biomarker: ', image.biom)
                pdf(image.biom, width = 9, height = 9)
                #output of region_ld.ld is in /SAN/biomed/biomed14/vyp-scratch/vincent/eQTLs/ ?
                # If INDEL ALLELES do not match exactly (for ex. are reversed from the reference EUR files in here /cluster/project8/vyp/vincent/toolsVarious/locuszoom/EUR/), skip for now:
                plotted = tryCatch(locuszoom.ugi(metal = pvalue_BF_file,
                  refSnp = pvalue_BF_df[pvalue_BF_df$results.snp==best.causal,"locus_snp"] , #rs10877835
                  title = 'A',
                  pvalCol='results.pvalues.df1',
                  legend = 'left',
                  markerCol='locus_snp',
                  ylab = '-log10( biomarker P-value )',
                  chrCol= 'chr',
                  posCol = 'pos',
                  chr = chr.name,
                  showGenes = TRUE,
                  show_xlab=FALSE,
                  temp.file.code = region_name,
                  start= pos.start ,
                  end = pos.end
                  ), error=function(e) NULL )
                dev.off()
                # If the plot is empty unlink:
                if (is.null(plotted)) unlink(image.biom)

        	#do.call(file.remove,list(list.files(wd, pattern="region_ld")))
	        #file.remove(paste(wd, region_name, "_df1", ".log", sep=""))
                message('Output pdf for eQTL: ', image.eqtl)
                pdf(image.eqtl, width = 9, height = 9)
                plotted = tryCatch(locuszoom.ugi(metal = pvalue_BF_file,
                  refSnp = pvalue_BF_df[pvalue_BF_df$results.snp==best.causal,"locus_snp"] , #rs10877835
                  title = 'B',
                  pvalCol='results.pvalues.df2',
                  legend = 'left',
                  markerCol='locus_snp',
                  ylab = '-log10( expression P-value )',
                  chrCol= 'chr',
                  posCol = 'pos',
                  chr = chr.name,
                  showGenes = TRUE,
                  show_xlab=FALSE,
                  temp.file.code = region_name,
                  start= pos.start ,
                  end = pos.end
                  ), error=function(e) NULL )
                dev.off()
                # If the plot is empty unlink:
                if (is.null(plotted)) unlink(image.biom)

	        #do.call(file.remove,list(list.files(wd, pattern="region_ld")))
	        #file.remove(paste(wd, region_name, "_df2", ".log", sep=""))

	        # Now merge the two pdfs??
        	#system(paste("convert ", region_name, "_df1.pdf", " ", region_name, "_df2.pdf", " ", region_name, ".pdf", sep=""))
	        #file.remove(paste(wd, region_name, "_df1.pdf", sep=""))
	        #file.remove(paste(wd, region_name, "_df2.pdf", sep=""))

                 res.temp$files= paste(as.character(paste(plot.fld, region_name, "_df1.pdf", sep="")), as.character(paste(plot.fld, region_name, "_df2.pdf", sep="")))
         }

         res.all <- rbind(res.all, res.temp)

      }
     }
   }

   res.all <- data.frame(res.all)
   res.all$ProbeID <- as.character(res.all$ProbeID)
   res.all$biom <- as.character(res.all$biom)
   res.all$Gene.name <- as.character(res.all$Gene.name)
   res.all$snp.eqtl <- as.character(res.all$snp.eqtl)
   res.all$best.causal <- as.character(res.all$best.causal)
   res.all$files <- as.character(res.all$files)



   res.all <- res.all[with(res.all, order(pp4, decreasing=T)),]
   outfname = paste(outfolder, prefix, '_summary.tab', sep='')
   write.table(x =  res.all , file = outfname, row.names = FALSE, quote = FALSE, sep = '\t')


}
