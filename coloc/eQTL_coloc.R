# R Script for running colocalisation analysis of two eQTL datasets 
# Written by Kitty Lo 
# On 9th June 2014 
##########################

coloc.eqtl.eqtl <- function(eqtl.f1, eqtl.f2, p12, doPlot, results.folder, match.by = "probe") { 

   ### HACK to make it work for the moment 
   #   For the dexa files, set eqtl.folder to empty 
   eqtl.folder = "/cluster/project8/vyp/eQTL_integration/"
   #eqtl.folder = ""
   suppressPackageStartupMessages({
     require(VPlib);
     require(coloc);
   });

   if ( doPlot ) { 
      plot.fld = paste(results.folder, "/plot/", sep="")
      pval.fld= paste(results.folder, "/pval/", sep="")
      if( !file.exists(plot.fld)) dir.create(file.path(plot.fld), showWarnings = TRUE)
      if( !file.exists(pval.fld)) dir.create(file.path(pval.fld), showWarnings = TRUE)
   } 

   data1 <- read.csv(eqtl.f1, sep =',', stringsAsFactors = FALSE) 
   data2 <- read.csv(eqtl.f2, sep =',', stringsAsFactors = FALSE)


   # Work out what is the shared list of probes in the two datasets 
   d1.probes <- unique(data1$ProbeID) 
   d2.probes <- unique(data2$ProbeID) 
   common.probes <- d1.probes[which(d1.probes %in% d2.probes)] 

   d1.genes <- unique(data1$Gene.name) 
   d2.genes <- unique(data2$Gene.name) 
   common.genes <- d1.genes[which(d1.genes %in% d2.genes)] 

   print(common.genes) 
   if (length(common.probes) == 0 & match.by == "probe") { 
      stop("No common probes found between the two datasets") 
   }
 
   if (length(common.genes) == 0 & match.by == "gene") { 
      stop("No common genes found between the two datasets") 
   } 

   # Assume that the number of samples is the same for every line in each dataset 
   N1 <- data1$N[1] 
   N2 <- data2$N[1] 

   res.all <- c() 

   if (match.by == "probe") { 
       n.matches <- length(common.probes) 
   }
 
   if (match.by == "gene") { 
       n.matches <- length(common.genes) 
   } 

   for (i in 1:n.matches)  { 
      if (match.by == "probe") { 
         match.1 <- which(data1$ProbeID == common.probes[i]) 
         match.2 <- which(data2$ProbeID == common.probes[i])
      }  
      if (match.by == "gene") { 
         match.1 <- which(data1$Gene.name == common.genes[i]) 
         match.2 <- which(data2$Gene.name == common.genes[i])
      }  


      chr.name <- as.character(data1[match.1[1], "CHR"])  
      region.1 <- read.csv(paste(eqtl.folder, as.character(data1[match.1[1], "output.file"]), sep ='')  ,
                         sep='\t', stringsAsFactors = FALSE) 
      region.2 <- read.csv(paste(eqtl.folder, as.character(data2[match.2[1], "output.file"]), sep='') ,
                         sep='\t', stringsAsFactors = FALSE) 

      pos.start <- min(min(region.1$POS), min(region.2$POS)) 
      pos.end <- max(max(region.1$POS), max(region.2$POS))  
      # Filter by MAF 
      region.1 <- subset(region.1, region.1$F > 0.03 & region.1$F < 0.97 )  
      region.2 <- subset(region.2, region.2$F > 0.03 & region.2$F < 0.97 )  

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

      pp3       <- coloc.res$summary[5]
      pp4       <- coloc.res$summary[6]
      snp.eqtl1 <- merged.data[which.min(merged.data$PVAL.x), "SNPID"]
      snp.eqtl2 <- merged.data[which.min(merged.data$PVAL.y), "SNPID"]
      min.pval.eqtl1 <- min(merged.data$PVAL.x)
      min.pval.eqtl2 <- min(merged.data$PVAL.y)
      best.causal = as.character(coloc.res$results$snp[which.max(coloc.res$results$SNP.PP.H4)])

      if (match.by == "probe") { 
          gene  <- region.1[which(region.1$ProbeID == common.probes[i]), "Gene.name"]
          gene.used <- gene[1] 
          probe.used <- common.probes[i] 
      } 
      if (match.by == "gene") { 
          gene  <- data1[match.1[1], "Gene.name"] 
          probe <- data1[match.1[1], "ProbeID"] 
          gene.used <- gene[1] 
          probe.used <- probe[i] 
      } 

      message("Gene: ", gene.used, ", Probe: ", probe.used) 
      res.temp = data.frame(ProbeID = probe.used, Gene.name = gene.used, pos.start=pos.start, 
                            pos.end=pos.end, snp.eqtl1 = snp.eqtl1, snp.eqtl2=snp.eqtl2, 
                            min.pval.eqtl1=min.pval.eqtl1, min.pval.eqtl2=min.pval.eqtl2, 
                            best.causal=best.causal, pp3, pp4, files = NA) 
      # Make a plot if the pp4 value is high enough  
      # Do it per probe for each of the two dataset 

      if (pp4 < 0.9 & doPlot) { 
           pvalue_BF_df = as.data.frame(coloc.res[2])
           region_name <- paste(gene.used, '.', probe.used, sep= '')
           pvalue_BF_file <- paste(pval.fld, 'pval_', region_name, '.txt', sep="")

###############################
## Plot for first eQTL dataset 
           ### LocusZoom arguments:
           pvalue_BF_df$chr = chr.name
           pvalue_BF_df$pos = merged.data[match(pvalue_BF_df$results.snp, merged.data$SNPID), "POS.x"]
           pvalue_BF_df$locus_snp <- paste(pvalue_BF_df$chr, pvalue_BF_df$pos, sep=":")
           pvalue_BF_df$locus_snp <- paste("chr", pvalue_BF_df$locus_snp, sep="")
           pvalue_BF_df$results.pvalues.df1 = merged.data[match(pvalue_BF_df$results.snp, merged.data$SNPID), "PVAL.x"]
           pvalue_BF_df$results.pvalues.df2 = merged.data[match(pvalue_BF_df$results.snp, merged.data$SNPID), "PVAL.y"]
           image.eqtl1 = paste(plot.fld, "/", region_name, '_df1.pdf', sep='')
           image.eqtl2 = paste(plot.fld, "/", region_name, '_df2.pdf', sep='')

           print(pvalue_BF_file)  
           print(pvalue_BF_df)  
           write.table(x =  pvalue_BF_df , file = pvalue_BF_file, row.names = FALSE, quote = FALSE, sep = '\t')

           message('Output pdf for eQTL 1: ', image.eqtl1)
           pdf(image.eqtl1, width = 9, height = 9)
           tryCatch({ 
                locuszoom.ugi(metal = pvalue_BF_file,
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
                  ) 
           }, error = function(e) {
                  print("something wrong with locuszoom")  

           }) 
           dev.off()

###############################
## Plot for second eQTL dataset 
           message('Output pdf for eQTL 2: ', image.eqtl2)
           pdf(image.eqtl2, width = 9, height = 9)
           tryCatch({ 
           locuszoom.ugi(metal = pvalue_BF_file,
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
                  )

           }, error = function(e) {
                  print("something wrong with locuszoom")  

           }) 
           dev.off()
           res.temp$files= paste(as.character(paste(plot.fld, region_name, "_df1.pdf", sep="")), 
                                 as.character(paste(plot.fld, region_name, "_df2.pdf", sep="")))
           

      } # End plotting 

      res.all <- rbind(res.all, res.temp) 
   } # End probe 

   res.all <- data.frame(res.all) 
   res.all <- res.all[with(res.all, order(pp4, decreasing=T)),]
   outfname = paste(results.folder, 'summary.tab', sep='')
   write.table(x =  res.all , file = outfname, row.names = FALSE, quote = FALSE, sep = '\t')
 
} 
