# /cluster/project8/vyp/kitty/scripts/eQTL_coloc.R 
wd = "/cluster/project8/vyp/eQTL_integration/data/SCZ_Ripke/"
setwd(wd)

# OUTPUT
 out_folder = "coloc/"
  dir.create(file.path(out_folder), showWarnings = FALSE)
 # Create summary data frame that will be saved in summary.file
 res_all <- data.frame()
 # Write the removed SNPs in a .txt file
 removed_list <- data.frame(Marker_removed = character(), reason = character())

########
# Dataset 1 is the input dataset for now.
# Dataset 2 is the eQTL dataset for now.
######## Dataset 2 options: eQTL
 eqtl_file = "/cluster/project8/vyp/eQTL_integration/data/brain_UKBEC/eQTLs/fgwas/fgwas_core_HIPP_summary_eQTLs.csv"
 type.dataset2 = "quant"
 n_sample_eqtl = 122
 s_dataset2 = 0.5

######## Dataset 1 options: biom
 biomarker_file = "/cluster/project8/vyp/eQTL_integration/data/SCZ_Ripke/summaryStats/scz.swe.ripke.2013.results.txt"

 no_header = FALSE  # set as TRUE if input file does not contain a header row
 n_sample_biom = 11244 
 ncases_sample_biom = 5001 # 6243 controls

 #Sample size of dataset  #This is hard to define for a meta-analysis: all SNPs must have same N
 marker_col = 1
 p_col = 8

 info = "Beta and SE" # info = "p-values" if only have p-values (no SE/beta)
 
 beta_col = 6 # Column of dataset with effect if effect is used (otherwise use p-values only). 
 se_col = 7 # StdErr  # Column of dataset with effect if effect is used (otherwise use p-values only)
 
 # If have chromosome and positions in hg19, already formatted: set as TRUE and give columns:
 have_chr_pos = TRUE
 if (have_chr_pos) {
   chr_col = 2
   pos_col = 3
   }

 pp4_filter = 0  # Minimum posterior probability for common signal 
 p1 = 1e-04 #1e-04  # prior probability a SNP is associated with trait 1  #Don't usually want to change this  ## Biomarker!!
 p2 = 1e-04  # prior probability a SNP is associated with trait 2  #Don't usually want to change this  ## Expression!!
 p12 = 1e-06 # prior probability a SNP is associated with both traits #Don't usually want to change this

 maf_filter = 0.001 # 0.05  #MAF filter applied to datasets
 rsq_filter = 0.3 #Imputation quality filter applied to datasets

 no_indels = FALSE # set as TRUE if don't want to include INDELs in coloc analysis

#########################################################################
########  COLOC SCRIPT PART1: FORMAT BIOMARKER FILE
#########################################################################
suppressPackageStartupMessages({
  require(GenomicRanges);
  require(coloc);
});

if (ncases_sample_biom ==0) { 
   type.dataset1 = "quant"
   s_dataset1=0.5 ## This will be ignored since the type is "quant"
   } else {
   type.dataset1 = "cc"
   s_dataset1 = ncases_sample_biom/n_sample_biom }  #  s = proportion of individuals that are cases (cases / N)


  ########### biomarker input data:
  cat("Uploading dataset ....\n")
  nrec = as.numeric(system(paste('cat ', biomarker_file, ' | wc -l', sep=''), intern=TRUE))

  if (no_header) {
     nrows=nrec
     header=FALSE } else {
     nrows=nrec -1 
     header=TRUE}

  if(info =="p-values" & !have_chr_pos) {
     biomarker.data <- tryCatch(read.table(pipe(paste("awk '{print $", marker_col, ",$", p_col, "}' ", biomarker_file, sep="")),  colClasses = c("character", "numeric"), nrows= nrows, comment.char="", header=header, stringsAsFactors =FALSE, col.names=c("SNP", "pvalues.df1")), error=function(e) NULL )
     }
  if(info =="p-values" & have_chr_pos) {      # Some numbers in position come out with scientific notation: import as character
     biomarker.data <- tryCatch(read.table(pipe(paste("awk '{print $", marker_col, ",$", p_col, ",$", chr_col, ",$", pos_col, "}' ", biomarker_file, sep="")),  colClasses = c("character", "numeric", "numeric", "character"), nrows= nrows, comment.char="", header=header, stringsAsFactors =FALSE, col.names=c("SNP", "pvalues.df1", "chr_name", "chrom_start")), error=function(e) NULL )
     }
  if(info =="Beta and SE" & !have_chr_pos) {
     biomarker.data <- tryCatch(read.table(pipe(paste("awk '{print $", marker_col, ",$", beta_col, ",$", se_col, ",$", p_col, "}' ", biomarker_file, sep="")),  colClasses = c("character", "numeric", "numeric", "numeric"), nrows= nrows, comment.char="", header=header, stringsAsFactors =FALSE, col.names=c("SNP", "beta.dataset1", "SE.df1", "pvalues.df1")), error=function(e) NULL )
     }
  if(info =="Beta and SE" & have_chr_pos) {
     biomarker.data <- tryCatch(read.table(pipe(paste("awk '{print $", marker_col, ",$", beta_col, ",$", se_col, ",$", chr_col, ",$", pos_col, ",$", p_col, "}' ", biomarker_file, sep="")),  colClasses = c("character", "numeric", "numeric", "numeric", "character", "numeric"), nrows= nrows, comment.char="", header=header, stringsAsFactors =FALSE, col.names=c("SNP", "beta.dataset1", "SE.df1", "chr_name", "chrom_start", "pvalues.df1")), error=function(e) NULL )
     }


  if ( is.null(biomarker.data ) )  {
     cat("Enter correct column with pvalues.\n") 
  }

  # Remove missing data
  biomarker.data = biomarker.data[complete.cases(biomarker.data),]

 if(info =="p-values") {
    if (class(biomarker.data[,2]) != "numeric" | max(biomarker.data[,2], na.rm=TRUE) >1)  {
     cat("Enter correct column with pvalues.\n") 
     }
 }
  if(info =="Beta and SE") {
    if (class(biomarker.data[,2]) != "numeric" & class(biomarker.data[,3]) != "numeric")  {
     cat("Enter correct column with pvalues.\n") 
     }
  }

 
  # If have a "chr" in front in marker column: take out
  if (length(grep("^chr", biomarker.data[,1]))>0) ( biomarker.data[,1] = gsub("^chr", "", biomarker.data[,1]) )


  if (!(any(grep("^rs", biomarker.data[,1])) | any(grep("^[0-9]{1,2}[:][1-9][0-9]*$", biomarker.data[,1])))) {  ### Change this to match exact {num num : numbers...} instead of ":"? But could have this :GGATT
     cat("Enter correct column with names of marker.\n") 
  }


  # If already have position
  if (have_chr_pos) { 
     biomarker.data$chrom_start <- as.numeric(as.character(biomarker.data$chrom_start))
     names(biomarker.data)[1] = "input_name"
     biomarker.data$chr_pos= paste(biomarker.data$chr_name, biomarker.data$chrom_start, sep=":")
     if(info =="Beta and SE") (biomarker.data = biomarker.data[,c("chr_pos", "chr_name", "chrom_start", "input_name", "beta.dataset1", "SE.df1", "pvalues.df1")])
     if(info =="p-values") (biomarker.data = biomarker.data[,c("chr_pos", "chr_name", "chrom_start", "input_name", "pvalues.df1")])
     }


  ## This takes a long time! But these could also be multiallelic SNPs?
  # Check if indels are all of this format "5:137454915:GAC_G", then I can use this:
  biomarker.data$indels = ifelse(biomarker.data[,1] %in% grep("[^0-9]$", biomarker.data[,1], value = TRUE, perl=TRUE), TRUE, FALSE)

  # If all the values in marker column are chr:pos.   
  if (length(grep("^[0-9]{1,2}[:][1-9][0-9]*$", biomarker.data[,1])) == length(biomarker.data[biomarker.data$indels==FALSE,1]) ) ( format_chr_pos <- TRUE ) else ( format_chr_pos <- FALSE   )

  if (format_chr_pos & !have_chr_pos) ( names(biomarker.data)[1] <- "chr_pos" )
  if (!format_chr_pos) {
  # Match rs to 1000 genomes
       #input.marker.name = names(biomarker.data) [1]
       rsid.data = biomarker.data[grep("^rs", biomarker.data[,1]),1]
       require(biomaRt)
       max <- 100000
       x <- seq_along(rsid.data)
       rsid_splits <- split(rsid.data, ceiling(x/max))
       pos<- data.frame()
       for (i in 1:length(rsid_splits)) {
          print(paste('Split rsids number ', i, sep=""))
          snp.db <- useMart("snp", dataset="hsapiens_snp")
          #pos = getBM(c('refsnp_id','chr_name','chrom_start'), filters = c('snp_filter'), values = rsid.data, mart = snp.db )
          pos_split = getBM(c('refsnp_id','chr_name','chrom_start'), filters = c('snp_filter'), values = rsid_splits[i], mart = snp.db )
       
          # pos = getBM(c('refsnp_id','allele', 'chr_name','chrom_start'), filters = c('snp_filter'), values = rsid.data, mart = snp.db )
          # pos$alleles <- gsub("/", "_", pos$allele)
          # Add the alleles to chr:pos:al1_al2   ex. 1:111864465:TC_T
          # pos$chr_pos_alleles <- paste(pos$chr_name, pos$chrom_start, pos$alleles, sep=":")
          pos <- rbind(pos, pos_split)
          print(dim(pos))
       }

       pos.na <- pos$refsnp_id[which(is.na(pos$chrom_start))]  # Some results from bomart come out as 'NA'
       pos <- subset(pos, !pos$refsnp_id %in% pos.na)
       no.match.genome<- rsid.data[!rsid.data %in% pos$refsnp_id]
       # Write the removed SNPs in a .txt file
       if (length(no.match.genome)>0) {
          removed_list <- rbind(removed_list, data.frame(Marker_removed = c(no.match.genome,pos.na), reason = "Does not map to genome"))
       }
       # For now remove SNPs that map to >1 region:
       multiple.match <- names(which(table(pos$refsnp_id) > 1))
       # Write the removed SNPs in .txt file
       if (length(multiple.match)>0) {
           removed_list <- rbind(removed_list, data.frame(Marker_removed = multiple.match, reason = "Map to >1 pos in genome"))
       }
       # Take out the multiple matches
       pos <- subset(pos, !pos$refsnp_id %in% multiple.match)
       pos$chr_pos <- paste(pos$chr_name, pos$chrom_start, sep=":")
       # Remove these SNPs from original dataset
       to.remove <-  biomarker.data[,1] %in% (c(multiple.match,no.match.genome))
       biomarker.data <- biomarker.data[!to.remove,]
       biomarker.data <- merge(biomarker.data, pos[,c(1,4)], by.x=1, by.y="refsnp_id", sort = FALSE)
       # biomarker.data$chr_pos <- ifelse(is.na(biomarker.data$chr_pos), biomarker.data[,1], biomarker.data$chr_pos) 
       # biomarker.data$chr_pos <- ifelse(to.replace, pos$chr_pos_alleles, biomarker.data[,1]) 

       # Reorder so you have chr_pos in first column and p-value in second column
       #find.col.df1 = c(which(names(biomarker.data)=="chr_pos"), 2, which(names(biomarker.data)==input.marker.name))
       #biomarker.data <- biomarker.data[,find.col.df1]
       #This is if want to keep all columns in dataframes: Must change function if want this!
       #ncol.df1 = 1:ncol(biomarker.data)
       #biomarker.data <- biomarker.data[,c(find.col.df1, ncol.df1[-find.col.df1])]

       if(info =="p-values") {
       biomarker.data <- biomarker.data[,c(4,2,1, 3)]   ## If take out INDEL column must change this!! ### *******
       # Must either change the name of original biomarker data or find a very unique name for marker--not 'SNP' b/c it messes up the merging
       names(biomarker.data)[3] <- "input_name"
       }
       if(info =="Beta and SE") {
       biomarker.data <- biomarker.data[,c(5,2,3,1,4)]   ## If take out INDEL column must change this!! ### *******
       # Must either change the name of original biomarker data or find a very unique name for marker--not 'SNP' b/c it messes up the merging
       names(biomarker.data)[4] <- "input_name"
       }

       #!!!!!!!!!!!!!!################### MUST ADD THIS TO WEBSITE!!!!!!!!!!!!!!!!!!!!!
       # Add the SNPs that were already in chr_pos format:
       biomarker.data$chr_pos <- ifelse(is.na(biomarker.data$chr_pos), biomarker.data$input_name, biomarker.data$chr_pos)
       #!!!!!!!!!!!!!!################### MUST ADD THIS TO WEBSITE!!!!!!!!!!!!!!!!!!!!!

  }


### I don't think this is necessary since above I have the condition that chr_pos is only numbers  
       #if (format_chr_pos)  {
           sex.chr<- biomarker.data$chr_pos[grep("^[XxYy][:]", biomarker.data$chr_pos)] ## For now remove sex chromosomes (since don't have in expression data)
           other <- biomarker.data$chr_pos[grep("^[A-Za-z]", biomarker.data$chr_pos)] ## For now remove sex chromosomes (since don't have in expression data)
           if (no_indels) ( indels<- biomarker.data$chr_pos[biomarker.data$indels==TRUE] ) else ( indels<-as.character()) ## For now remove indels
           to.remove <-  biomarker.data$chr_pos %in% c(sex.chr,other,indels) #}
####
       if (length(sex.chr)>0) {
           removed_list <- rbind(removed_list, data.frame(Marker_removed = sex.chr, reason = "Map to sex chromosome: no match for this at the moment"))
       }
       if (length(other)>0) {
           removed_list <- rbind(removed_list, data.frame(Marker_removed = other, reason = "Map to MT chromosome?"))
       }
       if (length(indels)>0) {
           removed_list <- rbind(removed_list, data.frame(Marker_removed = indels, reason = "Looks like an indel, SNP with multiple alleles, or position is NA"))
       }


       biomarker.data <- biomarker.data[!to.remove,]

 if (!have_chr_pos) {
       biomarker.data$chrom_start <- gsub("^[0-9]{1,2}[:]", "", biomarker.data$chr_pos)
       if (!no_indels) (biomarker.data$chrom_start = gsub("[:].*", "", biomarker.data$chrom_start))
       biomarker.data$chr_name <- gsub("[:].*", "", biomarker.data$chr_pos)
       biomarker.data$chrom_start <-as.numeric(biomarker.data$chrom_start)
       biomarker.data$chr_name <-as.numeric(biomarker.data$chr_name )
       #biomarker.data$chr_name <- as.numeric(lapply(strsplit(as.character(biomarker.data$chr_pos), "\\:", perl = TRUE), "[", 1))
       #biomarker.data$chrom_start <- as.numeric(lapply(strsplit(as.character(biomarker.data$chr_pos), "\\:"), "[", 2))
      }

      # Remove data that is not chr:pos
      # biomarker.data = biomarker.data[complete.cases(biomarker.data),]


#########################################################################
########  COLOC SCRIPT PART2:  IMPORT REGION FILE AND EQTL FILE
#########################################################################
# Make eqtl region small and take only regions with overlap with biomarker small

gap=200000
genes = read.csv("/cluster/project8/vyp/vincent/Software/pipeline/RNASeq/bundle/human/biomart/biomart_annotations_human.tab", sep = '\t')

eqtl_file = "/cluster/project8/vyp/eQTL_integration/data/brain_UKBEC/eQTLs/fgwas/fgwas_core_HIPP_summary_eQTLs.csv"
eqtl_summary <- read.csv(eqtl_file, sep =',')


  ###********** These are the line to filter by p-value to decrease time ******
  # Find overlapping regions only for signals with small p-val
  pval_filter_biom <- 10^(-3)
  biomarker.small <- subset(biomarker.data, biomarker.data$pvalues.df1 < pval_filter_biom)


  my.BM.GRanges <- GRanges(seqnames = biomarker.small$chr_name,
                              IRanges(start = biomarker.small$chrom_start, end= biomarker.small$chrom_start + 1))  


  my.QTL.GRanges <- GRanges(seqnames = genes$chromosome_name,
                              IRanges(start= genes$start_position -gap, end= genes$end_position +gap))

  my.overlap <- findOverlaps(query = my.BM.GRanges, subject = my.QTL.GRanges)
  regions.of.interest <- unique(my.overlap@subjectHits)
  biom.genes <- genes[ regions.of.interest, ]


  # Work out what is the shared list of genes in the two datasets 
   eqtl_summary$biom <- ifelse(eqtl_summary$Gene.name %in% biom.genes$external_gene_id, TRUE, FALSE)
   common.genes = eqtl_summary[eqtl_summary$biom,]
   #common.genes <- eqtl_summary[which(eqtl_summary$Gene.name %in% biom.genes$external_gene_id),]
   common.genes$range_start = biom.genes$start_position[match(common.genes$Gene.name, biom.genes$external_gene_id)] - gap
   common.genes$range_end = biom.genes$end_position[match(common.genes$Gene.name, biom.genes$external_gene_id)] + gap


  ### split biomarker data by chromosome
  message('For maximum speed, split the biomarker data by chromosomeom first')
  my_split_list <- list()
  for (chr in  unique(common.genes$CHR)) {
    message('Chr ', chr)
    my_split_list[[ as.character(chr) ]] <- subset( biomarker.data, chr_name == chr)
  }
  message('Done with the splitting')


#########################################################################
########  COLOC SCRIPT PART3: Go over regions that overlap between datasets
#########################################################################


  list.snps <- character(0)

  for (eqtl.region in 1:nrow(common.genes)) {

      my.chr <- common.genes$CHR[ eqtl.region ]
      probe.id <- common.genes$ProbeID[  eqtl.region]
      gene <- common.genes$Gene.name[  eqtl.region ]

      print(paste('Probe ', probe.id, ' ', gene, ' eqtl.region ', eqtl.region, sep=''))

      my.range <- c(as.numeric(common.genes$range_start[ eqtl.region ]), as.numeric(common.genes$range_end[ eqtl.region ]))

      # Import the correct p-value file (either conditional or not)
      eqtl.file <- as.character(common.genes$output.file[ eqtl.region ])
      # res2.eqtl <- get(load(paste(coloc_folder, eqtl.file, sep="")))
      eqtl.data <- read.table(eqtl.file, header=T)
      loc.biomarker <- subset(my_split_list[[as.character(my.chr)]], chrom_start >=  common.genes$range_start[ eqtl.region ] & chrom_start < common.genes$range_end[ eqtl.region ])
      # loc.biomarker <- subset(biomarker.data, chrom_start >=  common.genes$start[ eqtl.region ] & chrom_start < common.genes$stop[ eqtl.region ])

      if(info =="Beta and SE") {
        loc.biomarker$varbeta.dataset1 <-  loc.biomarker$SE.df1^2
     }

      # Want SNP to be chr:pos
      eqtl.data$chr_pos = ifelse(grepl("^rs", eqtl.data$SNPID), paste(eqtl.data$CHR, eqtl.data$POS, sep=":"), as.character(eqtl.data$SNPID))

      ####   Restrict genotype files to MAF and Rsq specified ?
      # For now use MAF from expression dataset    
      eqtl.data$MAF <- ifelse(eqtl.data$F<=0.5, eqtl.data$F, 1-eqtl.data$F)
      # eqtl.data = subset(eqtl.data, MAF > maf_filter & Rsq.eqtl > rsq_filter)

      #Generalise the column names (so unlikely to be repeated): 
      names(eqtl.data)[8] <- "beta.dataset2"
      names(eqtl.data)[9] <- "SE.df2"

      merged.df = tryCatch(merge(loc.biomarker,eqtl.data, by="chr_pos"), error=function(e) NULL )

      if (info =="p-values") {
             dataset1 = list(snp = merged.df$chr_pos, pvalues =  merged.df$pvalues.df1,
                N = n_sample_biom, type = s_dataset1)
         }
      if (info =="Beta and SE") {
             dataset1 = list(snp = merged.df$chr_pos, beta =  merged.df$beta.dataset1,
                 varbeta =  merged.df$varbeta.dataset1, N = n_sample_biom, type = s_dataset1)
         }

      dataset2 = list(snp = merged.df$chr_pos, beta = merged.df$beta.dataset2,
               varbeta = merged.df$SE.df2^2, N = n_sample_eqtl, type = s_dataset2)


      coloc.res <- coloc.abf(dataset1, dataset2,p1 = p1, p2 = p2, p12 = p12, MAF=maf1[match(eqtl.data$SNP, maf1$snp ) ,"maf"])


      # Write list of SNPs in input data that do not match eqtl data:
      if (format_chr_pos)  (list.snps.loc <- loc.biomarker$SNP[!loc.biomarker$SNP %in% eqtl.data$SNP])
      if (!format_chr_pos)  (list.snps.loc <- loc.biomarker$input_name[!loc.biomarker$SNP %in% eqtl.data$SNP])
      list.snps <- c(list.snps, list.snps.loc )
      list.snps <- list.snps[!duplicated(list.snps)]

    res.all <- rbind(res.all, data.frame(EnsemblID = as.character(common.genes$Gene.name[eqtl.region]), ProbeID = as.character(common.genes$ProbeID[eqtl.region]), PP3= as.numeric(coloc.res$summary[5]), PP4 = as.numeric(coloc.res$summary[6])))

}

   res.all <- res.all[with(res.all, order(PP4)),]

   write.table(x = res.all, file = "results_summary.txt", row.names = FALSE, quote = FALSE, sep = '\t')


#########################################################################
########  COLOC SCRIPT PART4: PLOT?
#########################################################################

      # Need to add p-value for df1 if don't have for plot and to add snp with lowest p-value:
       if(info =="Beta and SE") {
           pvalue_BF_df$results.pvalues.df1 = 2*pnorm(-abs(pvalue_BF_df$results.z.df1))
      }


