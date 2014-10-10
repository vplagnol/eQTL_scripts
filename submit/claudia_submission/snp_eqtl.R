# Claudia Giambartolomei 08/09/2014
# To make it faster I should load the genotype file once only and use data.frame in the function; I would have to know the chromosome beforehand in that case (since genotype is split by chr)
# Tried to make consistent with /cluster/project8/vyp/claudia/scripts/eQTL_scripts/Pickrell/create_Pickrell_input.R
  snp.assoc <- function (dataset, condition,
                         expression.file,
                         rsid=NA, SNP.chr=NA, SNP.position=NA, min.MAF = 0.03, imp.quality = 0.3, gap = 500000,
                         base.folder = '/cluster/project8/vyp/eQTL_integration/data/') {

      suppressPackageStartupMessages({
      require(snpStats);
      require(GenomicRanges);
      });

      # Find position of SNP if it is rsid
      if (is.na(SNP.position)) {
      library(biomaRt)
      snpmart = useMart("snp", dataset="hsapiens_snp")
      res = getBM(c("refsnp_id","chr_name","chrom_start"), values=rsid, filters="snp_filter", mart=snpmart)
      SNP.chr = res$chr_name
      SNP.position = res$chrom_start
      }

      # Load genotypes for the selected snp, if it exists
      geno = paste(base.folder, dataset, "/genotypes/chr", SNP.chr, sep="")
      message("Genotype file: ", geno)
      if(!file.exists(geno)) stop('Genotype file is not specified correctly')
      load(geno)
      if (nrow(genotypes$map[which(genotypes$map$chromosome==SNP.chr & genotypes$map$position==SNP.position),])==0) stop("SNP ", SNP.chr, ":", SNP.position, " not covered by genotype data for ", dataset)
      n.snp = which(genotypes$map$chromosome==SNP.chr & genotypes$map$position==SNP.position)
      genotypes$map <- genotypes$map[n.snp,]
      genotypes$genotypes <- genotypes$genotypes[,n.snp]

      # Load expression for the selected tissue 
      load(expression.file) 
      support = get(paste("support.", condition, sep="")) 
      # Subset by chromosome  
      # First need to find correct column names in expression data: ex. 'gene.chromosome' for liver, 'chromosome_name' for GTex, seqname (with chr in front) for probeset brain
      #support <- subset(support, chr.col.expr == SNP.chr)
      chr.col.expr = grep( names(support), pattern = 'chr')
      # UKBEC core has chr name as seqname
      if (dataset=="brain_UKBEC" && length(chr.col.expr)==0) (chr.col.expr = grep( names(support), pattern = 'seqname'))
      support[,chr.col.expr] = gsub("chr", "", support[,chr.col.expr])
      #Same for start/end of probes: gene.position.start/gene.position.end for liver; start_position/end_position for GTex; start/stop for probeset brain
      start.col.expr = grep( names(support), pattern = 'start')
      stop.col.expr = grep( names(support), pattern = 'stop') 
      if (length(stop.col.expr)==0) stop.col.expr = grep( names(support), pattern = 'end')
      support <- subset(support, support[,chr.col.expr] == SNP.chr)
      # Expression data   
      expression <- get(condition) 

      ####
      shared.samples <- intersect(dimnames(expression)[[2]], dimnames(genotypes$genotypes)[[1]])
      n.samples <- length(shared.samples)
      expression <- expression[, shared.samples]
      genotypes$genotypes <- genotypes$genotypes[ shared.samples, ]
      genotypes$fam <- genotypes$genotypes[ shared.samples, ]
      message('Nb of samples shared by expression and genotypes: ', n.samples)

      #### Filters applied to genetic data
      # Filter by imputation quality if column exists 
      my.sum <- col.summary(genotypes$genotypes)
      # v.good.SNPs <- (my.sum$MAF > min.MAF) & (my.sum$Certain.calls > min.certain.calls) &
      imp = ifelse(length(grep( names(genotypes$map), pattern = '^info\\.', value = TRUE))>0 && genotypes$map[,grep( names(genotypes$map), pattern = '^info\\.', value = TRUE)] > imp.quality, TRUE, FALSE)
      if (imp) (v.good.SNPs <- (my.sum$MAF > min.MAF) && imp) else (v.good.SNPs <- (my.sum$MAF > min.MAF))
      if (!v.good.SNPs) stop("SNP does not pass filters for good quality")
      #genotypes$genotypes <- genotypes$genotypes[, v.good.SNPs ]
      #genotypes$map<- genotypes$map[ v.good.SNPs, ]
      #my.sum <- my.sum[ v.good.SNPs, ]

      ####
      expression.frame <- data.frame(samples = shared.samples, exp = NA)
      row.names(expression.frame) <- shared.samples

      #### Now add covariates if there is such a file
      covariates.tab <- paste(base.folder, dataset, '/covariates/covariates_', condition, '.tab', sep = '')  ##look for condition specific covariates first
      if (!file.exists(covariates.tab)) {covariates.tab <- paste(base.folder, dataset, '/covariates/covariates.tab', sep = '')} ##otherwise the generic one

      if (file.exists(covariates.tab)) {
      covariates <- read.table(covariates.tab, header = TRUE, sep = '\t')
      row.names(covariates) <- covariates$id
      covariates <- covariates[ shared.samples, ]
      covar.labels <- subset ( names(covariates), names(covariates) != 'id')
      my.formula <- paste('exp ~ ', paste(covar.labels, collapse = ' + '))
      for (covar.loc in covar.labels) {expression.frame[, covar.loc] <- covariates[, covar.loc]}
      } else {
      my.formula <- 'exp ~ 1'
      }

      # Find probe regions (+/- gap) that cross this SNP
      my.SNP <- GRanges(seqnames = SNP.chr,
                              IRanges(start = SNP.position, end= SNP.position + 1))


      my.QTL.GRanges <- GRanges(seqnames = support[,chr.col.expr],
                              IRanges(start= support[,start.col.expr]-gap,end= support[,stop.col.expr]+gap))

      my.overlap <- findOverlaps(query = my.SNP, subject = my.QTL.GRanges)  ##this is where the efficiency happens
      regions.of.interest <- unique(my.overlap@subjectHits)
      support.loc <-support[ regions.of.interest, ]

      # Loop over each ProbeID
      if (nrow(support.loc)==0) stop("There are no Probes overlapping SNP of interest in expression data", dataset)

      res = data.frame()

      for (i in 1:nrow(support.loc)) {

       Gene.name <- support.loc[i,"Gene.name"]
       ProbeID <- support.loc[i,"ProbeID"]
       ensemblID <- support.loc[i, "ensemblID"]
       if (length(grep( names(support.loc), pattern = 'mean.expr', value = TRUE))>0) (mean.expr = support.loc[i,grep( names(support.loc), pattern = 'mean.expr', value = TRUE)])
 
       message(i, " ",  Gene.name, " ", ProbeID)

       #######
       loc.expression <- expression[ as.character(ProbeID), ] ### take the row of the expression data with the matching probe ID
       expression.frame$exp <- as.numeric(loc.expression)/sd(as.numeric(loc.expression))

       ####### Now we can properly compute the P-values
       message('Formula being used: ', my.formula)
       print(head(expression.frame))
       my.tests <- snp.rhs.estimates( snp.data = genotypes$genotypes, data = expression.frame, formula = as.formula (my.formula), family = 'gaussian')

       effect <- sapply(as(my.tests, 'list'), FUN = function(my.tests) {if (is.null(my.tests)) {return (NA);} else {return((my.tests)$beta)} } )
       SE <- sapply(as(my.tests, 'list'), FUN = function(my.tests) { if (is.null(my.tests)) {return (NA);} else {return(sqrt((my.tests)$Var.beta))} } )

       res.i <- data.frame (SNPID = dimnames(genotypes$genotypes)[[2]],
                     Z = signif(effect/SE, 3),
                     F = signif(my.sum$RAF, 3),
                     PVAL = signif(2*pnorm(-abs(effect/SE)), 3),
                     POS = genotypes$map$position,
                     CHR = SNP.chr,
                     N = sum(!is.na(expression.frame$exp)),
                     beta = signif(effect, 4),
                     se.beta = signif(SE, 4),
                     ensemblID = ensemblID,
                     Gene.name = Gene.name,
                     ProbeID = ProbeID,
                     dataset = dataset,
                     condition = condition)

       #### add mean.expression if it exists
       if (exists('mean.expr')) {
       res.i$mean.expr = mean.expr
       } else { 
       res.i$mean.expr = NA
       }

       res = rbind(res, res.i)
       }
       return(res)
     }

