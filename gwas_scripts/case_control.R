library(snpStats)
library(IRanges)

dataset <- 'celiac_vanHeel'
threshold <- 10^(-5)

GWAS.folder <- paste0('data/', dataset, '/GWAS')
if (!file.exists(GWAS.folder)) {dir.create(GWAS.folder)}

fgwas.folder <- paste0(GWAS.folder, '/fgwas_regions')
if (!file.exists(fgwas.folder)) {dir.create(fgwas.folder)}

###and now we clean the old files
file.remove(list.files( fgwas.folder, full.name = TRUE))

case.control.data <- paste0('data/', dataset, '/phenotypes/binary.RData')
load(case.control.data)
options(stringsAsFactors = FALSE)


final.fgwas.table <- data.frame()
first <- TRUE
true.region.number <- 1

for (chromosome in 22:1) {
  message('Chromosome ', chromosome)
  genotype.file <- paste0('data/celiac_vanHeel/genotypes/chr', chromosome)
  load(genotype.file)

  #### a key requirement is that the genotype file contains all the info so we first check that all the phenotypes row names are in teh genotype table
  missing.samples <-!(  row.names(binary.pheno) %in% row.names(genotypes$genotypes))
  if (sum(missing.samples) > 0) {stop('Some individuals in the phenotype table do not have a match in the genotype table')}


  case.control.pheno <- names(binary.pheno)

  for (pheno in case.control.pheno) {
    message('Now working with the column ', pheno)
    loc.geno <- genotypes$genotypes
    loc.map <- genotypes$map
    
    my.sum <- col.summary( loc.geno )
    
    if (first) {
      first <- FALSE
      my.table.pheno <- data.frame(samples = row.names(loc.geno) )
      row.names(my.table.pheno) <- my.table.pheno$samples
      #my.table.pheno[, 'cc' ] <- NA
      my.table.pheno[, 'cc' ] <- binary.pheno[ match( row.names(my.table.pheno), table = row.names(binary.pheno)) , pheno ]
    }
    
    my.formula <- 'cc ~ 1'
    my.tests <- snp.rhs.tests(snp.data = loc.geno, data = my.table.pheno, formula = as.formula (my.formula), score = TRUE, robust=TRUE, uncertain = TRUE, family = 'binomial')

    score= sapply(my.tests@score, "[[", 1)
    score.var= sapply(my.tests@score, "[[", 2)
    effect = score/score.var
    score.chisq = score^2/score.var
    pval <-  pchisq(score.chisq, df =1, lower.tail = FALSE)
    var.beta <- 1 / score.var
    SE <- sqrt(var.beta)

    res <- data.frame (SNPID = dimnames(loc.geno)[[2]],
                       Z = signif(effect/SE, 3),
                       F = signif(my.sum$RAF, 3),
                       PVAL = p.value ( my.tests),
                       POS = loc.map$position,
                       CHR = chromosome,
                       NCASE = sum( my.table.pheno[, 'cc' ] == 1, na.rm = TRUE),
                       NCONTROL = sum( my.table.pheno[, 'cc' ] == 2, na.rm = TRUE),
                       beta = signif(effect, 4),
                       se.beta = signif(SE, 4),
                       phenotype = pheno)
    


    ################### now split the whole thing into discrete regions
    min.p <- min(res$PVAL, na.rm = TRUE)
    best.SNP <- which.min( res$PVAL )
    best.SNP.name <- as.character(res$SNPID[ best.SNP ])
    region <- 1
    processed.file <- list()
    
    while (min.p < threshold) {
      
      r2.with.best.SNP <- ld (x = genotypes$genotypes, y = genotypes$genotypes[,best.SNP.name] , stats = 'R.squared')
      good.SNPs <- names(which(r2.with.best.SNP[,1] > 0.1))
      my.ld.range <- range(subset(res, SNPID %in% good.SNPs, POS, drop = TRUE), na.rm = TRUE)
      my.ld.range[1] <- my.ld.range[1] - 20000  ## extend a bit to be safe
      my.ld.range[2] <- my.ld.range[2] + 20000  ##extend on the other end

      ######## below some safety check to make sure that I do not go way too far
      my.ld.range[1] <- max( my.ld.range[1], res$POS[ best.SNP ] - 600000)  ##region no greater than 1.2 Mb basically
      my.ld.range[2] <- min( my.ld.range[2], res$POS[ best.SNP ] + 600000)

      message("Now analyzing region ", region, ' for chromosome ', chromosome, ' with min.p ', signif(min.p, 4), ', position ',  my.ld.range[1], '-',  my.ld.range[2])
      
      #####
      selected.SNPs <- res$POS > my.ld.range[ 1 ] & res$POS < my.ld.range[ 2 ]
      loc.res <- res[ selected.SNPs, ]
      processed.file[[ region ]] <- loc.res
      region <- region + 1

      
      ######### now we finalize the process
      res <- res[ ! selected.SNPs, ]
      min.p <- min(res$PVAL, na.rm = TRUE)
      best.SNP <- which.min( res$PVAL )
      best.SNP.name <- as.character(res$SNPID[ best.SNP ])      
    }
  }


  ######## now we need to merge regions that should be merged, otherwise it is messy
  n.regions.before <- length(processed.file)
  if (n.regions.before > 0) {  ## if there is at least one peak identified
    my.ranges <- data.frame(start = rep(NA, n.regions.before), end = rep(NA, n.regions.before))
    for (r in 1:n.regions.before) {
      my.range <- range ( processed.file[[ r ]]$POS )
      my.ranges$start[ r ] <- my.range[ 1 ]
      my.ranges$end[ r ] <- my.range[ 2 ]
    }
    
    irang <- IRanges(start = my.ranges$start - 50000, end = my.ranges$end + 50000)
    irang.red <- reduce(irang)
    overl <- findOverlaps (subject = irang.red, query = irang)
    
    for (r in 1:length(irang.red)) {
      regions.to.merge <- overl@queryHits[ overl@subjectHits == r  ]
      if (length(regions.to.merge) == 0) stop("Makes no sense")
      message("Number of regions to merge: ", length(regions.to.merge))
      combined.region <- data.frame()
      for (k in regions.to.merge) {combined.region <- rbind.data.frame( combined.region, processed.file[[ k ]]  )}
      combined.region <- combined.region[ order(combined.region$POS), ]

      
      ###############
      combined.region$SEGNUMBER <- paste0('chr', chromosome, '_region', true.region.number)
      final.fgwas.table <- rbind.data.frame( final.fgwas.table, combined.region)
      
      ### and now we print the merged region
      output.file <- paste0(fgwas.folder, '/GWAS_', pheno, '_chr', chromosome, '_region', true.region.number, '.tab')
      write.table(x = combined.region, file = output.file, sep = '\t', row.names = FALSE, quote = FALSE)
      true.region.number <- true.region.number + 1
    }
  }

  
  rm (list = c('res', 'best.SNP', 'min.p', 'output.file', 'region', 'loc.res', 'best.SNP.name', 'score', 'SE', 'my.range', 'my.ld.range'))
}

fgwas.file <- paste0( GWAS.folder, '/fgwas_', pheno, '_input.tab')
final.fgwas.table <- final.fgwas.table[ order(final.fgwas.table$CHR, final.fgwas.table$POS), ]
write.table(x = final.fgwas.table, file = fgwas.file, sep = '\t', row.names = FALSE, quote = FALSE)


my.cmd <- paste0('gzip -f ', fgwas.file)
system( my.cmd )
