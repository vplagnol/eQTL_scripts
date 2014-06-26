

stepwise.analysis <- function(expression.data, ##rows are probes, columns are samples, must have proper names
                              genotype.data,  ##snpStats object
                              snp.name,
                              base.formula = ' probe ~ 1 ',
                              verbose = FALSE,
                              covariates = NULL,
                              gene.names = NULL,
                              threshold = NA) {  ##genotype data

  library(snpStats)
  shared.samples <- intersect(dimnames(genotype.data)[[1]], dimnames(expression.data)[[2]])
  genotype.data <- genotype.data[ shared.samples, , drop = FALSE]
  expression.data <- expression.data[ , shared.samples]  ##we want non missing gene names and matched individual IDs
  
  my.probes <- dimnames(expression.data)[[1]]
  n.probes <- length(my.probes)
  if (is.na(threshold)) threshold <- 0.05/length(my.probes)

  message('Number of probes ', n.probes)
  
  dframe <- data.frame(samples = dimnames(expression.data)[[2]], stringsAsFactors = FALSE)
  if (!is.null(covariates)) {
    for (my.name in names(covar)) dframe[, my.name] <- covar[, my.name]
  }
  
  dframe$geno <- as.numeric((genotype.data[, snp.name])) - 1
  
  my.res <- data.frame( probe = my.probes, pval = NA)
  if (length(gene.names) == length(my.probes)) my.res$gene <- gene.names
  
  base.formula <- 'probe ~ 1'
  my.p <- 0
  step <- 0
  while (my.p < threshold) {
    if (verbose) print(base.formula)
    for (i in 1:n.probes) {
      probe <- my.probes[ i ]
      dframe$probe <- expression.data[ probe, ]
      
      mod.base <- lm (data = dframe, formula = base.formula, subset = !is.na(geno))
      mod.new <- lm (data = dframe, formula = paste(base.formula, ' + geno'))
      
      my.p <-  anova( mod.base, mod.new, test = 'F')[['Pr(>F)']][2]
      my.res$pval[ i ] <- coef(summary(mod.new))['geno',4]
      if (coef(mod.new)['geno'] < 0.00001)  my.res$pval[ i ] <- 1  ##dirty hack
      
    }
    if (verbose) print(my.res)
    my.p <- min( my.res$pval, na.rm = TRUE)
    message('Minimum p ', my.p)
    if (my.p < threshold) {
      step <- step + 1
      best.probe <- as.character(my.res$probe[ which.min(my.res$pval) ])
      if (verbose) {
        message('best probe is ', best.probe)
        print(my.res)
        #if (step == 8) browser()
      }
      
      dframe[, paste('X', best.probe, sep = '')] <- expression.data[ best.probe, ]
      base.formula <- paste(base.formula, ' +  X', best.probe, sep = '')
      base.model <- lm (data = dframe, formula = base.formula)
    }
  }
  return(list(formula = base.formula, nsteps = step))
}


correct.module.pvalues <- function (expression.data, genotype.data, covariates = NULL, gene.names, verbose = FALSE) {
  if (class(genotype.data) != 'matrix') stop('Class of genotype should be a matrix, even if there is a single column')
  if (!is.null(covariates)) {if (class(covariates) != 'data.frame') stop()}
  if (class(expression.data) != 'matrix') stop()

  if (verbose) {
    message('Running the module correction')
  }
  
  shared.samples <- intersect(dimnames(genotype.data)[[1]], dimnames(expression.data)[[2]])
  genotype.data.loc <- genotype.data[ shared.samples, , drop = FALSE]
  expression.data.loc <- expression.data[which (!is.na(gene.names)) , shared.samples]  ##we want non missing gene names and matched individual IDs
  gene.names <- gene.names [ which (!is.na(gene.names)) ]
  
  #pairwise.mat <- cor(t(expression.data.loc))
  my.unique.genes <- unique( gene.names )

  if (is.null(covariates)) {
    covariates.loc <- data.frame( sample = shared.samples, genotype = as.numeric(genotype.data.loc))
  } else {
    stop('Not implemented yet')
  }
  
  my.res <- data.frame ( Gene.name = my.unique.genes, uncorrected.pvalue = NA, corrected.pvalue = NA, stringsAsFactors = FALSE)
  for (i in 1:nrow(my.res)) {
    gene <- my.res$Gene.name [ i ]
    probes.for.that.gene <- dimnames(expression.data.loc)[[1]] [ which (gene.names == gene) ]
    covariates.loc$expression <- apply(MAR = 2, expression.data.loc[ probes.for.that.gene, , drop = FALSE], FUN = mean, na.rm = TRUE)
    my.res$uncorrected.pvalue[ i ] <- cor.test(covariates.loc$expression, covariates.loc$genotype)$p.value
  }

  #browser()

  most.correlated.gene.name <- my.res$Gene.name[ which.min(my.res$uncorrected.pvalue) ]
  probes.for.that.gene <- dimnames(expression.data.loc)[[1]] [ which (gene.names == most.correlated.gene.name) ]
  covariates.loc$most.cor <- apply(MAR = 2, expression.data.loc[ probes.for.that.gene, , drop = FALSE], FUN = mean, na.rm = TRUE)

  if (verbose) {
    message('Selected probes: ', probes.for.that.gene)
    print(sd(covariates.loc$most.cor))
  }
  print(head(my.res))
  for (i in 1:nrow(my.res)) {
    gene <- my.res$Gene.name [ i ]
    probes.for.that.gene <- dimnames(expression.data.loc)[[1]] [ which (gene.names == gene) ]
    covariates.loc$expression <- apply(MAR = 2, expression.data.loc[ probes.for.that.gene, , drop = FALSE], FUN = mean, na.rm = TRUE)

    #### Now the pairwise correlations to pick the most correlated probe
    #pairwise.loc <- pairwise.mat [probes.for.that.gene, ! dimnames(pairwise.mat)[[2]] %in% probes.for.that.gene, drop = FALSE]
    #gene.names.loc <- gene.names [ ! dimnames(pairwise.mat)[[2]] %in% probes.for.that.gene ]
    #pairwise.loc <- apply(pairwise.loc, MAR = 2, FUN = max)
    #most.correlated.gene.name <- gene.names.loc[ which.max( pairwise.loc ) ]
    #most.correlated <- names(which.max(pairwise.loc))
    #covariates.loc$most.cor <- expression.data.loc[ most.correlated , ]
    
    if (verbose) {
      message('Step 1:', gene, ' ', most.correlated.gene.name, ' ', signif(cor(covariates.loc$most.cor, covariates.loc$expression), 3))
      #if (gene == most.correlated.gene.name) browser()
    }
    
    mod1 <- lm (data = covariates.loc, formula = 'expression ~ most.cor + genotype')
    mod0 <- lm (data = covariates.loc, formula = 'expression ~ most.cor ', subset = !is.na(genotype))
    my.res$corrected.pvalue[ i ] <- anova(mod0, mod1)[[ 'Pr(>F)' ]][2]
  }

  return(my.res)
}

