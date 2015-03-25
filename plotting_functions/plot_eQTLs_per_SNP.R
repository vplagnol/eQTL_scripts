library(ggplot2)
library(snpStats)
source('/cluster/project8/vyp/eQTL_integration/scripts/eQTLs_scripts/tools.R')


######### this function plots one SNP/gene correlation across multiple studies

#snp.name <- 'rs2275888'
#chromosome <- '9'
#gene <- 'IFNB1'
#output.pdf <- 'IFNB1_steroids.pdf'

#setwd("/cluster/project8/vyp/eQTL_integration")


plot.eQTL <- function (choice.sets, snp.name, chromosome, gene.names, output.pdf, root = '/cluster/project8/vyp/eQTL_integration', free.scale = FALSE,
                       ylim = NULL,
                       hline = c()) {

  all.plots <- list()
  plot.nb <- 0

  n.sets <- length(choice.sets)
  message('Number of sets: ', n.sets)
  
  ## Load data
  all.data <- NULL
  for (i in 1:n.sets) {
    
    geno <- paste(root, '/data/', names(choice.sets)[i], '/genotypes/chr', chromosome, sep = '')
    message('Loading ', geno)
    load(geno)

    if (! snp.name %in% dimnames(genotypes$genotypes)[[2]] ) {stop('Missing ', snp.name, ' from ', geno)}

    geno.num <- as(genotypes$genotypes[, snp.name], 'numeric')
    A1 <- as.character(genotypes$map[ snp.name, 'allele.1'])
    A2 <- as.character(genotypes$map[ snp.name, 'allele.2'])
    alleles <- c( paste(c(A1, A1), collapse = ''), paste(c(A1, A2), collapse = ''), paste(c(A2, A2), collapse = ''))
    geno <- as.character(alleles[ round(geno.num + 1) ])
    geno.names <- dimnames(genotypes$genotypes)[[1]]
    
    for (condition in choice.sets[[ i ]]) {
      
      input.expression <- paste(root, '/data/', names(choice.sets)[i], '/expression_data/expression_', condition, '.RData', sep = '')
      message('Condition ', condition, ', file: ', input.expression)
      load(input.expression)
      
      support <- get(paste('support', condition, sep = '.'))
      expr.data <- get(condition)

      if (class(expr.data) == 'data.frame') expr.data <- as.matrix( expr.data )  ### a bit of a hack, will need fixing
      
      print(gene.names  %in% support$Gene.name | gene.names %in% support$ensemblID)
      for (loc.gene in gene.names) {
        message('Gene ', loc.gene)
        good.rows <-  which(row.names(support) == loc.gene | support$Gene.name == loc.gene | support$ensemblID == loc.gene)
        expr.data.loc <- expr.data[  good.rows, , drop = FALSE ]
        support.loc <- support[good.rows, ]
        print(support.loc)
        
        if (nrow(expr.data.loc) > 1) {
          message('Multiple matching probes ', nrow(expr.data.loc))
          print(dimnames(expr.data.loc)[[1]])
        }

        expr.data.loc <- expr.data.loc[which.max(apply(expr.data.loc, MAR  = 1, FUN = median)),]  ##take the highest expressing probe for now
        pretty.type=pretty.names(names(choice.sets)[i], condition)

        ##browser()
############### Now making sure the labels are pretty and involve the gene name if we selected by probes
        if (!is.na( support.loc$Gene.name[1]) && loc.gene == row.names(support.loc)[1]  && loc.gene != support.loc$Gene.name[1]) {  ## if the selected ID does not match the gene name
          my.label <- paste( row.names(support.loc)[ 1 ], ' / ', support.loc$Gene.name[1], sep = '')
          my.clean.gene <- support.loc$Gene.name[1]
        } else {

          if (is.na(support.loc$Gene.name[1])) {
            my.label <- loc.gene
            my.clean.gene <- loc.gene
          } else {
            my.label <- support.loc$Gene.name[1]
            my.clean.gene <- support.loc$Gene.name[1]
          }
        }
        
######
        
        my.data <- data.frame ( samples = geno.names, geno = factor(geno), dataset=pretty.type, gene = loc.gene, label = my.label, clean.gene = my.clean.gene )
        my.data$expr <- expr.data.loc[ as.character(my.data$samples) ]
        #my.data <- data.frame( samples = names(expr.data.loc), geno = factor(geno), expr = as.numeric(expr.data.loc), dataset=pretty.type, gene = loc.gene)
        all.data <- rbind(all.data, my.data)
      }
    }
  }

  if (length(gene.names) > 1) {
###if multiple genes, p-values must be per gene
    all.data$type <- all.data$label
  } else {
###if a single gene, p-values must be per dataset
    all.data$type <- all.data$dataset
  }

  ## Calc P.values
  P.value <- rep(NA, length(unique(all.data$type)))
  names(P.value) <- as.character(unique(all.data$type))
  for ( t in 1:length(names(P.value))){
    my.dset <- names(P.value)[t]
    message(my.dset)
    #browser()
    P.value[t] <- cor.test(as.numeric(all.data$geno[as.character(all.data$type) == my.dset]), all.data$expr[all.data$type == my.dset])$p.value
    message(my.dset, ' ', P.value[t])
  }

  message('Computing the nb of facets')
  #browser()
                                        #long cut way to find number of facets
  len <- length(levels(factor(all.data$type)))
  vars <- data.frame(expand.grid(levels(factor(all.data$type))))
  colnames(vars) <- c('type')


  nrows <- 1
  if ( len > 4 ) nrows <- 2
  if ( len > 8 ) nrows <- 3
  if ( len > 20 ) nrows <- 4
  if ( len > 30 ) nrows <- 5
  if ( len > 40 ) nrows <- 6
  if ( len > 50 ) nrows <- 7
  if (len > 60) nrows <- 8
  ncols <- ceiling(len / nrows)
  
  #browser()
  message('Now plotting the data')
  my.plot <- ggplot(data=all.data, aes(x = factor(geno), y = expr), xlab = snp.name) +
    geom_jitter(colour='grey50') +
      geom_boxplot(fill = "grey80", position=position_dodge(1), alpha=0.5) +
        theme_bw() +
          xlab(snp.name)  + 
            ylab (paste(ifelse (length(gene.names) == 1, gene.names, 'Gene'), 'expression'))


  if (free.scale) {
    #browser()
    dat <- data.frame(x = rep(2, len), y = as.numeric(tapply( all.data$expr, IND = all.data$type, FUN = min, na.rm = TRUE)), type=vars, labs=paste("P-value=", signif(P.value, 3), sep=""))
    my.plot <- my.plot + facet_grid(. ~ type, scales = 'free_y') +
      facet_wrap (facets = ~ type, nrow = nrows, ncol = ncols, scales = 'free_y') +
        geom_text(aes(x, y, label=labs, group=NULL),data=dat)
  } else {
    dat <- data.frame(x = rep(2, len), y = min(all.data$expr, na.rm = TRUE), type=vars, labs=paste("P-value=", signif(P.value, 3), sep=""))
    if (!is.null(ylim)) dat$y <- min(ylim)
    
    my.plot <- my.plot + facet_grid(. ~ type) +
      facet_wrap (facets = ~ type, nrow = nrows, ncol = ncols) + 
        geom_text(aes(x, y, label=labs, group=NULL),data=dat)

    if (!is.null(ylim)) my.plot <- my.plot + ylim(ylim)
  }

  if (length(hline) > 0) {
    message('Now adding a horizontal line')
    for( h in hline) {
      my.plot <- my.plot +  geom_hline(aes(yintercept= 0))
    }
  }
  
  ggsave (filename = output.pdf, plot = my.plot, width = ncols*4, height = 4*nrows, limitsize = FALSE)
  print(output.pdf)
  return (all.data)

}
