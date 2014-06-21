library(ggplot2)
library(snpStats)
library(gridExtra)

#snp.name <- 'rs2275888'
#chromosome <- '9'
#gene <- 'IFNB1'
#output.pdf <- 'IFNB1_steroids.pdf'

snp.name <- 'rs2846848'
chromosome <- '11'
gene <- 'BIRC3'
output.pdf <- 'BIRC3_steroids.pdf'

choice.sets <- list( LCL_dexamethasone_DiRienzo = c('treated', 'untreated'), WB_dexamethasone_DiRienzo = c('treated'))

pretty.names <- function(dataset, condition) {
  text <- NA
  
  if (dataset == 'LCL_dexamethasone_DiRienzo') {
    if (condition == 'treated') text <- 'LCL, Dex 8h'
    if (condition == 'untreated') text <- 'LCL, Resting'
  }

  if (dataset == 'WB_dexamethasone_DiRienzo') {
    if (condition == 'treated') text <- 'WB, Dex 8h'
    if (condition == 'untreated') text <- 'WB, Resting'
  }
  
  return(text)
}





all.plots <- list()
plot.nb <- 0

n.sets <- length(choice.sets)
for (i in 1:n.sets) {

  
  geno <- paste('data/', names(choice.sets)[i], '/genotypes/chr', chromosome, sep = '')
  message('Loading ', geno)
  load(geno)
  geno.num <- as(genotypes$genotypes[, snp.name], 'numeric')
  A1 <- genotypes$map[ snp.name, 'allele.1']
  A2 <- genotypes$map[ snp.name, 'allele.2']
  alleles <- c( paste(c(A1, A1), collapse = ''), paste(c(A1, A2), collapse = ''), paste(c(A2, A2), collapse = ''))
  geno <- alleles[ geno.num + 1 ]
  
  for (condition in choice.sets[[ i ]]) {

    input.expression <- paste('data/', names(choice.sets)[i], '/expression_data/expression_', condition, '.RData', sep = '')
    message('Condition ', condition, ', file: ', input.expression)
    load(input.expression)

    expr.data <- get(condition)
    expr.data <- expr.data[  which(support.treated$Gene.name == gene), ]

    if (class(expr.data) == 'matrix') {  ## if we have multiple probes
      expr.data <- expr.data[which.max(apply(expr.data, MAR  = 1, FUN = median)),]  ##take the highest expressing probe for now
    }

    
    my.data <- data.frame( geno = factor(geno), expr = as.numeric(expr.data))
    P.value <- cor.test(as.numeric(my.data$geno), my.data$expr)$p.value
    P.value <- signif(P.value, 3)
    
    message('P.value: ', P.value)
    
    my.plot <- qplot(data = my.data, x = geno, y = expr, xlab = snp.name, ylab = gene, geom=c("boxplot", "jitter")) +
      ggtitle ( pretty.names(names(choice.sets)[i], condition) )
    #my.plot <- my.plot + geom_text(data = NULL, aes(x = 2, y = min(my.data$expr) , label = paste('P =', P.value))) + theme_bw()
    my.plot <- my.plot + geom_text(data = NULL, aes(x = 2, y = min(my.data$expr) , label = 'is it pretty')) + theme_bw()
    
    plot.nb <- plot.nb + 1
    all.plots[[ plot.nb ]] <- my.plot 
 
  }
  
}

pdf(output.pdf, width = 9, height = 4)
if (plot.nb == 2)  {grid.arrange(all.plots[[1]], all.plots[[2]], nrow=1)}
if (plot.nb == 3)  {grid.arrange(all.plots[[1]], all.plots[[2]], all.plots[[3]], nrow=1)}
if (plot.nb == 4)  {grid.arrange(all.plots[[1]], all.plots[[2]], all.plots[[3]], all.plots[[4]], nrow=1)}
dev.off()
print(output.pdf)
