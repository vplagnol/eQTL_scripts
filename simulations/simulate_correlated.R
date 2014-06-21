
library(snpStats)
source('scripts/transeQTL_scripts/tools.R')

load('data/monocytes_Knight/genotypes/chr9')
load('data/monocytes_Knight/modules/LPS24_all_modules.RData')
load('data/monocytes_Knight/expression_data/expression_LPS24.RData')


shared.samples <- intersect(dimnames(genotypes$genotypes)[[1]], dimnames(LPS24)[[2]])
genotypes$genotypes <- genotypes$genotypes[ shared.samples, ]
LPS24 <- LPS24[, shared.samples]
my.SNP.chr9 <- as.numeric((genotypes$genotypes[, 'rs3898946'])) - 1


my.probes <- list.modules[[949]]$ProbeID
my.genes <- list.modules[[949]]$Gene.name
my.expression <- LPS24[ my.probes, ]



my.cov <- cov(t(my.expression))
n.samples <- ncol( my.expression )
n.probes <- nrow(my.expression)
my.means <- apply(my.expression, MAR = 1, FUN = mean)

library(mvtnorm)
my.like <- dmvnorm(x = t(my.expression), mean = my.means, sigma = my.cov, log=TRUE)

my.sum.like <- sum(my.like)


####### stepwise regression
print(summary(lm (my.expression['ILMN_1684576',] ~ geno)))
print(summary(lm (my.expression['ILMN_1684576',] ~ my.expression[1,] + geno))) ### surprising to see geno not contribute
geno <- my.SNP.chr9
dframe <- data.frame(samples = dimnames(my.expression)[[2]],
                     nID = as.numeric(gsub(pattern = '^Fairfax_', replacement = '', dimnames(my.expression)[[2]])))
dframe$covar.vector <- ifelse (dframe$nID %in% 1:288, 1, 2)


my.res <- stepwise.analysis (my.expression, gene.names = my.genes,
                             covar.vector = dframe$covar,
                             genotypes = geno)

stop()
sim.many.geno <- function() {
  pval <- 1
  while (pval > 10^-4) {
    geno <- rbinom (size = 2, n = .n.samples, prob = 0.3)
    mod <- lm( my.expression[1,] ~ geno + dframe$covar.vector)
    pval <-  coef(summary(mod))['geno',4]
  }
  return(geno)
}
my.sim.geno <- sim.many.geno()

my.res <- stepwise.analysis (my.expression,
                             covar.vector = dframe$covar,
                             genotypes = my.sim.geno)


print(cor.test(my.geno, data$x))


print(summary(lm (data$x ~ my.geno)))
print(summary(lm (data$y ~ my.geno)))
print(summary(lm (data$x ~ data$y + my.geno)))
print(summary(lm (data$y ~ data$x + my.geno)))
