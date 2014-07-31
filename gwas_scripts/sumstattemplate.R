#######################################################################################
#Written by Jonathan Griffiths 10 July 2014
#Produces summary stats for continuous metabolomics data & genotypes
#Run with all chromosomes via accompanying bash script sumstatsubmit.sh
#column format SNPID, CHR, POS, then repeat {BETA.biom, SE.biom, PVAL.biom, N.biom} for all metabolites.
################################################################################

#specify project in bash file

#load packages

source("http://bioconductor.org/biocLite.R")
biocLite("snpStats", lib='/home/zcqsjgr/R/x86_64-unknown-linux-gnu-library/3.0/snpStats/libs')##### Is this necessary? Running off cluster, seemed to give problems without this step every time.
library(snpStats, lib.loc='/home/zcqsjgr/R/x86_64-unknown-linux-gnu-library/3.0/snpStats/libs')

#load data
filemet=paste0('/cluster/project8/vyp/eQTL_integration/data/', project, '/phenotypes/continuous.RData')
load(filemet)
filegeno=paste0('/cluster/project8/vyp/eQTL_integration/data/', project, '/genotypes/chr', chr)
load(filegeno)
#########################################################
geno=genotypes[[1]]
map=genotypes[[2]]

#prepare initial columns of data frame
snp.loc=map$position
snp.chr=map$chromosome
snp.id=map$snp.name

biomarker=data.frame(snp.id,snp.chr,snp.loc)
names(biomarker)=c('SNPID', 'CHR', 'POS')

column=4

#########################################################

#metabolite loop to construct table
#column format SNPID, CHR, POS, then repeat {BETA.biom, SE.biom, PVAL.biom, N.biom} for all metabolites.
#standardise metabolite table for for(counter) generalisation
for(met in 1:dim(continuous.pheno)[2]){
metabolite=colnames(continuous.pheno)[met]

results=snp.rhs.tests(continuous.pheno[,met]~1, family='gaussian', snp.data=geno,data=continuous.pheno)
estimates=snp.rhs.estimates(continuous.pheno[,met]~1, family='gaussian', snp.data=geno,data=continuous.pheno)

BETA.biom=sapply(as(estimates, 'list'), FUN = function(estimates) {if (is.null(estimates)) {return (NA);} else {return((estimates)$beta)} } )
SE.biom=sapply(as(estimates, 'list'), FUN = function(estimates) { if (is.null(estimates)) {return (NA);} else {return(sqrt((estimates)$Var.beta))} } )
PVAL.biom=p.value(results)
N.biom=sample.size(results)

biomarker[,column]=BETA.biom
names(biomarker)[column]=paste0('BETA.',metabolite)
column=column+1

biomarker[,column]=SE.biom
names(biomarker)[column]=paste0('SE.',metabolite)
column=column+1

biomarker[,column]=PVAL.biom
names(biomarker)[column]=paste0('PVAL.',metabolite)
column=column+1

biomarker[,column]=N.biom
names(biomarker)[column]=paste0('N.',metabolite)
column=column+1
}

save(biomarker, file=paste0('/cluster/project8/vyp/eQTL_integration/data/', project, '/summaryStats/', oFile))
