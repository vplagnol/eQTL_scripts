#Written by Jonathan Griffiths 10 July 2014
#Produces summary stats for WHII metabolomics data & genotypes
#Looping for chromosomes provided by accompanying bash script sumstatsubmit.sh
################################################################################

source("http://bioconductor.org/biocLite.R")
biocLite("snpStats", lib='/home/zcqsjgr/R/x86_64-unknown-linux-gnu-library/3.0/snpStats/libs')
library(snpStats, lib.loc='/home/zcqsjgr/R/x86_64-unknown-linux-gnu-library/3.0/snpStats/libs')

metab=read.csv('/cluster/project8/vyp/eQTL_integration/data/WHII/rawData/metabolomics_phase5/dsg200k_9_metabolites.csv')
snpids=read.table('/cluster/project8/vyp/claudia/info_1000G_EUR.txt', header=T, nrows=17076866, colClasses=c('character','character','character','character','numeric'), sep='\t')

#remove duplicates (important for rownames)
dups=duplicated(metab$FID)
metab=metab[!dups,]

#adjust names for snpStats
rownames(metab)=metab$FID

#load genotype data
filename=paste0('/cluster/project8/vyp/eQTL_integration/data/WHII/rawData/genotypes/chr', chr, '_WHII_imputed.Rdata')
#load genotype
load(filename)
#########################################################

#fix genotype row names to match metab
rownames(genotypes)=strtrim(rownames(genotypes),10)

#prepare initial columns of data frame
snp.loc.temp=SNP.support$position.hg19
assign(paste('snp.locations','.',chr,sep=''),snp.loc.temp)
snp.chr.temp=rep(chr, length(snp.loc.temp))
assign(paste('snp.chromosome','.',chr,sep=''),snp.chr.temp)

snp.common=snpids$SNP%in%SNP.support$SNP
snp.id=snpids$rsid[snp.common]

biomarker=data.frame(snp.id,snp.chr.temp,snp.loc.temp)
names(biomarker)=c('SNPID', 'CHR', 'POS')

column=4

#########################################################

#metabolite loop to construct table

for(met in 4:236){
metabolite=colnames(metab)[met]

results=snp.rhs.tests(metab[,met]~1, family='gaussian', snp.data=genotypes,data=metab)
estimates=snp.rhs.estimates(metab[,met]~1, family='gaussian', snp.data=genotypes,data=metab)

BETA.biom=sapply(as(estimates, 'list'), FUN = function(estimates) {if (is.null(estimates)) {return (NA);} else {return((estimates)$beta)} } )
SE.biom=sapply(as(estimates, 'list'), FUN = function(estimates) { if (is.null(estimates)) {return (NA);} else {return(sqrt((estimates)$Var.beta))} } )
PVAL.biom=p.value(results)
N.biom=sample.size(results)

biomarker[,column]=BETA.biom
names(biomarker)[column]=paste('BETA.',metabolite, sep='')
column=column+1

biomarker[,column]=SE.biom
names(biomarker)[column]=paste('SE.',metabolite, sep='')
column=column+1

biomarker[,column]=PVAL.biom
names(biomarker)[column]=paste('PVAL.',metabolite, sep='')
column=column+1

biomarker[,column]=N.biom
names(biomarker)[column]=paste('N.',metabolite, sep='')
column=column+1
}

save(biomarker, file=paste0('/cluster/project8/vyp/eQTL_integration/data/WHII/summaryStats/', oFile))
