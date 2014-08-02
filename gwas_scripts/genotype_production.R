#produces rsid tagged genotypes from chr:location format
#also, creates the formalised list output
#nb will need to change script for titles of the base data e.g. WHII is chrXX_WHII_imputed_data.Rdata

source("http://bioconductor.org/biocLite.R")
biocLite("snpStats", lib=paste0('/home/', user, '/R/x86_64-unknown-linux-gnu-library/3.0/snpStats/libs'))##### Is this necessary? Running off cluster, seemed to give problems without this step every time.
library(snpStats, lib.loc=paste0('/home/', user, '/R/x86_64-unknown-linux-gnu-library/3.0/snpStats/libs'))

directory_to_chr_number='/cluster/project8/vyp/eQTL_integration/data/WHII/rawData/genotypes/chr'
rest_of_file_name='_WHII_imputed.Rdata'

#Load in genotype
load(paste0('/cluster/project8/vyp/eQTL_integration/data/WHII/rawData/genotypes/chr',chr,'_WHII_imputed.Rdata'))

#Load in Claudia's file
snpids=read.table('/cluster/project8/vyp/claudia/info_1000G_EUR.txt', header=T, nrows=17076866, colClasses=c('character','character','character','character','numeric'), sep='\t')

#Do 'genotypes'
snp.common=snpids$SNP%in%SNP.support$SNP
snp.id=snpids$rsid[snp.common]
loc=snpids$SNP[snp.common]
remove=which(duplicated(loc)) # had some issues with a chromosome in WHII, this was a fix. other problems might need alternate method
snp.id=snp.id[-loc]

row.names(genotypes)=strtrim(row.names(genotypes),10)

colnames(genotypes)=snp.id

#do 'map'
snp.name=snp.id
chromosome=SNP.support$chromosome
position=SNP.support$position.hg19
allele.1=SNP.support$Al1
allele.2=SNP.support$Al2

map=data.frame(snp.name, allele.1, allele.2, chromosome, position)

#put together
genotypes=list(genotypes,map)

save(genotypes, file=paste0('/cluster/project8/vyp/eQTL_integration/data/WHII/genotypes/chr',chr))
