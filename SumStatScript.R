#Script for production of matabolomics summary stats
#Written by Jonathan Griffiths
#Written to function on UGI machine - modifications needed for finalisation

metab=read.csv('/ugi/home/shared/jonathanGriffiths/data/WHII_metabolomics/WHII_dsg200k_9_metabolites.csv')

snpids=read.table('~/mnt/CS/claudia/info_1000G_EUR.txt', header=T, nrows=17076866, colClasses=c('character','character','character','character','numeric'), sep='\t')
ptm=proc.time()
#remove duplicates (important for rownames)
dups=duplicated(metab$FID)
metab=metab[!dups,]

#adjust names for snpStats
rownames(metab)=metab$FID

#loop for chromosomes
for(chr in 22:22){
filename=paste0('/ugi/home/shared/jonathanGriffiths/data/WHII_genotypes_imputed/CHR', chr, '/SNPStats/merged/chr', chr, '_WHII_imputed.Rdata')
#load genotype
load(filename)


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

#metabolite loop
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
save(biomarker, file=paste('/ugi/home/zcqsjgr/biom_chr',chr,'.Rdata', sep=''))
#save(biomarker, file=paste('/ugi/home/zcqsjgr/mnt/CS/eQTL_integration/data/WHII/summaryStats/biom_chr',chr,'.Rdata', sep=''))

}
proc.time()-ptm
