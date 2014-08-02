################## run on UGI machine (simple code); run AFTER submit_hits.sh has ENTIRELY completed
project='MRC1946'
mnt_to_vyp='~/mnt/CS/'
##################
signif.snps=data.frame()
for(chr in 1:22){
source_file=paste0(mnt_to_vyp, 'eQTL_integration/data/', project, '/summaryStats/signif/hits', chr, '.RData')
load(source_file)
signif.snps=rbind(signif.snps, hit.frame)
file.remove(source_file)
}

save(signif.snps, file=paste0(mnt_to_vyp, 'eQTL_integration/data/', project, '/summaryStats/signif/signif_1e-7.RData'))
