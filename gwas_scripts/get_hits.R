###################
p.cutoff=1e-7
###################

load(paste0('/cluster/project8/vyp/eQTL_integration/data/', project, '/summaryStats/biom_chr', chr, '.RData'))

hit.frame=data.frame()

columns.p=seq(from=6, to=dim(biomarker)[2], by=4)
	for(met in columns.p){
	hits=which(biomarker[,met]<p.cutoff)
	snp.ids=biomarker$SNPID[hits]
	snp.chr=biomarker$CHR[hits]
	snp.pos=biomarker$POS[hits]
		if(length(snp.ids)==0){
		metabolite=character(0)}
		else{
		length=nchar(names(biomarker)[met])
		metab=substring(first=6, last=length, text=names(biomarker)[met])
		metabolite=rep(metab, length(snp.ids))}
	snp.pval=biomarker[hits, met]
	snp.effect=biomarker[hits, met-2]
	snp.number=biomarker[hits, met+1]
	snp.se=biomarker[hits, met-1]
	temp.hit.frame=data.frame(snp.ids, snp.chr, snp.pos, metabolite, snp.pval, snp.effect, snp.number, snp.se)
	hit.frame=rbind(hit.frame, temp.hit.frame)
	}

save(hit.frame, file=paste0('/cluster/project8/vyp/eQTL_integration/data/', project, '/summaryStats/signif/', oFile))
