
require(GenomicRanges)
library(plyr)


biomarker_file = "/cluster/project8/vyp/eQTL_integration/data/SCZ_Ripke/summaryStats/scz.swe.ripke.2013.results.txt"
biomarker.data = read.table(biomarker_file, header=T)

min_pval = 5*10^(-8)
# Crate a region around min p-value and then find SNPs overlap that region and mark them with a "1", then do the same excluding SNPs in region==1..
all.sig <- length(which(biomarker.data$Pval < min_pval))


temp=biomarker.data 
region.nb = 1
x.with.region = data.frame()

for (i in 1:all.sig) {

  region<- c(temp$bp[which.min(temp$Pval)] - 150000, temp$bp[which.min(temp$Pval)] + 150000)
  chr<- temp$hg19chr[which.min(temp$Pval)]
  message('Min p-value for region.nb ', region.nb, ' in region chr', chr, ' ', region[1], ' ', region[2], ' is ', min(temp$Pval))

  my.BM.GRanges <- GRanges(seqnames = temp$hg19chr,
                              IRanges(start = temp$bp, end= temp$bp + 1))  


  my.QTL.GRanges <- GRanges(seqnames = chr,
                              IRanges(start=region[1],end= region[2]))

  my.overlap <- findOverlaps(query = my.BM.GRanges, subject = my.QTL.GRanges)

  snps.in.region <- temp$snpid[my.overlap@queryHits]
  #temp$region <- ifelse(temp$full.SNP %in% snps.in.region, region.nb, 'NA')


  x1<- temp[temp$snpid %in% snps.in.region,]
  x1$region<- region.nb
  x.with.region = rbind(x.with.region, x1)

  temp<- temp[!temp$snpid %in% snps.in.region,]
  if (min(temp$Pval) > min_pval) break;
  region.nb = region.nb +1

}

## We don't really care about this dataframe?
# x.with.region

  summary = ddply(x.with.region, .(region), summarise, chr=hg19chr[1], region.start=min(bp), region.stop=max(bp), min.p = min(Pval))

# write.table(x =  summary, file = paste('/ugi/scratch/claudia/coloc/', trait, '_summary.regions.sig.meta_analysis.txt', sep=''), row.names = FALSE, quote = FALSE, sep = '\t')



