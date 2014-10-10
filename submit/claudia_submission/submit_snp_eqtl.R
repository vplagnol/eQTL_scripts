# Claudia Giambartolomei 08/09/2014
source("/cluster/project8/vyp/claudia/scripts/eQTL_scripts/submit/claudia_submission/snp_eqtl.R")
##############################################
## Variables to set 
##############################################
try=c("--dataset=liver_Schadt", "--rsid=rs9982601", "--script=/cluster/project8/vyp/claudia/scripts/eQTL_scripts/submit/claudia_submission/submit_snp_eqtl.R")
#try=c("--dataset=GTex", "--tissue=brain", "--SNP.chr=21", "--SNP.position=35599128", "--script=/cluster/project8/vyp/claudia/scripts/eQTL_scripts/submit/claudia_submission/submit_snp_eqtl.R")
#$Rbin CMD BATCH --no-save --no-restore --rsid=$rsid --dataset=$dataset --SNP.chr=$chr.name --SNP.position=$pos $script
#data = c("liver_Schadt", "brain_UKBEC", "GTex")

getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}

myArgs <- getArgs()

if ('rsid' %in% names(myArgs)) rsid <- as.character( myArgs[[ 'rsid' ]])
if ('dataset' %in% names(myArgs)) dataset <- as.character( myArgs[[ 'dataset' ]])
if ('SNP.chr' %in% names(myArgs)) chr.name <- as.character( myArgs[[ 'SNP.chr' ]])
if ('SNP.position' %in% names(myArgs)) chr.name <- as.character( myArgs[[ 'SNP.position' ]])
if ('tissue' %in% names(myArgs)) tissue <- as.character( myArgs[[ 'tissue' ]])

base.folder  = "/cluster/project8/vyp/eQTL_integration/data/"

# If a tissue is specified, use that one, otherwise use all 
if (!exists("tissue")) ( cond.all= list.files(paste(base.folder, dataset, "/expression_data/", sep=""), pattern=".RData") )
if (exists("tissue")) ( cond.all = list.files(paste(base.folder, dataset, "/expression_data/", sep=""), pattern=tissue, ignore.case=T) )

res.all = data.frame()

for (fn in cond.all) {

  expression.file <- paste(base.folder, dataset, "/expression_data/", fn, sep="")
  condition = gsub(".RData", "", gsub("expression_", "", fn))

  res = snp.assoc(dataset=dataset, condition=condition, expression.file, rsid=rsid, min.MAF = 0.03, imp.quality = 0.3, gap = 500000, base.folder = '/cluster/project8/vyp/eQTL_integration/data/')
  # res = snp.assoc(dataset=dataset, condition=condition, expression.file, SNP.chr=SNP.chr, SNP.position=SNP.position, min.MAF = 0.03, imp.quality = 0.3, gap = 500000, base.folder = '/cluster/project8/vyp/eQTL_integration/data/')

  res.all = rbind(res.all, res)

}

write.table(x = res.all, file = "test.snp.tab", row.names = FALSE, quote = FALSE, sep = '\t')


