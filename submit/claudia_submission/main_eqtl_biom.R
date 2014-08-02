###################################################
## This script is called from a submit bash script 
args <- commandArgs(trailingOnly=TRUE)
biom.dataset=args[1]
eqtl.dataset=args[2]
biom.names=args[3]
condition=args[4]
chr.name=args[5]
type=args[6]


# If UCLEB, leave biom.names as NA and load later
if (biom.names=="ucleb") {
  biom.names.file <- "/cluster/project8/vyp/eQTL_integration/scripts/coloc/biom.names.RData"
  load(biom.names.file)
}

##############################################
## Variables to set 
##############################################

# Directory where the results will be stored 
results.folder = "/cluster/project8/vyp/eQTL_integration/coloc/"
# Directory on the cluster where the data are stored  
data.folder = "/cluster/project8/vyp/eQTL_integration/data/"

p12 <- 1e-6

##############################################
# Submission of /cluster/project8/vyp/eQTL_integration/scripts/coloc/eQTL_biom_coloc.R
source("/cluster/project8/vyp/eQTL_integration/scripts/coloc/eQTL_biom_coloc.R")

biom.rdata.fname <- paste(data.folder, biom.dataset, '/summaryStats/biom_chr', chr.name, '.RData', sep='')
eqtl.fname <- paste(data.folder, eqtl.dataset, '/eQTLs/fgwas/fgwas_', condition, '_summary_eQTLs.csv', sep='')
# ********** # Talk about this with Vincent
if (eqtl.dataset=="brain_UKBEC") (eqtl.fname=gsub("fgwas_", "fgwas_core_", eqtl.fname))

prefix = paste(condition, '_chr', chr.name, sep='')
out_main <- paste(results.folder, biom.dataset, '_', eqtl.dataset, '/', sep='')
outfolder <- paste(out_main, condition, '/', sep='')
if ( !file.exists(outfolder) ) dir.create(file.path(outfolder), showWarnings = FALSE, recursive=TRUE)

################# SCZ
## biom.rdata.fname <- paste(data.folder, biom.dataset, '/summaryStats/biom_chr', chr.name, '.RData', sep='')
## eqtl.fname <- paste(data.folder, eqtl.dataset, '/eQTLs/fgwas/fgwas_core_', condition, '_summary_eQTLs.csv', sep='')
##outfname <- paste(condition, 'chr', chr.name, '_summary.csv', sep='')
##if (file.exists(outfname)) (file.remove(outfname))
##biom.names = "scz"
##out_fld <- paste(out_main, condition, '/', sep='')
##dir.create(file.path(out_fld), showWarnings = FALSE)
##setwd(out_fld)


########################## Why can't I put these inside the function?
# For plotting with locuszoom
refFlat_path = "/cluster/project8/vyp/vincent/toolsVarious/locuszoom/refFlat.RData"
source('/cluster/project8/vyp/vincent/toolsVarious/locuszoom/call_locuszoom3_temp.R')
load(refFlat_path)
refFlatRaw <- refFlatRaw.VP


coloc.eqtl.biom(eqtl.fname, biom.rdata.fname, biom.names, p12, chr.name, type, plot=TRUE, outfolder, prefix)
