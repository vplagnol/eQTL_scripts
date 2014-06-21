# Submit using script : 'Submit_Phenotype_Associations_UCLEB.R'
args <- commandArgs(trailingOnly=TRUE)
biom.dataset=args[1]    
eqtl.dataset=args[2]    
condition=args[3]    
chr.name=args[4]   
type=args[5] 


# Submission of /cluster/project8/vyp/eQTL_integration/scripts/coloc/eQTL_biom_coloc.R
source("/cluster/project8/vyp/eQTL_integration/scripts/coloc/eQTL_biom_coloc.R")

main="/cluster/project8/vyp/eQTL_integration/"
out_main <- paste(main, 'coloc/', biom.dataset, '_', eqtl.dataset, '/', sep='')

biom.rdata.fname <- paste(main, 'data/', biom.dataset, '/summaryStats/biom_chr', chr.name, '.RData', sep='')
eqtl.fname <- paste(main, 'data/', eqtl.dataset, '/eQTLs/fgwas/fgwas_', condition, '_summary_eQTLs.csv', sep='')
p12 <- 1e-6
prefix = paste(condition, '_chr', chr.name, sep='')
#outfname <- paste(condition, 'chr', chr.name, '_summary.csv', sep='')
#if (file.exists(outfname)) (file.remove(outfname))
biom.names.file <- "/cluster/project8/vyp/eQTL_integration/scripts/coloc/biom.names.RData"
load(biom.names.file)
out_fld <- paste(out_main, condition, '/', sep='')
dir.create(file.path(out_fld), showWarnings = FALSE)
outfolder=out_fld


################# SCZ
## biom.rdata.fname <- paste(main, 'data/', biom.dataset, '/summaryStats/scz.swe.ripke.', chr.name, '.RData', sep='')
## eqtl.fname <- paste(main, 'data/', eqtl.dataset, '/eQTLs/fgwas/fgwas_core_', condition, '_summary_eQTLs.csv', sep='')
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


