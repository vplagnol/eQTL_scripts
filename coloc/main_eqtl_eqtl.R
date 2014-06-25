# Submit using script : 'Submit_Phenotype_Associations_UCLEB.R'
args <- commandArgs(trailingOnly=TRUE)
eqtl.1=args[1]    
eqtl.2=args[2]    
condition.1=args[3]    
condition.2=args[4]    
outfld=args[5] 
basedir=args[6] 

#eqtl.1="monocytes_Knight"    
#eqtl.2="Smith_macrophages"     
#condition.1="LPS24logFC"    
#condition.2="logFC"
 
# Submission of /cluster/project8/vyp/eQTL_integration/scripts/coloc/eQTL_biom_coloc.R
source(paste(basedir, "/coloc/eQTL_coloc.R", sep ='') )

main="/cluster/project8/vyp/eQTL_integration/"

eqtl.fname1 <- paste(main, 'data/', eqtl.1, '/eQTLs/fgwas/fgwas_', condition.1, '_summary_eQTLs.csv', sep='')
eqtl.fname2 <- paste(main, 'data/', eqtl.2, '/eQTLs/fgwas/fgwas_', condition.2, '_summary_eQTLs.csv', sep='')
p12 <- 1e-6
dir.create(file.path(outfld), showWarnings = FALSE)

########################## Why can't I put these inside the function?
# For plotting with locuszoom
refFlat_path = "/cluster/project8/vyp/vincent/toolsVarious/locuszoom/refFlat.RData"
source('/cluster/project8/vyp/vincent/toolsVarious/locuszoom/call_locuszoom3_temp.R')
load(refFlat_path)
refFlatRaw <- refFlatRaw.VP

message("eqtl 1:", eqtl.fname1) 
message("eqtl 2:", eqtl.fname1) 
coloc.eqtl.eqtl(eqtl.fname1, eqtl.fname2, p12, doPlot=TRUE, outfld) 


