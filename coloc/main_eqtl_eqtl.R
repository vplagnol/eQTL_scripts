###################################################
## This script is called from a submit bash script 


getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}

dataset1 <- 'default'
dataset2 <- 'default'
cond1 <- 'default'
cond2 <- 'default'

myArgs <- getArgs()

if ('dataset1' %in% names(myArgs)) dataset1 <- as.character( myArgs[[ 'dataset1' ]])
if ('dataset2' %in% names(myArgs)) dataset2 <- as.character( myArgs[[ 'dataset2' ]])
if ('cond1' %in% names(myArgs)) cond1 <- as.character( myArgs[[ 'cond1' ])]
if ('cond2' %in% names(myArgs)) cond2 <- as.character( myArgs[[ 'cond2' ])]
if ('match.by' %in% names(myArgs)) match.by <- as.character( myArgs[[ 'match.by' ])]

if (sum(c(dataset1, dataset2, cond1, cond2) == 'default') > 0) stop('All parameters must be set')


##############################################
## Variables to set 
##############################################

# Base directory of the github repository 
#scripts.folder = "/cluster/project8/vyp/kitty/eQTL_scripts/"
scripts.folder = "/cluster/project8/vyp/eQTL_integration/scripts/"
# Directory where the results will be stored 
#results.folder = "/cluster/project8/vyp/kitty/eQTL_scripts/coloc/results/"
results.folder = "/cluster/project8/vyp/eQTL_integration/coloc/"
# Directory on the cluster where the data are stored  
data.folder    = "/cluster/project8/vyp/eQTL_integration/data/"

dataset1 <- list ( LCL_dexamethasone_DiRienzo = 'logFC')
dataset2 <- list ( WB_dexamethasone_DiRienzo = 'logFC')
#dataset1 <- list ( monocytes_Knight = 'LPS24logFC')
#dataset2 <- list ( Smith_macrophages = 'logFC')
p12 <- 1e-6
match.by <- "gene" # either probe of gene 
##############################################

eqtl.dataset1 <- names(dataset1)
eqtl.dataset2 <- names(dataset2)
cond.1 <- as.character(dataset1[[1]])
cond.2 <- as.character(dataset2[[1]])

source(paste(scripts.folder, "/coloc/eQTL_coloc.R", sep ='') )
 
eqtl.fname1 <- paste(data.folder, eqtl.dataset1, '/eQTLs/fgwas/fgwas_', cond.1, '_summary_eQTLs.csv', sep='')
eqtl.fname2 <- paste(data.folder, eqtl.dataset2, '/eQTLs/fgwas/fgwas_', cond.2, '_summary_eQTLs.csv', sep='')
if ( !file.exists(results.folder) ) dir.create(file.path(results.folder), showWarnings = FALSE, recursive=TRUE)

# One folder for each eqtl-eqtl pair 
output.folder      =  paste(results.folder, "/", eqtl.dataset1, "_", eqtl.dataset2, "/", cond.1, "_", cond.2, sep = '')
if (!file.exists(output.folder)) dir.create( output.folder, recursive=TRUE )

########################## Why can't I put these inside the function?
# For plotting with locuszoom
refFlat_path = "/cluster/project8/vyp/vincent/toolsVarious/locuszoom/refFlat.RData"
source('/cluster/project8/vyp/vincent/toolsVarious/locuszoom/call_locuszoom3_temp.R')
load(refFlat_path)
refFlatRaw <- refFlatRaw.VP

coloc.eqtl.eqtl(eqtl.fname1, eqtl.fname2, p12, doPlot=TRUE, output.folder, match.by = match.by) 


