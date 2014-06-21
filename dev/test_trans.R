source('scripts/final_scripts/tools.R')

### define the dataset(s) of interest
dataset <- 'WB_dexamethasone_DiRienzo'
condition <- 'logFC'

### define the slice we are looking at
chromosome <- '22'
start <- 1
end <- .Machine$integer.max



##### basic option for matrixEQTL
covariates <- file <- name <- character()
errorCovariance = numeric();
useModel = modelLINEAR;

##### now define the output file that contains the eQTL data
region <- paste('chr', chromosome, '_', start, '_', end, sep = '')
if (start <= 1 && end >= 300*10^6) region <- paste('chr', chromosome, sep = '')
output <- file <- name <- paste('data/', dataset, '/eQTLs/matrixEQTL/', condition, '_', region, '_pval', -log10(pvOutputThreshold), '.tab', sep = '')

message('Output file in ', output <- file <- name)

output.file.GE <- paste('/scratch2/vyp-scratch2/vincent/eQTLs/', dataset, '_', condition, sep = '')
support.expression <- make.matEQTL.expression (dataset, condition, output.file.GE)

message('Output file in ', output <- file <- name)

output.file.GE <- paste('/scratch2/vyp-scratch2/vincent/eQTLs/', dataset, '_', condition, sep = '')
support.expression <- make.matEQTL.expression (dataset, condition, output.file.GE)

output.file <- paste('/scratch2/vyp-scratch2/vincent/eQTLs/', dataset, '_', chromosome, '_', start, '_', end, sep = '')
support.genotypes <- make.matEQTL.geno (dataset, chromosome, start, end, output.file)

