source('/cluster/project8/vyp/eQTL_integration/scripts/eQTLs_scripts/tools.R')

#force <- FALSE
#pvOutputThreshold = 1e-4

### define the dataset(s) of interest
#dataset <- 'WB_dexamethasone_DiRienzo'
#condition <- 'logFC'

### define the slice we are looking at
#chromosome <- '22'
#start <- 1
#end <- .Machine$integer.max




run.eQTL <- function( dataset,
                     condition,
                     chromosome, start = 1, end = 300*10^6,
                     pvOutputThreshold = 5, min.MAF = 0.03,
                     force = TRUE,
                     temp.folder = "/scratch2/vyp-scratch2/vincent/eQTLs",
                     base.folder = "/cluster/project8/vyp/eQTL_integration") {
  library(snpStats)
  library(MatrixEQTL)

  if (pvOutputThreshold < 1) {stop('Probably a misspecification of the pvalue output threshold, it should be given as -log10(p)')}
  
  ##### basic option for matrixEQTL
  covariates_file_name <- character()
  errorCovariance = numeric();
  useModel = modelLINEAR;

##### now define the output file that contains the eQTL data
  region <- paste('chr', chromosome, '_', start, '_', end, sep = '')
  if (start <= 1 && end >= 300*10^6) region <- paste('chr', chromosome, sep = '')

  oFolder <- paste(base.folder, '/data/', dataset, '/eQTLs/matrixEQTL/', condition, sep = '')
  if (!file.exists(oFolder)) dir.create(oFolder)
  message('Output folder is ', oFolder)
  output_file_name <- paste(oFolder, '/', condition, '_', region, '_pval', pvOutputThreshold, '.tab', sep = '')

######## Loading the genotype data
  message('Loading the genotype data')
  input.genotypes <- paste(base.folder, '/data/', dataset, '/genotypes/chr', chromosome, sep = '')
  load(input.genotypes)

######### Loading the expression data
  message('Loading the expression data')
  input.file <- paste(base.folder, '/data/', dataset, '/expression_data/expression_', condition, '.RData', sep = '')
  load(input.file)
  expression <- get(condition)
  support.expression <- get(paste('support', condition, sep = '.'))

#### now a check for row.names
  if ( sum(rownames(support.expression) != as.character(1:nrow(support.expression))) == 0) {
    message("There does not seem to be any rownames in the support file. Are you sure your data fit the guidelines?")
  }
  
  shared.samples <- intersect(dimnames(expression)[[2]], dimnames(genotypes$genotypes)[[1]])
  n.samples <- length(shared.samples)

  #### subset the samples so that we only look at shared expression/genotypes
  expression <- expression[, shared.samples]
  genotypes$genotypes <- genotypes$genotypes[ shared.samples, ]
  genotypes$fam <- genotypes$genotypes[ shared.samples, ]
  message('Nb of samples shared by expression and genotypes: ', n.samples)

####### apply the min.MAF threshold here
  if (min.MAF > 0) {
    message('Applying the MAF filter')
    my.sum <- col.summary( genotypes$genotypes )
    common.SNPs <- (my.sum$MAF >= min.MAF) & !is.na( my.sum$MAF )
    print(table(common.SNPs))
    genotypes$genotypes <- genotypes$genotypes[, common.SNPs ]
    genotypes$map <- genotypes$map[ common.SNPs, ]
  }
  
  
############
  message('Output file in ', output_file_name)
  
  output.file.GE <- paste(temp.folder, '/', dataset, '_', condition, '_chr', chromosome, sep = '')
  no.output <- make.matEQTL.expression (expression, output.file.GE)    

  output.file <- paste(temp.folder, '/', dataset, '_', condition, '_', chromosome, '_', start, '_', end, sep = '')
  support.genotypes <- make.matEQTL.geno (genotypes, chromosome, start, end, output.file)

  covar <- FALSE
  covar.file.name <- paste(base.folder, '/data/', dataset, '/covariates/covariates_', condition, '.tab', sep = '')
  if (!file.exists(covar.file.name)) {  covar.file.name <- paste(base.folder, '/data/', dataset, '/covariates/covariates.tab', sep = '')}

  
  if (file.exists(covar.file.name)) {
    message('Covariate file is ', covar.file.name)
    covariates <- read.table(covar.file.name, header = TRUE)
    covar <- TRUE
    row.names(covariates) <- covariates$id
    covariates <- covariates[ shared.samples, ]
    covariates <- covariates[, c('id', subset(names(covariates), names(covariates) != 'id')) ]  ##make sure that id is the first column
    matrixeQTL.covar.file <- paste(base.folder, '/data/', dataset, '/covariates/covariates_', condition, '_chr', chromosome, '.matrixeQTL', sep = '')
    write.table(x = t(as(covariates, 'matrix')), file = matrixeQTL.covar.file, row.names = TRUE, col.names = FALSE, quote = FALSE, sep = '\t')
  }
  
##### now some technical stuff for matrixEQTL
  if ( (!file.exists(output_file_name) || force)) {
    
    snps = SlicedData$new()
    snps$fileDelimiter = "\t"      # the TAB
    snps$fileOmitCharacters = "NA" # denote missing
    snps$fileSkipRows = 1          # one row of column
    snps$fileSkipColumns = 1       # one column of row
    snps$fileSliceSize = 2000      # read file in pieces of 2,000
    snps$LoadFile( output.file )
    
    gene = SlicedData$new()
    gene$fileDelimiter = "\t"  # the TAB character
    gene$fileOmitCharacters = "NA"  # denote missing values;
    gene$fileSkipRows = 1  # one row of column labels
    gene$fileSkipColumns = 1 # one column of row labels
    gene$fileSliceSize = 2000 # read file in slices of 2,000 rows
    gene$LoadFile(output.file.GE)
      
    cvrt = SlicedData$new();
    cvrt$fileDelimiter = "\t"; # the TAB character
    cvrt$fileOmitCharacters = "NA"; # denote missing values;
    cvrt$fileSkipRows = 1; # one row of column labels
    cvrt$fileSkipColumns = 1; # one column of row labels
    if( covar ) {
      message('There is a covariate file in ', matrixeQTL.covar.file)
      cvrt$LoadFile(matrixeQTL.covar.file);
    }
  
    me = Matrix_eQTL_engine(
      snps = snps,
      gene = gene,
      cvrt = cvrt,
      output_file_name = output_file_name,
      pvOutputThreshold = 10^-(pvOutputThreshold),
      useModel = useModel,
      errorCovariance = errorCovariance,
      verbose = TRUE,
      pvalue.hist = TRUE,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE);

    eQTL.data <- read.table(output_file_name, header = TRUE, stringsAsFactors = FALSE)
    eQTL.data$chromosome <- chromosome
    names(eQTL.data) <- replace(names(eQTL.data), names(eQTL.data) == 'gene', 'ProbeID') ###replace the poorly chosen column name gene with, instead, ProbeID
    
    ### add info from the SNP support file
    eQTL.data$position <- support.genotypes$position[ match(eQTL.data$SNP, row.names(support.genotypes)) ]
    eQTL.data$MAF <- support.genotypes$MAF[ match(eQTL.data$SNP, row.names(support.genotypes)) ]
    eQTL.data$Call.rate <- support.genotypes$Call.rate[ match(eQTL.data$SNP, row.names(support.genotypes)) ]
    
    eQTL.data$Gene.name <- support.expression$Gene.name[ match( eQTL.data$ProbeID, table = row.names(support.expression)) ]
    eQTL.data$ensemblID <- support.expression$ensemblID[ match( eQTL.data$ProbeID, table = row.names(support.expression)) ]

    if ( sum( c('gene.position.start', 'gene.position.end', 'gene.chromosome') %in% names(support.expression)) < 3 ) {
      message('Using standard annotations because the support file does not have all the key data')
      print(head(support.expression))
      annotations <- read.table('/cluster/project8/vyp/vincent/Software/pipeline/RNASeq/bundle/human/biomart/biomart_annotations_human.tab', header = TRUE)
      eQTL.data$gene.chromosome <- annotations$chromosome_name [ match( eQTL.data$ensemblID, table = annotations$EnsemblID) ]
      eQTL.data$gene.position.start <- annotations$start_position[ match( eQTL.data$ensemblID, table = annotations$EnsemblID) ]
      eQTL.data$gene.position.end <- annotations$end_position[ match( eQTL.data$ensemblID, table = annotations$EnsemblID) ]
    } else {  ### take the info from the support file if available
      eQTL.data$gene.chromosome <- support.expression$gene.chromosome[ match( eQTL.data$ProbeID, table = row.names(support.expression)) ]
      eQTL.data$gene.position.start <- support.expression$gene.position.start[ match( eQTL.data$ProbeID, table = row.names(support.expression)) ]
      eQTL.data$gene.position.end <- support.expression$gene.position.end[ match( eQTL.data$ProbeID, table = row.names(support.expression)) ]
    }

    
    eQTL.data$cis.eQTL <- (eQTL.data$chromosome == gsub(pattern = '^chr', replacement = '', eQTL.data$gene.chromosome)) & (abs(eQTL.data$position - eQTL.data$gene.position.start) < 0.5*10^6 | abs(eQTL.data$position - eQTL.data$gene.position.end) < 0.5*10^6)
    write.table(x = eQTL.data, file = output_file_name, sep = '\t', row.names = FALSE)
  } else {
    eQTL.data <- read.table(output_file_name, header = TRUE, stringsAsFactors = FALSE)
  }

  if (covar) clean <- file.remove(c(matrixeQTL.covar.file, output.file, output.file.GE))
  if (!covar) clean <- file.remove(c(output.file, output.file.GE))
  return(eQTL.data)
}






################## The function below essentially reproduces the run.eQTL one, but on a specified expression matrix
run.eQTL.selected <- function (genotypes, expression.matrix, expression.support = NULL, covariates = NULL, pvOutputThreshold = 2, min.MAF = 0.03) {
  library(fork)
  library(MatrixEQTL)
  errorCovariance = numeric();
  useModel = modelLINEAR;
  
  temp.folder <- '/SAN/biomed/biomed14/vyp-scratch/vincent/eQTLs'
  
  col.sum <- col.summary ( genotypes )
  genotypes <- genotypes[, col.sum$MAF > min.MAF ]

  jobid <- getpid()
  output.file <- paste(temp.folder, '/output_', jobid, '.tab', sep = '')
  output.file.geno <- paste(temp.folder, '/geno_', jobid, '.tab', sep = '')
  output.file.expression <- paste(temp.folder, '/expression_', jobid, '.tab', sep = '')

  
  shared.samples <- intersect(dimnames(expression.matrix)[[2]], dimnames(genotypes)[[1]])
  genotypes <- genotypes[ shared.samples, ]
  expression.matrix <- expression.matrix[, shared.samples ]
  
  no.output <- make.matEQTL.expression (expression.matrix, output.file.expression)
  no.output <- make.matEQTL.geno.basic (genotypes, output.file.geno)

  covar <- FALSE
  if (!is.null( covariates )) {
    covar <- TRUE
    covariates <- covariates[ shared.samples, ]
    output.file.covar <- paste(temp.folder, '/covar_', jobid, '.tab', sep = '')
    write.table(x = t(as(covariates, 'matrix')), file = matrixeQTL.covar.file, row.names = TRUE, col.names = FALSE, quote = FALSE, sep = '\t')
  }


  snps = SlicedData$new()
  snps$fileDelimiter = "\t"      # the TAB
  snps$fileOmitCharacters = "NA" # denote missing
  snps$fileSkipRows = 1          # one row of column
  snps$fileSkipColumns = 1       # one column of row
  snps$fileSliceSize = 2000      # read file in pieces of 2,000
  snps$LoadFile( output.file.geno )
  
  gene = SlicedData$new()
  gene$fileDelimiter = "\t"  # the TAB character
  gene$fileOmitCharacters = "NA"  # denote missing values;
  gene$fileSkipRows = 1  # one row of column labels
  gene$fileSkipColumns = 1 # one column of row labels
  gene$fileSliceSize = 2000 # read file in slices of 2,000 rows
  gene$LoadFile(output.file.expression)
  
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t"; # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1; # one row of column labels
  cvrt$fileSkipColumns = 1; # one column of row labels
  if( covar ) {
    message('There is a covariate file in ', output.file.covar)
    cvrt$LoadFile(output.file.covar)
  }
  
  me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output.file,
    pvOutputThreshold = 10^(-pvOutputThreshold),
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
  
  eQTL.data <- read.table(output.file, header = TRUE, stringsAsFactors = FALSE)
  names(eQTL.data) <- replace (names(eQTL.data), names(eQTL.data) == 'gene', 'ProbeID')
  
  if (!is.null(expression.support)) {
    loc.expression <- expression.support[ as.character(eQTL.data$ProbeID),]
    eQTL.data[, 'Gene.name'] <- loc.expression$Gene.name
  }
  
  rm <- file.remove(c(output.file.geno, output.file.expression))
  
  return(eQTL.data)
}
