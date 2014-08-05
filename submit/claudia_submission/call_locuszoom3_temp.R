source('/cluster/project8/vyp/vincent/toolsVarious/locuszoom/tools_locuszoom_temp.R')
plink_path = "/cluster/project8/vyp/vincent/Software/plink-1.07-x86_64/plink"
#refFlat_path = "/ugi/home/shared/vincent/locusZoom_UGI/data/refFlat.RData"
#load(refFlat_path)



transformation <- function(x) {-log10(x)}

#locuszoom.ugi <- function(title,
#                          metal,
#                          pvalCol,
#                          markerCol,
#                          chrCol,
#                          posCol,
#                          chr,
#                          start,
#                          end,
#                          prefix = 'default',
#                          legend = 'auto',
#                          showRefsnpAnnot = TRUE,
#                          xlab = 'Position',
#                          unit = 1000,
#                          showGenes = FALSE,
#                          ylab = '-log10(P-value)',
#                          show_xlab = TRUE,
#                          theme = 'publication') {

locuszoom.ugi <- function(...) {

  argList <- list(...)
  x <- paste(names(argList), unlist(argList), sep = "=")


  snpset <- NULL
  showRecomb <- FALSE
  showAnnot <- FALSE

  
  
  

  flags <- list(flank=FALSE,reloaded=FALSE);
  createdFiles <- list();
  refSnpPos <- empty.data.frame();
  
#
                                        # set program defaults -- may be overridden with command line arguments
                                        #
  default.args <- list(
  temp.file.code = 'ldfile',
  theme = NULL,                         # select a theme (collection of settings) for plot
  experimental = FALSE,                 # try some experimental features?
  pquery = FALSE,                       # is pquery available?
  format = "pdf",                       # file format (pdf or png or both)
  clean=TRUE,                           # remove temp files?
  build = "hg18",                       # build to use for position information
  metal = "metal.tbl",                  # metal output file
  alreadyTransformed=FALSE,             # are metal p-values already -log10() -transformed?
  pvalCol="P.value",                    # name for p-value column in metal file
  posCol="pos",                         # name for positions column in metal file
  chrCol = 'chr',                   
  markerCol="MarkerName",               # name for MarkerName column in metal file
  weightCol="Weight",                   # name for weights column in metal file
  weightRange=NULL,                     # use this instead of actual range of weights
  ymin=0,                               # min for p-value range (expanded to fit all p-vals if needed)
  ymax=10,                              # max for p-value range (expanded to fit all p-vals if needed)
  yat=NULL,                             # values for y-axis ticks
  xat=NULL,                             # values for x-axis ticks
show_xlab = FALSE,
  xnsmall=NULL,                         # number of digits after decimal point on x-axis labels
  chr = NULL,                           # chromosome
  start = NULL,                         # start of region (string, may include Mb, kb, etc.)
  end = NULL,                           # end of region (string, may include Mb, kb, etc.)
  flank = "300kb",                      # surround refsnp by this much
  xlabPos = -3.0,                       # position of xaxis label (in lines relative to bottom panel)
  ylabPos = -3.0,                       # position of yaxis label (in lines relative to left edge of panel)
  ylab = NULL,                            # override default label for y-axis
  recombPos = 3.0,                      # position of recomb label (in lines relative to right edge of panel)
  axisSize = 1,                         # sclaing factor for axes
  axisTextSize = 1,                     # sclaing factor for axis labels
  axisTextColor = "gray30",             # color of axis labels
  requiredGene = NULL,                  # gene name (string)
  refsnp = NULL,                        # snp name (string)
  refsnpName = NULL,                    # name given to refsnp on plot (usually same as refsnp)
  drawMarkerNames = TRUE,               # draw the rs# SNP names above them?
  denoteMarkersFile = NULL,              # file specifying marker names to highlight (along with brief label and/or color)
  refsnpTextColor = "black",            # color for ref snp label
  refsnpTextSize = 1,                   # sclaing factor for text size
  refsnpTextAlpha = 1,                  # alpha for ref snp label
  refsnpLineColor = "transparent",      # color for ref snp line (invisible by default)
  refsnpLineAlpha = 1,                 # alpha for ref snp line
  refsnpLineWidth = 1,                  # width of ref snp line
  refsnpLineType = 2,                   # type of ref snp line (2 = dashed, 1 = solid, etc.. R defaults)
  signifLine = NULL,                    # draw a horizontal line at significance threshold (specify in p-value scale)
pp4_value = "",
  cond_snps = NULL,                     # list of SNPs that remain significant after conditional analysis
  cond_pos = NULL,                      # list of conditional SNPs in chr:pos form
  title = "",                           # title for plot
  titleColor = "black",                 # color for title 
  titleFontFace = "plain",              # font face for title, use "italic" for genes
  titleCex = 2,                         # size change for title
  thresh = 1,                           # only get pvalues <= thresh   # this is now ignored.
  width = 9,                           # width of pdf (inches)
  height = 9,                           # height of pdf (inches)
  leftMarginLines = 5,                  # margin (in lines) on left
  rightMarginLines = 5,                 # margin (in lines) on right
  unit=1000,                         # bp per unit displayed in plot # Default was 1000000, changed it
  ldTable = "results.ld_point6",        # LD Table (for SQL)
  annot=NULL,                           # file for annotation 
  showAnnot=TRUE,                       # show annotation for each snp?
  showGenes=TRUE,                       # show genes?
  annotCol='annotation',                # column to use for annotation, if it exists
  annotPch='24,24,25,22,22,8,7,21,1',        # plot symbols for annotation
  condPch='4,16,17,15,25,8,7,13,12,9,10',    # plot symbols for groups of LD blocks of SNPs
  condRefsnpPch=23,                     # symbol to use for refsnps in conditional plot, NULL means use same symbol as other SNPs in group
  annotOrder=NULL,                      # ordering of annotation classes
  showRefsnpAnnot=TRUE,                 # show annotation for reference snp too?
  bigDiamond=FALSE,                     # put big diamond around refsnp?
  ld=NULL,                              # file for LD information for reference SNP
  cond_ld=NULL,                         # files for LD information for conditional SNPs
  ldCuts = "0,.2,.4,.6,.8,1",           # cut points for LD coloring
  ldThresh=NULL,                        # LD values below this threshold will be colored identically
  ldColors = "gray60,navy,lightskyblue,green,orange,red,purple3",                   # colors for LD on original LZ plots
  condLdColors = "gray60,#E41A1C,#377EB8,#4DAF4A,#984EA3,#FF7F00,#A65628,#F781BF",  # colors for conditional LD plots (red, blue, green, purple, orange)
  condLdLow = NULL,                     # color all SNPs with conditional LD in the lowest bin the same color (specify that color here)
  ldCol='rsquare',                      # name for LD column
  LDTitle=NULL,                         # title for LD legend
  smallDot = .4,                        # smallest p-value cex 
  largeDot = .8,                        # largest p-value cex 
  refDot = NULL,                        # largest p-value cex 
  rfrows = '10',                        # max number of rows for reflat genes
  fmrows = 10,                           # max number of rows for fine mapping regions
  gwrows = 3,                           # max number of rows for gwas hits
  warnMissingGenes = FALSE,             # should we warn about missing genese on the plot?
  warnMissingFineMap = TRUE,            # should we warn about missing fine mapping regions on the plot? 
  warnMissingGWAS = TRUE,               # should we wearn about missing gwas hits on the plots? 
  showPartialGenes = TRUE,              # should genes that don't fit completely be displayed?
  shiftGeneNames = TRUE,                # should genes that don't fit completely be displayed?
  geneFontSize = .8,                    # size for gene names
  geneColor = "navy",                   # color for genes
  snpset = "Affy500,Illu318,HapMap",    # SNP sets to show
  snpsetFile = NULL,                    # use this file for SNPset data (instead of pquery)
  rugColor = "gray30",                  # color for snpset rugs
  rugAlpha = 1,                         # alpha for snpset rugs
  metalRug = NULL,                      # if not null, use as label for rug of metal positions
  refFlat = NULL,                       # use this file with refFlat info (instead of pquery)
  fineMap = NULL,                       # give a file with fine mapping posterior probabilities
  gwasHits = NULL,                      # give a file with GWAS catalog hits (chr, pos, trait)
  showIso=FALSE,                        # show each isoform of gene separately
  showRecomb = FALSE,                    # show recombination rate?
  recomb=NULL,                          # rcombination rate file
  recombAxisColor=NULL,                 # color for reccomb rate axis labeing
  recombAxisAlpha=NULL,                 # color for reccomb rate axis labeing
  recombColor='blue',                   # color for reccomb rate on plot
  recombOver = FALSE,                   # overlay recombination rate? (else underlay it)
  recombFill = FALSE,                   # fill recombination rate? (else line only)
  recombFillAlpha=0.2,                  # recomb fill alpha
  recombLineAlpha=0.8,                  # recomb line/text alpha
  frameColor='gray30',                  # frame color for plots
  frameAlpha=1,                         # frame alpha for plots
  legendSize=.8,                        # scaling factor of legend
  legendAlpha=1,                        # transparency of legend background
  legendMissing=TRUE,                   # show 'missing' as category in legend?
  legend='auto',                        # legend? (auto, left, right, or none)
  hiStart=0,                            # start of hilite region
  hiEnd=0,                              # end of hilite region
  hiColor="blue",                       # hilite color
  hiAlpha=0.1,                          # hilite alpha
  clobber=TRUE,                         # overwrite files?
  reload=NULL,                          # .Rdata file to reload data from
  prelude=NULL,                         # code to execute after data is read but before plot is made (allows data modification)
  postlude=NULL,                        # code to execute after plot is made (allows annotation)
  prefix=NULL,                          # prefix for output files
  dryRun=FALSE                          # show a list of the arguments and then halt
  )

### default data
  
  refSnpPos <- data.frame()
  recrate.default <- data.frame(chr=NA, pos=NA, recomb=NA, chr=NA, pos=NA)[c(),,drop=FALSE]
  rug.default <- data.frame(snp=NA, chr=NA, pos=NA, snp_set=NA)[c(),,drop=FALSE]
  annot.default <- data.frame(snp=NA,annot_rank=NA) # [c(),,drop=FALSE]
  ld.default <- data.frame(snp1='rs0000', snp2='rs0001', build=NA, 
                           chr=0, pos1=0, pos2=2, midpoint=1, distance=2, 
                           rsquare=0, dprime=0, r2dp=0) # [c(),,drop=FALSE]
  
  refFlatRaw.default <- data.frame(geneName=NA, name=NA, chrom=NA, strand=NA, txStart=NA, txEnd=NA, 
                                   cdsStart=NA, cdsEnd=NA, exonCount=NA, exonStarts=NA, exonEnds=NA, status=NA)[c(),,drop=FALSE]
  
#########
### Must add this otherwise get an error about rug$snp_set
  recrate <- recrate.default
  annot <- annot.default
  rug <- rug.default
  ld <- ld.default
#  refFlatRaw <- refFlatRaw.default
                                        #args[['showRug']] = FALSE
#########
  
                                        #
                                        # read and process command line arguments
                                        #
 



  user.args <- ConformList(argv(x),names(default.args),message=TRUE)

  default.args <- ProcessThemes(default.args,user.args[['theme']])
  
  args <- ModifyList(default.args,user.args);
 
  args <- AdjustModesOfArgs(args);

 
  # VP
  #args[[ 'showGenes' ]] <- showGenes
  #args[[ 'legend' ]] <- legend
  ##args[[ 'unit' ]] <- unit
  #args[[ 'show_xlab' ]] <- show_xlab
  #args[[ 'title' ]] <- title
  #args[[ 'ylab' ]] <- ylab
  #args[[ 'pvalCol']] <- pvalCol


  if (args[[ 'showGenes' ]]) {args[[ 'rfrows' ]] <- 10} else {args[[ 'rfrows' ]] <- 0}
  message('Show Genes ', args[[ 'showGenes' ]], ' ', args[['rfrows']])

  

#  refFlatRaw <- refFlatRaw.VP
                                        


userFile <- list(
      recomb = !is.null(args[['recomb']]),
  snpsetFile = !is.null(args[['snpsetFile']]),
     refFlat = !is.null(args[['refFlat']]),
          ld = !is.null(args[['ld']]),
       annot = !is.null(args[['annot']])
  );

args <- MatchIfNull(args,'recombAxisAlpha','recombLineAlpha')
args <- MatchIfNull(args,'recombAxisColor','recombColor')


if ( args[['pquery']] ){
  GetData <- GetDataFromFileOrCommand
} else {
  GetData <- GetDataFromFileIgnoreCommand
}

args[['showRefsnpAnnot']] <- args[['showAnnot']] & args[['showRefsnpAnnot']];

args[['refsnpColor']] <- args[['ldColors']][length(args[['ldColors']])];

if ( args[['dryRun']] )  {
  message("Argument list:");
  message(paste("\t",names(args),'=', args, "\n"));
  q();
}

#
# read metal data or reload all.
#

if ( is.null(args[['reload']]) ) {
    if ( file.exists( args[['metal']]) ) {
      metal <- read.file(args[['metal']],sep="\t");
    } else {
      stop(paste('No such file: ', args[['metal']]));
    }
} else {
  if ( file.exists(args[['reload']]) ) {
     load( args[['reload']] );
     flags[['reloaded']] <- TRUE;
  } else {
     stop(paste("Stopping: Can't reload from", args[['reload']]));
  }
}
#
# column renaming in metal data.frame
#
if ( char2Rname(args[['pvalCol']]) %in% names(metal) ) {
  metal$P.value <- metal[ ,char2Rname(args[['pvalCol']]) ];
} else {
  stop(paste('No column named',args[['pvalCol']]));
}

#transformation <- SetTransformation( min(metal$P.value,na.rm=TRUE), max(metal$P.value,na.rm=TRUE), 
#                                    args[['alreadyTransformed']] );


args[['LDTitle']] <- SetLDTitle( args[['ldCol']],args[['LDTitle']] )

if ( args[['posCol']] %in% names(metal) ) {
  metal$pos <- metal[ ,args[['posCol']] ];
} else {
  stop(paste('No column named',args[['posCol']]));
}

if ( args[['chrCol']] %in% names(metal) ) {
  metal$chr <- as.character(metal[ ,args[['chrCol']] ]);
} else {
  stop(paste('No column named',args[['chrCol']]));
}

if ( char2Rname(args[['markerCol']]) %in% names(metal) ) {
  metal$MarkerName <- metal[ ,char2Rname(args[['markerCol']]) ];
} else {
  stop(paste('No column named',args[['markerCol']]));
}

# Cannot plot if p-val is exactly zero: take out:
metal <- metal[metal$P.value!=0.000000e+00,]

#
# if no region and no refsnp specified, choose best snp and range of data set:
#
if ( (is.null(args[['start']]) || is.null(args[['end']]) || is.null(args[['chr']]) ) && ( is.null(args[['refsnp']]) ) ) 
{
  args[['start']] <- min(metal$pos);
  args[['end']] <- max(metal$pos);
  args[['chr']] <- min(metal$chr);
  args[['refsnp']] <- as.character( metal$MarkerName[ order(metal$P.value)[1] ] );

  args <- ModifyList(list(prefix=paste('chr',
      args[['chr']],"_",args[['start']],"-",args[['end']],sep='')),
      args);

  args <- ModifyList(list(prefix='foo'),args);
  flags[['flank']] <- FALSE;
  

# if region but not refsnp, choose best snp as refsnp
} else if ( !is.null(args[['start']]) && !is.null(args[['end']]) && !is.null(args[['chr']]) && is.null(args[['refsnp']] ) ) 
{
  args <- ModifyList(
    list( refsnp = as.character( metal$MarkerName[ order(metal$P.value)[1] ] ) ),
    args
    );
  flags[['flank']] <- FALSE;

# if refsnp specifed but no region, select region flanking refsnp
} else if ( ( is.null(args[['start']]) || is.null(args[['end']]) || is.null(args[['chr']]) ) && (!is.null(args[['refsnp']]) ) ) 
{
  args <- ModifyList( args, list( flankBP=pos2bp(args[['flank']]) ) );

  refSnpPosFile <- paste(args[['refsnp']],"_pos.tbl",sep="");

  command <- paste("pquery snp_pos",
            " -defaults",
            " -sql",
            " Snp=", args[["refsnp"]],
            " Build=",args[["build"]],
            sep="");
  if ( is.null(refSnpPos) ) { args[['showRug']] = FALSE }
  refSnpPos <- GetData( refSnpPosFile, default=refSnpPos.default, command=command, clobber=TRUE);

  args[['refSnpPos']] <- as.character(refSnpPos$chrpos[1]);
  args[['refSnpBP']] <- pos2bp(refSnpPos$chrpos[1]);

  args <- ModifyList( args, list( start=args[['refSnpBP']] - args[['flankBP']] ) ) ;
  args <- ModifyList( args, list( end=args[['refSnpBP']] + args[['flankBP']] ) );
  args <- ModifyList( args, list( chr=refSnpPos$chr[1] ) );

  flags[['flank']] <- TRUE;

# else refsnp and region specified
} else {  
  flags[['flank']] <- FALSE;
}

# change refsnp to "none" if it was null, else leave as is
args <- ModifyList( list( refsnp = "none"), args);

args <- ModifyList( args, list( start=as.character(args[['start']]) ) );
args <- ModifyList( args, list( end=as.character(args[['end']]) ) );

# prefix
if (flags[['flank']]) {
  args <- ModifyList(
    list( prefix = paste(                   # #1
        args[['refsnp']],
        "_",   args[['flank']],
        sep="")
        ),
    args
    );
} else {
  args <- ModifyList(
    list( prefix = paste(                   # #2
        "chr", args[['chr']],
        "_",   args[['start']],
        "-",   args[['end']],
        sep="")
        ),
    args
    );
}

#log
args <- ModifyList(
  list( log = paste(args[['prefix']], ".log", sep="") ),
  args 
    );

#recomb
args <- ModifyList(
  list( recomb = paste(args[['prefix']], "_recomb", ".tbl", sep="") ),
  args 
    );

# annot
args <- ModifyList(
  list( annot = paste(args[['prefix']], "_annot", ".tbl", sep="") ),
  args 
    );

# ld
args <- ModifyList(
  list( ld = paste(args[['prefix']], "_ld", ".tbl", sep="") ),
  args 
    );

# snpsets
args <- ModifyList(
  list( snpsetFile = paste(args[['prefix']], "_snpsets", ".tbl", sep="") ),
  args 
    );

# pdf
args <- ModifyList(
  list( pdf = paste(args[['prefix']], ".pdf", sep="") ),
  args
  );

args <- ModifyList(
  list( png = paste(args[['prefix']], ".png", sep="") ),
  args
  );

args <- ModifyList(
  list( tiff = paste(args[['prefix']], ".tiff", sep="") ),
  args
  );

# rdata
args <- ModifyList(
  list( rdata = paste(args[['prefix']], ".Rdata", sep="") ),
  args
  );

# refFlat
args <- ModifyList(
  list( refFlat = paste(args[['prefix']], "_refFlat.txt", sep="") ),
  args
  );

args <- ModifyList(args, list( startBP=pos2bp(args[['start']]), endBP=pos2bp(args[['end']]) ));
args <- ModifyList(args, list( hiStartBP=pos2bp(args[['hiStart']]), hiEndBP=pos2bp(args[['hiEnd']]) ));

#######################################################
#
# now read other (non-metal) data
# This makes all the output go into log file, so take out for now
#sink(args[['log']]);

##if ( is.null(args[['reload']]) ) {

# maybe translate this python command (from m2zfast.py) into R to load pquery?? from pquery import *
# Or try to source the dbmeister.py file from the beginning?? system("/Applications/locuszoom/src/m2zfast.py")

##  # recombination rate
##  command <- paste("pquery recomb_in_region",
##      " -defaults",
##      " -sql",
##      " RecombTable=", args[["recombTable"]],
##      " Chr=",args[["chr"]],
##      " Start=",args[["start"]],
##      " End=",args[["end"]],
##      sep="");
##  if ( is.null(args[['recomb']]) && ! args[['pquery']] ) { args[['showRecomb']] <- FALSE }
##  tryCatch(
##    recrate <- GetData( args[['recomb']], default=recrate.default, 
##      command=command, clobber=!userFile[['recomb']] || args[['clobber']] ),
##    error = function(e) { warning(e) }
##    )
##  if ( prod(dim(recrate)) == 0 ) { args[['showRecomb']] <- FALSE }
##  cat("\n\n");
 
## Vincent added this line: check if it was already null?
 recrate <- NULL



##  # snpset positions

##  command <- paste("pquery snpset_in_region",
##      " -defaults",
##      " -sql",
##      ' "SnpSet=',args[["snpset"]],'"',
##      " Chr=",args[["chr"]],
##      " ChrStart=",args[["start"]],
##      " ChrEnd=",args[["end"]],
##      sep="");
##  rug <- GetData( args[['snpsetFile']], default=rug.default, command=command, 
##    clobber=!userFile[['snpsetFile']] || args[['clobber']] );

##  cat("\n\nsnpset summary:\n");
##  print(summary(rug));
##  cat("\n\n");

##  # annotation
##  if ( char2Rname(args[['annotCol']]) %in% names(metal) ) {  
##    if (is.null(args[['annotOrder']])) {
##      args[['annotOrder']] <- 
##        sort( unique( metal[,char2Rname(args[['annotCol']])] ) );
##    }

##    metal$annot <- MakeFactor(metal[,char2Rname(args[['annotCol']]) ], levels=args[['annotOrder']],
##            na.level='none')
##    pchVals <- rep(args[['annotPch']], length=length(levels(metal$annot)));
##    metal$pch <- pchVals[ as.numeric(metal$annot) ]
##    annot <- metal$annot
##  } 

##  cat("\nR-DEBUG: Loading annotation data...\n");
##  if( args[['showAnnot']] && ! 'pch'  %in% names(metal) ) { 
##    command <- paste("pquery snp_annot_in_region",
##          " -defaults",
##          " -sql",
##          " Chr=",args[["chr"]],
##          " Start=",args[["startBP"]],
##          " End=",args[["endBP"]],
##          sep="");
##    if ( is.null(args[['annot']]) && !args[['pquery']] ) { args[['showAnnot']] <- FALSE }
##    annot <- GetData( args[['annot']], annot.default, command=command, 
##      clobber=!userFile[['annot']] || args[['clobber']] )
##    if (prod(dim(annot)) == 0) { args[['showAnnot']] <- FALSE }
##    cat("\nR-DEBUG: Merging in annotation data...");
##    metal <- merge(metal, annot,  
##      by.x='MarkerName', by.y="snp",
##      all.x=TRUE, all.y=FALSE);
##    cat(" Done.\n");
##    print(head(metal));

##    metal$annot <- c('no annotation','framestop','splice','nonsyn','coding','utr','tfbscons','mcs44placental')[1+metal$annot_rank];
##    if ( is.null(args[['annotOrder']]) ) {
##      args[['annotOrder']] <- c('framestop','splice','nonsyn','coding','utr','tfbscons','mcs44placental','no annotation')
##    } 
##    metal$annot <- MakeFactor(metal$annot, levels=args[['annotOrder']],na.level='none') 
##    pchVals <- rep(args[['annotPch']], length=length(levels(metal$annot)));
##    metal$pch <- pchVals[ as.numeric(metal$annot) ]

##  }  else {

    if (! 'pch' %in% names(metal)) {
      metal$pch <- 21;
    }

    if (! 'annot' %in% names(metal) ) {
        metal$annot <- "none"
        metal$annot <- factor(metal$annot)
    }
##    annot <- data.frame();
##  }

  if (FALSE) {  # scraps from above
    cat('else: ');
      pchVals <- rep(args[['annotPch']], length=length(levels(metal$annot)));
      metal$pch <- pchVals[ as.numeric(metal$annot) ]
      annot <- metal$annot
      print(xtabs(~annot+pch,metal));
      print(metal[1:4,])
  }

# This makes all the output go into file, so take out for now
#  sink('annotationTally.txt')
#  print( args[['annotOrder']] )
#  print(args[['annotPch']])
#  print(args[['annotOrder']])
#  print(table(metal$annot))
#  print(table(metal$pch))
#  print(xtabs(~annot+pch,metal))
#This cancels the sink:
#  sink()
  # ld

##  command <- paste("pquery ld_in_region",
##      " -defaults",
##      " -sql",
##      " LDTable=", args[["ldTable"]],
##      " Chr=",args[["chr"]],
##      " Start=",args[["startBP"]],
##      " End=",args[["endBP"]],
##      sep="");
      


  if ( is.null(args[['ld']]) && ! args[['pquery']] ) { 
    args[['legend']] = 'none' 
  }
  
  # Load LD for reference SNP. 
##  ld <- GetData( args[['ld']], ld.default, command=command, clobber=!userFile[['ld']] || args[['clobber']] )
  
  cond_ld = NULL;
  for (ld_file in args[['cond_ld']]) {
    cond_ld = rbind(cond_ld,read.table(ld_file,header=T,sep="",comment.char="",stringsAsFactors=F));
  }
  
  cat("\n\n");

  if (! is.null(args[['metalRug']]) ) {
    metalRug <- data.frame(pos=metal$pos, snp_set=args[['metalRug']]);
    origRug <- data.frame(pos=rug$pos,snp_set=rug$snp_set)
    rug <- rbind(origRug,metalRug)
#    print(levels(rug))
  }
    
##  save(metal,annot,recrate,ld,args,rug,file='loaded.Rdata');


  if ( prod(dim(metal) ) < 1) { stop("No data read.\n"); }


  # Subset the data to the plotting region. 
  s <- metal$pos >= args[['startBP']] & metal$pos <= args[['endBP']] & metal$chr == args[['chr']] ;
  metal <- subset(metal, s);

  # merge LD info into metal data frame
  refSnp <- as.character(args[['refsnp']]);

  metal$group <- 1;
#  metal$LD <- NA;
  metal$ldcut <- NA;
  metal$group[metal$MarkerName == refSnp] <- length(args[['ldColors']]);
  

  #ld <- ld(x = genotypes[,as.character(metal$MarkerName)], y = genotypes[, refSnp], stats = 'R.squared')
  #X <- as(genotypes[,as.character(metal$MarkerName)], 'numeric')
  #ld <- apply(MAR = 2, X, FUN = cor, y = X[, refSnp], use = 'pairwise.complete.obs')
  #ld <- ld(x = genotypes[,as.character(metal$MarkerName)], y = genotypes[, refSnp], stats = 'R.squared')

#with plink:
temp.file <- '/SAN/biomed/biomed14/vyp-scratch/vincent/eQTLs/region_ld'
if ('temp.file.code' %in% names(args)) {
  temp.file.code <- as.character(args[['temp.file.code']]) ##VVV
  temp.file <- paste('/SAN/biomed/biomed14/vyp-scratch/vincent/eQTLs/', temp.file.code, sep = '')
}

message('Ref SNP is ', refSnp)
my.cmd <- paste(plink_path, " --noweb --silent --bfile /cluster/project8/vyp/vincent/toolsVarious/locuszoom/EUR/chr", args[['chr']], " --ld-snp ", refSnp, " --chr ", args[['chr']], " --from-bp ",  args[['startBP']], " --to-bp ", args[['endBP']], " --r2 --ld-window-r2 0 --ld-window 999999 --ld-window-kb 99999 --out ", temp.file,  sep="")
print(my.cmd)
ld_works=system(my.cmd, ignore.stdout=T, ignore.stderr=T)

if(ld_works==1) message("Reference SNP is not contained in the locusZoom reference files")
# system(paste(plink_path, " --noweb --silent --bfile /cluster/project8/vyp/vincent/toolsVarious/locuszoom/EUR/chr", args[['chr']], " --ld-snp ", refSnp, " --chr ", args[['chr']], " --from-bp ",  args[['startBP']], " --to-bp ", args[['endBP']], " --r2 --ld-window-r2 0 --ld-window 999999 --ld-window-kb 99999 --out /SAN/biomed/biomed14/vyp-scratch/vincent/eQTLs/region_ld", sep=""))

if(ld_works==0) {
result.file <- paste(temp.file, '.ld', sep = '')
ld<- read.table(result.file, header=TRUE)
#This will fill "NA" is there is no LD (Current Plink genotypes do not include INDELS!!)
#metal$LD <- ld$R2[ ld$SNP_B  %in% metal$MarkerName]
metal <- cbind(metal, ld$R2[match(metal$MarkerName, ld$SNP_B)])
names(metal)[ncol(metal)] <- "LD"


##metal$LD <- ld[ names(ld) %in% metal$MarkerName]



#   if (!is.null(ld)) {
#     message('LD is NULL')
#     # subset ld for reference SNP
#     snpCols <- which(apply(ld,2,Sniff,type="snp"))
#     if (length(snpCols) != 2) {
#       warning(paste("LD file doesn't smell right. (",length(snpCols)," SNP cols)",sep=""))
#       assign("warningMessages",c(warningMessages,"LD file doesn't smell right."), globalenv());
#       break;
#     }
    
#     w1 <- which ( ld[,snpCols[1]] == refSnp );
#     w2 <- which ( ld[,snpCols[2]] == refSnp );
#     c1 <- c(names(ld)[snpCols[1]],names(ld)[snpCols[2]],args[['ldCol']]); # "rsquare","dprime");
#     c2 <- c(names(ld)[snpCols[2]],names(ld)[snpCols[1]],args[['ldCol']]); # "rsquare","dprime");
#     ld1 <- ld[ w1, c1, drop=FALSE ]
#     ld2 <- ld[ w2, c2, drop=FALSE ]
#     names(ld1)[1:2] <- c("refSNP","otherSNP")
#     names(ld2)[1:2] <- c("refSNP","otherSNP")
#     lld <- rbind( ld1, ld2);


    
#     if (prod(dim(lld)) > 0) { 
#       metal <- merge(metal, lld,  
#         by.x='MarkerName', by.y="otherSNP",
#         all.x=TRUE, all.y=FALSE
#       );
      
#       if ( args[['ldCol']] %in% names(metal) ) {
#         metal$LD <- metal[ ,args[['ldCol']] ];
#      } else {
#        stop(paste('No column named',args[['ldCol']]));
#      }
      
      metal$ldcut <- cut(metal$LD,breaks=args[['ldCuts']],include.lowest=TRUE);
      metal$group <- 1 + as.numeric(metal$ldcut);
      metal$group[is.na(metal$group)] <- 1;
      metal$group[metal$MarkerName == refSnp] <- length(args[['ldColors']]) 
#     } else {
#       print('No usable LD information')
#       assign("warningMessages",c(warningMessages,'No usable LD information for reference SNP.'), globalenv());
#       warning("No usable LD information.");
#       args[['legend']] <- 'none';
#     }
#   }
  

  
  if (!is.null(cond_ld)) {
    # Pool all LD information together. 
    all_ld = rbind(ld,cond_ld);
    
    # Get the reference + conditional SNPs from the LD files. 
    ref_cond_snps_names = c(
      args[['refsnpName']],
      unlist(strsplit(args[['cond_snps']],","))
    );
    
    #ref_cond_snps = as.character(unique(all_ld$snp2));
    ref_cond_snps = c(
      args[['refsnp']],
      unlist(strsplit(args[['cond_pos']],","))
    );
    
    # For each SNP, find the (reference SNP, conditional SNP) with the highest LD value. 
    by_best_ld = by(all_ld,all_ld$snp1,function(x) { ind = which(x$rsquare == max(x$rsquare)); x[ind,]; })
    best_ld = Reduce(rbind,by_best_ld);
    best_ld = best_ld[,c("snp1","snp2",args[['ldCol']])];
    names(best_ld) = c("snp","best_ld_snp","best_ld");
    
    # Merge into metal. 
    metal = merge(metal,best_ld,by.x="MarkerName",by.y="snp",all.x=TRUE,all.y=FALSE);
    
    # Ref/conditional SNPs are in "best LD" with themselves. 
    metal[match(ref_cond_snps,metal$MarkerName),]$best_ld_snp = ref_cond_snps;
    metal[match(ref_cond_snps,metal$MarkerName),]$best_ld = 1;
    
    # Did the user choose to threshold LD? 
    threshold_cuts = function(ld_cuts,thresh) { 
      ld_cuts = as.numeric(ld_cuts);
      ld_cuts = ld_cuts[!ld_cuts <= thresh];
      c(0,thresh,ld_cuts);
    }
    
    if (!is.null(args[['ldThresh']])) {
      args[['ldCuts']] = threshold_cuts(args[['ldCuts']],args[['ldThresh']]);
    }
    
    # Break up best LD into bins. 
    num_ld_bins = length(args[['ldCuts']]) - 1;
    metal$best_ld_cut = cut(metal$best_ld,breaks=args[['ldCuts']],include.lowest=TRUE);
    metal$best_ld_cut = 1 + as.numeric(metal$best_ld_cut);
    metal$best_ld_cut[is.na(metal$best_ld_cut)] = 1;
   
    # Each SNP belongs to a group, which is the SNP it has highest LD with. 
    metal$best_ld_group = as.numeric(factor(metal$best_ld_snp,levels=ref_cond_snps)) + 1;
    metal$best_ld_group[is.na(metal$best_ld_group)] = 1;
        
    # Compute each SNP's color. 
    metal$ld_color_base = args[['condLdColors']][metal$best_ld_group];
    
    col_pick = function(x,num_bins,low_color=NULL) { 
      cols = c(
        args[['condLdColors']][1],
        tail(colorRampPalette(c("white",x))(num_bins + 1),-1) 
      );

      if (!is.null(low_color)) {
        cols[2] = low_color; # cols[1] is missing LD, cols[2] is lowest LD bin
      }

      cols;
    };
    
    # Assign colors to SNPs based on their best LD.  
    metal$ld_color = apply(metal,1,function(x) { 
      this_cut = as.numeric(x['best_ld_cut']);
      col_pick(x['ld_color_base'],num_ld_bins,args[['condLdLow']])[this_cut];
    });
    
    # Collect colors together for LD ribbon legend. 
    base_colors = args[['condLdColors']][seq(2,length(ref_cond_snps)+1)]
    cond_ld_colors = sapply(base_colors,function(x) tail(col_pick(x,num_ld_bins,args[['condLdLow']]),-1),simplify=F);
    names(cond_ld_colors) = ref_cond_snps_names;

    # Color the ref/cond SNPs a different color so they stand out. 
    metal$ld_color[metal$MarkerName %in% ref_cond_snps] = tail(args[['ldColors']],1);

    # Plotting symbols are based on the group, not the annotation.
    metal$pch = args[['condPch']][metal$best_ld_group];
    
    # Ref/cond SNPs have their own symbol if set by user. 
    if (!is.null(args[['condRefsnpPch']])) {
      metal$pch[metal$MarkerName %in% ref_cond_snps] = args[['condRefsnpPch']];
    }
  }
  
##  save_objs = c('metal','refSnp','args');
##  if (!is.null(cond_ld)) {
##    save_objs = c(save_objs,'ref_cond_snps','ref_cond_snps_names','all_ld','num_ld_bins','col_pick','base_colors','cond_ld_colors');
##  }
##  save(list=save_objs,file='temp.Rdata');

##  command <- paste("pquery refFlat_in_region",
##      " -defaults",
##      " -sql",
##      " Chrom=", chr2chrom(args[["chr"]]),  
##      " Start=",args[["start"]],
##      " End=",args[["end"]],
##      " Build=",args[["build"]],
##      sep="");

                                        #if (is.null(args[['refFlat']]) && ! args[['pquery']]) { args[['showGenes']] <- FALSE }
##  refFlatRaw <- GetData( args[['refFlat']], refFlatRaw.default, command=command, 
##    clobber = !userFile[['refFlat']] || args[['clobber']] );

  summary(refFlatRaw);

  # subset the refFlatdata
  s <- refFlatRaw$txEnd >= args[['startBP']] & 
     refFlatRaw$txStart <= args[['endBP']] & 
     refFlatRaw$chrom == chr2chrom(args[['chr']]
  );
  refFlatRaw <- subset(refFlatRaw, s);
##  save(refFlatRaw,args,file="refFlatRaw.Rdata");

  refFlat <- flatten.bed(refFlatRaw,multiplier=1/args[['unit']], args=args);
  summary(refFlat);


                                        # load fine mapping data
  fmregions = NULL;
  if (!is.null(args[['fineMap']])) {
    fmregions = LoadFineMap(args[['fineMap']]);
  }

  # subset fine mapping regions to plotting region
  if (!is.null(fmregions)) {
    fmregions = subset(fmregions,chr == args[['chr']]);
  }
  summary(fmregions);

  # load gwas hits
  gwas_hits = NULL;
  if (!is.null(args[['gwasHits']])) {
    gwas_hits = LoadGWASHits(args[['gwasHits']]);
  }

  # subset gwas hits to plotting region
  if (!is.null(gwas_hits)) {
    b_gwas = (gwas_hits$chr == args[['chr']]) & (gwas_hits$pos <= args[['endBP']] / 1E6) & (gwas_hits$pos >= args[['startBP']] / 1E6);
    gwas_hits = subset(gwas_hits,b_gwas);
  }
  summary(gwas_hits);
  
  # load marker denote, if available
  denote_markers = NULL;
  if (!is.null(args[['denoteMarkersFile']])) {
    denote_markers = read.table(args[['denoteMarkersFile']],header=T,sep="\t",comment.char="",stringsAsFactors=F);
    
    b_denote = (denote_markers$chr == args[['chr']]) & (denote_markers$pos <= args[['endBP']]) & (denote_markers$pos >= args[['startBP']]);
    denote_markers = denote_markers[b_denote,];
    denote_markers$pos = denote_markers$pos / 1E6;
    
    if (!"color" %in% names(denote_markers)) {
      denote_markers$color = "black";
    }
  }

  # adjust for position units
  metal$pos <- metal$pos / args[['unit']];
  recrate$pos <- recrate$pos / args[['unit']];
## This doesn't work anymore???
##  rug$pos <- rug$pos / args[['unit']];

#  cat("recrate summary:\n");
#  print(summary(recrate));
#  cat("\n\n");
#  cat("LD summary:\n");
#  print(summary(ld));
#  cat("\n\n");
#  cat("metal summary:\n");
#  print(summary(metal));
#  cat("\n\n");
##  save(metal,annot,recrate,refFlatRaw,refFlat,rug,fmregions,gwas_hits,file=args[['rdata']]);
##} else {
##  load(args[['rdata']]);
##}


if (is.character(args[['prelude']]) && file.exists(args[['prelude']])) {
  source(args[['prelude']]);
}

if ( prod(dim(rug)) == 0 || !("snp_set" %in% names(rug)) ) {
  nrugs <- 0;
} else {
  nrugs <- length(levels(rug$snp_set));
}

xRange <- range(metal$pos,na.rm=T);
xRange <- as.numeric(c(args[['start']],args[['end']])) / args[['unit']];
refFlat <- refFlat[ which( (refFlat$start <= xRange[2]) & (refFlat$stop >= xRange[1]) ), ]
yRange <- c(min(c(args[['ymin']],transformation(metal$P.value),na.rm=T)),
            max(c(args[['ymax']],transformation(metal$P.value)*1.1),na.rm=T));

recrateRange <- c(0,max(c(100,recrate$recomb),na.rm=T));
if (args[['experimental']]) { 
  recrate$recomb <- max(c(100,recrate$recomb),na.rm=T) - recrate$recomb;
  recrateRange <- c(0,max(c(100,recrate$recomb),na.rm=T));
}
recrateRange <- rev(recrateRange);
#print("recrateRange: ");
#print(recrateRange);

refSnp <- as.character(args[['refsnp']]);
refidx <- match(refSnp, metal$MarkerName);
if (!args[['showRefsnpAnnot']]) {
  metal$pch[refidx] <- 23;  # use a diamond for ref snp
}

if ( prod(dim(metal)) == 0 ) { 
  message ('No data to plot.'); 
} else {
  zplot(metal,ld,recrate,refidx,nrugs=nrugs,args=args,postlude=args[['postlude']], xRange = xRange, yRange = yRange, refFlat = refFlat, recrateRange = recrateRange);
}


                                        # This makes all the output go into file, so take out for now
#sink(args[['log']], append=TRUE);
#  grid.log(args,metal,ascii=TRUE);
#  cat('\n\n\n');
#  cat("List of genes in region\n");
#    cat("#######################\n");
  geneList <- make.gene.list(refFlat,
                             showIso = args[[ 'showIso' ]],
                             unit=args[['unit']],
                             args = args);

#print(geneList)

  if (! is.null(geneList)) {
    digits <- 7 + ceiling(log10(max(geneList$stop)));
    print(geneList,digits=digits);
  }
  cat('\n\n\n');

#This cancels the sink:
# sink();

# save(metal,refFlat,ld,recrate,refSnpPos,args,file='end.Rdata')
# CleanUp(args,refSnpPos,recrate,rug,ld,refFlatRaw);

date();
}
}
