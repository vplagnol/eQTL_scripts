


plot.eQTL <- function (chromosome, positions, pvalues, output.pdf, gene.list = NULL,
                       TF.list = c('/cluster/project8/vyp/FANTOM/data/fantom_enhancers.bed.gz', '/cluster/project8/vyp/FANTOM/data/fantom_promoters.bed.gz'),
                       gene.chromosome = NA, gene.position.start = NA, gene.position.end = NA, gene.name = NA, gene.context = TRUE) {

  library(Gviz)
  library(GenomicRanges)
  library(GenomicFeatures)


  
  ###### obtain the information about gene model
  geneModel.file <- paste('/cluster/project8/vyp/vincent/toolsVarious/ensemblAPI/data/by_chrom/canonical_all_exons_full_Human_hg19_chr', chromosome, '.bed', sep = '')
  message('Reading ', geneModel.file)
  geneModel <- read.table(file = geneModel.file, header = TRUE, sep = '\t')
  geneModel$strand <- ifelse (geneModel$strand == -1, '-', '+')
  geneModel$exon <- paste(geneModel$transcript, 1:nrow(geneModel), sep = '_')
  geneModel$chromosome <- paste('chr', geneModel$chromosome, sep = '')
  
  ######### Now we store the SNP data in a GRanges object
  data <- data.frame(positions = positions, pvalues = -log10(pvalues), col = 'black')
  
  my.range <- range(c(positions, gene.position.start, gene.position.end), na.rm = TRUE)
  SNP.pos <- GRanges(seqnames = paste('chr', chromosome, sep = ''),
                     IRanges (start = positions - 10, end = positions + 10 ))
  gen <- genome(SNP.pos)


  ########### The gene model track
  geneModel.loc <- subset(geneModel, symbol %in% gene.list)
  if (nrow(geneModel.loc) >= 1) {
    my.gene.names <- unique( as.character(geneModel.loc$symbol))
    grtrack <- GeneRegionTrack(geneModel.loc, genome = gen,  chromosome = chromosome, name = ifelse ( length(my.gene.names) == 1, my.gene.names, 'Genes'), showId = TRUE, transcriptAnnotation = "symbol")
    annot.list <- list( grtrack )
  } else {
    aTrack <- AnnotationTrack(start = gene.position.start,
                              width = gene.position.end - gene.position.start, chromosome = gene.chromosome, strand = c("*"), name = paste('Probe ', as.character(gene.name)), col = 'black', showId = FALSE)
 
    annot.list <- list( aTrack )
  }


###########  do we plot the surrounding genes?
  if (gene.context) {
    my.chrom <- paste('chr', chromosome, sep = '')
    my.local.genes <- subset( geneModel, chromosome == my.chrom & end > my.range[1] & start < my.range[2])
    my.gene.names <- unique( as.character(my.local.genes$symbol))
    geneModel.loc2 <- subset(geneModel, symbol %in% my.gene.names)
    if (nrow(geneModel.loc2) > 0) {
      grtrack2 <- GeneRegionTrack(geneModel.loc2, genome = gen,  chromosome = chromosome, symbol = geneModel.loc2$symbol, name = 'Genes', showId = TRUE, transcriptAnnotation = "symbol")
      annot.list <- append(annot.list, grtrack2)
    }
  }

  
############## Now the transcription factor tracks
  nTFs <- length(TF.list)

  if (nTFs >= 1) {
    for (i in 1:nTFs) {
      my.cols <- c('red', 'blue', 'green') [ 1: nTFs ]

      if (file.exists(TF.list[ i ])) {
        bed <- read.table( TF.list[ i ] )
        annot <- gsub(pattern = '.bed.gz', replacement = '', basename ( TF.list[ i ] ))
        message('Reading annotation file: ', TF.list[ i ], ' with ', nrow(bed), ' rows')
      } else {
        annot <- TF.list[ i ]
        input.file <- paste('/SAN/biomed/biomed14/vyp-scratch/ENCODE/TFpeaks/split/', annot, '.bed.gz', sep = '')
        print(input.file)
        bed <- read.table(input.file)
      }

      bed$V1 <- ifelse (grepl(pattern = '^chr', bed$V1), bed$V1, paste('chr', bed$V1, sep = ''))  ###add the "chr" if this tag is missing
      
      TF.bed <- subset(bed,  V1 %in% paste('chr', chromosome, sep = '') & V3 >= my.range[1] & V2 < my.range[2])

      message('Overlap size ', nrow(TF.bed))
      
      TF <- GRanges(seqnames = TF.bed[,1],
                    IRanges (start = TF.bed$V2, end = TF.bed$V3))
      
      atrack <- AnnotationTrack(reduce(TF), name = gsub(pattern = '_', replacement = ' ', annot), col = my.cols[ i ])
      annot.list <- append(annot.list, atrack)


########## Now we colour the SNPs in a smart way
      my.overlap <- findOverlaps(SNP.pos, TF)
      data$col <- ifelse ( 1:nrow(data) %in% my.overlap@queryHits, my.cols[i], as.character(data$col))
    }
  }


############## Now we plot the track with the P-values, mostly complicated by the choice of colours
  data.mat <- data.frame (black = -log10(pvalues))
  if (nTFs >= 1) {
    data.mat$black <- ifelse (data$col == 'black', data.mat$black, NA)
    my.lev <- c('black')
    for (i in 1:length(my.cols)) {
      loc.col <- my.cols[ i ]
      if (sum (data$col == loc.col) > 0) {
        my.lev <- c(my.lev, loc.col)
        data.mat[, loc.col ] <- ifelse (data$col == loc.col,  -log10(pvalues), NA)
      }
    }
  }

  dtrack <- DataTrack(data = data.mat, start = positions,
                      end = positions + 1, chromosome = chromosome, genome = gen,
                      col = names(data.mat),
                      name = "eQTL P-values",
                      groups = factor(my.lev, my.lev))
  annot.list <- append( annot.list, dtrack)

  ############ And now a basic axis track
  axisTrack <- GenomeAxisTrack()
  annot.list <- append( annot.list, axisTrack)  

  
  pdf(output.pdf)
  plotTracks(annot.list, from = my.range[1], to = my.range[2])
  dev.off()
  print(output.pdf)
  return(data)
}


#input.file <- 'data/WB_dexamethasone_DiRienzo/eQTLs/fgwas_individual_files/fgwas_logFC_rs274883_ENSG00000105499_PLA2G4C.tab'

#data <- read.table(file = input.file,
#                   sep = '\t', header = TRUE, stringsAsFactors = FALSE)

#plot.eQTL (chromosome = data$CHR[1], positions = data$POS, pvalues = data$PVAL, output.pdf = 'test.pdf', gene.list = 'PLA2G4C', TF.list = 'RELA')
