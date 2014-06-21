library(Gviz)
library(GenomicRanges)
library(GenomicFeatures)



prepare <- FALSE

if (prepare) {
  geneModels.VP <-read.table('/cluster/project8/vyp/vincent/toolsVarious/ensemblAPI/data/canonical_all_exons_full_Human_hg19.bed',header=T, sep = '\t', stringsAsFactors = F)

  geneModels.VP$strand <- ifelse (geneModels.VP$strand == -1, '-', '+')
  geneModels.VP$exon <- paste(geneModels.VP$transcript, 1:nrow(geneModels.VP), sep = '_')
  
  all.data <- read.table('data/eQTLs_no_annotations.tab', header = TRUE, stringsAsFactors = FALSE)
  gtrack <- GenomeAxisTrack()
    
  choice.TF <- c('RELA', 'SMARCA4')
  my.cols <- c('red', 'blue')
  bed.list <- list()
  nTFs <- length(choice.TF)

  
  for (i in 1:nTFs) {
    annot <- choice.TF[ i ]
    input.file <- paste('/SAN/biomed/biomed14/vyp-scratch/ENCODE/TFpeaks/split/', annot, '.bed.gz', sep = '')
    print(input.file)
    bed.list[[ i ]] <- read.table(input.file)
  }
}

region <- 'NSG1'
all.regions <- unique(all.data$REGION)
#all.regions <- 'PDGFRL'


for (region in all.regions) {
  
  data <- subset(all.data, REGION == region)
  data$col <- 'black'
  chr <- paste('chr', data$CHR[1], sep = '')
  my.range <- range(data$POS)

  SNP.pos <- GRanges(seqnames = chr,
                     IRanges (start = data$POS - 10, end = data$POS + 10 ))
  gen <- genome(SNP.pos)

  
  #biomartTrack <- BiomartGeneRegionTrack(genome="hg19", chromosome=chr, start=my.range[1], end=my.range[2], name="ENSEMBL")
  geneModels.loc <- subset(geneModels.VP, symbol == region)
  geneModels.loc$chromosome <- paste('chr', geneModels.loc$chromosome, sep = '')
  grtrack <- GeneRegionTrack(geneModels.loc, genome = gen,  chromosome = chr, name = paste("Gene Model", region))
  annot.list <- list( grtrack )
  
  
  for (i in 1:nTFs) {
    TF.bed <- subset(bed.list[[ i ]],  V1 == chr & V3 >= my.range[1] & V2 < my.range[2])
    TF <- GRanges(seqnames = TF.bed[,1],
                  IRanges (start = TF.bed$V2, end = TF.bed$V3))

    atrack <- AnnotationTrack(reduce(TF), name = choice.TF[ i ], col = my.cols[ i ])
    annot.list <- append(annot.list, atrack)

    my.overlap <- findOverlaps(SNP.pos, TF)
    data$col <- ifelse ( 1:nrow(data) %in% my.overlap@queryHits, my.cols[i], data$col)
  }
  gen <- genome(annot.list[[1]])
  
################

  data.mat <- data.frame (cblack = -log10(data$PVAL))
  data.mat$cblack <- ifelse (data$col == 'black', data.mat$cblack, NA)
  my.lev <- c('black')
  
  for (i in 1:length(my.cols)) {
    
    if (sum (data$col == my.cols[ i ]) > 0) {
      my.lev <- c(my.lev, my.cols[ i ])
      data.mat[, i+1] <- -log10(data$PVAL)
      data.mat[, i+1] <- ifelse (data$col == my.cols[ i ], data.mat[,i+1], NA)
    }
    
  }
  
  dtrack <- DataTrack(data = data.mat, start = data$POS,
                      end = data$POS + 1, chromosome = chr, genome = gen,
                      col = c('black', my.cols),
                      name = "eQTL P-values",
                      groups = factor(my.lev, my.lev))
  
############## Now we plot it all
  annot.list <- append( annot.list, dtrack)

  axisTrack <- GenomeAxisTrack()
  annot.list <- append( annot.list, axisTrack)
  
  output.pdf <- paste('fig/annotated/', region, '.pdf', sep = '')
  pdf(output.pdf)
  plotTracks(annot.list, from = my.range[1], to = my.range[2])
  dev.off()
  print(output.pdf)
}
