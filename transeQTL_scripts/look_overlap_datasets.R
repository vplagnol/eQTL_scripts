dataset1 <- 'WB_dexamethasone_DiRienzo'
dataset2 <- 'LCL_dexamethasone_DiRienzo'

condition <- 'logFC'
pval <- 5

for (chr in 22:1) {

  input.file1 <- paste('data/', dataset1, '/eQTLs/matrixEQTL/', condition, '_chr', chr, '_pval', pval, '.tab', sep = '')
  input.file2 <- paste('data/', dataset2, '/eQTLs/matrixEQTL/', condition, '_chr', chr, '_pval', pval, '.tab', sep = '')

  data1 <- read.table(input.file1, stringsAsFactors = FALSE, header = TRUE)
  data2 <- read.table(input.file2, stringsAsFactors = FALSE, header = TRUE)

  data1 <- subset(data1, ! is.na(ensemblID) & !cis.eQTL)
  data2 <- subset(data2, ! is.na(ensemblID) & !cis.eQTL)

  
  data1$sig <- paste(data1$SNP, data1$ensemblID)
  data2$sig <- paste(data2$SNP, data2$ensemblID)


  print(subset(data1, sig %in% intersect(data1$sig, data2$sig)))
  print(subset(data2, sig %in% intersect(data1$sig, data2$sig)))
}

