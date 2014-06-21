#dataset1 <- list(DC_MTB_Barreiro = 'logFC')
dataset2 <- list(LCL_dexamethasone_DiRienzo = 'logFC')
dataset1 <- list(WB_dexamethasone_DiRienzo = 'logFC')

gene <- 'OAS1'
#gene <- 'SERPINB10'

folder1 <- paste('data/', names(dataset1)[[1]], '/eQTLs/fgwas/fgwas_individual_files/', sep = '')
folder2 <- paste('data/', names(dataset2)[[1]], '/eQTLs/fgwas/fgwas_individual_files/', sep = '')

all.genes1 <- gsub(pattern = '.*_|.tab', replacement = '', list.files (folder1))
all.genes2 <- gsub(pattern = '.*_|.tab', replacement = '', list.files (folder2))

inter.genes <- intersect(all.genes2, all.genes1)

file1 <- list.files(folder1, pattern = gene, full.name = TRUE)
file2 <- list.files(folder2, pattern = gene, full.name = TRUE)

data1 <- read.table(file1, header = TRUE, sep = '\t')
data2 <- read.table(file2, header = TRUE, sep = '\t')



#print(subset(data1, SNPID %in% c('rs10774671', 'rs2285934')))
#print(subset(data2, SNPID %in% c('rs10774671', 'rs2285934')))

combined <- data1[, c('SNPID') , drop = FALSE]
combined$PVAL1 <- data1$PVAL

combined$PVAL2 <- data2$PVAL[ match(data1$SNPID, table = data2$SNPID) ]

