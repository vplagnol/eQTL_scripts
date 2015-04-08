find.trans.eQTLs <- function(base.folder = "/cluster/project8/vyp/eQTL_integration",
                             choice.sets,
                             chromosome.list = as.character(22:1),
                             min.MAF = 0.03,
                             pval.discovery = 10^(-12),
                             pval.validation = 10^(-5),
                             min.ngenes.per.module = 2,
                             pval.matrixeQTL = 5) {

  dataset <- names(choice.sets)[[ 1 ]]
  condition <- as.character(choice.sets[1])

  findings <- data.frame()
  
  message('Dataset: ', dataset, ' and condition: ', condition)


  for (chromosome in chromosome.list) {
    input.file <- paste(base.folder, '/data/', dataset, '/eQTLs/matrixEQTL/', condition, '/', condition, '_chr', chromosome, '_pval', pval.matrixeQTL, '.tab', sep = '')
    if (!file.exists(input.file)) {stop("File ", input.file, " does not exist")}
    message("Parsing ", input.file)
    data <- read.table(input.file, stringsAsFactors = FALSE, header = TRUE)

    data.discovery <- subset(data, p.value <  pval.discovery & MAF > min.MAF & !cis.eQTL & !gene.chromosome %in% c("X", "Y") & !is.na(gene.position.start))
    if (nrow(data.discovery) > 0) {
      message("Findings in step 1: ", nrow(data.discovery))
      validation.set <- subset(data, SNP %in% data.discovery$SNP & p.value < pval.validation & ! cis.eQTL )
      my.tab <- table(validation.set$SNP)
      good.SNPs <- names(subset(my.tab, my.tab >= min.ngenes.per.module))

      findings <- rbind.data.frame( findings, subset(data, SNP %in% good.SNPs))
      
    }
  }
  findings <- findings[ order(findings$SNP, findings$p.value), ]
  return(findings)
}


my.choice <- list("monocytes_TB_Nejentsev" = c('MTB'))
transeQTLs <- find.trans.eQTLs(choice.sets = my.choice)
transeQTLs <- subset(transeQTLs, !( chromosome == 17 & abs(44061023 - position) < 1000000)) ## pretty sure this is a silly artifact

write.csv(x = transeQTLs, file = "transeQTLs_Sergey_MTB.csv", row.names = FALSE)
