my.files <- list.files('coloc/UCLEB_GTex', recursive = TRUE, pattern = '.tab', full.name = TRUE)

final.frame <- data.frame()
for (file in my.files) {
  message(file)
  
  data <- read.table(file, sep = '\t', header = TRUE)
  data <- subset(data, pp4 > 0.5)

  if (nrow(data) > 0) {
    data$tissue <- basename(file)
    final.frame <- rbind.data.frame( final.frame, data )
  }
  
}


final.frame <- final.frame[ order(final.frame$pp4, decreasing = TRUE), ]
