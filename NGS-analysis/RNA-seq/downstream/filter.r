setwd('C:/Develop/pig/mimic_vs_m_nc')
library(dplyr)
data <- read.table("fpkm", header = TRUE, sep = "\t")
columns_to_filter <- data[,2:ncol(data)]

filtered_data <- data[apply(columns_to_filter, 1, function(x) all(x > 1)), ]
write.table(filtered_data, "filtered", sep = "\t", row.names = FALSE, col.names = TRUE)