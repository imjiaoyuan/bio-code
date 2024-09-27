# 数据示例：
# 	C0_1	C0_2	C0_3	C50_1	C50_2	C50_3	C500_1	C500_2	C500_3
# KYUSt_chr1.6353	96.0	96.0	48.0	83.0	37.0	52.0	99.0	66.0	63.0
# KYUSt_chr1.21862	10.0	7.0	11.0	7.0	9.0	6.0	10.0	7.0	14.0
# KYUSt_chr1.23596	0.0	16.0	8.0	5.0	2.0	0.0	1.0	3.0	10.0
# KYUSt_chr1.30554	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
# KYUSt_chr1.34605	23.0	19.0	10.0	101.0	30.0	52.0	77.0	47.0	16.0
# KYUSt_chr1.38890	37.0	65.0	75.0	45.0	91.0	129.0	42.0	69.0	76.0
# KYUSt_chr2.2084	3290.0	1101.0	1441.0	1018.0	3123.0	812.0	1477.0	500.0	2176.0
# KYUSt_chr2.4392	0.0	0.0	0.0	0.0	0.0	0.0	0.0	40.0	0.0
# KYUSt_chr2.4554	1096.0	1042.0	610.0	1001.0	1485.0	1732.0	1800.0	801.0	1059.0
# KYUSt_chr2.4720	1144.0	1483.0	1002.0	1065.0	1927.0	1537.0	1077.0	1206.0	1274.0
# KYUSt_chr2.5003	1.0	2.0	2.0	0.0	0.0	5.0	2.0	0.0	1.0

setwd('./RNA-seq/downstream')
rm(list = ls())  
Sys.setenv(LANGUAGE = "en")

library(pheatmap)
library(showtext)
font_add("Times_New_Roman", "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf")  # 添加新罗马字体
showtext_auto()  # 始终启用字体

# data <- read.table('./LP_Cd.txt', header = TRUE,  sep = '\t',  row.names = 1 ) # 读取数据，第一列为行名
# data <- as.data.frame(data) # 将数据转换为数据框格式
# data = data[rowSums(data)>10, ]
# data = data[,-c(1:5)]
# pheatmap(data, scale = "row", clustering_method = "ward.D", cutree_rows = 3, cutree_cols = 2 , main = "title")

# 进行标准化：scale = "row"
# 进行k-means 聚类：kmeans_k = 3
# 取消行、列的距离：cluster_rows = FALSE 和 cluster_cols = TRUE
# 为行、列指定不同的度量：clustering_distance_rows = "correlation", clustering_distance_cols = "manhattan"
# 使用 clustering_method 参数来指定不同的聚类方法，支持的方法：'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid'

# 图例
# 通过 legend_breaks 参数设置断点，legend_labels 参数设置断点处的标签：legend_breaks = c(-1, 0, 1), legend_labels = c("low", "median", "high"))
# 不显示图例：legend = F

# 设置边框颜色：border_color = "white"
# 删除边框：border = F

# 单元格的尺寸：cellwidth = 15, cellheight = 15

# 显示数值：display_numbers = TRUE
# 以科学计数法显示数值：number_format = "%.1e"
# 设置字体的颜色：number_color = "#4daf4a"
# 设置字体大小：fontsize_number = 10
# 为display_numbers传递一个参数来标记显示+还是-：display_numbers = matrix(ifelse(data > 100, "+", "-"), nrow = nrow(data))

# 分块
# 先取消行和列的聚类：cluster_cols = FALSE, cluster_rows = FALSE
# 指定几个单元格为一个块：gaps_col = 4, gaps_row = c(8,9))
# 根据层次聚类的结果进行分块：  cutree_rows = 3, cutree_cols = 2

# 标签
# 标题：main = "title"
# 不显示行和列的标签：  show_colnames = FALSE, show_rownames = FALSE
# 设置行和列的标签字体：  fontsize_row = 8, fontsize_col = 12
# 设置行和列标签的旋转角度：angle_col = 45, angle_row = 45 ，可选角度有270、0、45、90、315
# 统一设置标签字体大小：fontsize = 15
# 自定义行和列的标签：  labels_row = paste0("Gene", LETTERS[1:20]), labels_col = rep(c("Cancer", "Noraml"), each = 3)

# 注释
# 先分别构建分组信息
# annotation_row <- data.frame(
#   GeneClass = c(rep("nano", 200), rep("nano2", 41))
# )
# rownames(annotation_row) <- rownames(data)

# annotation_col <- data.frame(
#   type = rep(c("sample1", "sample2", "sample3"), each = 3)
# )
# rownames(annotation_col) <- colnames(data)  # 设置为行的注释

data <- read.table('./expression.txt', header = TRUE,  sep = '\t',  row.names = 1 )
data <- as.data.frame(data)
data = data[rowSums(data)>10, ]

# annotation_row <- data.frame(
#   GeneClass = c(rep("MT", 17))
# )
# rownames(annotation_row) <- rownames(data)

annotation_col <- data.frame(
  type = rep(c("C0", "C50", "C500"), each = 3)
)
rownames(annotation_col) <- colnames(data)

p <- pheatmap(data,
  scale = "row", 
  # show_rownames = FALSE,
  cluster_cols = FALSE, 
  # cluster_rows = FALSE,
  clustering_method = 'ward.D',
  cellwidth = 30, 
  cellheight = 10, 
  # color = rainbow(6), 
  # cellheight = 20, 
  # gaps_col = c(3, 3),
  fontfamily = "Times_New_Roman",
  angle_col = 45,
  # annotation_row = annotation_row, 
  annotation_col = annotation_col, 
  # annotation_legend = FALSE # 取消注释的图例
  )

save_pheatmap_pdf <- function(x, filename, width=15, height=15) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

save_pheatmap_pdf(p, "result.pdf")