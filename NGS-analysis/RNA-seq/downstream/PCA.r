# 数据示例：
# Geneid	Chr	Start	End	Strand	Length	C0_1	C0_2	C0_3	C50_1	C50_2	C50_3	C500_1	C500_2	C500_3
# KYUSt_chr1.5-E1	1	53340	53891	+	552	233	146	195	201	189	264	177	326	243
# KYUSt_chr1.6-E1	1	54648	55664	-	1017	2	4	8	12	6	3	20	47	17
# KYUSt_chr1.7-E1	1	59936	60283	+	348	3	0	6	7	0	0	1	0	0
# KYUSt_chr1.8-E1	1	61782	61841	+	60	0	0	0	0	0	0	2	0	0
# KYUSt_chr1.8-E2	1	62373	62621	+	249	9	21	81	34	17	8	57	12	9
# KYUSt_chr1.9-E10	1	66745	66819	+	75	0	0	0	0	0	0	0	0	0
# KYUSt_chr1.9-E11	1	67360	67443	+	84	0	0	1	0	0	0	0	0	0
# KYUSt_chr1.9-E12	1	69776	69819	+	44	0	0	0	0	0	0	0	0	0
# KYUSt_chr1.9-E13	1	71271	71760	+	490	0	0	0	0	0	0	0	0	0
# KYUSt_chr1.107-E3	1	677255	677949	-	695	4115	26070	15855	3440	11431	5642	6410	5576	11130
# KYUSt_chr1.107-E2	1	678069	679025	-	957	3339	11707	5568	4545	6474	6084	6897	4574	4191
# KYUSt_chr1.107-E1	1	679121	679679	-	559	3645	16809	6030	3822	9043	7565	6074	5887	2603

setwd('./RNA-seq/downstream')
rm(list = ls())  
Sys.setenv(LANGUAGE = "en")

library(clusterProfiler)
library(DESeq2)
library(stringr)
library(ggplot2)
library(ggpubr)
library(showtext)
font_add("Times_New_Roman", "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf")  # 添加新罗马字体
showtext_auto()

counts = read.csv(
    'gene.counts', 
    header = T,  
    sep = '\t', 
    row.names = "Geneid", 
    comment.char = '#', 
    check.name = F
)

# 删掉前五列
counts = counts[,-c(1:5)]

# 保留行相加大于10的数据
counts = counts[rowSums(counts)>10, ]

samples = data.frame(
    sampleID = c("C0_1", "C0_2", "C0_3", "C50_1", "C50_2", "C50_3", "C500_1", "C500_2", "C500_3"), 
    sample = c("sample1", "sample2", "sample3", "sample1", "sample2", "sample3", "sample1", "sample2", "sample3")
)

# 按照sampleID更改samples的行名
rownames(samples) = samples$sampleID

# 将因子型数据的默认排序设置为1=sample1, 2=sample2, 3=sample3
samples$sample = factor(samples$sample, levels = c('sample1', 'sample2', 'sample3'))

# 将数据转为matrix格式
counts = as.matrix(counts[rownames(samples)])

# 从矩阵中抽取DESeq可以运行的数据形式
dds = DESeqDataSetFromMatrix(countData = counts, colData = samples, design = ~sample)

# 使用DESeq做差异分析
dds = DESeq(dds)

# vst方差稳定转换
vsd = vst(dds, blind = F)

# 做PCA
# 使用sample来分组，并返回坐标轴用于下一步作图
plotPCA(vsd, intgroup = c('sample'), returnData = TRUE)

PCAreturn <- "
              PC1        PC2   group  sample   name
C0_1   -16.586637   7.444676 sample1 sample1   C0_1
C0_2   -12.125149 -18.965823 sample2 sample2   C0_2
C0_3   -26.368675  11.700192 sample3 sample3   C0_3
C50_1   26.720069   1.819072 sample1 sample1  C50_1
C50_2  -16.740632  -1.364539 sample2 sample2  C50_2
C50_3    8.106638 -13.586068 sample3 sample3  C50_3
C500_1  16.105152  28.532094 sample1 sample1 C500_1
C500_2   3.091405  -3.977518 sample2 sample2 C500_2
C500_3  17.797829 -11.602086 sample3 sample3 C500_3
"

# 将PCAreturn转换为表格形式
PCAreturn <- read.table(header = TRUE, text = PCAreturn)

# 将sample列改为因子型并且排序
PCAreturn$sample = factor(PCAreturn$sample, levels = c('sample1', 'sample2', 'sample3'))

# ggplot绘图
ggplot(PCAreturn) +
    geom_point(aes(PC1, PC2, color = sample), size = 2) +
    labs(x = "PC1", y = "PC2", face = "bold") +
    theme_light(base_size = 16) +
    theme(
        text = element_text(family = "Times_New_Roman", face = "bold", colour = "black", size = 30),  # 修改整体文本的字体大小
        axis.title = element_text(family = "Times_New_Roman", face = "bold", colour = "black", size = 30),  # 修改坐标轴标题的字体大小
        axis.text = element_text(family = "Times_New_Roman", face = "bold", colour = "black", size = 30),  # 修改坐标轴刻度标签的字体大小
        legend.text = element_text(family = "Times_New_Roman", face = "bold", colour = "black", size = 15),  # 修改图例文本的字体大小
        plot.title = element_text(family = "Times_New_Roman", face = "bold", colour = "black", size = 30)  # 修改图形标题的字体大小
    )
