# 先安装一些依赖的BiocManager的包，否则后续安装geneHapR容易出错
install.packages("BiocManager")
BiocManager::install(c("Biostrings", "GenomicRanges", "muscle", "IRanges", "rtracklayer", "trackViewer"))
# 安装geneHapR
install.packages("geneHapR")

# 首先把软件加载进来
library(geneHapR)

# 设置工作目录（windows的同学注意"\"和"/"的问题）
setwd("/home/yuanj/work/Haplotype-analysis/")

# 导入各种数据
gff <- import_gff("OsGHD7.gff3")                      # 导入GFF格式的注释数据
gff <- import_bed("12859_2023_5318_MOESM3_ESM.bed6")      # 导入BED格式的注释数据
pheno <- import_AccINFO("12859_2023_5318_MOESM4_ESM.tsv") # 导入表型数据
AccINFO <- import_AccINFO("12859_2023_5318_MOESM5_ESM.csv", 
                          sep = ",",                      # 分隔符号，默认为制表符"\t"
                          na.strings = "NA")              # 导入其他样本信息


# 导入VCF格式的基因型数据及变异信息的筛选
vcf <- import_vcf("OsGHD7.vcf.gz")
vcf <- filter_vcf(vcf, gff,
                  mode = "both",                       # both, POS, Type的其中一个
                  start = start, end = end, Chr = Chr, # 保留哪些位置的变异信息
                  type = "CDS", cusTyp = c("CDS"))     # 保留哪些类型的变异信息

# 导入表格形式的基因型数据
tbl <- read.csv("12859_2023_5318_MOESM2_ESM.geno")

# 导入其他格式的基因型数据请查看帮助
# 这里就不再一一展示了
help("geneHapR")
library(help = 'geneHapR')
browseVignettes('geneHapR')

# 一些基本参数的设置
geneID <- "OsGHD7"      # 基因ID
Chr <- "Chr7"           # 基因所处的染色体名称
start <- 9152403        # 基因的起始位置（染色体坐标）
end <- 9155185          # 基因的终止位置（染色体坐标）
hapPrefix <- "H"        # 单倍型名称的前缀

# 常用格式的单倍型鉴定，具体参数请查看帮助文档
# VCF:             vcf2hap()
# P.link(ped&map): plink.pedmap2hap()
# Fasta:           seqs2hap()
# HapMap:          hmp2hap()
# Table:           table2hap()

# 从VCF开始单倍型鉴定
hapResult <- vcf2hap(vcf, hapPrefix = hapPrefix,
                     hetero_remove = TRUE, # 移除包含杂合位点的样本
                     na_drop = TRUE) # 移除包含基因型缺失的样本
# 从表格形式的基因型数据开始单倍型分析
hapResult <- table2hap(vcf, hapPrefix = hapPrefix,
                       hetero_remove = TRUE, # 移除包含杂合位点的样本
                       na_drop = TRUE) # 移除包含基因型缺失的样本

# 对单倍型结果进行汇总整理
hapSummary <- hap_summary(hapResult, hapPrefix = hapPrefix)


# 将单倍型鉴定结果保存到硬盘
write.hap(hapResult, file = "GeneID.hapResult")
write.hap(hapSummary, file = "GeneID.hapSummary")

# 导入之前的单倍型分析结果
hapResult <- import_hap(file = "GeneID.hapResult")
hapSummary <- import_hap(file = "GeneID.hapSummary")

# 单倍型结果可视化分析
# 以表格形式展示各单倍型的基因型
plotHapTable(hapSummary,             # 单倍型结果
             hapPrefix = hapPrefix,  # 单倍型名称前缀
             angle = 45,             # 物理位置的角度
             displayIndelSize = 0,   # 图中展示最大的Indel大小
             title = geneID)         # 图片标题

# 在表格中添加注释信息
plotHapTable(hapSummary,             # 单倍型结果
             hapPrefix = hapPrefix,  # 单倍型名称前缀
             angle = 45,             # 物理位置的角度
             displayIndelSize = 4,   # 图中展示最大的Indel大小
             title = geneID,         # 图片标题
             INFO_tag = "ANN", tag_field = 11, geneName = geneID) # 添加的注释信息

# 在基因模式图上展示变异位点的信息
displayVarOnGeneModel(gff = gff, hapSummary = hapSummary,
                      startPOS = start-10,
                      endPOS = end+10,
                      CDS_h = 0.05, fiveUTR_h = 0.25, threeUTR_h = 0.25, # gene model parameters
                      cex = 0.8) # size of variants

# 单倍型网络分析
hapSummary[hapSummary == "DEL"] = "N"
hapnet <- get_hapNet(hapSummary,                  # 单倍型结果
                     AccINFO = AccINFO,           # 包含样本分类信息的数据框(data.frame)
                     groupName = "Subpopulation", # 含有样本分类信息的列名称
                     na.label = "Unknown")        # 未知分类样本的类别

plotHapNet(hapnet,                          # 单倍型网络
           scale = "log2",                  # 标准化方法"log10"或"log2"或"none"
           show.mutation = 2,               # 是否展示变异位点数量, 0,1,2,3
           col.link = 2, link.width = 2,    # 单倍型之间连线的颜色和宽度
           main = geneID,                   # 主标题
           pie.lim = c(0.5, 2),               # 圆圈的大小
           legend_version = 1,              # 图例形式（0或1）
           labels = T,                      # 是否在单倍型上添加label
           # legend = FALSE)                # 不添加图例
           # legend = TRUE)                 # 添加图例,但需要单击添加的位置
           legend = c(12,0),                # 图例的坐标
           cex.legend = 0.6)                # 图例中文字的大小

# 地理分布
AccINFO$Longitude <- as.numeric(AccINFO$Longitude)
AccINFO$Latitude <- as.numeric(AccINFO$Latitude)
hapDistribution(hapResult,             # 单倍型结果
                AccINFO = AccINFO,     # 含有地理坐标的数据框（data.frame）
                hapNames = c("H001", 
                             "H002", 
                             "H003"),  # 展示的单倍型名称建议不超过3个
                symbol.lim = c(3, 6),  # 圆圈的大小
                LON.col = "Longitude", # 经纬度所处的列名称
                LAT.col = "Latitude",  # 经纬度所处的列名称
                legend = "bottomleft", # 图例所处的位置
                cex.legend = 1,        # 图例大小
                lwd.pie = 0.2,         # 圆圈线条的粗细
                lwd = 1.5,             # 地图线条的粗细
                main = geneID)         # 主标题

# 连锁不平衡分析
plot_LDheatmap(hap = hapResult, # 单倍型结果
               add.map = TRUE,  # 是否添加基因模式图
               gff = gff,       # 注释信息
               Chr = Chr,       # 染色体名称
               start = start,   # 基因的起始位置
               end = end)       # 基因的终止位置（更多参数参见帮助文档）

# 表型关联分析
# 单个表型的分析
hapVsPheno(hap = hapResult,       # 单倍型分析结果
           pheno = pheno,         # 表型
           hapPrefix = hapPrefix, # 单倍型名称的前缀
           title = geneID,        # 主标题
           minAcc = 4,            # 参与p值计算所需的最小样本数
           symnum.args = list(    # 定义显著性标注方式
               cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
               symbols = c("***", "**", "*", "ns")),
           mergeFigs = TRUE)     # 结果包括两个图，是否融合成一张图

# 多个表型的分析
hapVsPhenos(hap = hapResult,
            pheno = pheno,
            hapPrefix = hapPrefix,
            title = geneID,
            compression = "lzw",                 # tiff文件的压缩方式
            res = 300, width = 12, height = 12,  # 图片大小的单位"inch"
            outPutSingleFile = TRUE,             # 只有pdf格式支持输出单个文件
            filename.surfix = "pdf",             # 文件格式: pdf, png, jpg, tiff, bmp
            filename.prefix = geneID)            # 文件名称为: prefix + pheno_name + surfix

# 位点效应估算（结果仅供参考）
# EFF <- siteEFF(hapResult, pheno)
# plotEFF(EFF, gff = gff,
#         Chr = Chr, start = start, end = end,
#         showType = c("five_prime_UTR", "CDS", "three_prime_UTR"), # see help(plotEFF)
#         y = "effect",                      # the means of y axis, one of effect or pvalue
#         ylab = "effect",                  # label of y axis
#         cex = 0.5,                         # Cex
#         legend.cex = 0.8,                  # legend size
#         main = geneID,                     # main title
#         CDS.height = 1,                    # controls the height of CDS, heights of others will be half of that
#         markMutants = TRUE,                # mark mutants by short lines
#         mutants.col = 1, mutants.type = 1, # parameters for appearance of mutants
#         pch = 20)                          # points type

# 逐位点比较变异效应
hapVsPhenoPerSite(hap = hapResult,              # 单倍型分析结果
                  pheno = pheno,                # 表型文件
                  phenoName = names(pheno)[10], # 表型名称
                  freq.min = 5)                 # 参与显著性计算的最小样本数
# 回车继续下一位点
# ESC退出当前命令

# 导入较大的VCF文件
filterLargeVCF()
filterLargeP.link()

# GUI操作
# 常规操作，加载geneHapR
library(geneHapR)
# 打开新世界的大门
startGUI.geneHapR()
# 剩下的交给你啦