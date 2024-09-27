rm(list = ls())  
Sys.setenv(LANGUAGE = "en")

library(geneHapR)

setwd("/home/yuanj/work/Haplotype-analysis/Os-SNP")

gff <- import_bed("/home/yuanj/work/Haplotype-analysis/Os-SNP/OsVIT1/OsVIT1.BED6")

tbl <- read.csv("/home/yuanj/work/Haplotype-analysis/Os-SNP/OsVIT1/Os04t0463400-01.CDS-1.geno")
# AccINFO <- import_AccINFO("/home/yuanj/work/Haplotype-analysis/Os-SNP/RiceVarMap2.csv",na.strings = "NA")

geneID <- "OsVIT1"      # 基因ID
Chr <- "Chr4"           # 基因所处的染色体名称
start <- 23134256      # 基因的起始位置（染色体坐标）
end <- 23137105     # 基因的终止位置（染色体坐标）
hapPrefix <- "H"        # 单倍型名称的前缀

hapResult <- table2hap(tbl, hapPrefix = hapPrefix,
            hetero_remove = TRUE, # 移除包含杂合位点的样本
            na_drop = TRUE) # 移除包含基因型缺失的样本

hapSummary <- hap_summary(hapResult, hapPrefix = hapPrefix)

png(file = "OsVIT1-1-1.png", width = 800, height = 600, units = "px", pointsize = 12, bg = "white")
plotHapTable(hapSummary,
            hapPrefix = hapPrefix,
            angle = 45,
            displayIndelSize = 0,
            title = geneID)
dev.off()


png(file = "OsVIT1-1-2.png", width = 800, height = 600, units = "px", pointsize = 12, bg = "white")
displayVarOnGeneModel(gff = gff, hapSummary = hapSummary,
                    startPOS = start-10,
                    endPOS = end+10,
                    CDS_h = 0.05, fiveUTR_h = 0.25, threeUTR_h = 0.25,
                    cex = 0.8)
dev.off()

# hapSummary[hapSummary == "DEL"] = "N"
# hapnet <- get_hapNet(hapSummary,                  # 单倍型结果
#                      AccINFO = AccINFO,           # 包含样本分类信息的数据框(data.frame)
#                      groupName = "Subpopulation", # 含有样本分类信息的列名称
#                      na.label = "Unknown")        # 未知分类样本的类别

# plotHapNet(hapnet,                          # 单倍型网络
#            scale = "log2",                  # 标准化方法"log10"或"log2"或"none"
#            show.mutation = 2,               # 是否展示变异位点数量, 0,1,2,3
#            col.link = 2, link.width = 2,    # 单倍型之间连线的颜色和宽度
#            main = geneID,                   # 主标题
#            pie.lim = c(0.5, 2),               # 圆圈的大小
#            legend_version = 1,              # 图例形式（0或1）
#            labels = T,                      # 是否在单倍型上添加label
#            # legend = FALSE)                # 不添加图例
#            # legend = TRUE)                 # 添加图例,但需要单击添加的位置
#            legend = c(12,0),                # 图例的坐标
#            cex.legend = 0.6)                # 图例中文字的大小