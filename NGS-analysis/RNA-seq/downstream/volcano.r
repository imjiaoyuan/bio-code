# 数据示例：
# 	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
# KYUSt_chr1.5-E1	213.082745455399	-0.159514113054599	0.32122397700057	-0.496582212025585	0.619483699498433	1
# KYUSt_chr1.6-E1	12.858979620238	-0.583646140965058	0.84821388423366	-0.688088407668978	0.491397110057002	1
# KYUSt_chr1.7-E1	1.93773722098783	0.307122539555924	2.44861358460283	0.125427115771613	0.900185422338187	NA
# KYUSt_chr1.8-E2	27.843256050414	0.92386601716311	1.12261950538863	0.822955607602139	0.41053323848971	0.993851867813233
# KYUSt_chr1.107-E3	9624.18154923748	1.2105461981764	0.753264345021032	1.60706690311036	0.108039692660408	0.826100628829101
# KYUSt_chr1.107-E2	5756.42276327944	0.282255247362519	0.52083207523245	0.541931384000394	0.587865775653251	1
# KYUSt_chr1.107-E1	6541.69156640286	0.396931212730972	0.672807064222869	0.589962908890457	0.555215516638021	1
# KYUSt_chr1.117-E1	87.1271890512162	0.356042837051628	0.492560961499353	0.722840145446841	0.469778099923088	1
# KYUSt_chr1.117-E2	569.916745538843	-0.32477765510939	0.363721126463563	-0.892930411458861	0.371894439924095	0.988626004998951
# KYUSt_chr1.117-E4	261.352157806559	-0.265783240477316	0.265421645074876	-1.00136234330979	0.316651662560004	0.975237273436443
# KYUSt_chr1.117-E6	58.3711351802128	-0.459967639584045	0.444678681825558	-1.03438203445175	0.300957599799505	0.968393795390404
# 对于新罗马字体的引用查看 https://yuanj.top/posts/m5q1u8y2/

setwd('./RNA-seq/downstream')
rm(list = ls())  
Sys.setenv(LANGUAGE = "en")

library(ggplot2)
library(showtext)
font_add("Times_New_Roman", "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf")
showtext_auto()

data <- read.table('./data.txt', header = TRUE)
data <- as.data.frame(data)

# 设置pvalue和logFC的阈值
cut_off_pvalue = 0.05
cut_off_logFC = 1

# 根据阈值分为上调基因down和下调基因down，无差异基因为No change，并保存到change列

data$change <- ifelse (
    data$pvalue < cut_off_pvalue & abs(data$log2FoldChange) >= cut_off_logFC, 
    ifelse(data$log2FoldChange > cut_off_logFC, 'Up', 'Down'),
    'No change'
)

p <- ggplot(
    data,
    aes(
        x = log2FoldChange,
        y = -log10(pvalue),
        colour = change,
    ))+
    geom_point(
        alpha = 0.4,
        size = 3.5,
    )+
    scale_color_manual(values = c("#546de5", "#d2dae2", "#ff4757"))+
    geom_vline(xintercept = c(-1, 1), lty = 4, col = "black", lwd = 0.8)+
    geom_hline(yintercept = -log10(cut_off_pvalue), lty = 4, col = "black", lwd = 0.8)+
    labs(
        x = "log2(fold change)",
        y = "-log10 (p-value)"
    )+
    theme_bw()+
    theme(
        axis.title.x = element_text(family = "Times_New_Roman", face = "bold", colour = "black", size = 30),  # 修改X轴标签字体大小
        axis.title.y = element_text(family = "Times_New_Roman", face = "bold", colour = "black", size = 30),  # 修改Y轴标签字体大小
        axis.text = element_text(family = "Times_New_Roman", face = "bold", colour = "black", size = 10),  # 修改轴刻度字体大小
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(family = "Times_New_Roman", face = "bold", colour = "black", size = 30)  # 修改图例文本字体大小
    )

plot <- p
plot

ggsave(plot, filename = "volcano.png", width = 10, height = 6, dpi = 300)