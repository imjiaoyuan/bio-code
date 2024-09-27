library(dplyr)
library(ggplot2)
# library(treedataverse)
library(ggtree)
library(ggtreeExtra)
library(showtext)

font_add("Times_New_Roman", "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf")  # 添加新罗马字体
showtext_auto()

library(treeio)
data <- read.tree("tree.nwk")
tree = fortify(data)

main_p = ggtree(tree,
    mapping = NULL,
    layout = "circular", # 常用circular，daylight
    open.angle = 0, # 部分支持，如fan
    mrsd = NULL, # 时间轴
    as.Date = FALSE,
    yscale = "none",
    yscale_mapping = NULL,
    ladderize = TRUE, # 阶梯状排列树
    right = FALSE,
    branch.length = "branch.length", # "none"就会让branch末端都对齐
    root.position = 0,
    xlim = NULL,
    #线段风格
    color="black", size=0.5, linetype=1)

main_p+
    geom_tiplab(align = T, linetype = 3, linesize = 0.5, family = "Times_New_Roman", face = "italic", size=4, color="black")  # 显示基因名称
    # geom_nodelab(mapping = aes(label=node))
    # geom_cladelabel(node=136, label = "Euarchotoglires", align = T, offset = 0.32, color = "red", barsize = 2) +
    # geom_cladelabel(node=120, label = "Afrotheria", align = T, offset = 0.32, color = "blue",fontsize = 6) +
    # geom_cladelabel(node=122, label = "Marsupialia", align = T, offset = 0.32, color = "green",angle = 45)