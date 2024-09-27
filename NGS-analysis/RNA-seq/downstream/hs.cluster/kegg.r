library(clusterProfiler)
library(org.Hs.eg.db)

genes <- read.delim(
    'human_genes.txt',  # 需要Entrez ID
    header = TRUE,
    stringsAsFactors = FALSE
)[[1]]

kegg <- enrichKEGG(
    gene = genes,
    organism = 'hsa',
    keyType = 'kegg',
    pAdjustMethod = 'fdr',
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
)

write.table(
    kegg,
    'entich.kegg.txt',
    sep = '\t',
    row.names = FALSE,
    quote = FALSE
)

barplot(kegg)  #富集柱形图
dotplot(kegg)  #富集气泡图