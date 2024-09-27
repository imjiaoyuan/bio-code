# 数据示例：
# "Geneid"	"Normal-1-1_FPKM"	"Normal-1-2_FPKM"	"Normal-1-3_FPKM"	"HS-1-1_FPKM"	"HS-1-2_FPKM"	"HS-1-3_FPKM"	"Normal-2-1_FPKM"	"Normal-2-2_FPKM"	"Normal-2-3_FPKM"	"HS-2-1_FPKM"	"HS-2-2_FPKM"	"HS-2-3_FPKM"
# "Os01g0100100"	254855195.911414	242248722.316865	273253833.049404	587052810.902896	594548551.959114	539011925.042589	186030664.39523	145826235.093697	125383304.940375	103918228.279387	110391822.827939	100851788.756388
# "Os01g0100200"	1	887311.446317657	887311.446317657	1774622.89263531	3549245.78527063	5323868.67790595	1	887311.446317657	1	1	1	1
# "Os01g0100300"	1	1	1	1	1	1	1	1	1	1	1	1
# "Os01g0100400"	48588616.3813049	32855159.648311	51365108.7459509	136973623.322536	130032392.410921	124479407.681629	75890791.3003239	51827857.473392	59694585.8398889	9254974.54881999	16658954.187876	15270708.005553
# "Os01g0100466"	1	1	1	1	1	1	1	1	1	1	1	1
# "Os01g0100500"	610261026.10261	571107110.711071	669216921.692169	919441944.19442	886138613.861386	851485148.514851	698019801.980198	548154815.481548	536903690.369037	458595859.585959	527452745.274527	504500450.045004
# "Os01g0100600"	188775510.204082	179591836.734694	202040816.326531	398469387.755102	369387755.102041	339795918.367347	218367346.938776	152551020.408163	142346938.77551	1.75e+08	186224489.795918	137755102.040816
# 使用fpkm标准化之后的数据

# setwd('./pub_share/')
rm(list = ls())  
Sys.setenv(LANGUAGE = "en")

library(DESeq2)

fpkm = read.csv(
    'fpkm_output', 
    header = T,  
    sep = '\t', 
    row.names = "Geneid", 
    comment.char = '#', 
    check.name = F
)

# 删掉前五列
# fpkm = fpkm[,-c(1:4)]

# 将小数四舍五入为整数
fpkm <- round(fpkm)

# 将空值全部写为1
numeric_mask <- sapply(fpkm, is.numeric)
fpkm[numeric_mask] <- lapply(fpkm[numeric_mask], function(x) ifelse(is.numeric(x) & x < 1, x + 100, x))

# 保留行相加大于10的数据
# fpkm <- fpkm[rowSums(fpkm)>10, ]

# 保证所有小数都转化为整数
fpkm[-1, ] <- apply(fpkm[-1, ], 2, as.integer)

# 如果有缺失值，删除这些行/列
missing_values <- sum(is.na(fpkm))
if (missing_values > 0) {
  fpkm <- na.omit(fpkm)
}

samples = data.frame(
    sampleID = c("Normal-1-1_FPKM", "Normal-1-2_FPKM", "Normal-1-3_FPKM", "HS-1-1_FPKM", "HS-1-2_FPKM", "HS-1-3_FPKM", "Normal-2-1_FPKM", "Normal-2-2_FPKM", "Normal-2-3_FPKM", "HS-2-1_FPKM", "HS-2-2_FPKM", "HS-2-3_FPKM"), 
    sample = c("sample1", "sample1", "sample1", "sample2", "sample2", "sample2", "sample3", "sample3", "sample3", "sample4", "sample4", "sample4")
)

# 按照sampleID更改samples的行名
rownames(samples) = samples$sampleID

# 将因子型数据的默认排序设置为1=sample1, 2=sample2, 3=sample3
samples$sample = factor(samples$sample, levels = c('sample1', 'sample2', 'sample3', 'sample4'))

# 构建 DESeqDataSet 对象
dds = DESeqDataSetFromMatrix(countData = fpkm, colData = samples, design = ~sample)

# 计算差异倍数并获得 p 值
# parallel = TRUE 可以多线程运行，在数据量较大时建议开启
dds_count <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)

# 查看样本上调还是下调
sampl1_vs_sample2 <- results(dds_count, contrast = c('sample', 'sample1', 'sample2'))

sampl3_vs_sample4 <- results(dds_count, contrast = c('sample', 'sample3', 'sample4'))

result1 <- data.frame(sampl1_vs_sample2, stringsAsFactors = FALSE, check.names = FALSE)
result2 <- data.frame(sampl3_vs_sample4, stringsAsFactors = FALSE, check.names = FALSE)

write.table(result1, 'sampl1_vs_sample2.DESeq2.txt', col.names = NA, sep = '\t', quote = FALSE)
write.table(result2, 'sampl3_vs_sample3.DESeq2.txt', col.names = NA, sep = '\t', quote = FALSE)