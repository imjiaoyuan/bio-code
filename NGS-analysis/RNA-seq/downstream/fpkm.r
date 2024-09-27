# 数据实例：
# # Program:featureCounts v2.0.6; Command:"featureCounts" "-T" "5" "-t" "exon" "-g" "gene_id" "-a" "IRGSP-1.0_representative_transcript_exon_2024-01-11.gtf" "-o" "gene.counts" "-p" "IRGSP_genome_73" "IRGSP_genome_74" "IRGSP_genome_75" "IRGSP_genome_76" "IRGSP_genome_77" "IRGSP_genome_78" "IRGSP_genome_79" "IRGSP_genome_80" "IRGSP_genome_81" "IRGSP_genome_82" "IRGSP_genome_83" "IRGSP_genome_84" 																	
# Geneid	Chr	Start	End	Strand	Length	Normal-1-1	Normal-1-2	Normal-1-3	HS-1-1	HS-1-2	HS-1-3	Normal-2-1	Normal-2-2	Normal-2-3	HS-2-1	HS-2-2	HS-2-3
# Os01g0100100	chr01;chr01;chr01;chr01;chr01;chr01;chr01;chr01;chr01;chr01;chr01;chr01	2983;3354;4357;5457;7136;8028;8232;8408;9210;10102;10274;10504	3268;3616;4455;5560;7944;8150;8320;8608;9615;10187;10430;10815	+;+;+;+;+;+;+;+;+;+;+;+	2935	748	711	802	1723	1745	1582	546	428	368	305	324	296
# Os01g0100200	chr01;chr01	11218;12152	12060;12435	+;+	1127	0	1	1	2	4	6	0	1	0	0	0	0
# Os01g0100300	chr01;chr01	11372;12146	12042;12284	-;-	810	0	0	0	0	0	0	0	0	0	0	0	0
# Os01g0100400	chr01;chr01;chr01;chr01;chr01	12721;13906;14359;14969;15266	13813;14271;14437;15171;15685	+;+;+;+;+	2161	105	71	111	296	281	269	164	112	129	20	36	33
# Os01g0100466	chr01;chr01	12808;13880	13782;13978	-;-	1074	0	0	0	0	0	0	0	0	0	0	0	0
# Os01g0100500	chr01;chr01;chr01;chr01;chr01;chr01;chr01;chr01	16399;17383;17558;18501;18968;19142;19531;19734	16976;17474;18258;18571;19057;19321;19629;20144	+;+;+;+;+;+;+;+	2222	1356	1269	1487	2043	1969	1892	1551	1218	1193	1019	1172	1121
# Os01g0100600	chr01;chr01;chr01;chr01;chr01;chr01	22841;23572;23962;24492;25445;25883	23281;23847;24033;24577;25519;26892	+;+;+;+;+;+	1960	370	352	396	781	724	666	428	299	279	343	365	270
# Os01g0100650	chr01	25861	26424	-	564	0	0	0	0	0	0	0	0	0	0	1	0
# 即counts矩阵，将counts矩阵转换为fpkm标准化矩阵

rm(list = ls())
# setwd('C:/Develop/pig/mimic_vs_m_nc')
counts <- read.csv(
    'counts',
    header = TRUE,
    sep = '\t',
    comment.char = '#',
    check.names = FALSE
)

depths <- colSums(counts[, 6:ncol(counts)])
for (clm in colnames(counts)[6:ncol(counts)]) {
    col_fpkm <- paste0(clm, "_FPKM")
    counts[col_fpkm] <- (counts[, clm]) / (counts[, "Length"] * 10^6  / depths[clm])
}
counts = counts[,-c(2:13)]
numeric_mask <- sapply(counts, is.numeric)
# counts[numeric_mask] <- lapply(counts[numeric_mask], function(x) ifelse(x < 1, x + 1, x))
write.table(counts, file = 'fpkm', sep = '\t', row.names = FALSE)