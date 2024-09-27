library("org.Hs.eg.db") 

g_symbol <- read.delim(
    'symbol.txt',
    header = TRUE,
    stringsAsFactors = FALSE
)[[1]]

entrez_id = mapIds(x = org.Hs.eg.db,
              keys = g_symbol,
              keytype = "SYMBOL",
              column = "ENTREZID")

write.table(
    entrez_id,
    'entrezid.txt',
    sep = '\t',
    row.names = FALSE,
    quote = FALSE
)