library(phyloseq)
library(ggjoy)
library(dplyr)
library(ggtree)
library(ggplot2)
data("GlobalPatterns")

GP <- GlobalPatterns
GP <- prune_taxa(taxa_sums(GP) > 600, GP)
sample_data(GP)$human <- get_variable(GP, "SampleType") %in%
c("Feces", "Skin")

mergedGP <- merge_samples(GP, "SampleType")
mergedGP <- rarefy_even_depth(mergedGP,rngseed=394582)
mergedGP <- tax_glom(mergedGP,"Order")

melt_simple <- psmelt(mergedGP) %>%
filter(Abundance < 120) %>%
select(OTU, val=Abundance)
p <- ggtree(mergedGP) +
geom_tippoint(aes(color=Phylum), size=1.5)
facet_plot(p, panel="Abundance", data=melt_simple,
geom = geom_joy, mapping = aes(x=val,group=label,
fill=Phylum),
color='grey80', lwd=.3,
stat = "binline")