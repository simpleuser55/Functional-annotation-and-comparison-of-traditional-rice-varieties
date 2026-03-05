# Load KEGG Enrichment File:
library(readr)
library(dplyr)
library(ggplot2)

kegg <- read.delim("TBtools.KEGG.enrichment.result.xls", 
                   header = TRUE,
                   sep = "\t",
                   quote = "",
                   comment.char = "")


# Map columns and convert numeric columns:
colnames(kegg) <- c("Pathway",
                    "MainClass",
                    "GeneCount",
                    "AllSelected",
                    "BgGeneHits",
                    "AllBackground",
                    "Pvalue",
                    "RichFactor",
                    "GeneList",
                    "Qvalue")

# Then:
kegg$Qvalue <- as.numeric(kegg$Qvalue)
kegg$RichFactor <- as.numeric(kegg$RichFactor)
kegg$GeneCount <- as.numeric(kegg$GeneCount)

# Then filter significant:
kegg_sig <- kegg %>%
  filter(Qvalue < 0.05)

# Focus on Stress + Metabolic Pathways:
kegg_focus <- kegg_sig %>%
  filter(grepl("stress|MAPK|hormone|zeatin|glutathione|lipid|fatty|terpenoid|secondary|signal|oxidative",
                Pathway,
                ignore.case = TRUE))

# Prepare Bubble Data:
kegg_focus$logFDR <- -log10(kegg_focus$Qvalue)

# Check Near-Significant KEGG Pathways:
kegg %>%
  arrange(Pvalue) %>%
  head(10)

# Then:
kegg_plot <- kegg %>%
  arrange(Pvalue) %>%
  head(6)

kegg_plot$logP <- -log10(kegg_plot$Pvalue)

# Plot with raw p-value:
ggplot(kegg_plot, aes(x = RichFactor,
                      y = reorder(Pathway, RichFactor),
                      size = GeneCount,
                      color = logP)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_classic(base_size = 13) +
  labs(title = "KEGG Pathway Enrichment (Top Pathways by Raw P-value)",
       x = "Rich Factor",
       y = "",
       color = "-log10(P-value)",
       size = "Gene Count")

