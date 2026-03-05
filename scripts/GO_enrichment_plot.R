# In R:
install.packages("ggplot2")
install.packages("dplyr")
install.packages("readr")

# Load your GO file:
library(ggplot2)
library(dplyr)
library(readr)

go <- read.delim("Biologcal GO.txt", header=TRUE)

# Clean columns:
lines <- readLines("Biologcal GO.txt")

go <- read.delim("Biologcal GO.txt",
                 skip = 11,     # 12 - 1
                 header = TRUE,
                 sep = "\t",
                 quote = "",
                 comment.char = "")

# Then: 
colnames(go) <- c("Term",
                  "Ref",
                  "GeneList",
                  "Expected",
                  "Direction",
                  "FoldEnrichment",
                  "Pvalue",
                  "FDR")

# Prepare Significant GO Terms:
library(dplyr)
library(ggplot2)

go$FDR <- as.numeric(go$FDR)

go_sig <- go %>%
  filter(FDR < 0.05) %>%
  arrange(desc(FoldEnrichment)) %>%
  head(15)

go_sig$logFDR <- -log10(go_sig$FDR)

# Filter for stress keywords:
go_stress <- go %>%
  filter(FDR < 0.05) %>%
  filter(grepl("stress|defense|immune|salt|heat|cold|oxidative|detox|abiotic|biotic", 
                Term, ignore.case = TRUE)) %>%
  arrange(desc(FoldEnrichment)) %>%
  head(15)

go_stress$logFDR <- -log10(as.numeric(go_stress$FDR))

# Plot:
ggplot(go_stress, aes(x = reorder(Term, logFDR),
                      y = logFDR)) +
  geom_bar(stat = "identity", fill = "#8B0000") +
  coord_flip() +
  theme_classic(base_size = 13) +
  labs(title = "Stress-Related GO Enrichment in Seeraga Samba (HIGH + MODERATE Variants)",
       x = "",
       y = "-log10(FDR)") +
  theme(plot.title = element_text(face = "bold"))
