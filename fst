# Load data
file.choose()
pi_data <- read.table("/Users/carlaperscky/Desktop/wd/diversity_100kb.windowed.pi", header = TRUE)
# Plot nucleotide diversity
library(ggplot2)
ggplot(pi_data, aes(x = BIN_START, y = PI)) +
  geom_line() +
  labs(title = "Nucleotide Diversity (π) Across Chromosome NC_087445.1",
       x = "Genomic Position (bp)",
       y = "Nucleotide Diversity (π)") +
  theme_minimal()
