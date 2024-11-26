
#repeat for each species vcf file
vcftools --vcf /Users/carlaperscky/Desktop/wd/bactrianus.recode.vcf \
         --chr NC_087445.1 \
         --maf 0.05 \
         --window-pi 100000 \
         --window-pi-step 100000 \
         --out diversity_bactrianus

# Load data in RStudio
file.choose()
pi_data <- read.table("/Users/carlaperscky/Desktop/wd/diversity_bactrianus.windowed.pi", header = TRUE)
# Plot nucleotide diversity
library(ggplot2)
ggplot(pi_data, aes(x = BIN_START, y = PI)) +
  geom_line() +
  labs(title = "Nucleotide Diversity (π) Across Chromosome NC_087445.1",
       x = "Genomic Position (bp)",
       y = "Nucleotide Diversity (π)") +
  theme_minimal()
