
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
  labs(title = "C. bactrianus Nucleotide Diversity (π) Across Chr 10",
       x = "Genomic Position (bp)",
       y = "Nucleotide Diversity (π)") +
  theme_minimal()
file.choose()
pi_data2 <- read.table("/Users/carlaperscky/Desktop/wd/diversity_dromedarius.windowed.pi", header = TRUE)
ggplot(pi_data2, aes(x = BIN_START, y = PI)) +
  geom_line() +
  labs(title = "C. dromedarius Nucleotide Diversity (π) Across Chr 10",
       x = "Genomic Position (bp)",
       y = "Nucleotide Diversity (π)") +
  theme_minimal()
file.choose()
pi_data3 <- read.table("/Users/carlaperscky/Desktop/wd/diversity_dromedarius_chr1.windowed.pi", header = TRUE)
ggplot(pi_data3, aes(x = BIN_START, y = PI)) +
  geom_line() +
  labs(title = "C. dromedarius Nucleotide Diversity (π) Across Chr 1",
       x = "Genomic Position (bp)",
       y = "Nucleotide Diversity (π)") +
  theme_minimal()
file.choose()
pi_data4 <- read.table("/Users/carlaperscky/Desktop/wd/diversity_bactrianus_chr1.windowed.pi",  header = TRUE)
ggplot(pi_data4, aes(x = BIN_START, y = PI)) +
  geom_line() +
  labs(title = "C. bactrianus Nucleotide Diversity (π) Across Chr 1",
       x = "Genomic Position (bp)",
       y = "Nucleotide Diversity (π)") +
  theme_minimal()
file.choose()
pi_data5 <- read.table("/Users/carlaperscky/Desktop/wd/diversity_dromedarius_chr4.windowed.pi", header = TRUE)
ggplot(pi_data5, aes(x = BIN_START, y = PI)) +
  geom_line() +
  labs(title = "C. dromedarius Nucleotide Diversity (π) Across Chr 4",
       x = "Genomic Position (bp)",
       y = "Nucleotide Diversity (π)") +
  theme_minimal()
file.choose()
pi_data6 <- read.table("/Users/carlaperscky/Desktop/wd/diversity_bactrianus_chr4.windowed.pi", header = TRUE)
ggplot(pi_data6, aes(x = BIN_START, y = PI)) +
  geom_line() +
  labs(title = "C. bactrianus Nucleotide Diversity (π) Across Chr 4",
       x = "Genomic Position (bp)",
       y = "Nucleotide Diversity (π)") +
  theme_minimal()
