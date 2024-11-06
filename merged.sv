##imported vcf file into Rstudio, total observations 57062 and tried many things
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("browseVCF")
install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install(c("GenomicRanges", "Organism.dplyr"))
BiocManager::install()
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("browseVCF")
library(browseVCF)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("browseVCF")
library(browseVCF)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("VariantAnnotation")
Yes
library(VariantAnnotation)
dataset_path <- file.choose()
print(dataset_path)
vcf_file <- "/Users/carlaperscky/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SEAGRASS RSTUDIO/merged.sv copy.vcf"
vcf_data <- readVcf(vcf_file, genome = NA)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("VariantAnnotation")
library(VariantAnnotation)
vcf_file <- "/Users/carlaperscky/Library/Mobile Documents/com~apple~CloudDocs/Desktop/SEAGRASS RSTUDIO/merged.sv.copy.vcf"
vcf_data <- readVcf(vcf_file, genome = "")
dataset_path <- file.choose()
print(dataset_path)
vcf_file <- "/Users/carlaperscky/Desktop/merged.sv.vcf"
vcf_data <- readVcf(vcf_file, genome = "mCamDro1.fa")
summary(vcf_data)
high_quality_variants <- vcf_data[info(vcf_data)$X44 >= 27]
summary(high_quality_variants)
high_quality_variants <- vcf_data[which(info(vcf_data)$QUAL > 30)]
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("StructuralVariantAnnotation")
library(StructuralVariantAnnotation)
 ##then i found the write pdf document that explained vcfR
library(vcfR)
library(dplyr)
library(ggplot2)
install.packages("vcfR")
install.packages("devtools")
devtools::install_github("knausb/vcfR")
library(vcfR)
dataset_path <- file.choose()
print(dataset_path)
vcf <- read.vcfR("/Users/carlaperscky/Desktop/RStudio/merged.sv.vcf")
  ##i was having problems with the columns, so I tried colnames
colnames(merged.sv) <- c("chrom", "pos", "id", "ref", "alt", "qual", "filter", "info")
vcf_df <- as.data.frame(vcf@fix)
  ##once i turned the vcf file into a df it became easier to manipulate
precise_variants <- vcf_df %>%
  filter(grepl("PRECISE", INFO) & QUAL >= 27)
 ##after filtering theres a total of 30858 SNPs 
sv_types <- table(gsub(".*SVTYPE=([A-Z]+).*", "\\1", precise_variants$INFO))
print(sv_types)
  BND   DEL   DUP   INS   INV 
 4331 18223  2707  4272  1325 
