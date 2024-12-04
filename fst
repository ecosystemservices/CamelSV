
vcftools --vcf /Users/carlaperscky/Desktop/wd/clean.vcf.recode.vcf
--remove-filtered-all --weir-fst-pop
/Users/carlaperscky/Desktop/wd/dromedarius.txt --weir-fst-pop
/Users/carlaperscky/Desktop/wd/bactrianus.txt --maf 0.05 --out
fst_output

install.packages("vcfR")
library(vcfR)
library(dplyr)
library(ggplot2)
fst_data <- read.table("fst_output.weir.fst", header = TRUE)
head(fst_data)
fst_data <- fst_data %>%
  mutate(POS = POS / 1e6) %>%
  group_by(CHROM) %>%
  mutate(cum_POS = POS + lag(cumsum(ifelse(CHROM != lag(CHROM, default = CHROM[1]), max(POS), 0)), default = 0)) %>%
  ungroup()
fst_data$chr_color <- as.numeric(factor(fst_data$CHROM)) %% 2
#establish chr labels
fst_data$CHR <- as.factor(fst_data$CHROM)
chr_labels <- fst_data %>%
  group_by(POS) %>%
  summarize(midpoint = mean(cum_POS))
ggplot(fst_data, aes(x = cum_POS / 1e6, y = WEIR_AND_COCKERHAM_FST, color = as.factor(CHROM))) +
  geom_point(size = 0.5) +
  scale_color_manual(values = rep(c("blue", "red"), length(unique(fst_data$CHROM)))) +
  scale_x_continuous(breaks = chr_labels$midpoint / 1e6, labels = chr_labels$POS) +
  labs(
    title = "Genome-wide FST",
    x = "Chromosome",
    y = "FST"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
#remove NA data
fst_data <- fst_data %>%
  filter(!is.na(WEIR_AND_COCKERHAM_FST))
ggplot(fst_data, aes(x = cum_POS / 1e6, y = WEIR_AND_COCKERHAM_FST, color = as.factor(CHROM))) +
  geom_point(size = 0.5) +
  scale_color_manual(values = rep(c("blue", "red"), length(unique(fst_data$CHROM)))) +
  scale_x_continuous(breaks = chr_labels$midpoint / 1e6, labels = chr_labels$POS) +
  scale_y_continuous(expand = c(0, 0)) +  # Start y-axis at 0
  labs(
    title = "Genome-wide FST",
    x = "Chromosome Position in Mb",
    y = expression(F[ST])  # Use LaTeX-style formatting for FST
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()  # Remove grid lines
  )
ggplot(fst_data, aes(x = cum_POS / 1e6, y = WEIR_AND_COCKERHAM_FST, color = as.factor(CHROM))) +
  geom_point(size = 0.5) +
  scale_color_manual(values = rep(c("blue", "red"), length(unique(fst_data$CHROM)))) +
  scale_x_continuous(
    breaks = seq(0, max(fst_data$cum_POS) / 1e6, by = 20),  # Tick marks every 5 Mb
    labels = scales::comma  # Format x-axis labels as numbers with commas
  ) +
  scale_y_continuous(
    limits = c(0, 1),  # Start y-axis at 0 and auto-adjust maximum
    expand = c(0, 0)    # Remove extra padding around the y-axis
  ) +
  labs(
    title = "Genome-wide FST",
    x = "Genomic Position (Mb)",
    y = expression(F[ST])  # Use LaTeX-style formatting for FST
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()  # Remove grid lines
  )
fst_data <- fst_data %>%
  group_by(CHROM) %>%
  mutate(
    chr_start = min(POS),   # Start position of each chromosome
    chr_end = max(POS),     # End position of each chromosome
    chr_offset = lag(cumsum(max(POS)), default = 0),  # Offset for cumulative positions
    cum_POS = POS + chr_offset  # Cumulative position for each SNP
  ) %>%
  ungroup()
ggplot(fst_data, aes(x = cum_POS, y = WEIR_AND_COCKERHAM_FST, color = as.factor(CHROM))) +
  geom_point(size = 0.5) +
  scale_color_manual(values = rep(c("blue", "red"), length(unique(fst_data$CHROM)))) +
  scale_x_continuous(
    breaks = chr_labels$midpoint,  # Use chromosome midpoints for breaks
    labels = chr_labels$CHROM,    # Label chromosomes in ascending order
    expand = c(0, 0)              # Remove padding on x-axis
  ) +
  scale_y_continuous(
    limits = c(0, 1),            # Start y-axis at 0
    expand = c(0, 0)              # Remove padding on y-axis
  ) +
  labs(
    title = "Genome-wide FST Manhattan Plot",
    x = "Chromosome Position",
    y = expression(F[ST])  # Format FST using LaTeX-style
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",         # Remove legend
    panel.grid = element_blank(),     # Remove grid lines
    axis.text.x = element_text(angle = 0, hjust = 1),  # Rotate x-axis labels
    axis.ticks.x = element_line()    # Ensure x-axis ticks are visible
  )
# Summary statistics
summary(fst_data$WEIR_AND_COCKERHAM_FST)

# Visualize distribution
ggplot(fst_data, aes(x = WEIR_AND_COCKERHAM_FST)) +
  geom_histogram(bins = 30, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of FST Values", x = expression(F[ST]), y = "Count") +
  theme_minimal()
fst_threshold <- quantile(fst_data$WEIR_AND_COCKERHAM_FST, 0.97, na.rm = TRUE)

## FST outliers

fst_outliers <- fst_data %>%
  filter(WEIR_AND_COCKERHAM_FST >= fst_threshold)
head(fst_outliers)
ggplot(fst_data, aes(x = cum_POS, y = WEIR_AND_COCKERHAM_FST, color = as.factor(chr_color))) +
  geom_point(size = 0.5) +
  geom_point(data = fst_outliers, aes(x = cum_POS, y = WEIR_AND_COCKERHAM_FST), color = "red", size = 1) +
  labs(
    title = "Genome-wide FST with Outliers",
    x = "Genomic Position (Mb)",
    y = expression(F[ST])
  ) +
    scale_color_manual(values = c("black", "gray")) +
  scale_y_continuous(
    limits = c(0, 1),            # Start y-axis at 0
    expand = c(0, 0)              # Remove padding on y-axis
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
outlier_summary <- fst_outliers %>%
  group_by(CHROM) %>%
  summarize(
    num_outliers = n(),
    max_FST = max(WEIR_AND_COCKERHAM_FST),
    peak_position = POS[which.max(WEIR_AND_COCKERHAM_FST)]
  ) %>%
  arrange(desc(max_FST))
print(outlier_summary)

## outlier peaks
window_size <- 1
fst_data <- fst_data %>%
  mutate(POS_bp = POS * 1e6)
regions_of_interest <- fst_data %>%
  filter(
    CHROM %in% unique(outlier_summary$CHROM),
    POS >= (outlier_summary$peak_position - window_size),
    POS <= (outlier_summary$peak_position + window_size)
  )

# FST Windowed

ggplot(fst_data, aes(x = POS, y = WEIR_AND_COCKERHAM_FST, color = as.factor(chr_color))) +
  geom_point(alpha = 0.7, size = 0.8) +
  facet_wrap(~CHROM, scales = "free_x", nrow = 1) +
  geom_hline(yintercept = 0.6, color = "red", linetype = "dashed") +
  scale_color_manual(values = c("black", "gray")) +
  labs(
    title = "FST per variant",
    x = "Genomic Position (Mb)",
    y = expression(F[ST])
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
## USE THIS PLOT FOR FST_DATA

ggplot(fst_data, aes(x = cum_POS, y = WEIR_AND_COCKERHAM_FST, color = as.factor(chr_color))) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("black", "gray")) +
  scale_y_continuous(
    limits = c(0, 1),            # Start y-axis at 0
    expand = c(0, 0)              # Remove padding on y-axis
  ) +
  labs(
    title = "Genome-wide FST Manhattan Plot",
    x = "Chromosome Position (in Mb)",
    y = expression(F[ST])  # Format FST using LaTeX-style
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",         # Remove legend
    panel.grid = element_blank(),     # Remove grid lines
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
  )

## FST Windowed
library(ggplot2)
file.choose()
fst_window <- read.table("/Users/carlaperscky/Desktop/wd/Rstudio2024/fst_sliding_window.windowed.weir.fst", header = TRUE)
fst_window$CHROM <- as.factor(fst_window$CHROM)
fst_window$chr_color <- as.numeric(factor(fst_window$CHROM)) %% 2
ggplot(fst_window, aes(x = BIN_START, y = MEAN_FST, color = as.factor(chr_color))) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("black", "gray")) +
  scale_x_continuous(
    breaks = chr_labels$midpoint,  # Use chromosome midpoints for breaks
    labels = chr_labels$CHROM,    # Label chromosomes in ascending order
    expand = c(0, 0)              # Remove padding on x-axis
  ) +
  scale_y_continuous(
    limits = c(0, 1),            # Start y-axis at 0
    expand = c(0, 0)              # Remove padding on y-axis
  ) +
  labs(
    title = "Mean FST Manhattan Plot",
    x = "Chromosome Position (in Mb)",
    y = expression(F[ST])  # Format FST using LaTeX-style
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",         # Remove legend
    panel.grid = element_blank(),     # Remove grid lines
    axis.text.x = element_text(angle = 0, hjust = 1),  # Rotate x-axis labels
    axis.ticks.x = element_line()    # Ensure x-axis ticks are visible
  )



