library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)  # For any pivoting if needed

# Replace with your file path
file_path <- 'modkit_pileup/A549_rep1_all_mods_aligned_pileup.bed'

# Column names
column_names <- c(
  'chrom', 'start', 'end', 'modified_base_code', 'score', 'strand',
  'thick_start', 'thick_end', 'color', 'valid_coverage',
  'percent_modified', 'n_mod', 'n_canonical', 'n_other_mod',
  'n_delete', 'n_fail', 'n_diff', 'n_nocall'
)

# Load the data (assuming tab-separated, skip comments if any)
df <- read_tsv(file_path, col_names = column_names, comment = '#')

# Display the first few rows
head(df)

# Basic info (similar to df.info() in Python)
str(df)

# Summary statistics
summary(df)

# Unique modification codes
unique_mod_codes <- unique(df$modified_base_code)
print(paste("Unique modification codes:", paste(unique_mod_codes, collapse = ", ")))

# Unique chromosomes
unique_chroms <- unique(df$chrom)
print(paste("Unique chromosomes:", paste(unique_chroms, collapse = ", ")))

# Check for missing values
missing_per_col <- colSums(is.na(df))
print("Missing values per column:")
print(missing_per_col)

# Filter for minimum coverage
min_coverage <- 10
filtered_df <- df %>% filter(valid_coverage >= min_coverage)

print(paste("Original rows:", nrow(df), ", Filtered rows:", nrow(filtered_df)))

# Histogram of Percent Modified
ggplot(filtered_df, aes(x = percent_modified)) +
  geom_histogram(bins = 50, fill = "blue", color = "black") +
  labs(title = "Distribution of Percent Modified", x = "Percent Modified", y = "Frequency") +
  theme_minimal()

# Boxplot of Percent Modified by Modification Code
ggplot(filtered_df, aes(x = modified_base_code, y = percent_modified)) +
  geom_boxplot(fill = "lightblue") +
  labs(title = "Percent Modified by Modification Code", x = "Modification Code", y = "Percent Modified") +
  theme_minimal()

# Boxplot of Percent Modified by Chromosome
ggplot(filtered_df, aes(x = chrom, y = percent_modified)) +
  geom_boxplot(fill = "lightgreen") +
  labs(title = "Percent Modified by Chromosome", x = "Chromosome", y = "Percent Modified") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Select a chromosome (change as needed, taking the first unique one)
chrom_to_plot <- unique(filtered_df$chrom)[1]
chrom_df <- filtered_df %>% filter(chrom == chrom_to_plot)

# Scatter Plot of Modification Along a Chromosome
ggplot(chrom_df, aes(x = start, y = percent_modified)) +
  geom_point(alpha = 0.5, size = 1) +
  labs(title = paste("Percent Modified Along", chrom_to_plot), x = "Position", y = "Percent Modified") +
  theme_minimal()

# Average percent modified per chromosome
avg_per_chrom <- filtered_df %>%
  group_by(chrom) %>%
  summarise(mean_percent_modified = mean(percent_modified, na.rm = TRUE)) %>%
  ungroup()

print("Average Percent Modified per Chromosome:")
print(avg_per_chrom)

# Plot it
ggplot(avg_per_chrom, aes(x = chrom, y = mean_percent_modified)) +
  geom_bar(stat = "identity", fill = "purple") +
  labs(title = "Average Percent Modified per Chromosome", x = "Chromosome", y = "Average Percent Modified") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Average per modification code
avg_per_code <- filtered_df %>%
  group_by(modified_base_code) %>%
  summarise(mean_percent_modified = mean(percent_modified, na.rm = TRUE)) %>%
  ungroup()

print("Average Percent Modified per Code:")
print(avg_per_code)


# explore pseudouridine only
pseudouridine_df <- filtered_df %>% filter(modified_base_code == "17802")

# histogram of pseudouridine
ggplot(pseudouridine_df, aes(x = percent_modified)) +
  geom_histogram(bins = 50, fill = "blue", color = "black") +
  labs(title = "Distribution of Percent Modified for Pseudouridine", x = "Percent Modified", y = "Frequency") +
  theme_minimal()

# boxplot of pseudouridine by chromosome
# simplify chorm names, use fivth part sep by "|"
pseudouridine_df$chrom <- sapply(strsplit(pseudouridine_df$chrom, "\\|"), function(x) x[5])

# remove the third part of the chrom name
ggplot(pseudouridine_df, aes(x = chrom, y = percent_modified)) +
  geom_boxplot(fill = "lightgreen") +
  labs(title = "Percent Modified by snRNA variant for Pseudouridine", x = "snRNA variant", y = "Percent Modified") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/pseudouridine_by_variant.png", width = 15, height = 10)
