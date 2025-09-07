
# Step 1: Get all sv_overlap files
sv_files <- list.files("~/tre1_eqtl/alleqtls_analysis/testset_07.30/SV_overlap_test_all", pattern = ".*\\_allsvs.txt$", full.names = TRUE)

# Step 2: Collect all allele–SNP–SV–genotype records
allele_sv_records <- list()

for (file in sv_files) {
  sv_data <- fread(file, data.table = FALSE)
  
  if (nrow(sv_data) == 0) next
  
  # Extract genotype from filename
  genotype <- sub("(.*)_allsvs\\.txt$", "\\1", basename(file))
  
  
  # Define SV identity as combination of start, end, and feature
  sv_data$sv_id <- paste0(sv_data$sv_start, "_", sv_data$sv_end, "_", sv_data$sv_feature)
  
  allele_sv_records[[length(allele_sv_records) + 1]] <- data.frame(
    top_SNP = sv_data$SNP,
    genotype = genotype,
    allele = sv_data$Allele,
    sv_id = sv_data$sv_id,
    sv_feature = sv_data$sv_feature,
    gene = sv_data$gene,
    stringsAsFactors = FALSE
  )
}
n_distinct(allele_sv_records$gene)

head(allele_sv_records)
##assigning a unique name to unique sv id for our understanding

library(dplyr)
library(stringr)

allele_sv_records <- do.call(rbind, allele_sv_records)

allele_sv_records <- allele_sv_records %>%
  mutate(
    chr = sub("_.*", "", top_SNP),
    sv_start = as.numeric(sub("_.*", "", sv_id)),
    sv_end = as.numeric(str_extract(sv_id, "(?<=_)[0-9]+(?=_)"))
  )

# Create unique list of svs and sort by chr number and start
unique_svs <- allele_sv_records %>%
  select(chr, sv_start, sv_end) %>%
  distinct() %>%
  mutate(chr_num = as.numeric(str_extract(chr, "\\d+"))) %>%
  arrange(chr_num, sv_start) %>%
  mutate(sv_name = paste0("sv", row_number()))

allele_sv_records <- allele_sv_records %>%
  left_join(unique_svs, by = c("chr", "sv_start", "sv_end"))



library(dplyr)

setwd("~/tre1_eqtl/alleqtls_analysis/noneandtrans_control/SV_overlap_control_all")

colnames(merged_data_none) = c("genotype","allele", "top_SNP" ,"gene")

# Left join: keep all rows from merged_data
full_data <- merged_data_none %>%
  left_join(allele_sv_records, by = c("genotype", "allele", "top_SNP", "gene"))

n_distinct(full_data$gene)
write.csv(full_data, "full_datasvstestallsvs.csv")

full_data_test = fread("~/tre1_eqtl/alleqtls_analysis/test_set/full_datasvsallgeno.csv")

full_data2 <- full_data %>% filter(!genotype %in% c("B104", "B73", "Mo17"))

n_distinct(full_data2$gene)

head(full_snp_sv_df)
library(dplyr)
#here we are filtering to keep only those snps that have ref and alternate allele in our dataset because in some cases all the genotypes would have the same allele
# Find top_SNPs with more than one unique allele
snps_with_multiple_alleles <- full_data2 %>%
  distinct(top_SNP, allele) %>%
  group_by(top_SNP) %>%
  summarise(n_alleles = n(), .groups = "drop") %>%
  filter(n_alleles > 1)

# filter to keep only those snps and then filter for heterozygosity
filtered_df <- allele_sv_records %>%
  filter(top_SNP %in% snps_with_multiple_alleles$top_SNP)

hets = filtered_df %>%
  filter(str_detect(allele, ":"))

n_distinct(hets$top_SNP)
n_distinct(hets$gene)

filtered_df <- filtered_df %>%
  filter(!str_detect(allele, ":"))

n_distinct(filtered_df$top_SNP)
n_distinct(filtered_df$gene)

n_distinct(full_data2$top_SNP)
n_distinct(full_data2$gene)
library(dplyr)

library(dplyr)

##next step is again filter for maf
# Count genotypes per SNP × allele
allele_counts <- filtered_df %>%
  distinct(top_SNP, genotype, allele) %>%
  group_by(top_SNP, allele) %>%
  summarise(n_genotypes = n(), .groups = "drop")

# Keep only allele groups with ≥3 genotypes
valid_alleles <- allele_counts %>%
  filter(n_genotypes >= 3)

# Count how many such alleles each SNP has
valid_snps <- valid_alleles %>%
  group_by(top_SNP) %>%
  summarise(n_alleles = n(), .groups = "drop") %>%
  filter(n_alleles > 1)  # Must have exactly two alleles meeting the ≥3 threshold


filtered_df2 <- filtered_df %>%
  semi_join(valid_alleles, by = c("top_SNP", "allele")) %>%
  filter(top_SNP %in% valid_snps$top_SNP)

n_distinct(filtered_df2$top_SNP)

head(filtered_df2)
head(filtered_df2)
hets_check2 = filtered_df2 %>%
  filter(str_detect(allele, ":"))

library(dplyr)
library(tidyr)

##usually all the nones get filtered by now but if some allele for a gene has na assign none


# Replace NA sv_name with "none" so it can be treated as a category
df <- filtered_df2 %>%
  mutate(sv_name = ifelse(is.na(sv_name), "none", sv_name))

# Get list of SVs per allele per SNP
allele_sv_summary <- df %>%
  group_by(top_SNP, allele) %>%
  summarise(
    svs = list(sort(unique(sv_name))),
    n_svs = sum(sv_name != "none"),
    .groups = "drop"
  )

# Pivot to compare both alleles
allele_sv_wide <- allele_sv_summary %>%
  group_by(top_SNP) %>%
  mutate(allele_index = paste0("allele_", row_number())) %>%
  ungroup() %>%
  pivot_wider(
    names_from = allele_index,
    values_from = c(allele, svs, n_svs),
    names_glue = "{.value}_{allele_index}"
  )

# Compare SV profiles across alleles
compare_svs <- allele_sv_wide %>%
  rowwise() %>%
  mutate(
    shared_svs = list(intersect(svs_allele_1, svs_allele_2)),
    unique_to_allele_1 = list(setdiff(svs_allele_1, svs_allele_2)),
    unique_to_allele_2 = list(setdiff(svs_allele_2, svs_allele_1))
  ) %>%
  mutate(
    sv_diff_type = case_when(
      identical(svs_allele_1, svs_allele_2) ~ "same_svs",
      length(unlist(shared_svs)) == 0 ~ "completely_different_svs",
      TRUE ~ "partial_overlap"
    )
  ) %>%
  ungroup()

n_distinct(compare_svs$top_SNP)
head(compare_svs, 10)
head(allele_sv_summary)

library(ggplot2)


# Calculate % for each category
plot_data <- compare_svs %>%
  count(sv_diff_type) %>%
  mutate(percentage = round(100 * n / sum(n), 1))

# plot
ggplot(plot_data, aes(x = sv_diff_type, y = percentage, fill = sv_diff_type)) +
  geom_col() +
  geom_text(aes(label = paste0(percentage, "%")), vjust = -0.5, size = 4.5) +
  labs(
    title = "SV Profile Differences Between Alleles",
    x = "SV Difference Type",
    y = "Percentage of SNPs",
    fill = "SV Difference"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 14)



ggplot(compare_svs, aes(x = sv_diff_type, fill = sv_diff_type)) +
  geom_bar() +
  labs(
    title = "SV Profile Differences Between Alleles",
    x = "SV Difference Type",
    y = "Number of SNPs",
    fill = "SV Difference"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_brewer(palette = "Set2")





write.csv(full_data, "full_data_control_allsvs.csv")

library(dplyr)

compare_svs_export <- compare_svs %>%
  mutate(across(where(is.list), ~ sapply(., function(x) paste(x, collapse = ";"))))

# Now write to CSV
write.csv(compare_svs_export, "sv_comparisonresults_all_test.csv", row.names = FALSE)


