
library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)

# Load data
filtered_test <- fread("~/Documents/GitHub/New_results/Filteredforeverything_JAMES08.28.csv", 
                       data.table = FALSE, header = TRUE)
n_distinct(filtered_test$top_SNP)
# Remove first two columns and normalize NA/"" -> "none"
filtered_test <- filtered_test[, -(1:2)]
df <- filtered_test %>%
  mutate(sv_name = coalesce(na_if(sv_name, ""), "none"))

# Get total unique SNPs for percentage calculation
total_snps <- n_distinct(df$top_SNP)

# 1) Create genotype-level SV lists
exclude_none <- FALSE
geno_sets <- df %>%
  distinct(top_SNP, allele, genotype, sv_name) %>%
  { if (exclude_none) filter(., sv_name != "none") else . } %>%
  group_by(top_SNP, allele, genotype) %>%
  summarise(svs = list(sort(unique(sv_name))), .groups = "drop")

# 2) Calculate SV frequency across genotypes (per SNP)
low_thr <- 3   # Minimum genotypes an SV must appear in
high_thr <- 11  # Maximum genotypes an SV can appear in (excludes those in 12-14)

sv_counts <- geno_sets %>%
  select(top_SNP, genotype, svs) %>%
  unnest_longer(svs, values_to = "sv_name") %>%
  mutate(sv_name = na_if(sv_name, ""),
         sv_name = ifelse(sv_name == "none", NA_character_, sv_name)) %>%
  filter(!is.na(sv_name)) %>%
  distinct(top_SNP, genotype, sv_name) %>%
  count(top_SNP, sv_name, name = "n_genotypes")

# FIXED: Use AND instead of OR for mid-frequency filtering
sv_mid <- sv_counts %>%
  filter(n_genotypes >= low_thr & n_genotypes <= high_thr)

# 3) Keep ONLY mid-frequency SVs
geno_sets_mid <- geno_sets %>%
  unnest_longer(svs, values_to = "sv_name") %>%
  mutate(sv_name = na_if(sv_name, ""),
         sv_name = ifelse(sv_name == "none", NA_character_, sv_name)) %>%
  semi_join(sv_mid %>% select(top_SNP, sv_name), 
            by = c("top_SNP", "sv_name")) %>%
  group_by(top_SNP, allele, genotype) %>%
  summarise(svs = list(sort(unique(sv_name))), .groups = "drop") %>%
  mutate(svs = lapply(svs, function(x) if (length(x) == 0) "none" else x))

# 4) Calculate allele sizes
allele_sizes <- geno_sets %>%
  distinct(top_SNP, allele, genotype) %>%
  count(top_SNP, allele, name = "allele_size")

# Identify major allele for each SNP
majors <- allele_sizes %>%
  group_by(top_SNP) %>%
  slice_max(order_by = allele_size, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(top_SNP, major_allele = allele)

# 5) Perfect-on-allele test using mid-frequency set
sv_long <- geno_sets_mid %>%
  unnest_longer(svs, values_to = "sv_name") %>%
  mutate(sv_name = na_if(sv_name, ""),
         sv_name = ifelse(sv_name == "none", NA_character_, sv_name)) %>%
  filter(!is.na(sv_name)) %>%
  distinct(top_SNP, allele, genotype, sv_name)

# Count carriers within each allele
sv_by_allele <- sv_long %>%
  count(top_SNP, sv_name, allele, name = "n_present") %>%
  left_join(allele_sizes, by = c("top_SNP", "allele"))

# Define perfect matches
sv_per_allele <- sv_by_allele %>%
  group_by(top_SNP, sv_name) %>%
  mutate(
    total_present = sum(n_present),
    present_in_all = n_present == allele_size & n_present > 0,
    present_in_others = total_present - n_present,
    absent_in_others = present_in_others == 0,
    perfect_on_this_allele = present_in_all & absent_in_others
  ) %>%
  ungroup()

# Keep only perfect matches
sv_perfect_any_allele <- sv_per_allele %>%
  filter(perfect_on_this_allele) %>%
  left_join(majors, by = "top_SNP") %>%
  mutate(matches_major = allele == major_allele) %>%
  select(
    top_SNP, sv_name, matched_allele = allele,
    matches_major, carriers = n_present,
    allele_size, n_other_carriers = present_in_others
  ) %>%
  arrange(top_SNP, sv_name, matched_allele)

# ========== NEW AGGREGATION AT SNP LEVEL ==========

# Aggregate perfect matches at SNP level
snp_level_summary <- sv_perfect_any_allele %>%
  group_by(top_SNP) %>%
  summarise(
    n_perfect_svs = n_distinct(sv_name),
    n_major_allele_svs = sum(matches_major),
    n_minor_allele_svs = sum(!matches_major),
    has_major_perfect = any(matches_major),
    has_minor_perfect = any(!matches_major),
    # List all perfect SVs for this SNP
    perfect_svs_list = list(unique(sv_name)),
    # Detailed info for each SV-allele pair
    sv_allele_details = list(
      tibble(
        sv = sv_name,
        allele = matched_allele,
        is_major = matches_major
      )
    )
  ) %>%
  # Classify SNP based on which allele has perfect matches
  mutate(
    perfect_type = case_when(
      has_major_perfect & has_minor_perfect ~ "both_alleles",
      has_major_perfect ~ "major_only",
      has_minor_perfect ~ "minor_only",
      TRUE ~ "none"
    )
  )

# Calculate overall summary statistics based on UNIQUE SNPs
n_snps_with_perfect <- nrow(snp_level_summary)
n_snps_major_only <- sum(snp_level_summary$perfect_type == "major_only")
n_snps_minor_only <- sum(snp_level_summary$perfect_type == "minor_only")
n_snps_both <- sum(snp_level_summary$perfect_type == "both_alleles")

overall_summary <- tibble(
  total_snps_in_dataset = total_snps,
  n_snps_with_perfect_sv_association = n_snps_with_perfect,
  pct_snps_with_perfect = round(100 * n_snps_with_perfect / total_snps, 2),
  n_snps_perfect_major_only = n_snps_major_only,
  n_snps_perfect_minor_only = n_snps_minor_only,
  n_snps_perfect_both_alleles = n_snps_both,
  pct_major_only = round(100 * n_snps_major_only / total_snps, 2),
  pct_minor_only = round(100 * n_snps_minor_only / total_snps, 2),
  pct_both_alleles = round(100 * n_snps_both / total_snps, 2),
  total_sv_allele_pairs = nrow(sv_perfect_any_allele),
  avg_svs_per_snp = round(mean(snp_level_summary$n_perfect_svs), 2)
)

# Create a detailed report for each SNP with its perfect SV associations
snp_sv_catalog <- snp_level_summary %>%
  select(top_SNP, perfect_type, n_perfect_svs, sv_allele_details) %>%
  unnest(sv_allele_details) %>%
  arrange(top_SNP, sv, allele)

# Print summaries
cat("\n======== OVERALL SUMMARY ========\n")
print(overall_summary)

cat("\n======== SNP TYPE DISTRIBUTION ========\n")
print(table(snp_level_summary$perfect_type))

cat("\n======== TOP 10 SNPs BY NUMBER OF PERFECT SVs ========\n")
print(snp_level_summary %>% 
        arrange(desc(n_perfect_svs)) %>% 
        select(top_SNP, n_perfect_svs, perfect_type) %>%
        head(10))

setwd("~/Documents/GitHub/New_results")

# Save results
write.csv(overall_summary, "overall_summary_snp_level_test10.16.csv", row.names = FALSE)
write.csv(snp_level_summary, "snp_level_perfect_associations_test10.16.csv", row.names = FALSE)
write.csv(snp_sv_catalog, "detailed_snp_sv_catalog_test10.16.csv", row.names = FALSE)
write.csv(sv_perfect_any_allele, "all_perfect_sv_allele_pairs_test10.16.csv", row.names = FALSE)

# Optional: Create a simple summary showing one row per SNP with concatenated SV names
snp_simple_summary <- snp_level_summary %>%
  mutate(
    perfect_svs_concat = sapply(perfect_svs_list, paste, collapse = ";")
  ) %>%
  select(top_SNP, perfect_type, n_perfect_svs, perfect_svs_concat) %>%
  arrange(top_SNP)

write.csv(snp_simple_summary, "snp_simple_summary.csv", row.names = FALSE)
