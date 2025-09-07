library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)   

filtered_test <- fread("~/Documents/GitHub/Genetic-linkage-of-structural-variations/Step3_mergingSVsandprefiltering/Test_set/output_files/filtered_final_test/FilteredforMAF08.28.csv",
                                         data.table = FALSE, header = TRUE)
head(filtered_test)
n_distinct(filtered_test$genotype)
filtered_test <- filtered_test[, -(1:2)] 

df <- filtered_test %>%
  mutate(sv_name = coalesce(na_if(sv_name, ""), "none"))  # normalize NA/"" -> "none"
n_distinct(df$sv_name)
exclude_none <- FALSE



# 1) Genotype-level SV lists (unfiltered)

geno_sets <- df %>%
  distinct(top_SNP, allele, genotype, sv_name) %>%
  { if (exclude_none) filter(., sv_name != "none") else . } %>%
  group_by(top_SNP, allele, genotype) %>%
  summarise(svs = list(sort(unique(sv_name))), .groups = "drop")

n_distinct(geno_sets$svs)
n_distinct(df$top_SNP)
n_distinct(geno_sets$top_SNP)

head(geno_sets)

n_unique_ids_in_genosets <- geno_sets %>%
  unnest_longer(svs, values_to = "sv_name") %>%
  filter(!is.na(sv_name), sv_name != "none") %>%
  distinct(sv_name) %>%
  nrow()
n_unique_ids_in_genosets

library(dplyr)
library(tidyr)



# 2) SV frequency across genotypes (per SNP)

low_thr  <- 3
high_thr <- 11

sv_counts <- geno_sets %>%
  select(top_SNP, genotype, svs) %>%
  unnest_longer(svs, values_to = "sv_name") %>%
  mutate(sv_name = na_if(sv_name, ""),
         sv_name = ifelse(sv_name == "none", NA_character_, sv_name)) %>%
  filter(!is.na(sv_name)) %>%
  distinct(top_SNP, genotype, sv_name) %>%     # each (SNP, genotype, SV) once
  count(top_SNP, sv_name, name = "n_genotypes")

sv_extremes <- sv_counts %>% filter(n_genotypes < low_thr | n_genotypes > high_thr)
sv_mid      <- sv_counts %>% filter(n_genotypes >= low_thr, n_genotypes <= high_thr)


# 3 Keep ONLY mid-frequency SVs
geno_sets_mid <- geno_sets %>%
  unnest_longer(svs, values_to = "sv_name") %>%
  mutate(sv_name = na_if(sv_name, ""),
         sv_name = ifelse(sv_name == "none", NA_character_, sv_name)) %>%
  semi_join(sv_mid %>% select(top_SNP, sv_name), by = c("top_SNP", "sv_name")) %>%
  group_by(top_SNP, allele, genotype) %>%
  summarise(svs = list(sort(unique(sv_name))), .groups = "drop") %>%
  mutate(svs = lapply(svs, function(x) if (length(x) == 0) "none" else x))



# 4) Allele sizes (from the UNFILTERED geno_sets!)

allele_sizes <- geno_sets %>%
  distinct(top_SNP, allele, genotype) %>%
  count(top_SNP, allele, name = "allele_size")

allele_sizes %>% arrange(top_SNP, desc(allele_size), allele)

majors <- allele_sizes %>%
  group_by(top_SNP) %>%
  slice_max(order_by = allele_size, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(top_SNP, major_allele = allele)

# 5) Perfect-on-allele test using the mid-frequency set

sv_long <- geno_sets_mid %>%
  unnest_longer(svs, values_to = "sv_name") %>%
  mutate(sv_name = na_if(sv_name, ""),
         sv_name = ifelse(sv_name == "none", NA_character_, sv_name)) %>%
  filter(!is.na(sv_name)) %>%
  distinct(top_SNP, allele, genotype, sv_name)

sv_space <- sv_long %>% distinct(top_SNP, sv_name)

# Carriers within each allele
sv_by_allele <- sv_long %>%
  count(top_SNP, sv_name, allele, name = "n_present") %>%
  left_join(allele_sizes, by = c("top_SNP", "allele"))   # allele_size from unfiltered set

# Perfect definition
sv_per_allele <- sv_by_allele %>%
  group_by(top_SNP, sv_name) %>%
  mutate(
    total_present        = sum(n_present),
    present_in_all       = n_present == allele_size & n_present > 0,
    present_in_others    = total_present - n_present,
    absent_in_others     = present_in_others == 0,
    perfect_on_this_allele = present_in_all & absent_in_others
  ) %>%
  ungroup()

# Keep only perfects; label major/minor
sv_perfect_any_allele <- sv_per_allele %>%
  filter(perfect_on_this_allele) %>%
  left_join(majors, by = "top_SNP") %>%
  mutate(matches_major = allele == major_allele) %>%
  select(
    top_SNP, sv_name,
    matched_allele = allele,
    matches_major,
    carriers = n_present,
    allele_size,
    n_other_carriers = present_in_others
  ) %>%
  arrange(top_SNP, sv_name, matched_allele)

# Optional boolean per (SNP, SV)
sv_any_flag <- sv_per_allele %>%
  group_by(top_SNP, sv_name) %>%
  summarise(
    perfect_any_allele = any(perfect_on_this_allele),
    matched_allele = paste(allele[perfect_on_this_allele], collapse = ";"),
    .groups = "drop"
  )


overall_total_pairs   <- nrow(sv_space)
overall_perfect_pairs <- nrow(distinct(sv_perfect_any_allele, top_SNP, sv_name))
overall_major_pairs   <- sv_perfect_any_allele %>% summarise(n = sum(matches_major)) %>% pull(n)
overall_minor_pairs   <- overall_perfect_pairs - overall_major_pairs

overall_summary <- tibble(
  total_snp_sv_pairs_considered = overall_total_pairs,
  n_perfect_pairs_any_allele    = overall_perfect_pairs,
  pct_perfect_any_allele        = round(100 * overall_perfect_pairs / pmax(overall_total_pairs, 1L), 2),
  n_perfect_pairs_match_major   = overall_major_pairs,
  n_perfect_pairs_match_minor   = overall_minor_pairs
)

per_snp_counts <- sv_perfect_any_allele %>%
  group_by(top_SNP) %>%
  summarise(
    n_perfect_total = n_distinct(sv_name),
    n_perfect_major = sum(matches_major),
    n_perfect_minor = n_perfect_total - n_perfect_major,
    .groups = "drop"
  ) %>%
  left_join(majors, by = "top_SNP")

sv_catalog <- sv_perfect_any_allele %>%
  group_by(top_SNP, matched_allele) %>%
  summarise(
    n_sv = n_distinct(sv_name),
    svs  = list(sort(unique(sv_name))),
    .groups = "drop"
  ) %>%
  arrange(top_SNP, matched_allele)

# Quick peeks
overall_summary
head(per_snp_counts, 10)
sv_catalog %>%
  group_by(top_SNP) %>%
  slice_max(n_sv, n = 5, with_ties = TRUE) %>%
  ungroup()


write.csv(overall_summary, "overall_summary_test07.30.csv", row.names = FALSE)















