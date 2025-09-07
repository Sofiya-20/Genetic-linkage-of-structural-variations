library(dplyr)
library(tidyr)


##final code till 07.27

filtered_test <- fread("~/tre1_eqtl/alleqtls_analysis/analysis_08.25/testset_07.30/fulldatatest_filteredforMAF.csv", data.table = F, header = T)
head(filtered_test)
filtered_test = filtered_test[, -(1:2)]

df <- filtered_test %>%
  mutate(sv_name = ifelse(is.na(sv_name), "none", sv_name))

exclude_none <- FALSE
fix_cutoff <- 1.0   # use 1.0 for strict fixation; try 0.95 or 0.90 for near-fixation
chr10_10031025

geno_sets <- df %>%
  distinct(top_SNP, allele, genotype, sv_name) %>%
  { if (exclude_none) filter(., sv_name != "none") else . } %>%
  group_by(top_SNP, allele, genotype) %>%
  summarise(svs = list(sort(unique(sv_name))), .groups = "drop")


geno_sets <- df %>%
  mutate(sv_name = coalesce(na_if(sv_name, ""), "none")) %>%  # NA/"" -> "none"
  distinct(top_SNP, allele, genotype, sv_name) %>%
  group_by(top_SNP, allele, genotype) %>%
  summarise(svs = list(sort(unique(sv_name))), .groups = "drop")

sum(geno_sets$svs == "none" )


head(geno_sets)

geno_sets %>% head()

head(geno_sets)

low_thr  <- 3
high_thr <- 11

# 1) Count, for each SNP, how many genotypes carry each SV
sv_counts <- geno_sets %>%
  select(top_SNP, genotype, svs) %>%
  tidyr::unnest_longer(svs, values_to = "sv_name") %>%   # if your tidyr is older, use tidyr::unnest(svs)
  mutate(sv_name = na_if(sv_name, ""),
         sv_name = ifelse(sv_name == "none", NA_character_, sv_name)) %>%
  filter(!is.na(sv_name)) %>%
  distinct(top_SNP, genotype, sv_name) %>%  # each (SNP, genotype, SV) counted once
  count(top_SNP, sv_name, name = "n_genotypes")

# 2a) SVs that are present in <3 or >11 genotypes (extremes)
sv_extremes <- sv_counts %>%
  filter(n_genotypes < low_thr | n_genotypes > high_thr)

# 2b) SVs present in 3..11 genotypes (mid-frequency)
sv_mid <- sv_counts %>%
  filter(n_genotypes >= low_thr, n_genotypes <= high_thr)

# 3a) Rebuild geno_sets keeping ONLY extreme SVs
geno_sets_extremes <- geno_sets %>%
  tidyr::unnest_longer(svs, values_to = "sv_name") %>%
  mutate(sv_name = na_if(sv_name, ""),
         sv_name = ifelse(sv_name == "none", NA_character_, sv_name)) %>%
  left_join(sv_extremes, by = c("top_SNP", "sv_name")) %>%
  filter(!is.na(n_genotypes)) %>%  # keep extremes
  group_by(top_SNP, allele, genotype) %>%
  summarise(svs = list(sort(unique(sv_name))), .groups = "drop") %>%
  mutate(svs = lapply(svs, function(x) if (length(x) == 0) "none" else x))

# 3b) Rebuild geno_sets keeping ONLY mid-frequency SVs
geno_sets_mid <- geno_sets %>%
  tidyr::unnest_longer(svs, values_to = "sv_name") %>%
  mutate(sv_name = na_if(sv_name, ""),
         sv_name = ifelse(sv_name == "none", NA_character_, sv_name)) %>%
  left_join(sv_mid, by = c("top_SNP", "sv_name")) %>%
  filter(!is.na(n_genotypes)) %>%  # keep mid-frequency
  group_by(top_SNP, allele, genotype) %>%
  summarise(svs = list(sort(unique(sv_name))), .groups = "drop") %>%
  mutate(svs = lapply(svs, function(x) if (length(x) == 0) "none" else x))

# Optional: a quick table of extremes per SNP
sv_extremes_summary <- sv_extremes %>%
  arrange(top_SNP, desc(n_genotypes))


n_distinct(geno_sets_mid$top_SNP)

library(dplyr)
library(tidyr)

# 0) Allele sizes and (optional) major label
allele_sizes <- geno_sets_mid %>%
  distinct(top_SNP, allele, genotype) %>%
  count(top_SNP, allele, name = "allele_size")

allele_sizes %>% arrange()

## why do I still have allele sizes as 1 when I have filtered everything??

majors <- allele_sizes %>%
  group_by(top_SNP) %>%
  slice_max(order_by = allele_size, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(top_SNP, major_allele = allele)

# 1) Long form SV presence per genotype (drop "none"/blanks)
sv_long <- geno_sets_mid %>%
  unnest_longer(svs, values_to = "sv_name") %>%
  mutate(sv_name = na_if(sv_name, ""),
         sv_name = ifelse(sv_name == "none", NA_character_, sv_name)) %>%
  filter(!is.na(sv_name)) %>%
  distinct(top_SNP, allele, genotype, sv_name)

sv_space = sv_long %>% distinct(top_SNP, sv_name)

# 2) Count carriers of each SV within each allele
sv_by_allele <- sv_long %>%
  count(top_SNP, sv_name, allele, name = "n_present") %>%
  left_join(allele_sizes, by = c("top_SNP","allele"))

# 3) Test "perfect on this allele": present in ALL of its genotypes, NONE of others
sv_per_allele <- sv_by_allele %>%
  group_by(top_SNP, sv_name) %>%
  mutate(
    total_present     = sum(n_present),
    present_in_all    = n_present == allele_size & n_present > 0,
    present_in_others = total_present - n_present,
    absent_in_others  = present_in_others == 0,
    perfect_on_this_allele = present_in_all & absent_in_others
  ) %>%
  ungroup()

# 4) Keep only perfect matches (works for either major or minor)
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
  arrange 


# Optional: a per-(SNP, SV) boolean
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
  pct_perfect_any_allele        = round(100 * overall_perfect_pairs / overall_total_pairs, 2),
  n_perfect_pairs_match_major   = overall_major_pairs,
  n_perfect_pairs_match_minor   = overall_minor_pairs
)

# Per-SNP: how many perfect SVs, split by allele and by major/minor
per_snp_counts <- sv_perfect_any_allele %>%
  group_by(top_SNP) %>%
  summarise(
    n_perfect_total   = n_distinct(sv_name),
    n_perfect_major   = sum(matches_major),
    n_perfect_minor   = n_perfect_total - n_perfect_major,
    .groups = "drop"
  ) %>%
  left_join(majors, by = "top_SNP")

# Compact catalog: list the perfectly matching SVs for each allele at each SNP
sv_catalog <- sv_perfect_any_allele %>%
  group_by(top_SNP, matched_allele) %>%
  summarise(
    n_sv  = n_distinct(sv_name),
    svs   = list(sort(unique(sv_name))),
    .groups = "drop"
  ) %>%
  arrange(top_SNP, matched_allele)

# Quick peeks
overall_summary
head(per_snp_counts, 10)

sv_catalog %>% group_by(top_SNP) %>% slice_max(n_sv, n = 5, with_ties = TRUE) %>% ungroup()


write.csv(overall_summary, "overall_summary_control07.30.csv")
write.csv(geno_sets_mid, "geno_sets_mid_test07.30.csv")

utils::write.table(geno_sets_mid, "geno_sets_mid_test07.30.csv")
