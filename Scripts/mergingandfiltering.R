library(dplyr)
library(stringr)

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


colnames(merged_data_none) = c("genotype","allele", "top_SNP" ,"gene")

# Left join: keep all rows from merged_data
full_data <- merged_data_none %>%
  left_join(allele_sv_records, by = c("genotype", "allele", "top_SNP", "gene"))

n_distinct(full_data$gene)
write.csv(full_data, "full_datasvstestallsvs.csv")

full_data_test = fread("~/tre1_eqtl/alleqtls_analysis/test_set/full_datasvsallgeno.csv")

# 0) Start from your full data
# full_data2: columns include top_SNP, allele, genotype, sv_name, gene, ...

full_data_test = fread("~/tre1_eqtl/alleqtls_analysis/analysis_08.25/control_selection/control_snps_quantiles/full_data_control_allsvs.csv")

full_data2 <- full_data_test %>% filter(!genotype %in% c("B104", "B73", "Mo17"))

n_distinct(full_data2$genotype)
# 1) Remove heterozygous calls first
no_hets <- full_data2 %>%
  filter(!str_detect(allele, ":"))

# 2) Keep SNPs that are biallelic after het removal
biallelic_snps <- no_hets %>%
  distinct(top_SNP, allele) %>%
  count(top_SNP, name = "n_alleles") %>%
  filter(n_alleles == 2) %>%
  pull(top_SNP)

filtered_df <- no_hets %>%
  filter(top_SNP %in% biallelic_snps)

# 3) Count genotypes per SNP × allele and total per SNP
allele_counts <- filtered_df %>%
  distinct(top_SNP, genotype, allele) %>%
  count(top_SNP, allele, name = "n_genotypes") %>%
  group_by(top_SNP) %>%
  mutate(total = sum(n_genotypes)) %>%
  ungroup()

# 4) Mid-frequency filter: keep alleles with ≥3 and ≤ (total - 3)
#    (When total == 14, this is exactly 3..11)
valid_alleles <- allele_counts %>%
  filter(n_genotypes >= 3, n_genotypes <= pmax(total - 3, 0))

# 5) Keep only SNPs where *exactly two* alleles pass the mid-frequency filter
valid_snps <- valid_alleles %>%
  count(top_SNP, name = "n_valid") %>%
  filter(n_valid == 2) %>%
  pull(top_SNP)

filtered_df2 <- filtered_df %>%
  semi_join(valid_alleles, by = c("top_SNP", "allele")) %>%
  filter(top_SNP %in% valid_snps)

# quick checks
n_distinct(filtered_df2$top_SNP)
n_distinct(filtered_df2$gene)


df <- filtered_test %>%
  mutate(sv_name = ifelse(is.na(sv_name), "none", sv_name))

exclude_none <- FALSE


prana = fread("~/tre1_eqtl/alleqtls_analysis/analysis_08.25/control_selection/control_snps_quantiles/fulldatacontrol_filteredforMAF.csv")


fix_cutoff <- 1.0   # use 1.0 for strict fixation; try 0.95 or 0.90 for near-fixation

geno_sets <- df %>%
  distinct(top_SNP, allele, genotype, sv_name) %>%
  { if (exclude_none) filter(., sv_name != "none") else . } %>%
  group_by(top_SNP, allele, genotype) %>%
  summarise(svs = list(sort(unique(sv_name))), .groups = "drop")

geno_sets <- geno_sets %>%
  mutate(
    svs = lapply(svs, function(x) {
      x <- sort(unique(x))
      x <- x[!is.na(x) & nzchar(x)]   # drop NAs and ""
      if (length(x) == 0) "none" else x
    })
  )

geno_sets <- df %>%
  mutate(sv_name = coalesce(na_if(sv_name, ""), "none")) %>%  # NA/"" -> "none"
  distinct(top_SNP, allele, genotype, sv_name) %>%
  group_by(top_SNP, allele, genotype) %>%
  summarise(svs = list(sort(unique(sv_name))), .groups = "drop")

sum(geno_sets$svs == "none" )


head(geno_sets)

geno_sets %>% head()


head(geno_sets)


## second part testing for perfect match

write.csv(filtered_df2, "Filteredforeverything_JAMES08.28.csv")


