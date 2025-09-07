library(data.table)
library(dplyr)
library(readr)
library(ggplot2)
setwd("~/tre1_eqtl/alleqtls_analysis/analysis_08.25")
write.csv(map_dt, "map_dt.csv")

control_positions = fread("~/tre1_eqtl/alleqtls_analysis/finalsets/noneandtrans_control/control_positions.csv", data.table = F)

control_positions = control_positions %>% mutate(ext_start = pmax(0, start - 1e6), ext_end = end + 1e6 )

map_dt = fread("~/tre1_eqtl/alleqtls_analysis/analysis_08.25/map_dt.csv", data.table = F)

cand <- copy(control_positions) %>% as.data.table()

cand = cand %>% select(gene, chr, ext_start, ext_end)
#filtering snps in the extended overlap region
snps_iv <- copy(map_dt)[, `:=`(ext_start = pos, ext_end = pos)]
str(cand)
str(snps_iv)
snps_iv$chr = as.numeric(snps_iv$chr)

setkey(cand, chr, ext_start, ext_end)
setkey(snps_iv, chr,ext_start, ext_end)

cand_pool <- foverlaps(snps_iv, cand, nomatch = 0L)  # SNPs that fall in any gene window

#getting the MAF of both from plink
##plink –bfile  "bugeater_693geno_filter2_maf_het"--exclude “cissnps.txt”--recode A--out “all_snps_outputexceptcis"
## plink --bfile bugeater_693geno_filter2_maf_het  --extract control_
##selected_quantiles.txt --freq --out control_selected_quantiles

cis_maf <- read_table("~/tre1_eqtl/alleqtls_analysis/finalsets/testset_07.30/cis_SNP_freq.frq") %>%
  select(SNP, MAF) %>% mutate(SNP = trimws(SNP))

all_maf <- read_table("~/tre1_eqtl/alleqtls_analysis/analysis_08.25/all_snps_freq_exceptcis.frq") %>%  # or freq for all SNPs in your map
  select(SNP, MAF) %>% mutate(SNP = trimws(SNP))

# Attach MAF to control candidate pool
cand_pool <- cand_pool %>%
  as.data.frame() %>%
  left_join(all_maf, by = "SNP") %>%
  filter(!is.na(MAF))

# Define MAF bins and compute target & pool proportions
breaks <- quantile(cis_maf$MAF, probs = seq(0, 1, by = 0.05), na.rm = TRUE)
breaks[1] <- 0; breaks[length(breaks)] <- 0.5   # clamp

cis_bins <- cis_maf %>%
  mutate(bin = cut(MAF, breaks, include.lowest = TRUE, right = FALSE))

pool_bins <- cand_pool %>%
  mutate(bin = cut(MAF, breaks, include.lowest = TRUE, right = FALSE))

# proportions
target_p <- cis_bins %>%
  count(bin, name = "n") %>%
  mutate(p = n / sum(n)) %>%
  select(bin, p)

pool_p <- pool_bins %>%
  count(bin, name = "n") %>%
  mutate(p = n / sum(n)) %>%
  select(bin, p)

# merge to get weights per bin: w_bin = target_p / pool_p
weights_df <- full_join(target_p, pool_p, by = "bin", suffix = c("_target", "_pool")) %>%
  mutate(w = ifelse(is.na(p_target) | is.na(p_pool) | p_pool == 0, 0, p_target / p_pool)) %>%
  select(bin, w)

# Assign per-SNP sampling weights and sample 1 SNP per gene
cand_pool <- pool_bins %>%
  left_join(weights_df, by = "bin") %>%
  mutate(w = ifelse(is.na(w) | !is.finite(w), 0, w))

# Some bins may have zero weight if absent in cis set; add a tiny floor to avoid zero-probability genes
w_floor <- 1e-6
cand_pool <- cand_pool %>%
  mutate(w = pmax(w, w_floor))

# Sample 1 SNP per gene, probability ∝ weight
set.seed(42)
control_selected <- cand_pool %>%
  group_by(gene) %>%
  # normalize weights within gene
  mutate(pw = w / sum(w)) %>%
  slice_sample(n = 1, weight_by = pw, replace = FALSE) %>%
  ungroup() %>%
  select(gene, SNP, chr, pos, MAF)

#compare distributions
# KS test
ks <- ks.test(cis_maf$MAF, control_selected$MAF)

# density plot
plot_df <- bind_rows(
  cis_maf %>% mutate(dataset = "cis"),
  control_selected %>% select(MAF) %>% mutate(dataset = "control")
)



write.csv(control_selected, 'control_selected_quantiles.csv')
quantiles = fread("~/tre1_eqtl/alleqtls_analysis/analysis_08.25/control_selected_quantiles.csv", data.table = F)

plot_df2 = bind_rows(cis_maf %>% mutate(dataset = "cis"),
                     quantiles %>% select(MAF) %>% mutate(dataset = "control")
)

#plot to see the MAF density distribution
ggplot(plot_df2, aes(MAF, color = dataset, fill = dataset)) +
  geom_density(alpha = 0.3) +
  labs(x = "Minor Allele Frequency", y = "Density",
       subtitle = paste0("KS p = ", signif(ks$p.value, 3))) +
  theme_minimal(base_size = 12)














