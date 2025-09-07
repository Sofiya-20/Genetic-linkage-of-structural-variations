full_keys <- geno_sets %>% distinct(top_SNP, allele, genotype)

# MID-FREQUENCY
geno_sets_mid <- geno_sets %>%
  unnest_longer(svs, values_to = "sv_name") %>%
  mutate(sv_name = na_if(sv_name, ""),
         sv_name = ifelse(sv_name == "none", NA_character_, sv_name)) %>%
  semi_join(sv_mid %>% select(top_SNP, sv_name),
            by = c("top_SNP", "sv_name")) %>%
  group_by(top_SNP, allele, genotype) %>%
  summarise(svs = list(sort(unique(sv_name))), .groups = "drop") %>%
  right_join(full_keys, by = c("top_SNP","allele","genotype")) %>%
  mutate(svs = lapply(svs, function(x) if (is.null(x) || length(x)==0 || all(is.na(x))) "none" else x))


perfect_on_allele <- function(geno_sets_filtered,
                              allele_sizes_unfiltered,
                              majors,
                              fix_cutoff = 1.0  # 1.0 = strict; try 0.95 for near-fix
){
  sv_long <- geno_sets_filtered %>%
    unnest_longer(svs, values_to = "sv_name") %>%
    mutate(sv_name = na_if(sv_name, ""),
           sv_name = ifelse(sv_name == "none", NA_character_, sv_name)) %>%
    filter(!is.na(sv_name)) %>%
    distinct(top_SNP, allele, genotype, sv_name)
  
  sv_space <- sv_long %>% distinct(top_SNP, sv_name)
  
  sv_by_allele <- sv_long %>%
    count(top_SNP, sv_name, allele, name = "n_present") %>%
    left_join(allele_sizes_unfiltered, by = c("top_SNP","allele")) %>%
    mutate(prop_present = n_present / pmax(allele_size, 1L))
  
  sv_per_allele <- sv_by_allele %>%
    group_by(top_SNP, sv_name) %>%
    mutate(
      total_present      = sum(n_present),
      present_in_all     = prop_present >= fix_cutoff & n_present > 0,
      present_in_others  = total_present - n_present,
      absent_in_others   = present_in_others == 0,
      perfect_on_this_allele = present_in_all & absent_in_others
    ) %>% ungroup()
  
  sv_perfect_any_allele <- sv_per_allele %>%
    filter(perfect_on_this_allele) %>%
    left_join(majors, by = "top_SNP") %>%
    mutate(matches_major = allele == major_allele) %>%
    transmute(
      top_SNP, sv_name,
      matched_allele = allele,
      matches_major,
      carriers = n_present,
      allele_size,
      n_other_carriers = present_in_others
    ) %>%
    arrange(top_SNP, sv_name, matched_allele)
  
  list(
    sv_space = sv_space,
    sv_per_allele = sv_per_allele,
    sv_perfect_any_allele = sv_perfect_any_allele
  )
}
mid_res  <- perfect_on_allele(geno_sets_mid, allele_sizes, majors, fix_cutoff = 1.00)


summarize_perfects <- function(res){
  overall_total_pairs   <- nrow(res$sv_space)
  overall_perfect_pairs <- res$sv_perfect_any_allele %>% distinct(top_SNP, sv_name) %>% nrow()
  overall_major_pairs   <- sum(res$sv_perfect_any_allele$matches_major)
  tibble(
    total_snp_sv_pairs_considered = overall_total_pairs,
    n_perfect_pairs_any_allele    = overall_perfect_pairs,
    pct_perfect_any_allele        = round(100 * overall_perfect_pairs / pmax(overall_total_pairs, 1L), 2),
    n_perfect_pairs_match_major   = overall_major_pairs,
    n_perfect_pairs_match_minor   = overall_perfect_pairs - overall_major_pairs
  )
}

overall_summary_mid   <- summarize_perfects(mid_res)
overall_summary_ext   <- summarize_perfects(ext_res)
overall_summary_near  <- summarize_perfects(near_res)






write.csv(overall_summary_mid, "overall_summary_mid_test08.30.csv")






