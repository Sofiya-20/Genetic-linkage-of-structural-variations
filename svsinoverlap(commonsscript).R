
control_SNP_alleles = fread("~/tre1_eqtl/alleqtls_analysis/controlandtrans_control/control_snps_alleles.raw", data.table = F)
control_SNP_alleles_frq <- fread("~/tre1_eqtl/alleqtls_analysis/controlandtrans_control/control_SNP_freq.frq", data.table = F)
write.csv( top_SNP_alleles, "control_snp_alleles.csv")


control_SNP_alleles <- control_SNP_alleles[,-c(2:6)]
head(control_SNP_alleles)

# 0 = homozygous REF  
# 1 = heterozygous  
# 2 = homozygous ALT

colnames(conctrol_SNP_alleles) <- sub("_$", "", sub("^(([^_]*_){2}).*", "\\1", colnames(control_SNP_alleles)))


# A1 is the ALT allele = 2
# A2 is the REF allele = 0

colnames(control_SNP_alleles)[2]
for (i in 2:ncol(control_SNP_alleles)) {
  ALT <- control_SNP_alleles_frq$A1[which(control_SNP_alleles_frq$SNP == colnames(control_SNP_alleles)[i])]
  REF <- control_SNP_alleles_frq$A2[which(control_SNP_alleles_frq$SNP == colnames(control_SNP_alleles)[i])]
  HET <- paste0(REF, ":", ALT)
  
  control_SNP_alleles[which(control_SNP_alleles[, i] == 0), i] <- REF
  control_SNP_alleles[which(control_SNP_alleles[, i] == 1), i] <- HET
  control_SNP_alleles[which(control_SNP_alleles[, i] == 2), i] <- ALT
}

write.csv(control_SNP_alleles, "control_snps_alleles.csv")

library(dplyr)
sampled_genenames_control = unique_snps %>% select(top_SNP, gene)

## I have the alleles associate with the snps filtering to keep only the NAM parents

common = fread("~/tre1_eqtl/genotypescommon.txt", header = FALSE)

subset_control_SNP_alleles <- list()

snp_names <- colnames(top_SNP_alleles)[-1]

for (snp in snp_names) {

  df <- top_SNP_alleles %>% select(FID, all_of(snp))
  
  # Subset by taxa
  subset_df <- df %>% filter(FID %in% common$V1)
  colnames(subset_df)[1:2] <- c("Taxa", "Allele")
  
  
  if (nrow(subset_df) > 0) {
    subset_control_SNP_alleles[[snp]] <- subset_df
  }
}


## getting the gene name associated with each SNP

merged_data_list <- list()

for (snp in names(subset_control_SNP_alleles)) {
  df <- subset_control_SNP_alleles[[snp]]
  
  # Get corresponding gene name
  gene <- unique_snps$gene[unique_snps$SNP == snp]
  

  gene_name <- if (length(gene) == 1) gene else snp
  
 
  df$SNP <- snp
  df$gene <- gene_name
  
 
  merged_data_list[[snp]] <- df
}

merged_data_control <- do.call(rbind, merged_data_list)


merged_data_list <- list()

for (snp in names(subset_control_SNP_alleles)) {
  df <- subset_control_SNP_alleles[[snp]]
  
  # Get all genes linked to this SNP
  genes <- unique_snps$gene[unique_snps$top_SNP == snp]
  
  # For each gene, create a copy of df with that gene name
  for (gene_name in genes) {
    df_copy <- df
    df_copy$SNP <- snp
    df_copy$gene <- gene_name
    merged_data_list[[paste0(snp, "_", gene_name)]] <- df_copy
  }
}
merged_data_control <- do.call(rbind, merged_data_list)

n_distinct(merged_data_control$gene)

##now saving this information separately for each genotype since we have sv files per genotype
B73_geneid =fread("~/Downloads/B73_geneModels_v5_v2.csv", data.table = F, header = T)
B73_geneid = B73_geneid[, (1:4)]


B73_geneid$extended_start = B73_geneid$start - 2000
str(B73_geneid)
B73_geneid$start = as.numeric(B73_geneid$start)


positions = merge(merged_data, B73_geneid, by.x = "gene", by.y = "ID")
n_distinct(merged_data$top_SNP)
setdiff(merged_data$gene, positions$gene)

control_positions = merge(merged_data_control, B73_geneid, by.x = "gene", by.y = "ID")

unique_genotypes <- unique(control_positions[["Taxa"]])

# Loop through each genotype and save a file
for (gen in unique_genotypes) {
  df_gen <- control_positions[control_positions[["Taxa"]] == gen, ]
  
  # Clean genotype name for file safety
  safe_gen <- gsub("[^A-Za-z0-9_\\-]", "_", gen)
  
  filename <- paste0("subset_", safe_gen, ".txt")
  write.table(df_gen, filename, sep = "\t", quote = FALSE, row.names = FALSE)
}

##now we are perfoming overlap between the subset files for each genotype and the SV files while keeping only deletion SVs


library(GenomicRanges)
library(data.table)
library(dplyr)

subset_files <- list.files("~/tre1_eqtl/alleqtls_analysis/testset_07.30", pattern = "^subset_.*\\.txt$", full.names = TRUE)


for (subset_file in subset_files) {
  
  # Extract genotype name
  genotype_name <- gsub("subset_|\\.txt", "", basename(subset_file))
  
  # Define corresponding SV file path
  sv_file <- paste0("~/SVs/Zm-", genotype_name, "-REFERENCE-NAM-1.0_SV_knobs_centromeres_vs_B73_coordinates.bed")
  
  # Skip if SV file doesn't exist
  if (!file.exists(sv_file)) {
    message("SV file not found for genotype: ", genotype_name)
    next
  }
  
  # Load subset and SV data
  genotype_data <- fread(subset_file, data.table = FALSE)
  sv_data <- fread(sv_file, data.table = FALSE, header = FALSE)
  colnames(sv_data) <- c("chr", "start", "end", "length", "feature", "method")
  
  # Keep only DEL-type SVs
  sv_data$feature_type <- sub("feature_type=", "", sv_data$feature)
  sv_data <- sv_data[sv_data$feature_type == "DEL", ]
  if (nrow(sv_data) == 0) {
    message("No DEL features found for genotype: ", genotype_name)
    next
  }
  
  # Clean up chromosome names and adjust coordinates
  sv_data$chr <- gsub("^chr", "", sv_data$chr)
  sv_data$start <- sv_data$start + 1
  
  # Create GRanges for genes with required metadata
  gr_genes <- GRanges(
    seqnames = genotype_data$chr,
    ranges = IRanges(start = genotype_data$extended_start,
                     end = genotype_data$end)
  )
  mcols(gr_genes)$gene <- genotype_data$gene
  mcols(gr_genes)$Taxa <- genotype_data$Taxa
  mcols(gr_genes)$Allele <- genotype_data$Allele
  mcols(gr_genes)$SNP <- genotype_data$SNP
  
  # Create GRanges for SVs
  gr_svs <- GRanges(
    seqnames = sv_data$chr,
    ranges = IRanges(start = sv_data$start, end = sv_data$end),
    sv_feature = sv_data$feature_type,
    method = sv_data$method
  )
  
  # Find overlaps
  hits <- findOverlaps(gr_genes, gr_svs)
  if (length(hits) == 0) {
    message("No overlaps found for genotype: ", genotype_name)
    next
  }
  
  # Construct matched output
  matched_df <- data.frame(
    gene = mcols(gr_genes)$gene[queryHits(hits)],
    Taxa = mcols(gr_genes)$Taxa[queryHits(hits)],
    Allele = mcols(gr_genes)$Allele[queryHits(hits)],
    SNP = mcols(gr_genes)$SNP[queryHits(hits)],
    chr = as.character(seqnames(gr_genes)[queryHits(hits)]),
    extended_start = start(gr_genes)[queryHits(hits)],
    gene_end = end(gr_genes)[queryHits(hits)],
    sv_start = start(gr_svs)[subjectHits(hits)],
    sv_end = end(gr_svs)[subjectHits(hits)],
    sv_feature = mcols(gr_svs)$sv_feature[subjectHits(hits)],
    sv_method = mcols(gr_svs)$method[subjectHits(hits)]
  )
  
  # Save output
  output_file <- paste0("~/tre1_eqtl/alleqtls_analysis/testset_07.30/SV_overlap_test_del/", genotype_name, "_delsvs.txt")
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  write.table(matched_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  message("Processed: ", genotype_name)
}

## have provided all the above results in the folder, rerun for confirmation

for (subset_file in subset_files) {
  
  # Extract genotype name
  genotype_name <- gsub("subset_|\\.txt", "", basename(subset_file))
  
  # Define corresponding SV file path
  sv_file <- paste0("~/SVs/Zm-", genotype_name, "-REFERENCE-NAM-1.0_SV_knobs_centromeres_vs_B73_coordinates.bed")
  
  # Skip if SV file doesn't exist
  if (!file.exists(sv_file)) {
    message("SV file not found for genotype: ", genotype_name)
    next
  }
  
  # Load subset and SV data
  genotype_data <- fread(subset_file, data.table = FALSE)
  sv_data <- fread(sv_file, data.table = FALSE, header = FALSE)
  colnames(sv_data) <- c("chr", "start", "end", "length", "feature", "method")
  
  # Clean up chromosome names and adjust coordinates
  sv_data$chr <- gsub("^chr", "", sv_data$chr)
  sv_data$start <- sv_data$start + 1
  
  # Extract feature type from attribute column
  sv_data$feature_type <- sub("feature_type=", "", sv_data$feature)
  
  # Create GRanges for genes with required metadata
  gr_genes <- GRanges(
    seqnames = genotype_data$chr,
    ranges = IRanges(start = genotype_data$extended_start,
                     end = genotype_data$end)
  )
  mcols(gr_genes)$gene <- genotype_data$gene
  mcols(gr_genes)$Taxa <- genotype_data$Taxa
  mcols(gr_genes)$Allele <- genotype_data$Allele
  mcols(gr_genes)$SNP <- genotype_data$SNP
  
  # Create GRanges for SVs (all types)
  gr_svs <- GRanges(
    seqnames = sv_data$chr,
    ranges = IRanges(start = sv_data$start, end = sv_data$end),
    sv_feature = sv_data$feature_type,
    method = sv_data$method
  )
  
  # Find overlaps
  hits <- findOverlaps(gr_genes, gr_svs)
  if (length(hits) == 0) {
    message("No overlaps found for genotype: ", genotype_name)
    next
  }
  
  # Construct matched output
  matched_df <- data.frame(
    gene = mcols(gr_genes)$gene[queryHits(hits)],
    Taxa = mcols(gr_genes)$Taxa[queryHits(hits)],
    Allele = mcols(gr_genes)$Allele[queryHits(hits)],
    SNP = mcols(gr_genes)$SNP[queryHits(hits)],
    chr = as.character(seqnames(gr_genes)[queryHits(hits)]),
    extended_start = start(gr_genes)[queryHits(hits)],
    gene_end = end(gr_genes)[queryHits(hits)],
    sv_start = start(gr_svs)[subjectHits(hits)],
    sv_end = end(gr_svs)[subjectHits(hits)],
    sv_feature = mcols(gr_svs)$sv_feature[subjectHits(hits)],
    sv_method = mcols(gr_svs)$method[subjectHits(hits)]
  )
  
  # Save output
  output_file <- paste0("~/tre1_eqtl/alleqtls_analysis/testset_07.30/SV_overlap_test_all/", genotype_name, "_allsvs.txt")
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  write.table(matched_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  message("Processed: ", genotype_name)
}




























