library(genio)
SNPs = read_plink("dataset_ld_pruned")
all_mat = t(SNPs$X)
snp_matrix = all_mat

# Function to calculate MAF for a single SNP column
calculate_maf <- function(snp_column) {
  # Count alleles: homozygous reference (1) counts as 2 reference alleles,
  # heterozygous (2) counts as 1 reference and 1 alternative allele,
  # homozygous alternative (3) counts as 2 alternative alleles.
  ref_allele_count <- sum(snp_column == 0) * 2 + sum(snp_column == 1)
  alt_allele_count <- sum(snp_column == 2) * 2 + sum(snp_column == 1)

  # Total alleles
  total_alleles <- ref_allele_count + alt_allele_count

  # Calculate frequencies
  ref_freq <- ref_allele_count / total_alleles
  alt_freq <- alt_allele_count / total_alleles

  # MAF is the minimum of the two frequencies
  maf <- min(ref_freq, alt_freq)
  return(maf)
}

# Apply the MAF calculation across all SNPs (columns) and get a MAF vector
maf_vector <- apply(snp_matrix, 2, calculate_maf)


# Filter SNPs with MAF > 0.2
snp_matrix_filtered <- snp_matrix[, maf_vector > 0.2]

# Remove dup cols
snp_matrix_filtered_uni <- t(unique(t(snp_matrix_filtered)))

#keep the first 2000 SNPs for genetic example for geeVerse
snp_data <- snp_matrix_filtered_uni[,sample(1:NCOL(snp_matrix_filtered_uni),2000)]
saveRDS(snp_data,"chr22_hapgen2_data1.RDS")
