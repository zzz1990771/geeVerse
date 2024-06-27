#' A Simulated Genetic Data from HapGen2
#'
#' The `simuGene` dataset contains 500 SNPs simulated  data from
#' a commonly used tool for genetic data, HapGen2. We re-sampled existing genotype
#' data to create this simulated data. The genotype data we resample from is the
#' publicly available 1000 Genomes Project data. More specifically, we use resampled
#' from chromosome 14.
#'
#' @format A data frame with 1000 rows (subjects) and 500 columns (SNPs).
#'
#' @examples
#' data(simuGene)
#' head(simuGene)
"simuGene"
