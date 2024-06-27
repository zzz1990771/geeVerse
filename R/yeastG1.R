#' A Subset of Yeast Cell Cycle Gene Expression Data (G1 Phase)
#'
#' The `yeastG1` dataset contains gene expression data from the yeast cell
#' cycle during the G1 phase.
#' The original dataset (Spellman et al. 1998) includes expression levels for 6178 genes measured at 18 time points.
#' And this is a subset of 283 cell-cycled-regularized genes observed over 4 time
#' points at G1 stage and the standardized binding probabilities of a total of 96 TFs
#' obtained from
#'
#' @format A data frame with 1132 rows and 99 columns.
#'
#' The dataset contains gene expression levels for the following transcription factors:
#' ABF1, ACE2, ADR1, ARG80, ARG81, ARO80, ASH1, BAS1, CAD1, CBF1, CIN5, CRZ1,
#' CUP9, DAL81, DAL82, DIG1, DOT6, FHL1, FKH1, FKH2, FZF1, GAL4, GAT1, GAT3,
#' GCN4, GCR1, GCR2, GLN3, GRF10.Pho2., GTS1, HAL9, HAP2, HAP3, HAP4, HAP5,
#' HIR1, HIR2, HMS1, HSF1, IME4, INO2, INO4, IXR1, LEU3, MAC1, MAL13, MATa1,
#' MBP1, MCM1, MET31, MET4, MIG1, MOT3, MSN1, MSN4, MSS11, MTH1, NDD1, NRG1,
#' PDR1, PHD1, PHO4, PUT3, RAP1, RCS1, REB1, RFX1, RGM1, RLM1, RME1, ROX1,
#' RPH1, RTG1, RTG3, SFP1, SIG1, SIP4, SKN7, SMP1, SOK2, SRD1, STB1, STE12,
#' STP1, STP2, SUM1, SWI4, SWI5, SWI6, YAP1, YAP5, YAP6, YFL044C, YJL206C,
#' ZAP1, ZMS1
#'
#' @source Spellman, P. T., Sherlock, G., Zhang, M. Q., Iyer, V. R., Anders, K., Eisen, M. B., ... & Futcher, B. (1998).
#'   Comprehensive identification of cell cycle-regulated genes of the yeast Saccharomyces cerevisiae by microarray hybridization.
#'   Molecular biology of the cell, 9(12), 3273-3297.
#' @source Wang, L., Zhou, J., and Qu, A. (2012). Penalized generalized estimating equations for high-dimensional longitudinal data anaysis. \emph{Biometrics}, \bold{68}, 353--360.
#'
#' @examples
#' data(yeastG1)
#' head(yeastG1)
"yeastG1"
