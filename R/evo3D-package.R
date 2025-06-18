#' evo3D: Patch-Level Selection Analysis Using Structure-Informed MSA Subsets
#'
#' Defines sliding windows in 3D space based on protein structure, extracts
#' corresponding subsets from a codon-aligned MSA, and computes selection statistics
#' such as nucleotide diversity, Tajima's D, and haplotype diversity.
#'
#' @docType package
#' @keywords internal
#' @importFrom Rcpp sourceCpp
#' @useDynLib evo3D, .registration = TRUE
"_PACKAGE"