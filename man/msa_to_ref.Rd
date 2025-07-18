% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MSA_MODULE.R
\name{msa_to_ref}
\alias{msa_to_ref}
\title{Extract Reference and Peptide Sequence from MSA}
\usage{
msa_to_ref(
  msa,
  ref_method = 1,
  force_seqtype = NULL,
  genetic_code = 1,
  verbose = 0
)
}
\arguments{
\item{msa}{A character file path to a FASTA file, a matrix, or an object returned by \code{bio3d::read.fasta()}.}

\item{ref_method}{Reference extraction method: one of \code{"most_complete"}, \code{"consensus"}, or a numeric row index.}

\item{force_seqtype}{Optional sequence type: \code{"protein"}, \code{"nucleotide"}, or \code{NULL} to auto-detect.}
}
\value{
A list with the following elements:
\itemize{
  \item \code{msa_mat}: The standardized alignment matrix.
  \item \code{ref}: The reference sequence (DNA or protein).
  \item \code{pep}: The translated peptide sequence (if nucleotide input).
  \item \code{seq_type}: The detected or specified sequence type.
}
}
\description{
Wrapper function that standardizes input, extracts a reference sequence from an MSA, detects or applies a sequence type, and translates the reference to peptide if needed.
}
