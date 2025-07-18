% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ALN_MODULE.R
\name{.map_aln_to_positions}
\alias{.map_aln_to_positions}
\title{Map Aligned Amino Acid Positions to Nucleotide Positions}
\usage{
.map_aln_to_positions(aln_mat, residue_df, chain = NA)
}
\arguments{
\item{aln_mat}{A character matrix from \code{.align_sequences()}, with two rows: reference and structure-derived.}

\item{residue_df}{Data frame of patch or interface residues, typically from \code{identify_patches()}.}

\item{chain}{Optional chain ID string (or vector) to restrict mapping to a subset of chains.}
}
\value{
Modified \code{aln_mat} with alignment positions replaced by codon and PDB residue indices.
}
\description{
Overlays alignment columns with residue-level and codon-level position indices for downstream
mapping of structural patches to nucleotide-level coordinates.
}
\details{
Typically used after \code{.align_sequences()} and \code{identify_patches()} to annotate alignment matrix rows.
}
\keyword{internal}
