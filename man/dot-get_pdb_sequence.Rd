% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PDB_MODULE.R
\name{.get_pdb_sequence}
\alias{.get_pdb_sequence}
\title{Extract Sequence from PDB}
\usage{
.get_pdb_sequence(pdb, chain = NA, in_module = F)
}
\arguments{
\item{pdb}{A PDB object from \code{bio3d}.}

\item{chain}{Optional chain ID(s).}

\item{pdb_path}{Path to a PDB file.}
}
\value{
A named character vector of sequences, one per chain.
}
\description{
Retrieves the amino acid sequence from a PDB file or object.
}
