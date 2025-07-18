% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EXTENSION_UTILS.R
\name{identify_epitopes}
\alias{identify_epitopes}
\title{Identify Antibody–Antigen Contacts}
\usage{
identify_epitopes(
  pdb,
  ag_chain = NULL,
  h_chain = NULL,
  l_chain = NULL,
  dist_cutoff = 5
)
}
\arguments{
\item{pdb}{A \code{bio3d} PDB object.}

\item{ag_chain}{Chain ID of the antigen.}

\item{h_chain}{Chain ID of the antibody heavy chain.}

\item{l_chain}{Chain ID of the antibody light chain.}

\item{dist_cutoff}{Maximum distance (in Å) to define a contact. Default is 5.}
}
\value{
A list with:
\item{epitope}{Residue IDs (e.g., \code{"35_A_+42_A_"}) on the antigen within contact range.}
\item{paratope_h}{Heavy chain residues contacting the antigen.}
\item{paratope_l}{Light chain residues contacting the antigen.}
\item{contacts}{Matrix of all atom-level contacts (row = antigen residue, col = antibody residue).}
}
\description{
Identifies epitope and paratope residues from a PDB structure based on inter-chain atom distances.
}
