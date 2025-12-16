# evo3D

**Structure-informed evolutionary analysis in R**

evo3D is a newly released package.  
If you encounter issues, have suggestions, or would like to request new features, please open an issue - community feedback is appreciated.

<img src="man/figures/evo3d_hex_b.png" width="200"/>

An R package for structure-informed evolutionary analysis that links multiple-sequence alignments to 3D protein structures and returns spatially defined sequence windows (“patches”) for downstream population-genetic and evolutionary statistics. evo3D supports surface-, interface-, and multimer-aware analyses and is designed for flexible, end-to-end workflows via a single wrapper, with modular functions exposed for advanced use.

---
## Key Features

- External-dependency-free R package with cross-platform support
- Links nucleotide or protein MSAs to 3D structures and extracts spatially defined sequence windows (“patches”)
- Patch definitions support surface exposure, protein–protein interfaces, multimeric assemblies, and homomultimer merging
- Returns patch-level MSAs (“spatial haplotypes”) for arbitrary downstream evolutionary or population-genetic analyses
- Includes built-in calculations for site- and patch-level statistics (e.g., entropy, nucleotide diversity, haplotype diversity, Tajima’s D)
- Explicit and inspectable MSA–structure alignment, enabling user validation and correction
- Supports integration across multiple structural models and multiple MSAs
- Exports results to tabular formats and maps statistics to PDB B- and Q-factor fields for direct structural visualization
---

## Installation

To install the latest version from GitHub:

```r
# (Optional) Remove any previously installed version
if ("evo3D" %in% rownames(installed.packages())) {
  remove.packages("evo3D")
}

# 1. Install evo3D (and devtools if needed)
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("bbroyle/evo3D")

# 2. Install msa (used for MSA–structure alignment)
# Note: this package can take some time to install
# and may be replaced by DECIPHER, Biostrings, or simple pairwise alignment in the future
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("msa", update = FALSE)
```

## Quick Example

```r
library(evo3D)

# Example data shipped with evo3D (replace with your own)
msa_path <- system.file("extdata", "rh5_pfalc.fasta", package = "evo3D")
pdb_path <- system.file("extdata", "rh5_6mpv_AB.pdb", package = "evo3D")

# run_evo3d() is designed for single analysis runs
# chain = "auto" (default) handles MSA–PDB chain mapping
results <- run_evo3d(msa_path, pdb_path, detail_level = 2)

# Codon-indexed windows corresponding to structure-informed patches
res_df <- results$evo3d_df

# With detail_level = 2, patch-level MSA subsets ("spatial haplotypes")
# are returned directly in R and can be used for any downstream analysis
msa_subs <- results$final_msa_subsets

# Optionally write patch-level MSAs to FASTA
# Files are named by their anchoring (central) codon
write_patch_fastas(msa_subs, output_dir = "testing")
```

Let's quikly cover run_evo3d() results. Results are in a structured list -- with the following entries:<br/>
$evo3d_df -- dataframe holding msa to pdb alignment information, 3D codon patch information, and any calculated statistics<br/>
$final_msa_subsets -- list of msa subsets named by their ancoring codon <br/>
$msa_info_sets -- outputs of module 1 msa_to_ref()<br/>
$pdb_info_sets -- outputs of module 2 pdb_to_patch()<br/>
$aln_info_sets -- outputs of module 3 aln_msa_to_pdb()<br/>
$call_info -- meta data of the analysis run inlcuding 3D sliding window paramters and msa and pdb file paths

## Tunable workflow

Patch definitions and statistical analyses can be customized directly through the
`run_evo3d()` wrapper. Below, we analyze the Hepatitis C virus E1/E2 complex,
construct surface-defined patches, and compute Tajima’s D for each patch.

```r
library(evo3D)

# Define HCV inputs
msa_path1 <- system.file("extdata", "e1_hepc_sorted.aln", package = "evo3D")
msa_path2 <- system.file("extdata", "e2_hepc_sorted.aln", package = "evo3D")
pdb_path  <- system.file("extdata", "e1e2_8fsj.pdb",  package = "evo3D")

# Inspect default control parameters
show_evo3d_defaults("stat")
show_evo3d_defaults("pdb")

# Enable Tajima's D calculation
stat_controls <- list(calc_tajima = TRUE)

# Define fixed-size surface patches
# Here: patches of 15 residues, restricted to residues with SASA ≥ 10 Å²
pdb_controls <- list(
  dist_cutoff = NA,
  max_patch   = 15,
  sasa_cutoff = 10
)

# Run evo3D on the E1/E2 complex
results <- run_evo3d(
  list(msa_path1, msa_path2),
  pdb_path,
  stat_controls = stat_controls,
  pdb_controls  = pdb_controls
)

# Inspect patch-level results
res_df <- results$evo3d_df

# Write Tajima’s D to the PDB B-factor field for visualization
write_to_pdb(
  results,
  stat_name = "tajima",
  outfile   = "testing/hcv_tajima.pdb"
)

# In PyMOL (example):
# select stat, b > -99
# spectrum b, selection=stat

```

## License

This package is released under the MIT License.  
The structural solvent accessibility logic is adapted from the
[DSSP project](https://github.com/PDB-REDO/dssp),
which is licensed under the BSD 2-Clause License.
See `inst/LICENSE.note` for details.

## Preprint

A preprint describing the evo3D framework and methodology is available at:  
https://doi.org/10.64898/2025.12.12.694041

## Cite

```bibtex
@article{broyles2025evo3d,
  title   = {evo3D: a generalised framework for structure-informed evolutionary analysis, implemented in R},
  author  = {Broyles, Bradley Keith and He, Qixin},
  year    = {2025},
  doi     = {10.64898/2025.12.12.694041},
  url     = {https://doi.org/10.64898/2025.12.12.694041},
  note    = {Preprint}
}
```

## Contact

Brad Broyles  
PhD Candidate, Computational and Structural Biology  
Purdue University  
bbroyle@purdue.edu

