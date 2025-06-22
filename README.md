# evo3D
<img src="man/figures/evo3d_hex_b.png" width="200"/>

An R package for structure-aware population genetics, enabling selection and diversity analysis over sliding "patch" windows of protein surfaces (or structures in general). A wide range of analysis can be completed through the simple wrapper run_evo3d(), with the underlying modules ( msa_to_ref(), pdb_to_patch(), and aln_msa_to_pdb() ) also available for tailored use.

---

## Key Features

- Runs fully within R -- no need to link external dependencies
- Maps structure-defined windows (including across mutli-chain complexes and protein-protein interfaces) to MSAs for patch-level MSA subsetting
- Computes selection and diversity statistics (nucleotide diversity, haplotype diversity, Tajima’s D, etc.) or exports patch-level MSA subsets for downstream analysis
- Includes tunable parameters for patch size, surface filtering, and more
- Transparent workflow exposing internal MSA-to-PDB alignments for user inspection
- Supports multiple structural models to improve MSA coverage
- Optionally writes computed statistics to the B-factor column for PDB visualization

---

## Installation

To install from GitHub:

```r
# 1. Install evo3D (and devtools if you don’t have it)

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("bbroyle/evo3D")

# 2. Install msa (necessary for msa to pdb alignments)
# this package can take a while to install,
# we may switch in future to DECIPHER or BioStrings

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("msa", update = FALSE)

```

## Quick Example

```r
library(evo3D)

# getting file paths -- you can replace with your own data #
msa_path = system.file("extdata", "rh5_pfalc.fasta", package = "evo3D")
pdb_path = system.file("extdata", "rh5_6mpv_AB.pdb", package = "evo3D")

# run_evo3D is designed for single analysis runs # --- for batch runs run_evo3D_batch will be provided soon #
# chain = 'auto' by default or set to 'B' for this example #
results = run_evo3d(msa_path, pdb_path, chain = 'B') 

write_stat_to_bfactor(results, stat_name = "hap", outfile = "rh5_hap_div.pdb")

# !!! scores are not selection at that amino acid -- but at the patch centered on that amino acid !!!
# also some residues don't have a score for this statistic - they have been given arbitrarily low score of -10 #
# to color structure try these commands in pymol #
# select stat, b > -10
# spectrum b, selection=stat

```

Let's quikly cover run_evo3d() results. Results are in a structured list -- with the following entries:<br/>
$evo3d_df -- dataframe holding msa to pdb alignment information, 3D codon patch information, and calculated statistics<br/>
$final_msa_subsets -- list of msa subsets named on codon at center of patch or interface id<br/>
$msa_info_sets -- outputs of module 1 msa_to_ref()<br/>
$pdb_info_sets -- outputs of module 2 pdb_to_patch()<br/>
$aln_info_sets -- outputs of module 3 aln_msa_to_pdb()<br/>
$call_info -- meta data of the analysis run inlcuding 3D sliding window paramters and msa and pdb file paths

## Tunable patch parameters
Here we will set patch size to 10 Angstrom radius (defulat 15 Angstrom) and include all residues (defualt is residues with relative solvent accessibility of at least 0.1).
Many other parameters can be adjusted through msa_controls, pdb_controls, and stat_controls to run_evo3d(). 
Additionally we will save our msa subsets and dataframe to file with argumnets write_patch_fastas and write_evo3d_df.
Let's also turn off statistics calculation sense we are saving msa subsets for downstream analysis

```r
library(evo3D)

# load the example data again
msa_path = system.file("extdata", "rh5_pfalc.fasta", package = "evo3D")
pdb_path = system.file("extdata", "rh5_6mpv_AB.pdb", package = "evo3D")

# lets look at default parameters #
show_evo3d_defaults()

# we want to change patch.dist.cutoff and patch.rsa.cutoff in pdb_controls #
results2 = run_evo3d(msa_path, pdb_path, pdb_controls = list(patch.dist.cutoff = 10, patch.rsa.cutoff = 0),
                     write_patch_fastas = TRUE, write_evo3d_df = TRUE, output_dir = 'rh5_10ang_0rsa',
                     run_selection = FALSE)

```

## Running on protein complexes (multiple MSAs), multiple proteins (complementary MSA coverage)

```r

# protein complexes with multiple MSAs
msa1 = system.file("extdata", "e1_hepc_sorted.aln", package = "evo3D")
msa2 = system.file("extdata", "e2_hepc_sorted.aln", package = "evo3D")
pdb = system.file("extdata", "e1e2_8fsj.pdb", package = "evo3D")
results1 = run_evo3d(list(msa1, msa2), pdb)

# multiple pdbs that increase MSA coverage
msa = system.file("extdata", "rh5_pfalc.fasta", package = "evo3D")
pdb1 = system.file("extdata", "rh5_6mpv_AB.pdb", package = "evo3D")
pdb2 = system.file("extdata", "rh5_8q5d_A.pdb", package = "evo3D")
results2 = run_evo3d(msa, list(pdb1, pdb2))

# homomultimers (overload chain argument)
# data not inlcuded in package #
results3 = run_evo3d(msa, pdb, chain = 'ABC')

```

## Running step-wise evo3D modules

```r
library(evo3D)

# same example files as above #
msa_path = system.file("extdata", "rh5_pfalc.fasta", package = "evo3D")
pdb_path = system.file("extdata", "rh5_6mpv_AB.pdb", package = "evo3D")

# read in msa, get a reference sequence, and translate that sequence #
msa_info = msa_to_ref(msa_path)

# Some helpful functions for detecting chains of interest #
# .auto_detect_chain shows % of kmers for pdb chains matching peptide reference sequence #
evo3D:::.auto_detect_chain(msa_info$pep, pdb_path)
# -or- #
# .plot_chain_map() gives quick viewing access without leaving R #
evo3D:::.plot_chain_map(pdb_path)

# read in pdb, calculate distance matrix, calculate solvent accessibiliity, calculate patches #
pdb_info = pdb_to_patch(pdb_path, chain = 'B')
  
# generate alignment between msa and pdb and create msa_subsets #
aln_info = aln_msa_to_pdb(msa_info, pdb_info, chain = 'B')
  
# calculate selection (3 stats from pegas package - tajima's D, nucleotide diversity, haplotype diversity #
selection_df = run_pegas_three(aln_info$msa_subsets, aln_info$aln_df)

```

## License

This package is released under the MIT License.  
The structural solvent accessibility logic is adapted from the [DSSP project](https://github.com/PDB-REDO/dssp),  
licensed under the BSD 2-Clause License. See `inst/LICENSE.note` for details.

## Conatct

Brad Broyles  
PhD Candidate, Computational and Structural Biology  
Purdue University  
bbroyle@purdue.edu

