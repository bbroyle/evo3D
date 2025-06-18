# evo3D
<img src="man/figures/evo3d_hex_b.png" width="200"/>

**evopatchr** is an R package for structure-aware population genetics, enabling patch-level evolutionary analysis of protein surfaces. It integrates selection metrics with 3D structural data to identify spatially clustered signals of diversity and selection.

Designed with applications in immunology, virology, and comparative genomics, `evo3D` supports multi-chain proteins, multi-model structures, and antibody–antigen complexes.

---

## Key Features

- Computes selection metrics (π, haplotype diversity, Tajima’s D) at surface-defined patches
- Detects and integrates antibody epitopes from PDB structures
- Supports multi-chain and multi-model PDBs for robust structural inference
- Maps structure-defined residues to codon-aligned MSA windows
- Generates B-factor-encoded PDBs and per-residue selection tables
- Fully self-contained: no reliance on DSSP, MAFFT, or external alignment tools

---

## Installation

To install from GitHub:

```r
# install.packages("devtools")
devtools::install_github("bbroyle/evopatchr")

# msa and alignment functionalities require 'msa' from Bioconducter
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("msa")
```

## Quick Example

```r
library(evopatchr)

msa_path <- system.file("extdata", "rh5_pfalc.fasta", package = "evopatchr")
pdb_path <- system.file("extdata", "rh5_4wat.pdb", package = "evopatchr")

# run_patchr is designed for single msa and single pdb runs #
result <- run_patchr_single(msa_path = msa_path,
                             pdb_path = pdb_path,
                             chain = 'A')

write_stat_to_bfactor(result$selection_df,
                      result$pdb_info$pdb,
                      stat_name = "hap",
                      outfile = "rh5_hap_div.pdb")
```

## Running step-wise (more control)
### more examples of step-wise runs at end of README  

```r
library(evopatchr)

# single MSA and PDB (or mmCIF) // essentially run_patchr_single() ----

msa_path <- system.file("extdata", "rh5_pfalc.fasta", package = "evopatchr")
pdb_path <- system.file("extdata", "rh5_4wat.pdb", package = "evopatchr")

# read in msa #
msa_info <- WRAPPER_msa_to_ref(msa_path = msa_path)

# read in pdb #
pdb_info <- WRAPPER_pdb_to_patch(pdb_path = pdb_path,
                                chain = c('A'))
  
# generate alignment between msa and pdb / and create msa_subsets #
aln_info <- WRAPPER_align_msa_pdb(msa_info = msa_info,
                                 pdb_info = pdb_info, 
                                 chain = 'A', coverage_plot = T)
  
# calculate selection #
selection_df <- run_pegas_three(aln_info$msa_subsets, pdb_info$residue_df)

```



## License

This package is released under the MIT License.  
Portions of the structural accessibility logic are adapted from the [DSSP project](https://github.com/PDB-REDO/dssp),  
licensed under the BSD 2-Clause License. See `inst/LICENSE.note` for details.

## Conatct

Brad Broyles
PhD Candidate, Computational and Structural Biology,
Purdue University
bbroyle@purdue.edu

## Running multi-chain (heterodimer)

```r
library(evopatchr)

msa_path1 <- system.file("extdata", "e1_hepc.aln", package = "evopatchr")
msa_path2 <- system.file("extdata", "e2_hepc.aln", package = "evopatchr")
pdb_path <- system.file("extdata", "e1e2_8fsj.pdb", package = "evopatchr")

# read in two MSA #
msa1 <- WRAPPER_msa_to_ref(msa_path = msa_path1)
msa2 <- WRAPPER_msa_to_ref(msa_path = msa_path2)

# read in pdb #
pdb_info <- WRAPPER_pdb_to_patch(pdb_path = pdb_path,
                                chain = c('A', 'E'))

# each chain gets a alignment to msa #
aln_info1 <- WRAPPER_align_msa_pdb(msa = msa1, 
                                  pdb_info = pdb_info, 
                                  chain = 'A')

aln_info2 <- WRAPPER_align_msa_pdb(msa = msa2,
                                  pdb_info = pdb_info, 
                                  chain = 'E')

# build cross chain msa subsets #
msas <- extend_msa(aln_info1$msa_subsets, 
                  aln_info2$msa_subsets)

# calculate selection #
selection_df <- run_pegas_three(msas, pdb_info$residue_df)

```

## Running multi-chain (homotrimer)

```r
library(evopatchr)

msa_path <- system.file("extdata", "spike_sarscov2.fa", package = "evopatchr")
pdb_path <- system.file("extdata", "spike_7fb0.pdb", package = "evopatchr")

# read in one msa #
msa1 <- WRAPPER_msa_to_ref(msa_path = msa_path)

# read in pdb ** note distance method set to 'ca' ~ faster (but not required) ** #
pdb_info <- WRAPPER_pdb_to_patch(pdb_path = pdb_path,
                                chain = c('A', 'B', 'C'),
                                distance_method = 'ca')

# each chain gets a alignment to msa #
aln_info1 <- WRAPPER_align_msa_pdb(msa_info = msa1, 
                                  pdb_info = pdb_info, 
                                  chain = 'A')

aln_info2 <- WRAPPER_align_msa_pdb(msa_info = msa1,
                                  pdb_info = pdb_info, 
                                  chain = 'B')

aln_info3 <- WRAPPER_align_msa_pdb(msa_info = msa1,
                                  pdb_info = pdb_info, 
                                  chain = 'C')

# build cross chain msa subsets #
msas <- extend_msa(aln_info1$msa_subsets, 
                  aln_info2$msa_subsets)

msas <- extend_msa(msas,
                  aln_info3$msa_subsets)
  
# calculate selection #
selection_df <- run_pegas_three(msas, pdb_info$residue_df)

```

## Running multi model (complementary resolved regions in PDB)

```r
library(evopatchr)

msa_path <- system.file("extdata", "spike_sarscov2.fa", package = "evopatchr")
pdb_path1 <- system.file("extdata", "spike_7fb0.pdb", package = "evopatchr")
pdb_path2 <- system.file("extdata", "spike_7fb1.cif", package = "evopatchr")

# read in MSA #
msa_info <- WRAPPER_msa_to_ref(msa_path = msa_path)

# read in pdb #
pdb_info1 <- WRAPPER_pdb_to_patch(pdb_path = pdb_path1,
                                chain = c('A'))

pdb_info2 <- WRAPPER_pdb_to_patch(pdb_path = pdb_path2,
                                 chain = c('A'))

# each chain gets a alignment to msa #
aln_info1 <- WRAPPER_align_msa_pdb(msa_info = msa_info,
                                 pdb_info = pdb_info1, 
                                 chain = 'A')

aln_info2 <- WRAPPER_align_msa_pdb(msa_info = msa_info,
                                  pdb_info = pdb_info2, 
                                  chain = 'A')

# extend capture windows based on complimentary structures #
merged_nuc_windows <- extend_nuc_windows(aln_info1$nuc_patches, 
                                        aln_info2$nuc_patches)

# generate msa subsets using combined window information #
msas <- extract_msa_subsets(msa_info$msa_mat, merged_nuc_windows)
 
# calculate selection ~ dataframe is based on codon position in msa #
selection_df <- run_pegas_three(msas)

```

## Running selection (including epitope) ** NEEDS WORK **

```r
library(evopatchr)

msa_path <- system.file("extdata", "rh5_pfalc.fasta", package = "evopatchr")
pdb_path <- system.file("extdata", "rh5_6rcu.pdb", package = "evopatchr")

# read in MSA #
msa_info <- WRAPPER_msa_to_ref(msa_path = msa_path)

# read in pdb #
pdb_info <- WRAPPER_pdb_to_patch(pdb_path = pdb_path, chain = 'A')

# grab epitope information #
pdb <- .standardize_pdb_input(pdb_path = pdb_path, chain = 'all')
epitope_info1 <- identify_epitopes(pdb, ag_chain = 'A', h_chain = 'B', l_chain = 'C')
epitope_info2 <- identify_epitopes(pdb, ag_chain = 'A', h_chain = 'D', l_chain = 'E')

# add to residue_df
pdb_info$residue_df[nrow(pdb_info$residue_df) + 1,] <- NA
pdb_info$residue_df[nrow(pdb_info$residue_df), 8] <- 'epi:BC>A'
pdb_info$residue_df[nrow(pdb_info$residue_df), 9] <- epitope_info1$epitope

pdb_info$residue_df[nrow(pdb_info$residue_df) + 1,] <- NA
pdb_info$residue_df[nrow(pdb_info$residue_df), 8] <- 'epi:DE>A'
pdb_info$residue_df[nrow(pdb_info$residue_df), 9] <- epitope_info2$epitope

# generate alignment between msa and pdb / and create msa_subsets #
aln_info <- WRAPPER_align_msa_pdb(msa_info = msa_info,
                                 pdb_info = pdb_info, 
                                 chain = 'A', coverage_plot = T)

# calculate selection #
selection_df <- run_pegas_three(aln_info$msa_subsets, pdb_info$residue_df)

```
