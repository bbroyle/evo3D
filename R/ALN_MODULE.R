# ------------------------------------------ #
# ALN_MODULE.R
# utilities and module for aligning msa to pdb
# module (wrapper) is at the end of this script
# Brad Broyles
# ------------------------------------------ #

# .calculate_coverage() ----

#' Calculate Aligned Coverage Ranges
#'
#' Identifies contiguous non-gap regions in each sequence of an alignment matrix and formats them as range strings.
#' Optionally includes ranges of mismatched positions between aligned sequences.
#'
#' Called by \code{.align_sequences()}.
#'
#' @param aln_mat A character matrix representing a sequence alignment (rows = sequences, columns = positions).
#' @param mismatch Optional numeric vector of mismatch positions to summarize (e.g., from sequence comparison).
#'
#' @return A named list. Each element is a character vector of ranges for a sequence (e.g., \code{"5:25"}), plus an optional \code{mismatch} element.
#' @keyword internal
.calculate_coverage = function(aln_mat, mismatch) {
  covered_regions = list()

  # convert alignments into ranges
  for (i in 1:nrow(aln_mat)) {
    seq_name = rownames(aln_mat)[i]
    # ignore gap positions
    positions = which(aln_mat[i, ] != "-")

    if (length(positions) > 0) {
      # find breaks in sequence
      breaks = c(0, which(diff(positions) > 1), length(positions))
      ranges = character()

      for (j in 1:(length(breaks) - 1)) {
        start_idx = breaks[j] + 1
        end_idx = breaks[j + 1]
        start_pos = positions[start_idx]
        end_pos = positions[end_idx]
        ranges = c(ranges, paste0(start_pos, ":", end_pos))
      }

      covered_regions[[seq_name]] = ranges
    } else {
      covered_regions[[seq_name]] = character(0)
    }
  }

  # convert mismatches into ranges
  if (length(mismatch) > 0) {
    # find breaks in sequence
    breaks = c(0, which(diff(mismatch) > 1), length(mismatch))
    ranges = character()

    for (j in 1:(length(breaks) - 1)) {
      start_idx = breaks[j] + 1
      end_idx = breaks[j + 1]
      start_pos = mismatch[start_idx]
      end_pos = mismatch[end_idx]
      ranges = c(ranges, paste0(start_pos, ":", end_pos))
    }

    covered_regions$mismatch = ranges
  } else {
    covered_regions$mismatch = NA
  }


  return(covered_regions)
}

# .plot_coverage() ----

# ** remove this function for now ** #

# #' Plot Sequence Alignment Coverage
# #'
# #' Visualizes aligned coverage regions as horizontal bars for each sequence.
# #' Currently only supports range-style input (e.g., \code{"5:25"}), with placeholder handling for mismatches.
# #'
# #' @param coverage A named list of ranges (e.g., from \code{.calculate_coverage()}), with each element representing a sequence.
# #'
# #' @return A \code{ggplot2} plot object for visual inspection of coverage.
# #' @examples
# #' \dontrun{
# #' coverage <- .calculate_coverage(aln_matrix)
# #' plot <- .plot_coverage(coverage)
# #' print(plot)
# #' }
# #' @keyword internal
# .plot_coverage = function(coverage) {
#   # convert coverage to data frame for plotting
#   plot_data <- data.frame(
#     sequence = character(),
#     start = numeric(),
#     end = numeric(),
#     stringsAsFactors = FALSE
#   )
#
#   # see if mismatch has ranges or if it is NA #
#   if(is.na(coverage$mismatch[1])){
#     coverage$mismatch = '0:0'
#   }
#
#   for (seq_name in names(coverage)) {
#     ranges <- coverage[[seq_name]]
#     for (range in ranges) {
#       range_parts <- as.numeric(strsplit(range, ":")[[1]])
#       plot_data <- rbind(plot_data, data.frame(
#         sequence = seq_name,
#         start = range_parts[1],
#         end = range_parts[2],
#         stringsAsFactors = FALSE
#       ))
#     }
#   }
#
#   # shrink name if too long (max 20 char then ...) #
#   plot_data$sequence = apply(plot_data, 1, function(x) {
#     if (nchar(x[1]) > 20) {
#       return(paste0(substr(x[1], 1, 20), " ..."))
#     } else {
#       return(x[1])
#     }
#   })
#
#   plot_data$sequence = factor(plot_data$sequence, levels = unique(plot_data$sequence))
#
#   # *** STILL NEED TO INCORPORATE MISMATCH POINTS *** #
#   plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = start, xend = end, y = sequence, yend = sequence)) +
#     ggplot2::geom_segment(linewidth = 3) +
#     ggplot2::theme_bw() +
#     ggplot2::xlab("Alignment Position") +
#     ggplot2::ylab("Sequence")
#
#   return(plot)
#
# }


# .align_sequences() ----

#' Align Reference and PDB Sequences
#'
#' Aligns a pair of amino acid sequences using ClustalOmega (via the \pkg{msa} package),
#' computes positional mismatches, and summarizes coverage across the alignment.
#'
#' Used internally to align MSA or PDB-derived sequences to a reference
#'
#' @param sequences A named character vector of two protein sequences: the reference first, the structure-derived second.
#' @param user_supplied_alignment Placeholder for pre-aligned input (currently unused).
#'
#' @return A list with elements:
#' \describe{
#'   \item{aln_mat}{Character matrix of aligned sequences.}
#'   \item{coverage}{Named list of coverage ranges for each sequence.}
#' }
#' @export
.align_sequences = function(sequences, user_supplied_alignment = NULL){

  # user_supplied_alignment not supported yet -- soon #

  # need to track original gap (from msa)
  # then introduced gaps (from alignment)
  # also track matching positions
  if (!requireNamespace("msa", quietly = TRUE)) {
    stop(
      paste(
        "The 'msa' package is required for alignment.",
        "",
        "To install, run:",
        "  if (!requireNamespace(\"BiocManager\", quietly = TRUE))",
        "    install.packages(\"BiocManager\")",
        "  BiocManager::install(\"msa\")",
        sep = "\n"
      ),
      call. = FALSE
    )
  }

  # use msa::msa() ~ with baked in clustal omega (defualt to GONNET sub matrix)
    # could make invisible()
  aln = suppressMessages(msa::msa(sequences, method = 'ClustalOmega', type = 'protein', order = 'input'))

  aln_set = as.character(aln@unmasked)
  aln_chars = strsplit(aln_set, '')

  # which position are identical
  mismatch = which(aln_chars[[1]] != aln_chars[[2]])

  # only keep mismatch if they are aa to aa (not - or X)
  mismatch = mismatch[!aln_chars[[1]][mismatch] %in% c('-', 'X')]
  mismatch = mismatch[!aln_chars[[2]][mismatch] %in% c('-')]

  # unpack alignment into matrix
  aln_mat = do.call(rbind, aln_chars)
  rownames(aln_mat) = names(sequences)

  # calculate coverage and plot #
  coverage = .calculate_coverage(aln_mat, mismatch)

  # return alignment matrix and coverage plot
  return(list(
    aln_mat = aln_mat,
    coverage = coverage
  ))
}

# 6/17/25 -- try wrapping msa::msa with invisible() to hide message "gonnet"

# .map_aln_to_positions() ----

#' Map Aligned Amino Acid Positions to Nucleotide Positions
#'
#' Overlays alignment columns with residue-level and codon-level position indices for downstream
#' mapping of structural patches to nucleotide-level coordinates.
#'
#' Typically used after \code{.align_sequences()} and \code{identify_patches()} to annotate alignment matrix rows.
#'
#' @param aln_mat A character matrix from \code{.align_sequences()}, with two rows: reference and structure-derived.
#' @param residue_df Data frame of patch or interface residues, typically from \code{identify_patches()}.
#' @param chain Optional chain ID string (or vector) to restrict mapping to a subset of chains.
#'
#' @return Modified \code{aln_mat} with alignment positions replaced by codon and PDB residue indices.
#' @export
.map_aln_to_positions = function(aln_mat, residue_df, chain = NA){
  # aln_mat is from align_sequences() aln_mat$aln_mat
  # residue_df is from identify_windows()

  # filter for chain of interest #
  if(!any(is.na(chain))){
    residue_df = residue_df[residue_df$orig_chain %in% chain,]
  }

  # replace aa residues with codon position #
  # assuming first is ref sequence #
  ref_pos = which(aln_mat[1,]!='-')
  aln_mat[1,ref_pos] = 1:length(ref_pos)

  # for each set in residue_df update that row in aln_mat
  pos = residue_df$residue_id

  # fill position in pdb seq with position #
  ref_pos = which(aln_mat[2,]!='-')
  aln_mat[2,ref_pos] = pos

  return(aln_mat)
}

# .map_patches_to_codons() ----

#' Map Protein Patch Regions to Codon Coordinates
#'
#' Converts patch-level residue groupings (e.g., structural surface patches) into codon-aligned nucleotide coordinates
#' based on a dual-layer alignment matrix produced by \code{.map_aln_to_positions()}.
#'
#' This function creates \code{codon_patch} strings for each residue, encoding the codons that make up its structural context.
#' Also assigns a unique \code{msa_subset_id} used for downstream patch extraction.
#'
#' @param pos_aln A 2-row alignment matrix with codon indices (row 1) and residue indices (row 2), from \code{.map_aln_to_positions()}.
#' @param residue_df Data frame of patch-level residues, containing at minimum \code{residue_id} and \code{patch} columns.
#'
#' @return A data frame matching residues to codon-based patch identifiers, with columns:
#' \describe{
#'   \item{residue_id}{PDB residue identifier}
#'   \item{codon}{Codon position in the alignment}
#'   \item{codon_patch}{Codon-based patch identifier (e.g., "15+22+27")}
#'   \item{msa_subset_id}{Unique ID used to extract the corresponding MSA subset}
#' }
#' @export
.map_patches_to_codons = function(pos_aln, residue_df){

  # 6/14 update -- now we will add a msa_subset_id column #
  # msa are not built yet but the column is present when it is needed #
  # keeps submodular functions behaving well outside of modules and wrappers #

  # Create base df from ALL codons in pos_mat
  all_codons = pos_aln[1,]
  all_residues = pos_aln[2,]

  # Start with codon-based rows
  codon_df = data.frame(
    residue_id = all_residues,
    codon_patch = NA,
    codon = all_codons,
    stringsAsFactors = FALSE
  )

  # Add non-codon rows (interfaces, unmapped chains)
  non_codon_residues = residue_df$residue_id[!residue_df$residue_id %in% all_residues]
  if(length(non_codon_residues) > 0){
    extra_df = data.frame(
      residue_id = non_codon_residues,
      codon_patch = NA,
      codon = NA,
      stringsAsFactors = FALSE
    )
    codon_df = rbind(codon_df, extra_df)
  }

  # Now map patches to this complete df
  for(i in 1:nrow(residue_df)){
    patch = residue_df$patch[i]
    if(is.na(patch)) next

    patch_residues = unlist(strsplit(patch, '\\+'))
    codons = which(pos_aln[2,] %in% patch_residues)

    if(length(codons) > 0){
      codon_nums = pos_aln[1, codons]
      codon_nums = codon_nums[codon_nums != '-']

      if(length(codon_nums) > 0){
        # Use actual codon numbers without filling gaps
        codon_nums_sorted = sort(as.numeric(codon_nums))
        patch_str = paste0(codon_nums_sorted, collapse = '+')

        # Assign to matching residue
        codon_df$codon_patch[codon_df$residue_id == residue_df$residue_id[i]] = patch_str
      }
    }
  }

  # ADD msa_subset_id to table #
  codon_df$msa_subset_id = ifelse(
    !is.na(codon_df$codon),
    paste0('codon_', codon_df$codon),
    ifelse(
      grepl('^interface', codon_df$residue_id),
      codon_df$residue_id,  # interfaces already have valid names
      paste0('residue_', codon_df$residue_id)  # numeric residues need prefix
    )
  )

  return(codon_df)
}

# .extract_msa_subsets() ----

#.extract_msa_subsets() ----

#' Extract Codon-Aligned Nucleotide MSA Windows
#'
#' Subsets a nucleotide multiple sequence alignment based on codon-level patches
#' derived from 3D structural neighborhoods. Each subset corresponds to a structural patch.
#'
#' Assumes codon numbering starts at 1, with codon 1 = positions 1:3 in the MSA, codon 2 = 4:6, etc.
#' Used internally to extract windows for diversity and selection analysis.
#'
#' @param msa A nucleotide multiple sequence alignment in matrix form (e.g., from \code{ape::as.DNAbin()}).
#' @param codon_patches A data frame with \code{codon_patch} and \code{msa_subset_id} columns (from \code{.map_patches_to_codons()}).
#'
#' @return A named list of nucleotide MSA subsets, each a matrix corresponding to one patch.
#' @export
.extract_msa_subsets = function(msa, codon_patches){
  # msa is the alignment matrix #
  # codon_patches is the patches to extract #

  # drop patches that dont have codon positions #
  patches = codon_patches[!is.na(codon_patches$codon_patch),]

  msa_subset = list()
  for(i in 1:nrow(patches)){
    codons = patches$codon_patch[i]
    codon_pos = unlist(strsplit(codons, '\\+'))

    # convert codon positions to nucleotide positions #
    nuc_pos = c()
    for(codon in codon_pos){
      codon_num = as.numeric(codon)
      # each codon spans 3 nucleotides: codon 1 = nucs 1:3, codon 2 = nucs 4:6, etc.
      nuc_range = ((codon_num - 1) * 3 + 1):(codon_num * 3)
      nuc_pos = c(nuc_pos, nuc_range)
    }

    # subset #
    fasta = msa[, nuc_pos]
    msa_subset[[i]] = fasta

    # use the pre-computed msa_subset_id
    names(msa_subset)[i] = patches$msa_subset_id[i]
  }

  # return list of subsets #
  return(msa_subset)
}

# aln_msa_to_pdb() ----

#' Align Protein Structure to MSA and Extract Codon-Level Windows
#'
#' Module for aligning a structure-derived protein sequence to the reference sequence in a multiple
#' sequence alignment (MSA), mapping 3D patches to codon-level positions, and extracting corresponding
#' nucleotide windows. This forms the bridge between structure and sequence in the evo3D pipeline.
#'
#' This function is a core module in the evo3D workflow and is called internally by \code{run_evo3d()}.
#' It chains together multiple submodules to produce codon-aligned MSA subsets for each structural patch.
#'
#' @param msa_info A named list output from \code{WRAPPER_msa_to_ref()}, containing the MSA matrix and reference peptide sequence.
#' @param pdb_info A named list output from \code{WRAPPER_pdb_to_patch()}, including PDB-derived sequences and residue/patch annotations.
#' @param chain Character string indicating the chain to analyze. Required.
#' @param coverage_plot Logical; if \code{TRUE}, alignment coverage plot will be shown.
#' @param drop_unused_residues Logical; if \code{TRUE}, drops PDB residues not linked to codon-level info.
#'
#' @return A list with:
#' \describe{
#'   \item{aln_coverage}{Named list of alignment coverage ranges.}
#'   \item{aln_df}{Data frame mapping residues to codons, amino acids, and patch identifiers.}
#'   \item{msa_subsets}{Named list of nucleotide MSA matrices, one per patch.}
#' }
#' @export
aln_msa_to_pdb = function(msa_info, pdb_info, chain = NA, coverage_plot = FALSE, drop_unused_residues = TRUE){

  # msa_info ~ must be list object from WRAPPER_msa #
  # pdb_info ~ must be list object from WRAPPER_pdb #

  # chain must be specified
  if(missing(chain) || any(is.na(chain))) {
    stop("Chain must be specified. Use .auto_detect_chain() for automatic detection.")
  }

  # step 0: prep data #
  pep = msa_info$pep
  seq = pdb_info$seq_set[chain]

  residue_df = pdb_info$residue_df

  # step 1: align sequences # #NOTE LETS MOVE COVERAGE PLOT OUTSIDE ALIGNMENTS #
  aln = .align_sequences(c(pep, seq), generate_plot = coverage_plot)

  # step 2: map alignment to positions #
  pos_mat = .map_aln_to_positions(aln$aln_mat,
                                 residue_df,
                                 chain = chain)

  # step 3: map patches to nucleotides #
  codon_patches = .map_patches_to_codons(pos_mat, residue_df)

  # FILTER: drop unused PDB residues if requested
  if(drop_unused_residues) {
    keep_rows = !is.na(codon_patches$codon) |
      !is.na(codon_patches$codon_patch) |
      grepl('^interface_', codon_patches$residue_id)

    codon_patches = codon_patches[keep_rows, ]
  }

  # step 3.5 (adding aa info to codon_patches) #

  # match codon positions to get the column indices
  matches = match(codon_patches$codon, pos_mat[1,])

  # pull amino acids directly from alignment matrix
  codon_patches$ref_aa = aln$aln_mat[1, matches]
  codon_patches$pdb_aa = aln$aln_mat[2, matches]

  # step 4: grab subsets of MSA #
  msa_subsets = .extract_msa_subsets(msa_info$msa_mat, codon_patches)

  # return data #
  return(list(
    aln_coverage = aln$coverage,
    aln_df = codon_patches,
    msa_subsets = msa_subsets)
  )

}
