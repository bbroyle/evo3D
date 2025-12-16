# --------------------------------------------------------------- #
# ALN MODULE — align msa to pdb and derive codon-level windows
# internal .functions() are wrapped by exported aln_msa_to_pdb()
#
# NOTES:
# 1. alignment uses msa::msa(…, method = "ClustalOmega", type = "protein")
#    default sub matrix BLOSUM65; order = "input"
# 2. mismatches are computed on aa-to-aa positions only (ignores '-' and 'X')
# 3. .map_aln_to_positions() overlays indices:
#    ref codon positions and pdb residue_ids are added as extra columns
# 4. patches are converted to codon ids via .map_patches_to_codons()
#    codon_id format: "<codon>_<msaLabel>"; patch strings are '+'-delimited
# 5. gap-map handling:
#    - if a residue aligns to '-', it is considered gap-mapped
#    - drop_gap_map = TRUE removes gap-centered patches and prunes gap members
# 6. cleanup/rebuild step:
#    .fix_gap_map_and_dedup_codons() deduplicates codons and, when max_patch is set,
#    rebuilds patches from nearest valid residues using pdb_info$residue_dist
#    (optional dist_cutoff and only_exposed_in_patch respected)
# 7. patch modes:
#    - "codon": size/uniqueness tracked at codon level
#    - "residue": size built on residues, codons assigned afterward
# 8. interfaces:
#    residue_ids beginning with "interface" are held out during rebuild
#    and reattached unchanged at the end
# 9. multi-MSA support:
#    patches may span multiple MSAs (e.g., "1_msa1+2_msa1+1_msa2");
#    subsets are extracted per MSA and cbind()-merged
# 10. outputs:
#    - coverage: per-run alignment metadata (no aln matrix)
#    - aln_df: residue ↔ codon mapping with patch/codon strings and metadata
#    - msa_subsets: nucleotide windows per patch (samples × nucleotides)
#
# email: bbroyle@purdue.edu
# --------------------------------------------------------------- #

# .auto_detect_chain() ----

#' Automatically detect matching PDB chains by k-mer coverage
#'
#' Compare a peptide sequence against each chain in a PDB structure and return
#' similarity scores based on k-mer coverage. Higher scores indicate closer
#' correspondence between the query peptide and a given PDB chain.
#'
#' This function is intended for internal use in automatic chain-mapping
#' workflows.
#'
#' @param pep Amino acid sequence as a single character string.
#' @param pdb PDB structure input; may be a file path, a \code{bio3d} PDB object,
#'   or a standardized PDB object. A PDB identifier alone is not supported.
#' @param k Integer; k-mer size used for coverage calculation (default 4).
#' @param in_module Logical; internal flag. If \code{TRUE}, skips PDB validation
#'   steps handled upstream (default \code{FALSE}).
#'
#' @return Named numeric vector of k-mer coverage scores for each PDB chain,
#'   sorted in descending order. Names correspond to chain identifiers.
#'
#' @details
#' Coverage is defined as the fraction of k-mers in a PDB chain sequence that
#' are also present in the query peptide sequence. This measure is asymmetric
#' and reflects how well a PDB chain is represented within the query.
#'
#' @keywords internal
.auto_detect_chain = function(pep, pdb, k = 4, in_module = F){

  # Changed to coverage instead of jaccard
  kmer_coverage = function(pdb_seq, msa_seq) {
    # seq is too short for kmer - just return 0
    if (nchar(pdb_seq) < k || nchar(msa_seq) < k) return(0)

    pdb_kmers = substring(pdb_seq, 1:(nchar(pdb_seq) - k + 1), k:(nchar(pdb_seq)))
    msa_kmers = substring(msa_seq, 1:(nchar(msa_seq) - k + 1), k:(nchar(msa_seq)))

    # What fraction of PDB kmers are found in MSA?
    return(length(intersect(pdb_kmers, msa_kmers)) / length(pdb_kmers))
  }

  # if not in module validate pdb (handled in .get_pdb_sequences()) #
  seq_set = .get_pdb_sequence(pdb, in_module = in_module)

  dist = sapply(seq_set, function(x) kmer_coverage(x, pep))

  # sort by descending order and return
  dist = sort(dist, decreasing = T)

  return(dist)
}


# .calculate_coverage() ----

#' Calculate aligned coverage ranges
#'
#' Identify contiguous non-gap regions in aligned reference and structure
#' sequences and return them as start:end ranges. Also returns the alignment
#' indices of mismatched positions (excluding gaps and 'X').
#'
#' Typically used downstream of \code{.map_aln_to_positions()} to summarize
#' aligned coverage and mismatch locations.
#'
#' @param pos_sets A single position-mapped alignment matrix or a named list of
#'   such matrices, as returned by \code{.map_aln_to_positions()}. Each matrix
#'   must contain columns \code{ref_aa} and \code{pdb_aa}.
#'
#' @return A list (same length as \code{pos_sets}). Each element is a list with:
#' \describe{
#'   \item{ref_aa}{Character vector of non-gap ranges in the reference sequence
#'     (e.g., \code{"5:25"}).}
#'   \item{pdb_aa}{Character vector of non-gap ranges in the structure sequence.}
#'   \item{mismatch}{Integer vector of alignment indices with mismatched residues.}
#' }
#' @keywords internal

.calculate_coverage = function(pos_sets) {

  # normalize input to list
  if (is.matrix(pos_sets)) {
    pos_sets = list(pos_sets)
  }

  collapse_ranges = function(idx) {
    if (length(idx) == 0L) return(character(0))
    brks = c(0L, which(diff(idx) > 1L), length(idx))
    sapply(seq_len(length(brks) - 1L), function(i) {
      a = idx[brks[i] + 1L]
      b = idx[brks[i + 1L]]
      if (a == b) as.character(a) else paste0(a, ":", b)
    })
  }

  out = vector("list", length(pos_sets))
  names(out) = names(pos_sets)

  for (k in seq_along(pos_sets)) {

    pos_mat = pos_sets[[k]]

    ref_idx = which(pos_mat[, "ref_aa"] != "-")
    pdb_idx = which(pos_mat[, "pdb_aa"] != "-")

    mismatch_idx = which(
      pos_mat[, "ref_aa"] != "-" &
        pos_mat[, "pdb_aa"] != "-" &
        pos_mat[, "ref_aa"] != "X" &
        pos_mat[, "pdb_aa"] != "X" &
        pos_mat[, "ref_aa"] != pos_mat[, "pdb_aa"]
    )

    out[[k]] = list(
      ref_aa   = collapse_ranges(ref_idx),
      pdb_aa   = collapse_ranges(pdb_idx),
      mismatch = mismatch_idx
    )
  }

  out
}

# .align_sequences() ----

#' Align reference and structure sequences
#'
#' Align a reference amino acid sequence against a structure-derived sequence
#' using Clustal Omega via the \pkg{msa} package.
#'
#' Intended for internal use when mapping reference sequences to PDB-derived
#' sequences prior to downstream patch analyses.
#'
#' @param sequences Named character vector of length two containing the reference
#'   sequence and the structure-derived sequence. Names are used only to preserve
#'   input order (e.g., \code{c(ref = "MKT...", pdb = "MK...")}).
#' @param submat Character; substitution matrix used for alignment. One of
#'   \code{"BLOSUM30"}, \code{"BLOSUM40"}, \code{"BLOSUM50"}, \code{"BLOSUM65"},
#'   \code{"BLOSUM80"}, or \code{"Gonnet"} (default \code{"BLOSUM65"}).
#'
#' @return List with:
#' \describe{
#'   \item{aln_mat}{Two-column character matrix with aligned sequences. Columns are
#'   \code{ref_aa} and \code{pdb_aa}.}
#' }
#'
#' @keywords internal

.align_sequences = function(sequences, submat = 'BLOSUM65'){

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
  aln = suppressMessages(msa::msa(sequences, method = 'ClustalOmega', type = 'protein',
                                  order = 'input', substitutionMatrix = submat,
                                  verbose = FALSE))

  aln_set = as.character(aln@unmasked)
  aln_chars = strsplit(aln_set, '')

  # unpack alignment into matrix
  aln_mat = do.call(cbind, aln_chars)
  colnames(aln_mat) = c('ref_aa', 'pdb_aa')

  # return alignment matrix and coverage plot
  return(list(
    aln_mat = aln_mat
  ))
}

# .map_aln_to_positions() ----

#' Map Alignment Columns to Codon and Residue Indices
#'
#' Annotates an alignment matrix by overlaying reference codon positions and
#' PDB residue indices onto aligned amino acids. This preserves the aligned
#' sequence characters while adding positional context for downstream mapping.
#'
#' Typically used after \code{.align_sequences()} and \code{identify_patches()}
#' to connect reference codons with structure-derived residues.
#'
#' @param aln_mat A character matrix returned by \code{.align_sequences()},
#'   with two rows (reference and structure-derived sequences).
#' @param residue_df Data frame of patch or interface residues, such as that
#'   returned by \code{identify_patches()}, containing at least
#'   \code{residue_id} and \code{orig_chain}.
#' @param chain Optional chain identifier (string or vector). If supplied,
#'   restricts mapping to residues belonging to the given chain(s).
#'
#' @return A modified character matrix with additional columns:
#' \describe{
#'   \item{codon}{Reference codon indices for aligned positions.}
#'   \item{residue_id}{PDB residue indices aligned to the structure sequence.}
#' }
#' @keywords internal

.map_aln_to_positions = function(aln_mat, residue_df, chain = NA){

  # filter for chain of interest #
  if(!any(is.na(chain))){
    residue_df = residue_df[residue_df$orig_chain %in% chain,]
  }

  # replace ref aa residues with codon position #
  # assuming first is ref sequence #

  # lets keep aa info in tack -- otherwise we just merge later #
  aln_mat = cbind(aln_mat, '-')

  ref_pos = which(aln_mat[,1]!='-')
  aln_mat[ref_pos,5] = 1:length(ref_pos)

  # for each set in residue_df update that row in aln_mat
  pos = residue_df$residue_id

  aln_mat = cbind(aln_mat, '-')
  # fill position in pdb seq with position #
  aa_pos = which(aln_mat[,2]!='-')
  aln_mat[aa_pos,6] = pos

  # rename
  colnames(aln_mat)[5:6] = c('codon', 'residue_id')

  return(aln_mat)
}

# .map_patches_to_codons() ----

#' Convert Residue Patches to Codon-Based Coordinates
#'
#' Maps structural patch definitions from residue space into codon-aligned
#' coordinates, using a MSA-PDB alignment table from
#' \code{.map_aln_to_positions()}. This produces codon-based patch identifiers
#' that can be used for downstream nucleotide-level extraction.
#'
#' Each residue inherits codon patch membership, with options to remove
#' gap-mapped residues. When gaps are dropped, both gap-seeded patches and
#' membership of gap residues in other patches are removed. If a residue
#' distance matrix is supplied, \code{max_dist} is recalculated after
#' gap-removal.
#'
#' @param aln_table Alignment table from \code{.map_aln_to_positions()},
#'   containing at least columns \code{codon}, \code{msa}, \code{pdb},
#'   \code{residue_id}, \code{ref_aa}, and \code{pdb_aa}.
#' @param residue_df Data frame of patch-level residues (e.g. from
#'   \code{pdb_to_patch()}), with columns including \code{residue_id},
#'   \code{patch}, \code{patch_len}, \code{exposed}, and \code{max_dist}.
#' @param drop_gap_map Logical; whether to drop gap-mapped residues. If
#'   \code{TRUE} (default), patches centered on gap residues are removed,
#'   and gap residues are excluded from other patches.
#' @param dist_mat Optional residue–residue distance matrix (from
#'   \code{pdb_to_patch()}). If supplied and \code{drop_gap_map = TRUE},
#'   \code{max_dist} is updated after gap removal.
#'
#' @return A data frame aligning residues with codon-based patch information,
#'   including:
#'   \describe{
#'     \item{codon_id}{Codon identifier string (e.g. \code{"15_msa1"}).}
#'     \item{codon_patch}{“+”-delimited string of codon IDs making up the patch.}
#'     \item{codon_len}{Number of codons in the patch.}
#'     \item{unique_codon}{Number of unique codons represented.}
#'     \item{gap_map_count}{Number of gap-mapped residues removed from the patch.}
#'     \item{patch, patch_len}{Residue-based patch membership (cleaned).}
#'     \item{exposed, max_dist}{Exposure and distance metadata (updated if applicable).}
#'   }
#' @keywords internal


.map_patches_to_codons = function(aln_table, residue_df, drop_gap_map = TRUE, dist_mat = NA){

  # in aln_msa_to_pdb() always drop gap map #
  # in run_evo3d() always drop gap map, but dont touch dist_mat - later function #
  # as manual function drop_gap_map can be turned off #

  # copy over original residue patch, patch len, exposure, and max_dist #
  aln_table$orig_patch = residue_df$patch[match(aln_table$residue_id, residue_df$residue_id)]
  aln_table$orig_patch_len = residue_df$patch_len[match(aln_table$residue_id, residue_df$residue_id)]
  aln_table$exposed = residue_df$exposed[match(aln_table$residue_id, residue_df$residue_id)]
  aln_table$max_dist = residue_df$max_dist[match(aln_table$residue_id, residue_df$residue_id)]

  # Open new columns
  aln_table$codon_patch = NA
  aln_table$codon_len = NA
  aln_table$unique_codon = NA
  aln_table$gap_map_count = NA

  aln_table$patch = NA
  aln_table$patch_len = NA

  # Set up codon_id
  aln_table$codon_id = paste0(aln_table$codon, '_', aln_table$msa)
  aln_table$codon_id[is.na(aln_table$codon)] = NA
  aln_table$codon_id[aln_table$codon == '-'] = NA

  # check for gap_map residues #
  gap_ids = aln_table$residue_id[which(aln_table$codon == '-')]

  # if gap_map is dropped remove gap_map seeds #
  if(drop_gap_map){
    aln_table$orig_patch[aln_table$residue_id %in% gap_ids] = NA
  }

  # for each row - copy patch, *remove gaps*, return codons #
  for(i in seq_len(nrow(aln_table))){
    if(is.na(aln_table$orig_patch[i])) next

    p = strsplit(aln_table$orig_patch[i], '\\+')[[1]]

    # count gap map #
    g = which(p %in% gap_ids)

    if(drop_gap_map){
      p = setdiff(p, gap_ids)
    }

    # find codon positions #
    c = aln_table$codon_id[match(p, aln_table$residue_id)]

    # fill table #
    aln_table$patch[i] = paste0(p, collapse = '+')
    aln_table$patch_len[i] = length(p)

    aln_table$codon_patch[i] = paste0(c, collapse = '+')
    aln_table$codon_len[i] = length(c)
    aln_table$unique_codon[i] = length(unique(c))

    aln_table$gap_map_count[i] = length(g)

    # update max_dist if available #
    if(length(g) > 0 && is.matrix(dist_mat) && drop_gap_map){
      # replace gap map dist #
      id = aln_table$residue_id[i]

      if(grepl('^interface', id)) next

      aln_table$max_dist[i] = max(dist_mat[id, p], na.rm = T)
    }

  }

  # clean up table #

  # drop original patch info #
  aln_table = aln_table[ , !(names(aln_table) %in% c("orig_patch", "orig_patch_len")) ]

  # reorder rest of columns #
  # codon_id, codon, msa, pdb, residue_id, ref_aa, pdb_aa, codon_patch,
  # codon_len, unique_codon, gap_map_count, exposed, max_dist, patch, patch_len #

  # and anything else #
  # your desired order
  wanted_order = c(
    "codon_id", "codon", "msa", "pdb", "residue_id", "ref_aa", "pdb_aa",
    "codon_patch", "codon_len", "unique_codon", "gap_map_count",
    "exposed", "max_dist", "patch", "patch_len"
  )

  # columns in data that match your desired order
  present = intersect(wanted_order, names(aln_table))

  # everything else (not in your desired order)
  extra = setdiff(names(aln_table), wanted_order)

  # final reordering
  aln_table = aln_table[ , c(present, extra)]


  # last clean up if codon_patch is '' set all codon info to NA #
  # possible they are rebuilt - so save them #
  #aln_table$codon_patch[aln_table$codon_patch == ''] = NA
  #aln_table$codon_len[aln_table$codon_patch == ''] = NA
  #aln_table$unique_codon[aln_table$codon_patch == ''] = NA

  return(aln_table)
}


# .extract_msa_subsets_single() ----

#' Extract a Codon-Aligned Subset from a single MSA
#'
#' Internal helper to \code{.extract_msa_subsets()}. Converts codon indices
#' (without MSA labels) into nucleotide positions (3 bases per codon) and
#' subsets a single MSA accordingly.
#'
#' @param msa A nucleotide alignment matrix (samples × sites).
#' @param codon_patches Data frame with \code{msa_subset_id} and
#'   \code{codon_patch} (a “+”-delimited string of codon indices).
#'
#' @return A list of one or more nucleotide MSA subsets, each labeled by
#'   \code{msa_subset_id}.
#' @keywords internal

.extract_msa_subsets_single = function(msa, codon_patches, seq_type = 'nucleotide'){
  # works across MSA's #
  # drop patches that dont have codon positions #
  # now just works on one row at a time? #
  patches = codon_patches[!is.na(codon_patches$codon_patch),]

  msa_subset = list()
  for(i in 1:nrow(patches)){
    codons = patches$codon_patch[i]
    codon_pos = unlist(strsplit(codons, '\\+'))

    # convert codon positions to nucleotide positions #
    nuc_pos = c()
    if(seq_type == 'protein'){
      nuc_pos = as.numeric(codon_pos)
    } else {
      for(codon in codon_pos){
        codon_num = as.numeric(codon)
        # each codon spans 3 nucleotides: codon 1 = nucs 1:3, codon 2 = nucs 4:6, etc.
        nuc_range = ((codon_num - 1) * 3 + 1):(codon_num * 3)
        nuc_pos = c(nuc_pos, nuc_range)
      }

    }

    # subset #
    fasta = msa[, nuc_pos, drop = FALSE]
    msa_subset[[i]] = fasta

    # use the pre-computed msa_subset_id
    names(msa_subset)[i] = patches$msa_subset_id[i]
  }

  # return list of subsets #
  return(msa_subset)
}

# .extract_msa_subsets() ----

#' Extract Codon-Aligned Nucleotide Subsets Across Multiple MSAs
#'
#' Subsets one or more nucleotide multiple sequence alignments into
#' codon-aligned windows corresponding to structural patches. Codon patches
#' may reference multiple MSAs (e.g., \code{"1_msa1+2_msa1+1_msa2"}), in which
#' case codons are extracted from each MSA separately and then merged.
#'
#' Typically used after \code{.map_patches_to_codons()} to build
#' nucleotide-level MSAs for patch-based diversity and selection analyses.
#'
#' @param msa_info A named list of MSA entries. Each entry must contain
#'   \code{msa_mat}, a nucleotide alignment matrix (samples × sites), usually
#'   from \code{ape::as.DNAbin()} or similar.
#' @param codon_patches Data frame with at least \code{codon_patch} and
#'   \code{msa_subset_id} columns, as produced by
#'   \code{.map_patches_to_codons()}. Codon patches are represented as
#'   “+”-delimited strings of codon indices, each suffixed by its MSA label
#'   (e.g., \code{"15_msa1+22_msa1+7_msa2"}).
#' @param use_sample_names Logical; whether to collapse cross-MSA or cross-chain
#'   patches so that resulting subsets share the same FASTA identifiers.
#'   Default is \code{TRUE}.
#'
#' @return A named list of nucleotide alignment subsets, one per patch. Each
#'   element is a merged matrix (samples × nucleotides) containing the codons
#'   drawn from one or more MSAs according to the patch definition.
#' @keywords internal

.extract_msa_subsets = function(msa_info, codon_patches, use_sample_names = TRUE,
                                genetic_code = 1, fill_char = '-'){
  # works across MSA's #
  # schedules extract_msa_subsets_single #

  # still need to incorporate use_sample_names #

  # two options --
  # 1. pull apart patches into msa1, msa2, msaX columns and schedule blocks to extract_msa_subsets #
  # 2. send row by row -- still pulling apart #
  # -- doing option 2

  # add path for protein MSA #
  msas = unique(codon_patches$msa)
  seqtypes = unlist(lapply(msas, function(x){
    msa_info[[x]]$seq_type
  }))

  global_seqtype = ifelse(any(seqtypes == 'protein'), 'protein', 'nucleotide')

  # add new column for each identifier #
  msa_subsets = list()
  for(i in seq_len(nrow(codon_patches))){
    msa_subsets[[i]] = NULL
    p = codon_patches$codon_patch[i]
    if(is.na(p)) next

    p = strsplit(p, '\\+')[[1]]
    id = codon_patches$msa_subset_id[i]

    # separate codons from msa_ids # (assumes id always present)
    codons = gsub('_.+', '', p)
    msas = gsub('.+_', '', p)

    msa_set = unique(msas)
    msa_list = list()
    for(j in seq_len(length(msa_set))){
      msa_name = msa_set[j]
      codons2 = codons[msas == msa_name]
      seqtype = msa_info[[msa_name]]$seq_type
      hold = .extract_msa_subsets_single(msa_info[[msa_name]]$msa_mat,
                                                  data.frame(msa_subset_id = id, codon_patch = paste(codons2, collapse = '+')),
                                                  seq_type = seqtype)[[1]]


      # if seqtype nuc and global is protein - then translate these #
      if(seqtype == 'nucleotide' && global_seqtype == 'protein'){
        # should worry about '---' in frame being being '-' instead of 'X'
        len = ncol(hold)
        gap_grid = matrix(0, nrow = nrow(hold), ncol = len / 3)

        for(ii in seq_len(ncol(gap_grid))){
          st = ii * 3 - 2
          en = ii * 3

          ro = which(apply(hold[,st:en, drop = FALSE], 1, function(x){all(x == '-')}))
          gap_grid[ro,ii] = 1

        }

        hold = t(apply(hold, 1, function(x){
          y = seqinr::translate(x, numcode = genetic_code)
        }))

        # if ncol was 3 need to reformat #
        if(len == 3){
          hold = t(hold)
        }

        # 12.9 what was this fixing ? #
        #if(length(hold) == 1){
        #  rn = colnames(hold)
        #  rownames(hold) = rn
        #  colnames(hold) = NULL
        #}

        pos = which(gap_grid == 1)
        hold[pos] = '-'

      }

      msa_list[[j]] = hold

    }

    # right here -- there needs to be a checking for sample_names
    msa_subset = .join_msa_subsets(msa_list, fill_char = fill_char, use_sample_names = use_sample_names)

    msa_subsets[[i]] = msa_subset
    names(msa_subsets)[i] = id

  }

  # clean up set to only the ones with pulls #
  keep = !vapply(msa_subsets, is.null, logical(1))
  msa_subsets = msa_subsets[keep]

  # just return
  return(msa_subsets)
}

# .join_msa_subsets() ----

#' Join MSA subset matrices
#'
#' Combine a list of MSA subset matrices column-wise into a single matrix. When
#' \code{use_sample_names = TRUE}, only samples present in all subsets (common
#' rownames) are retained prior to joining. When \code{FALSE}, subsets are joined
#' by row order and shorter matrices are padded with \code{fill_char}.
#'
#' @param msa_list List of MSA subset matrices to join.
#' @param fill_char Character used to pad missing rows when
#'   \code{use_sample_names = FALSE} (default \code{"-"}).
#' @param use_sample_names Logical; if \code{TRUE} (default), join on common
#'   rownames. If \code{FALSE}, join by row order with padding.
#'
#' @return A character matrix formed by column-binding the input subsets.
#'
#' @keywords internal

.join_msa_subsets = function(msa_list, fill_char = '-', use_sample_names = TRUE){


  # if joining on use_sample_names just take whats in common #
  if(use_sample_names){
    name_list = lapply(msa_list, rownames)
    common_names = Reduce(intersect, name_list)

    # ADD STOP MESSAGE IF no common_names #
    if(!length(common_names)){
      stop('NO COMMON NAMES - try turning off use_sample_names in aln_controls')
    }

    # cbind these common names #
    len = length(msa_list)
    sub_set = list()
    for(i in seq_len(len)){
      sub_set[[i]] = msa_list[[i]][match(common_names, rownames(msa_list[[i]])),,drop=FALSE]
    }

    msa_set = do.call(cbind, sub_set)
  } else {
    # so can join different lengths but need to fill shorter msas #
    # replace with rowname combo? or generic row name #

    # what is the row lengths of these sets
    lens = unlist(lapply(msa_list, nrow))
    mlen = max(lens)

    for(i in seq_along(lens)){
      if(nrow(msa_list[[i]]) == mlen) next

      # fill in rows #
      clen = nrow(msa_list[[i]])
      flen = mlen - clen

      nmat = matrix(fill_char, nrow = flen, ncol = ncol(msa_list[[i]]))
      rownames(nmat) = paste0('artificial_fill', seq_len(flen))

      msa_list[[i]] = rbind(msa_list[[i]], nmat)

    }

    # all equal rows #
    msa_set = do.call(cbind, msa_list)

  }

  return(msa_set)
}


# .fix_gap_map_and_dedup_codons ----

#' Fix Gap-Mapped Residues and Deduplicate Codons in Patches
#'
#' Cleans and rebuilds patch definitions by removing gap-mapped residues and
#' ensuring codons are not duplicated within codon-based patches. When a maximum
#' patch size (\code{max_patch}) is specified, patches are rebuilt from the
#' residue distance matrix to satisfy size constraints.
#'
#' Behavior depends on \code{patch_mode}:
#' \itemize{
#'   \item In \code{"codon"} mode, codons within each \code{codon_patch} are
#'   deduplicated. If \code{max_patch} is set, patches smaller than the target
#'   are rebuilt using nearest valid codons until the limit is reached.
#'   \item In \code{"residue"} mode, patches are rebuilt at the residue level,
#'   with codons reassigned from residue membership.
#' }
#'
#' Interface residues (with IDs beginning \code{"interface"}) are held out and
#' reattached after processing.
#'
#' @param residue_df Data frame of patch-level residues, typically from
#'   \code{.map_patches_to_codons()}, with columns \code{residue_id},
#'   \code{patch}, \code{codon_patch}, \code{codon}, \code{pdb_aa},
#'   \code{patch_len}, \code{codon_len}, \code{unique_codon}, and
#'   \code{exposed}.
#' @param patch_mode Either \code{"codon"} or \code{"residue"}, indicating how
#'   patches are rebuilt and deduplicated.
#' @param dist_cutoff Numeric; optional maximum residue–residue distance used
#'   when rebuilding patches. Distances greater than this are excluded.
#' @param max_patch Integer; maximum allowed patch size. If \code{NA}, no patch
#'   rebuilding is performed.
#' @param pdb_info PDB metadata list containing at least \code{residue_dist}, a
#'   residue–residue distance matrix.
#' @param only_exposed_in_patch Logical; whether only exposed residues can be
#'   included when rebuilding patches.
#'
#' @return The input \code{residue_df} with updated patch definitions. Columns
#'   \code{patch}, \code{patch_len}, \code{codon_patch}, \code{codon_len}, and
#'   \code{unique_codon} are rebuilt as needed. Interface residues are preserved
#'   and reattached unchanged.
#' @keywords internal

.fix_gap_map_and_dedup_codons = function(residue_df, patch_mode, dist_cutoff,
                                    max_patch, pdb_info, only_exposed_in_patch) {

  # if patch_mode = codon, dedup whole table #
  if (patch_mode == "codon") {
    # Deduplicate codons within each codon_patch
    for (i in seq_len(nrow(residue_df))) {
      p = residue_df$codon_patch[i]
      if (is.na(p)) next
      p = unique(strsplit(p, "\\+")[[1]])
      residue_df$unique_codon[i] = length(p)
      residue_df$codon_patch[i] = paste(p, collapse = "+")
    }
  }


  # ---- Variable patch scheme ----
  if (is.na(max_patch)) {
    # nothing to do
    return(residue_df)
  }

  # ---- Fixed patch scheme (max_patch is set) ----

  # Step 0: separate interfaces
  ro = grep("^interface", residue_df$residue_id)
  interface_df = NULL
  if (length(ro)) {
    interface_df = residue_df[ro, ]
    residue_df = residue_df[-ro, ]
  }

  # Step 1: rebuild patches depending on mode
  if (patch_mode == "codon") {
    ro = which(!is.na(residue_df$codon_patch) & residue_df$unique_codon < max_patch)

    if(length(ro)){
      # grab valid replacements #
      valid = residue_df[residue_df$codon != "-" &
                            residue_df$pdb_aa != "-" &
                            (if (only_exposed_in_patch) residue_df$exposed else TRUE),
                          c("codon_id", "residue_id"), drop = FALSE]

      # loop through each row and rebuild with valid sets #
      if(nrow(valid)){

        # trim dist_mat
        d = pdb_info$residue_dist
        d = d[, valid$residue_id, drop = FALSE]

        # make mapping
        res2cod = stats::setNames(valid$codon_id, valid$residue_id)

        # loop through ro #
        for (i in ro) {
          resi = residue_df$residue_id[i]
          if (!resi %in% rownames(d)) next

          d2 = sort(d[resi, , drop = TRUE])
          if (!is.na(dist_cutoff)) d2 = d2[d2 <= dist_cutoff]
          if (!length(d2)) next

          cod_seq = res2cod[names(d2)]
          uniq_cod = cod_seq[!duplicated(cod_seq)]
          set = uniq_cod[seq_len(min(max_patch, length(uniq_cod)))]

          residue_df$patch_len[i] = length(set)
          residue_df$patch[i] = paste(names(set), collapse = "+")
          residue_df$codon_patch[i] = paste(set, collapse = "+")
          residue_df$codon_len[i] = length(set)
          residue_df$unique_codon[i]= length(unique(set)) # codon_len and unique_codon should be the same #
        }
      }


      }
    } else if (patch_mode == "residue") {
      ro = which(!is.na(residue_df$codon_patch) & residue_df$codon_len < max_patch)

      if(length(ro)){
        # grab valid replacements #
        valid = residue_df[residue_df$codon != "-" &
                              residue_df$pdb_aa != "-" &
                              (if (only_exposed_in_patch) residue_df$exposed else TRUE),
                            c("codon_id", "residue_id"), drop = FALSE]

        # loop through each row and rebuild with valid sets #
        if(nrow(valid)){

          # trim dist_mat
          d = pdb_info$residue_dist
          d = d[, valid$residue_id, drop = FALSE]

          # make mapping
          res2cod = stats::setNames(valid$codon_id, valid$residue_id)

          # loop through ro #
          for (i in ro) {
            resi = residue_df$residue_id[i]
            if (!resi %in% rownames(d)) next

            d2 = sort(d[resi, , drop = TRUE])
            if (!is.na(dist_cutoff)) d2 = d2[d2 <= dist_cutoff]
            if (!length(d2)) next

            take_k = min(length(d2), max_patch)
            set = names(d2)[seq_len(take_k)]

            residue_df$patch_len[i] = length(set)
            residue_df$patch[i] = paste(set, collapse = "+")
            residue_df$codon_patch[i] = paste(res2cod[set], collapse = "+")
            residue_df$codon_len[i] = length(res2cod[set])
            residue_df$unique_codon[i] = length(unique(res2cod[set]))
          }
        }
      }
    }

  # Step 2: reattach interfaces
  if (!is.null(interface_df)) {
    residue_df = rbind(residue_df, interface_df)
  }

  # last step set exposed and max_dist to NA for gap mapped positions #
  ro = grepl('-', residue_df$codon)
  residue_df$exposed[ro] = NA
  residue_df$max_dist[ro] = NA

  return(residue_df)
}
# aln_msa_to_pdb() ----

#' Align Structure to Reference MSA and Generate Codon-Level Windows
#'
#' Aligns a structure-derived protein sequence to the reference sequence in one
#' or more nucleotide MSAs, maps structural patches to codon positions, and
#' extracts corresponding nucleotide windows. This function links 3D structure
#' to sequence diversity by chaining together multiple submodules:
#' \enumerate{
#'   \item Sequence alignment of PDB vs. MSA reference (\code{.align_sequences()}).
#'   \item Mapping aligned positions to codon and residue indices
#'         (\code{.map_aln_to_positions()}).
#'   \item Converting structural patches to codon-level definitions
#'         (\code{.map_patches_to_codons()}).
#'   \item Optional patch cleanup: gap-map removal and codon deduplication
#'         (\code{.fix_gap_map_and_dedup_codons()}).
#'   \item Extraction of nucleotide MSA subsets for each patch
#'         (\code{.extract_msa_subsets()}).
#' }
#'
#' This is the core alignment module in the evo3D workflow and is called
#' internally by \code{run_evo3d()}.
#'
#' @param msa_info A named list from \code{msa_to_ref()}, containing
#'   at least \code{msa_mat} (nucleotide MSA matrix) and \code{pep}
#'   (reference peptide sequence).
#' @param pdb_info A named list from \code{WRAPPER_pdb_to_patch()}, including
#'   \code{seq_set} (PDB-derived peptide sequences), \code{residue_df}
#'   (residue/patch annotations), and \code{residue_dist} (residue distance
#'   matrix).
#' @param chain Character string specifying the PDB chain to analyze. Ignored if
#'   \code{run_grid} is provided.
#' @param subset_msa Logical; if \code{TRUE}, nucleotide MSA subsets are
#'   extracted for each patch. Default \code{TRUE}.
#' @param verbose Integer; level of console reporting. Default \code{1}.
#' @param run_grid Optional data frame specifying combinations of MSA, PDB, and
#'   chain to process. If supplied, overrides \code{chain}.
#' @param drop_gap_map Logical; if \code{TRUE}, drops patches centered on
#'   gap-mapped residues. Default \code{TRUE}.
#' @param fix_gap_map_and_dedup_codons Logical; if \code{TRUE}, calls
#'   \code{.fix_gap_map_and_dedup_codons()} to remove gap residues from patch
#'   membership and ensure codon deduplication. Default \code{TRUE}.
#' @param patch_mode Either \code{"codon"} or \code{"residue"}, controlling
#'   how patches are deduplicated and rebuilt.
#' @param max_patch Integer; maximum patch size. If \code{NA}, patches are
#'   variable-length.
#' @param dist_cutoff Numeric; maximum residue–residue distance when expanding
#'   patches. Default \code{NULL} (no cutoff).
#' @param only_exposed_in_patch Logical; if \code{TRUE}, restricts patch
#'   rebuilding to exposed residues only.
#' @param user_aln Optional alignment object previously returned by
#'   \code{aln_msa_to_pdb()} that has been manually modified (e.g., via
#'   \code{adjust_aln()}). If provided, this alignment is reused instead of
#'   recomputing a new alignment.
#'
#' @return A list with three components:
#' \describe{
#'   \item{coverage}{List of alignment metadata for each MSA–PDB–chain run.}
#'   \item{aln_df}{Data frame mapping residues to codons, amino acids, and patch
#'         identifiers, after gap/deduplication cleanup.}
#'   \item{msa_subsets}{Named list of nucleotide MSA matrices, one per patch
#'         (if \code{subset_msa = TRUE}).}
#' }
#' @export

aln_msa_to_pdb = function(msa_info, pdb_info, chain = 'auto',
                          subset_msa = TRUE,
                          verbose = 1,
                          run_grid = NULL,
                          drop_gap_map = TRUE, # aln will drop patches centered on gap map residues
                          fix_gap_map_and_dedup_codons = TRUE, # aln will drop gap map residues from patch membership
                          patch_mode, max_patch, dist_cutoff, only_exposed_in_patch,
                          user_aln = NULL){

  # msa_info ~ must be list object from msa_to_ref() # ~ can be list of msa_infos
  # pdb_info ~ must be list object from ONE pdb_to_patch() #

  if(is.null(user_aln)){
  # step 0: prep data ----
  # in module always a run grid -- outside module can be chain specified #

  # If run grid supplied - ignore chain info #
  if(!is.null(run_grid)){
    chain = 'run_grid'
  }

  # at the end lets make sure to have a run_grid format
  if(F){
    # probably dont need chain mapping here 9.30.25 #
    # chain must be specified or 'auto'
    if(chain == 'auto') {
      chain = .auto_detect_chain(msa_info$pep, pdb_info$pdb)
      cat('\tDetected chain:', names(chain)[1], '\n')
      cat('\tAt', round(chain[1], 3) * 100, '% identity', '\n')
      chain = names(chain)[1]  # use the first chain detected
    }
  }

  # step 1: align sequences ----
  # run alignment for each row of run_grid
  if(verbose > 0){
    cat('\taln_msa_to_pdb: aligning msa(s) to pdb ')
  }

  aln_sets = list()
  for(i in 1:nrow(run_grid)){
    pep = msa_info[[run_grid$msa[i]]]$pep
    seq = pdb_info$seq_set[run_grid$chain[i]]

    if(any(is.na(seq))){
      stop(paste0('No sequence found for chain ', chain, '. Check PDB input.'))
    }

    aln_sets[[i]] = .align_sequences(c(pep, seq))

    # add msa and pdb info here (pdb_id and msa_id can pull from wrapper level -- as module it just is pdb1, and msa1-n)
    aln_sets[[i]]$aln_mat = cbind(aln_sets[[i]]$aln_mat, msa = run_grid$msa[i], pdb = run_grid$pdb[i])

    names(aln_sets)[i] = paste0(run_grid$msa[i], '_', run_grid$pdb[i], '_', run_grid$chain[i])
  }

  # could stack same chain now -- or later #

  # step 2: map alignment to positions ----
  if(verbose > 0){
    cat('\taln_msa_to_pdb: mapping alignment to codon and pdb positions\n')
  }

  # grab residue df #
  residue_df = pdb_info$residue_df
  pos_sets = list()
  for(i in 1:nrow(run_grid)){
    pos_sets[[i]] = .map_aln_to_positions(aln_sets[[i]]$aln_mat, residue_df, chain = run_grid$chain[i])
    names(pos_sets)[i] = paste0(run_grid$msa[i], '_', run_grid$pdb[i], '_', run_grid$chain[i])
  }

  } else {
    # use user adjusted pos_mat #
    residue_df = pdb_info$residue_df
    pos_sets = user_aln$pos_mat
  }

  # calculate coverage on pos_mat #
  coverage = .calculate_coverage(pos_sets)

  # step 2.5: create large rbind() list -- every residue_id has row ----
  aln_table = do.call(rbind, pos_sets)
  aln_table = as.data.frame(aln_table, stringsAsFactors = FALSE)

  # at this point need to recover "interfaces" from residue_df and place them as residue_id #
  # any interfaces or drug binding sites (can be coded as interface) #
  ro = grep('^interface', residue_df$residue_id)
  if(length(ro)>0){
    interface_table = data.frame(
      residue_id = residue_df$residue_id[ro],
      ref_aa = NA,
      pdb_aa = NA,
      msa = NA,
      pdb = run_grid$pdb[1],
      codon = NA
    )

    aln_table = rbind(aln_table, interface_table)
  }


  # step 3: map patches to codons ----
  if(verbose > 0){
    cat('\taln_msa_to_pdb: converting pdb patches to codon\n')
  }

  codon_patches = .map_patches_to_codons(aln_table, residue_df,
                                         drop_gap_map = TRUE)

  # STEP 4: add MSA subset ids ----
  # here tied to residue_id #
  codon_patches$msa_subset_id = paste0(codon_patches$residue_id, '_', codon_patches$pdb)
  codon_patches$msa_subset_id[is.na(codon_patches$codon_patch)] = NA

  # STEP 5: fix gap map and dedup codons ----
  if(fix_gap_map_and_dedup_codons){
    codon_patches = .fix_gap_map_and_dedup_codons(residue_df = codon_patches,
                                                  patch_mode = patch_mode,
                                                  max_patch = max_patch,
                                                  dist_cutoff = dist_cutoff,
                                                  only_exposed_in_patch = only_exposed_in_patch,
                                                  pdb_info = pdb_info)
  }

  #table(codon_patches$codon_len)
  #table(codon_patches$unique_codon)

  # step 4: grab subsets of MSA ----
  msa_subsets = NULL
  if(subset_msa){
    msa_subsets = .extract_msa_subsets(msa_info, codon_patches)
  }

  # step 5: gather coverage data for individual aliments and return ----

  # also fixed exposure to FALSE if NA
  #codon_patches$exposed[is.na(codon_patches$exposed)] = FALSE

  # return data #
  results = list(
    coverage = coverage,
    aln_df = codon_patches,
    pos_mat = pos_sets,
    msa_subsets = msa_subsets)

  class(results) = 'evo3D_aln_info'

  return(results)

}



# adjust_aln() ----

#' Adjust alignment rows by shifting residue assignments
#'
#' Utility to manually correct alignment mappings by moving residue identifiers
#' and aligned structure characters from one row to another. This is primarily
#' used to fix small alignment offsets introduced during automated alignment or
#' post-processing.
#'
#' If \code{pos_mat} is a list of matrices (e.g., for multimers), the adjustment
#' is applied independently to each element.
#'
#' @param pos_mat A character matrix containing at least columns
#'   \code{"residue_id"} and \code{"pdb_aa"}, or a list of such matrices.
#' @param from Integer vector of source row indices.
#' @param to Integer vector of destination row indices. Each destination row
#'   must contain a gap (\code{fill_char}) in \code{residue_id}.
#' @param fill_char Character used to represent gaps (default \code{"-"}).
#'
#' @return A modified alignment matrix (or list of matrices) with rows adjusted.
#'
#' @export

adjust_aln = function(pos_mat, from = NULL, to = NULL, fill_char = "-") {

  # if multimers
  if (is.list(pos_mat)) {
    return(lapply(
      pos_mat,
      adjust_aln,
      from = from,
      to   = to,
      fill_char = fill_char
    ))
  }

  # --- single-matrix case below ---

  if (length(from) != length(to)) {
    stop("arguments from and to need the same length")
  }

  for (i in seq_along(from)) {
    if (pos_mat[to[i], "residue_id"] != "-") {
      msg = paste0(
        "Cannot replace row ", to[i],
        " because it is not a gap: ",
        pos_mat[to[i], "residue_id"]
      )
      stop(msg)
    }

    # copy
    pos_mat[to[i], "residue_id"] = pos_mat[from[i], "residue_id"]
    pos_mat[to[i], "pdb_aa"]     = pos_mat[from[i], "pdb_aa"]

    # blank old
    pos_mat[from[i], "residue_id"] = fill_char
    pos_mat[from[i], "pdb_aa"]     = fill_char
  }

  pos_mat
}
