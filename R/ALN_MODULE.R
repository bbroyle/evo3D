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
#' @keywords internal
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
  aln = suppressMessages(msa::msa(sequences, method = 'ClustalOmega', type = 'protein', order = 'input', substitutionMatrix = 'BLOSUM65'))

  aln_set = as.character(aln@unmasked)
  aln_chars = strsplit(aln_set, '')

  # which position are identical
  mismatch = which(aln_chars[[1]] != aln_chars[[2]])

  # only keep mismatch if they are aa to aa (not - or X)
  mismatch = mismatch[!aln_chars[[1]][mismatch] %in% c('-', 'X')]
  mismatch = mismatch[!aln_chars[[2]][mismatch] %in% c('-')]

  # unpack alignment into matrix
  aln_mat = do.call(cbind, aln_chars)
  colnames(aln_mat) = c('ref_aa', 'pdb_aa')

  # calculate coverage and plot #
  coverage = .calculate_coverage(t(aln_mat), mismatch)

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
#' @keywords internal

.map_aln_to_positions = function(aln_mat, residue_df, chain = NA){

  # aln_mat is from .align_sequences() aln_mat$aln_mat
  # residue_df is from pdb_to_patch()

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
  wanted_order <- c(
    "codon_id", "codon", "msa", "pdb", "residue_id", "ref_aa", "pdb_aa",
    "codon_patch", "codon_len", "unique_codon", "gap_map_count",
    "exposed", "max_dist", "patch", "patch_len"
  )

  # columns in data that match your desired order
  present <- intersect(wanted_order, names(aln_table))

  # everything else (not in your desired order)
  extra <- setdiff(names(aln_table), wanted_order)

  # final reordering
  aln_table <- aln_table[ , c(present, extra)]


  # last clean up if codon_patch is '' set all codon info to NA #
  # possible they are rebuilt - so save them #
  #aln_table$codon_patch[aln_table$codon_patch == ''] = NA
  #aln_table$codon_len[aln_table$codon_patch == ''] = NA
  #aln_table$unique_codon[aln_table$codon_patch == ''] = NA

  return(aln_table)
}

# .extract_msa_subsets() ----

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
#' @keywords internal
.extract_msa_subsets_single = function(msa, codon_patches){
  # works across MSA's #
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

.extract_msa_subsets = function(msa_info, codon_patches, use_sample_names = TRUE){
  # works across MSA's #
  # schedules extract_msa_subsets_single #

  # still need to incorporate use_sample_names #

  # two options --
  # 1. pull apart patches into msa1, msa2, msaX columns and schedule blocks to extract_msa_subsets #
  # 2. send row by row -- still pulling apart #
  # -- doing option 2

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
      msa_list[[j]] = .extract_msa_subsets_single(msa_info[[msa_name]]$msa_mat, data.frame(msa_subset_id = id, codon_patch = paste(codons2, collapse = '+')))[[1]]
    }

    # right here -- there needs to be a checking for sample_names
    msa_subset = do.call(cbind, msa_list)

    msa_subsets[[i]] = msa_subset
    names(msa_subsets)[i] = id

  }

  # clean up set to only the ones with pulls #
  keep = !vapply(msa_subsets, is.null, logical(1))
  msa_subsets = msa_subsets[keep]

  # just return
  return(msa_subsets)
}


# .fix_gap_map_and_dedup_codons ----

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
  ro <- grep("^interface", residue_df$residue_id)
  interface_df <- NULL
  if (length(ro)) {
    interface_df <- residue_df[ro, ]
    residue_df <- residue_df[-ro, ]
  }

  # Step 1: rebuild patches depending on mode
  if (patch_mode == "codon") {
    ro <- which(!is.na(residue_df$codon_patch) & residue_df$unique_codon < max_patch)

    if(length(ro)){
      # grab valid replacements #
      valid <- residue_df[residue_df$codon != "-" &
                            residue_df$pdb_aa != "-" &
                            (if (only_exposed_in_patch) residue_df$exposed else TRUE),
                          c("codon_id", "residue_id"), drop = FALSE]

      # loop through each row and rebuild with valid sets #
      if(nrow(valid)){

        # trim dist_mat
        d <- pdb_info$residue_dist
        d <- d[, valid$residue_id, drop = FALSE]

        # make mapping
        res2cod <- setNames(valid$codon_id, valid$residue_id)

        # loop through ro #
        for (i in ro) {
          resi <- residue_df$residue_id[i]
          if (!resi %in% rownames(d)) next

          d2 <- sort(d[resi, , drop = TRUE])
          if (!is.na(dist_cutoff)) d2 <- d2[d2 <= dist_cutoff]
          if (!length(d2)) next

          cod_seq <- res2cod[names(d2)]
          uniq_cod <- cod_seq[!duplicated(cod_seq)]
          set <- uniq_cod[seq_len(min(max_patch, length(uniq_cod)))]

          residue_df$patch_len[i] = length(set)
          residue_df$patch[i]       <- paste(names(set), collapse = "+")
          residue_df$codon_patch[i] <- paste(set, collapse = "+")
          residue_df$codon_len[i]   <- length(set)
          residue_df$unique_codon[i]<- length(unique(set)) # codon_len and unique_codon should be the same #
        }
      }


      }
    } else if (patch_mode == "residue") {
      ro <- which(!is.na(residue_df$codon_patch) & residue_df$codon_len < max_patch)

      if(length(ro)){
        # grab valid replacements #
        valid <- residue_df[residue_df$codon != "-" &
                              residue_df$pdb_aa != "-" &
                              (if (only_exposed_in_patch) residue_df$exposed else TRUE),
                            c("codon_id", "residue_id"), drop = FALSE]

        # loop through each row and rebuild with valid sets #
        if(nrow(valid)){

          # trim dist_mat
          d <- pdb_info$residue_dist
          d <- d[, valid$residue_id, drop = FALSE]

          # make mapping
          res2cod <- setNames(valid$codon_id, valid$residue_id)

          # loop through ro #
          for (i in ro) {
            resi <- residue_df$residue_id[i]
            if (!resi %in% rownames(d)) next

            d2 <- sort(d[resi, , drop = TRUE])
            if (!is.na(dist_cutoff)) d2 <- d2[d2 <= dist_cutoff]
            if (!length(d2)) next

            take_k <- min(length(d2), max_patch)
            set <- names(d2)[seq_len(take_k)]

            residue_df$patch_len[i]   = length(set)
            residue_df$patch[i]       <- paste(set, collapse = "+")
            residue_df$codon_patch[i] <- paste(res2cod[set], collapse = "+")
            residue_df$codon_len[i]   <- length(res2cod[set])
            residue_df$unique_codon[i]<- length(unique(res2cod[set]))
          }
        }
      }
    }

  # Step 2: reattach interfaces
  if (!is.null(interface_df)) {
    residue_df <- rbind(residue_df, interface_df)
  }

  return(residue_df)
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
#' @param drop_unused_residues Logical; if \code{TRUE}, drops PDB residues not linked to codon-level info.
#'
#' @return A list with:
#' \describe{
#'   \item{aln_coverage}{Named list of alignment coverage ranges.}
#'   \item{aln_df}{Data frame mapping residues to codons, amino acids, and patch identifiers.}
#'   \item{msa_subsets}{Named list of nucleotide MSA matrices, one per patch.}
#' }
#' @export
aln_msa_to_pdb = function(msa_info, pdb_info, chain = 'auto',
                           subset_msa = TRUE,
                           verbose = 1,
                           run_grid = NULL,
                           drop_gap_map = TRUE, # aln will drop patches centered on gap map residues
                           fix_gap_map_and_dedup_codons = TRUE, # aln will drop gap map residues from patch membership
                           patch_mode, max_patch, dist_cutoff, only_exposed_in_patch,
                           merge_type = 'hold', merge_exposure = 'hold'){

  # msa_info ~ must be list object from msa_to_ref() # ~ can be list of msa_infos
  # pdb_info ~ must be list object from ONE pdb_to_patch() #

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
    aln_sets[[i]]$aln_mat <- cbind(aln_sets[[i]]$aln_mat, msa = run_grid$msa[i], pdb = run_grid$pdb[i])

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

  table(codon_patches$codon_len)
  table(codon_patches$unique_codon)

  # step 4: grab subsets of MSA ----
  msa_subsets = NULL
  if(subset_msa){
    msa_subsets = .extract_msa_subsets(msa_info, codon_patches)
  }

  # step 5: gather coverage data for individual aliments and return ----
  aln_sets <- lapply(aln_sets, function(x) {
    x$aln_mat <- NULL
    x
  })

  # also fixed exposure to FALSE if NA
  #codon_patches$exposed[is.na(codon_patches$exposed)] = FALSE

  # return data #
  return(list(
    coverage = aln_sets,
    aln_df = codon_patches,
    msa_subsets = msa_subsets)
  )

}


