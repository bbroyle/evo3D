# --------------------------------------------------------------- #
# PDB MODULE
# Utilities for processing PDB and mmCIF structure data and returning
# standardized pdb_info objects for downstream evo3D analysis.
#
# Internal helper functions (.functions) are wrapped by the
# user-facing pdb_to_patch() interface.
#
# NOTES:
# 1. Hydrogens and HETATM records are removed by default.
# 2. Both RSA and SASA are supported; SASA is recommended where available.
# 3. Occlusion and interface chains are supported for accessibility and
#    contact-based analyses.
#
# Contact: bbroyle@purdue.edu
# --------------------------------------------------------------- #

# .split_chain() ----
#' Split chain identifiers into single-character IDs
#'
#' Internal helper to normalize the \code{chain} argument. If \code{chain} is
#' \code{NA}, it is returned unchanged. Otherwise, each element of \code{chain}
#' is split into single-character chain identifiers and flattened into a
#' character vector.
#'
#' @param chain Character vector of chain identifiers, or \code{NA}.
#'
#' @return A character vector of single-character chain IDs, or \code{NA}.
#'
#' @keywords internal
#' @noRd
.split_chain = function(chain) {
  if (length(chain) == 1L && is.na(chain)) {
    return(chain)
  }
  unlist(strsplit(chain, "", fixed = TRUE), use.names = FALSE)
}

# .plot_chain_map () ----

#' Plot PCA map of PDB chains (internal)
#'
#' Internal utility to visualize the spatial layout of PDB chains using a
#' 2D PCA projection of Cα coordinates.
#'
#' @param pdb PDB input (file path or bio3d PDB object).
#' @param chain Optional chain identifier(s); \code{NA} uses all chains.
#'
#' @keywords internal
.plot_chain_map = function(pdb, chain = NA) {

  # local helper to compute coords
  pca_pdb = function(pdb, chain) {
    pdb = .standardize_pdb_input(pdb)
    ca = pdb$atom[pdb$atom$elety == "CA", ]

    if (!(length(chain) == 1L && is.na(chain))) {
      ca = ca[ca$chain %in% chain, ]
    }

    coords = ca[, c("x","y","z")]
    pc = stats::prcomp(coords, center = TRUE, scale. = FALSE)
    mds_coords = pc$x[, 1:2]

    data.frame(
      residue_id = ca$residue_id,
      chain = as.character(ca$chain),
      x = mds_coords[, 1],
      y = mds_coords[, 2]
    )

  }

  # if chain present split into mulit_chain #
  chain = .split_chain(chain)

  # compute coords
  plot_df = pca_pdb(pdb, chain)

  # median coords for labels
  label_df = stats::aggregate(cbind(x, y) ~ chain, data = plot_df, FUN = stats::median)

  # make sure chain is character
  chain_ids = unique(plot_df$chain)

  if (length(chain_ids) == 0) {
    stop("No chains found for plotting")
  }

  # if only one chain, still assign a color
  cols = if (length(chain_ids) == 1) "black" else grDevices::rainbow(length(chain_ids))
  col_map = stats::setNames(cols, chain_ids)

  # plot points
  plot(
    plot_df$x, plot_df$y,
    col = col_map[plot_df$chain],
    pch = 19, cex = 0.7,
    xlab = "", ylab = "",
    axes = FALSE
  )

  # add labels at medians
  graphics::text(label_df$x, label_df$y, labels = label_df$chain, col = "black", font = 2, cex = 1)

  # add legend if more than one chain
  if (length(chain_ids) > 1) {
    graphics::legend("bottomleft", legend = chain_ids, col = cols, pch = 19, bty = "n", cex = 1)
  }



}


# .standardize_pdb_input() ----

#' Standardize PDB/mmCIF input
#'
#' Reads a PDB or mmCIF file (or accepts an existing `bio3d` PDB object) and
#' returns a cleaned PDB object ready for downstream evo3D analysis. This
#' standardization step enforces evo3D assumptions about chain identifiers and
#' internal delimiters, and appends a unique per-atom residue identifier.
#'
#' @param pdb PDB input: a file path, or an object returned by
#'   \code{bio3d::read.pdb()} or \code{bio3d::read.cif()}.
#' @param force_file_type Optional character to override auto-detection; one of
#'   \code{"cif"} or \code{"pdb"}. Default \code{NULL} uses the file extension.
#'
#' @details
#' The following preprocessing steps are applied:
#' \itemize{
#'   \item Insert codes with \code{NA} are replaced by \code{""} (empty string).
#'   \item Chain identifiers must be single characters; multi-character chain IDs
#'     result in an error (current evo3D limitation).
#'   \item The characters \code{+}, \code{_}, \code{(}, and \code{)} are reserved
#'     internally by evo3D. If present in \code{atom$resno}, \code{atom$chain},
#'     or \code{atom$insert}, the function stops with an error.
#'   \item A unique residue identifier is added as \code{atom$residue_id},
#'     constructed as \code{"resno_chain_insert"}.
#' }
#'
#' @return A \code{"pdb"} object (from \pkg{bio3d}) with an additional
#'   \code{atom$residue_id} column.
#'
#' @export

.standardize_pdb_input = function(pdb, force_file_type = NULL){

  # validate pdb (file or pdb object)
  if (length(pdb) == 1L && is.character(pdb)) {

    if (!file.exists(pdb)) stop("PDB file does not exist: ", pdb)

    if (!is.null(force_file_type)) {
      if (force_file_type == "cif") pdb = bio3d::read.cif(pdb, multi = FALSE, rm.insert = FALSE, rm.alt = TRUE)
      else if (force_file_type == "pdb") pdb = bio3d::read.pdb(pdb, multi = FALSE, rm.insert = FALSE, rm.alt = TRUE)
      else stop("Invalid force_file_type. Use 'cif' or 'pdb'.")
    } else {
      if (grepl("cif$", pdb, ignore.case = TRUE)) pdb = bio3d::read.cif(pdb, multi = FALSE, rm.insert = FALSE, rm.alt = TRUE)
      else pdb = bio3d::read.pdb(pdb, multi = FALSE, rm.insert = FALSE, rm.alt = TRUE)
    }

  } else if (!inherits(pdb, "pdb")) {
    stop("pdb provided is not a recognized format.\n Please provide a file path or an object from bio3d::read.pdb()/read.cif().")
  }

  # check for required info #
  if (is.null(pdb$atom) || !is.data.frame(pdb$atom)) {
    stop("Input contains no atom table (pdb$atom).")
  }

  # fix insert NA to '' #
  pdb$atom$insert = as.character(pdb$atom$insert)
  pdb$atom$insert[is.na(pdb$atom$insert)] = ""

  # fix chain NA and check for multi char #
  pdb$atom$chain = as.character(pdb$atom$chain)
  pdb$atom$chain[is.na(pdb$atom$chain)] = ""

  # if multi char chains - stop #
  # internally multi char are used chains signify multimer chains #
  lens = nchar(pdb$atom$chain)
  if (any(lens != 1L)) {
    bad = unique(pdb$atom$chain[lens != 1L])
    stop(
      "Invalid chain IDs detected (must be single character):\n ",
      paste(bad, collapse = ", "),
      "\n evo3D currently supports only single-character chain identifiers."
    )
  }

  # characters reserved internally
  reserved_chars = c("\\+", "_", "\\(", "\\)")

  for (pat in reserved_chars) {
    in_resno = any(grepl(pat, pdb$atom$resno))
    in_chain = any(grepl(pat, pdb$atom$chain))
    in_insert = any(grepl(pat, pdb$atom$insert))

    if (any(c(in_resno, in_chain, in_insert))) {
      stop(
        sprintf(
          "Reserved character '%s' present in pdb info (resno, chain, insert). evo3D reserves this character internally.",
          gsub("\\\\", "", pat)
        )
      )
    }
  }

  # add residue id information #
  pdb$atom$residue_id = paste0(pdb$atom$resno, '_', pdb$atom$chain, '_', pdb$atom$insert)

  # pass along pdb #
  return(pdb)
}

# .get_pdb_sequence() ----

#' Extract PDB chain sequences
#'
#' Retrieves amino acid sequences for one or more chains from a PDB structure.
#' This function is used internally by PDB alignment and chain-detection
#' utilities.
#'
#' @param pdb A PDB object from \pkg{bio3d}, or a file path passed to
#'   \code{.standardize_pdb_input()} if \code{in_module = FALSE}.
#' @param chain Optional character vector of chain identifiers. If \code{NA},
#'   sequences for all protein chains are extracted.
#' @param in_module Logical. If \code{TRUE}, assumes \code{pdb} has already been
#'   standardized (default \code{FALSE}).
#'
#' @return A named character vector of amino acid sequences, one per chain.
#'
#' @keywords internal

.get_pdb_sequence = function(pdb, chain = NA, in_module = FALSE){

  #if running in module (.standardize_pdb_input() is already run) #
  if(!in_module){
    pdb = .standardize_pdb_input(pdb)
    chain = .split_chain(chain)
  }

  if(length(chain) == 1L && is.na(chain)) {
    chain = unique(pdb$atom$chain)
  }

  # get pdb sequences with names as chains
  aa_seq = sapply(chain, function(x){
    paste0(
      bio3d::pdbseq(bio3d::trim.pdb(pdb, chain = x, 'protein')),
      collapse = '')
  }, USE.NAMES = TRUE)

  return(aa_seq)
}

# .calculate_residue_distance() ----

#' Residue-wise distance matrix
#'
#' Computes a pairwise Euclidean distance matrix between residues of a PDB
#' structure. Distances are calculated from atomic coordinates and aggregated
#' at the residue level using one of several supported methods.
#'
#' This function is used internally by patch-definition, exposure, and
#' structure-based analysis routines in evo3D.
#'
#' @param pdb A PDB object from \pkg{bio3d}, or a file path passed to
#'   \code{.standardize_pdb_input()} if \code{in_module = FALSE}.
#' @param chain Optional character vector of chain identifiers. If \code{NA},
#'   all protein chains are retained.
#' @param distance_method Character string specifying the residue distance
#'   definition. One of:
#'   \itemize{
#'     \item \code{"ca"}: Cα–Cα distances.
#'     \item \code{"all"}: minimum inter-atomic distance between residues.
#'     \item \code{"centroid"}: Euclidean distance between residue centroids.
#'   }
#' @param in_module Logical. If \code{TRUE}, assumes \code{pdb} has already been
#'   standardized (default \code{FALSE}).
#'
#' @return A symmetric numeric matrix of residue–residue distances, with
#'   residue identifiers as row and column names.
#'
#' @keywords internal

.calculate_residue_distance = function(pdb, chain = NA, distance_method = 'ca',
                                       in_module = FALSE){

  # distance_method can be 'ca', 'all', 'centroid'
  # hydrogens always removed from calculation #
  if (!distance_method %in% c("ca","all","centroid")) {
    stop("distance_method must be one of: 'ca', 'all', 'centroid'")
  }

  #if running in module (.standardize_pdb_input() is already run) #
  if(!in_module){
    pdb = .standardize_pdb_input(pdb)
    chain = .split_chain(chain)
  }

  # trim chains
  if(!(length(chain) == 1L && is.na(chain))) {
    pdb = bio3d::trim.pdb(pdb, chain = chain)
  }

  # pdb should be protein only and no H (what about glycan distance?) #
  # if chain was NA - compute distances between all protein chains #
  pdb = bio3d::trim.pdb(pdb, 'protein')
  pdb = bio3d::trim.pdb(pdb, 'noh')

  # set distance mat
  if(distance_method == 'ca'){

    pdb = bio3d::trim.pdb(pdb, 'calpha')
    # get atom distance matrix (could output symmetrical)
    res_dist = bio3d::dm.xyz(pdb$xyz, grpby = pdb$atom$residue_id, mask.lower = FALSE)

    # set row and column names
    colnames(res_dist) = rownames(res_dist) = unique(pdb$atom$residue_id)

    # set diag to 0
    diag(res_dist) = 0

  } else if (distance_method == 'centroid'){
    uni_res = unique(pdb$atom$residue_id)
    centroids = data.frame(residue_id = uni_res,
                           x = NA_real_, y = NA_real_, z = NA_real_,
                           stringsAsFactors = FALSE)

    for(i in 1:nrow(centroids)){
      ro = which(pdb$atom$residue_id == centroids$residue_id[i])
      centroids$x[i] = mean(pdb$atom$x[ro], na.rm = TRUE)
      centroids$y[i] = mean(pdb$atom$y[ro], na.rm = TRUE)
      centroids$z[i] = mean(pdb$atom$z[ro], na.rm = TRUE)
    }

    # this one quick enough to supply to dist() - other methods go through bio3d::dm.xyz()
    res_dist = stats::dist(matrix(as.numeric(unlist(centroids[,2:4])), nrow = nrow(centroids)))
    res_dist = as.matrix(res_dist)

    # set row and column names
    colnames(res_dist) = rownames(res_dist) = centroids$residue_id

  } else {
    # all atom method
    # get atom distance matrix (could output symmetrical)
    res_dist = bio3d::dm.xyz(pdb$xyz, grpby = pdb$atom$residue_id, mask.lower = FALSE)

    # set row and column names
    colnames(res_dist) = rownames(res_dist) = unique(pdb$atom$residue_id)

    # set diag to 0
    diag(res_dist) = 0
  }

  # return residue wise distance matrix
  return(res_dist)

}

# .calculate_accessibility() ----

#' Calculate solvent accessibility
#'
#' Computes per-residue solvent accessible surface area (SASA) using a
#' DSSP-style algorithm reimplemented in R and C++.
#'
#' Hydrogen atoms are removed prior to calculation. Residues are temporarily
#' renumbered sequentially to ensure contiguous indexing during SASA
#' computation, then mapped back to their original PDB residue numbers, chain
#' identifiers, and insert codes. Accessibility is reported as both absolute
#' SASA (Å²) and normalized relative solvent accessibility (RSA).
#'
#' This function is used internally by surface-exposure and patch-definition
#' routines, but is also exposed for direct use.
#'
#' @param pdb A PDB object from \pkg{bio3d}, or a file path passed to
#'   \code{.standardize_pdb_input()} if \code{in_module = FALSE}.
#' @param chain Optional character vector of chain identifiers to retain.
#'   Default \code{NA} retains all protein chains.
#' @param method Character string specifying the normalization scheme; one of
#'   \code{"rose"}, \code{"miller"}, \code{"theoretical_tien"}, or
#'   \code{"empirical_tien"} (default \code{"rose"}).
#' @param drop_incomplete Logical. If \code{TRUE} (default), residues missing
#'   backbone atoms (N, CA, C, O) are excluded, matching DSSP behavior.
#' @param in_module Logical. If \code{TRUE}, assumes \code{pdb} has already been
#'   standardized (default \code{FALSE}).
#'
#' @return A data frame with one row per residue, containing:
#' \itemize{
#'   \item \code{residue_index}: sequential residue index used internally
#'   \item \code{aa}: one-letter amino acid code
#'   \item \code{orig_resno}, \code{orig_chain}, \code{orig_insert}: original PDB identifiers
#'   \item \code{sasa}: absolute solvent accessible surface area (Å²)
#'   \item \code{rsa}: relative solvent accessibility (0–1, capped at 1)
#'   \item \code{residue_id}: concatenated unique identifier
#' }
#'
#' @export

.calculate_accessibility = function(pdb, chain = NA, method = 'rose', drop_incomplete = TRUE, in_module = FALSE){
  # return residue wise solvent accessibility

  #if running in module (.standardize_pdb_input() is already run) #
  if(!in_module){
    pdb = .standardize_pdb_input(pdb)
    chain = .split_chain(chain)
  }

  # check method is valid #
  if (!method %in% c("rose","miller","theoretical_tien","empirical_tien")) {
    stop("method must be one of: rose, miller, theoretical_tien, empirical_tien")
  }

  # trim chains
  if(!(length(chain) == 1L && is.na(chain))) {
    pdb = bio3d::trim.pdb(pdb, chain = chain)
  }

  # otherwise using all protein chains #

  # MAX ACC VALUE TABLE ------------------------------------------- #
  max_acc = data.frame(
    residue = c('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
                'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'),
    theoretical_tien = c(129.0, 274.0, 195.0, 193.0, 167.0, 223.0, 225.0, 104.0, 224.0, 197.0,
                         201.0, 236.0, 224.0, 240.0, 159.0, 155.0, 172.0, 285.0, 263.0, 174.0),
    empirical_tien = c(121.0, 265.0, 187.0, 187.0, 148.0, 214.0, 214.0, 97.0, 216.0, 195.0,
                       191.0, 230.0, 203.0, 228.0, 154.0, 143.0, 163.0, 264.0, 255.0, 165.0),
    miller = c(113.0, 241.0, 158.0, 151.0, 140.0, 183.0, 189.0, 85.0, 194.0, 182.0,
               180.0, 211.0, 204.0, 218.0, 143.0, 122.0, 146.0, 259.0, 229.0, 160.0),
    rose = c(118.1, 256.0, 165.5, 158.7, 146.1, 186.2, 193.2, 88.1, 202.5, 181.0,
             193.1, 225.8, 203.4, 222.8, 146.8, 129.8, 152.5, 266.3, 236.8, 164.5),
    stringsAsFactors = FALSE
  )

  # PREP PDB ------------------------------------------------------ #
  # Extract atom data
  atoms = pdb$atom

  # Drop hydrogen atoms and non-ATOM records.
  # Filter using both element symbol (elesy) and atom name (elety),
  # as element annotations can be unreliable in some PDB files.
  atoms = atoms[
    atoms$type == "ATOM" &
      atoms$elesy != "H" &
      !grepl("^H", atoms$elety),
  ]

  # Handle residue numbering (each aa needs unique)
  atoms$orig_resno = atoms$resno
  atoms$resno = 1
  atoms$insert[is.na(atoms$insert)] = ''
  atoms$chain = paste0(atoms$chain, '_', atoms$insert)

  # renumber sequentially (keep atoms from same residue together)
  prev = atoms[1,]
  for(i in 2:nrow(atoms)){
    if(atoms$orig_resno[i] == prev$orig_resno &&
       atoms$chain[i] == prev$chain &&
       atoms$insert[i] == prev$insert){
      atoms$resno[i] = prev$resno
    } else {
      next_res = prev$resno + 1
      atoms$resno[i] = next_res
      prev = atoms[i,]
    }
  }

  # drop incomplete residues (same as original DSSP)
  if(drop_incomplete){
    # function to check completeness
    complete_check = function(resno){
      types = atoms$elety[atoms$resno == resno]
      all(c('N', 'CA', 'C', 'O') %in% types)
    }

    # check for complete residues
    resno_list = unique(atoms$resno)
    incomplete = sapply(resno_list, function(x) !complete_check(x))

    if(any(incomplete)){
      atoms = atoms[!atoms$resno %in% resno_list[incomplete],]
    }
  }

  # SETTING DATASETS FOR DSSP - rewrite  --------------------------- #
  resno_list = unique(atoms$resno)
  residue_df = data.frame(
    residue_index = resno_list,
    chain_id = sapply(resno_list, function(i) atoms$chain[atoms$resno == i][1]),
    residue_name = sapply(resno_list, function(i) bio3d::aa321(atoms$resid[atoms$resno == i][1])),
    orig_resno = sapply(resno_list, function(i) atoms$orig_resno[atoms$resno == i][1]),
    stringsAsFactors = FALSE
  )

  atom_df = data.frame(
    residue_index = atoms$resno,
    atom_type = atoms$elety,
    x = atoms$x,
    y = atoms$y,
    z = atoms$z,
    stringsAsFactors = FALSE
  )

  # RUNNING DSSP REWRITE ------------------------------------------- #
  accessibility = calculateDSSPAccessibility(atom_df, residue_df)
  residue_df$sasa = round(accessibility, 3)

  # add rsa
  residue_df$max_acc = max_acc[match(residue_df$residue_name, max_acc$residue), method]
  residue_df$rsa = residue_df$sasa / residue_df$max_acc
  residue_df$rsa[residue_df$rsa > 1] = 1 # if rsa above 1 (set at 1)
  residue_df$rsa = round(residue_df$rsa, 2)

  # return with residue id
  residue_df$residue_id = paste0(residue_df$orig_resno, '_', residue_df$chain_id)

  # clean up a little bit #
  residue_df$aa = residue_df$residue_name
  residue_df$orig_chain = gsub('_.*', '', residue_df$chain_id)
  residue_df$orig_insert = gsub('.+_', '', residue_df$chain_id)

  # simplify output !! might want to have orig_resno, chain, insert seperate #
  residue_df = residue_df[, c('residue_index', 'aa', 'orig_resno', 'orig_chain', 'orig_insert', 'sasa', 'rsa', 'residue_id')]

  return(residue_df)
}


# .is_exposed() ----

#' Classify exposed residues
#'
#' Marks residues as exposed based on solvent accessibility thresholds.
#' Exposure can be defined using relative solvent accessibility (RSA),
#' absolute solvent accessible surface area (SASA), or both. When both
#' criteria are provided, their combination is controlled by
#' \code{use_rsa_sasa}.
#'
#' This function is called internally after
#' \code{.calculate_accessibility()} to append an \code{exposed} column
#' to the residue-level data frame.
#'
#' @param residue_df Data frame of residues, typically the output of
#'   \code{.calculate_accessibility()}, containing \code{rsa} and/or
#'   \code{sasa} columns.
#' @param rsa_cutoff Numeric threshold for RSA-based exposure classification.
#'   May be a single value (interpreted as \code{>=}) or a length-2 numeric
#'   range. Default \code{NA} disables RSA filtering.
#' @param sasa_cutoff Numeric threshold for SASA-based exposure classification.
#'   May be a single value (interpreted as \code{>=}) or a length-2 numeric
#'   range. Default \code{NA} disables SASA filtering.
#' @param use_rsa_sasa Character string, either \code{"and"} or \code{"or"},
#'   specifying how RSA and SASA criteria are combined when both are provided.
#'   Default \code{"and"}. If either cutoff is \code{NA}, this is
#'   automatically set to \code{"and"}.
#'
#' @return The input \code{residue_df} with an additional logical column
#'   \code{exposed}, indicating whether each residue meets the exposure
#'   criteria.
#'
#' @keywords internal

.is_exposed = function(residue_df, rsa_cutoff = NA,
                       sasa_cutoff = NA, use_rsa_sasa = 'and') {

  # local helper to handle rsa and sasa ranges #
  .in_range <- function(x, cutoff) {
    if (length(cutoff) == 1L) {
      x >= cutoff
    } else {
      x >= cutoff[1] & x <= cutoff[2]
    }
  }

  # checking input types #
  if (length(rsa_cutoff) > 2L || length(sasa_cutoff) > 2L) {
    stop("rsa and sasa cutoffs must be length 1 or 2")
  }

  if (length(rsa_cutoff) == 2L && anyNA(rsa_cutoff)) {
    stop("rsa_cutoff ranges cannot contain NA")
  }

  if (length(sasa_cutoff) == 2L && anyNA(sasa_cutoff)) {
    stop("sasa_cutoff ranges cannot contain NA")
  }

  # if either or both are NA, set to AND #
  if (anyNA(rsa_cutoff) || anyNA(sasa_cutoff)) {
    use_rsa_sasa = "and"
  }

  # build exposure under rsa (if available)
  if (!is.na(rsa_cutoff[1])) {
    pass_rsa = .in_range(residue_df$rsa, rsa_cutoff)
  } else {
    pass_rsa = TRUE
  }

  # build exposure under sasa (if available)
  if (!is.na(sasa_cutoff[1])) {
    pass_sasa = .in_range(residue_df$sasa, sasa_cutoff)
  } else {
    pass_sasa = TRUE
  }

  # combine
  if (use_rsa_sasa == "or") {
    residue_df$exposed = (pass_rsa | pass_sasa)
  } else {
    # AND works here for no filters (both true),
    # one filter (one is always true, other is test)
    # and both filters (explicitly looking for and)
    residue_df$exposed = (pass_rsa & pass_sasa)
  }

  return(residue_df)
}


# .identify_patches() ----

#' Identify surface patches
#'
#' Defines residue-level surface patches by growing spatial neighborhoods
#' around exposed seed residues. Patches are constructed using residue–residue
#' distances, with optional constraints on maximum distance and patch size.
#'
#' Seed residues are always those marked as \code{exposed} in
#' \code{residue_df}. If exposure thresholds upstream are permissive (e.g.,
#' cutoffs set to zero), all residues may qualify as seeds.
#'
#' @param dist_mat Numeric residue–residue distance matrix, with residue IDs
#'   as row and column names.
#' @param residue_df Data frame of residues containing exposure information
#'   (including \code{exposed}, \code{rsa}, and/or \code{sasa}).
#' @param dist_cutoff Numeric. Maximum distance (in Å) for including neighboring
#'   residues in a patch (default \code{15}).
#' @param max_patch Numeric. Optional maximum number of residues per patch
#'   (default \code{NA}, unlimited).
#' @param only_exposed_in_patch Logical. If \code{TRUE} (default), restricts
#'   patch members to residues marked as exposed.
#'
#' @return The input \code{residue_df} with additional columns:
#' \itemize{
#'   \item \code{patch}: concatenated residue IDs for each patch, joined by \code{+}
#'   \item \code{patch_len}: number of residues in the patch
#'   \item \code{max_dist}: maximum distance (Å) from the seed residue to any
#'     patch member
#' }
#'
#' @keywords internal

.identify_patches = function(dist_mat, residue_df,
                             dist_cutoff = 15, max_patch = NA,
                             only_exposed_in_patch = TRUE){

  # grab dist mat subset #
  res_ids = intersect(rownames(dist_mat), residue_df$residue_id)
  if (!length(res_ids)) stop("No overlap between dist_mat rownames and residue_df$residue_id")
  dist_mat = dist_mat[res_ids, res_ids, drop = FALSE]

  # first is to get seed residues #
  # only exposed residues are seeds #
  seeds = residue_df[!is.na(residue_df$exposed) & residue_df$exposed,]

  # If only exposed can be in patch -- shrink to seeds #
  if(only_exposed_in_patch){
    in_set = intersect(seeds$residue_id, rownames(dist_mat))
    dist_mat = dist_mat[in_set, in_set, drop = FALSE]
  }

  # get neighbors (could be dist_cutoff, or could be max_patch, or both)
  seeds$patch = NA_character_
  seeds$patch_len = NA_integer_
  seeds$max_dist = NA_real_

  for(i in seq_len(nrow(seeds))){
    center = seeds$residue_id[i]

    # get neighbors (sorted by distance) #
    neighbors = sort(dist_mat[center,])

    # filtered by dist_cut #
    if(!is.na(dist_cutoff)){
      neighbors = neighbors[neighbors <= dist_cutoff]
    }

    # filtered by max_patch #
    if(!is.na(max_patch)){
      neighbors = neighbors[seq_len(min(length(neighbors), max_patch))]
    }

    # add to seeds
    if (length(neighbors)) {
      seeds$patch[i] = paste0(names(neighbors), collapse = "+")
      seeds$patch_len[i] = length(neighbors)
      seeds$max_dist[i] = round(max(neighbors), 2)
    }

  }

  # add back to residue_df
  idx = match(seeds$residue_id, residue_df$residue_id)
  ok  = !is.na(idx)

  residue_df$patch = NA_character_
  residue_df$patch_len = NA_integer_
  residue_df$max_dist = NA_real_

  residue_df$patch[idx] = seeds$patch
  residue_df$patch_len[idx] = seeds$patch_len
  residue_df$max_dist[idx] = seeds$max_dist

  return(residue_df)
}

# .identify_interfaces() ----

#' Identify interface contacts
#'
#' Identifies residue contacts between two chain sets in a structure using
#' atom–atom distance criteria (via \code{bio3d::binding.site()}).
#'
#' Chains are first trimmed to protein atoms and hydrogens are removed. Contacting
#' residues from \code{chain} are returned as a \code{+}-concatenated string of
#' residue identifiers.
#'
#' @param pdb A \pkg{bio3d} \code{"pdb"} object, or a file path if
#'   \code{in_module = FALSE} (standardized via \code{.standardize_pdb_input()}).
#' @param chain Chain identifier(s) defining the primary set of residues.
#'   May be a vector and/or a compact string (e.g. \code{"ABC"}).
#' @param interface_chain Chain identifier(s) defining the interacting set.
#'   May be a vector and/or a compact string. If \code{NA}, returns
#'   \code{name = NA} and \code{interf = NA}.
#' @param dist_cutoff Numeric. Maximum distance (Å) for a contact (default \code{5}).
#' @param in_module Logical. If \code{TRUE}, assumes \code{pdb} and chain inputs
#'   have already been standardized/split (default \code{FALSE}).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{name}: interface label of the form
#'     \code{"interface_[chain]_[interface_chain]"}
#'   \item \code{interf}: \code{+}-concatenated residue identifiers for contacting residues
#' }
#'
#' @keywords internal

.identify_interface = function(pdb, chain = NA, interface_chain = NA,
                               dist_cutoff = 5, in_module = FALSE){

  if (length(interface_chain) == 1L && is.na(interface_chain)){
    message('No interface chain given')
    return(list(name = NA_character_, interf = NA_character_))
  }

  if(!in_module){
    pdb = .standardize_pdb_input(pdb)
    chain = .split_chain(chain)
    interface_chain = .split_chain(interface_chain)
  }

  # remove H and HETATM
  pdb = bio3d::trim.pdb(pdb, 'protein')
  pdb = bio3d::trim.pdb(pdb, 'noh')

  # Grab chains of interest #
  pdb1 = bio3d::trim.pdb(pdb, chain = c(chain))
  pdb2 = bio3d::trim.pdb(pdb, chain = c(interface_chain))

  # add residue_id to chain so i can recover
  pdb1$atom$chain = pdb1$atom$residue_id

  dist_mat = bio3d::binding.site(pdb1, pdb2, cutoff = dist_cutoff)

  interf = dist_mat$resnames
  interf = gsub('.+\\(', '', interf)
  interf = gsub('\\)', '', interf)
  interf = paste0(interf, collapse = '+')

  # return interface window #
  # format interf_chains_interacting #

  interf_name = paste0('interface_', paste0(chain, collapse = ''), '_', paste0(interface_chain, collapse = ''))

  return(list(
    name = interf_name,
    interf = interf
  ))
}

# pdb_to_patch() ----

#' Extract surface patches from a PDB structure
#'
#' High-level wrapper that computes chain sequences, residue–residue distances,
#' solvent accessibility, exposure status, and surface patches from a protein
#' structure.
#'
#' @param pdb PDB input: file path, \code{bio3d} object, or standardized PDB.
#' @param chain Chain identifier(s) to analyze; \code{NA} uses all chains.
#' @param interface_chain Optional chain identifier(s) defining interface contacts.
#' @param occlusion_chain Optional chain identifier(s) used to occlude the surface
#'   during accessibility calculations.
#' @param distance_method Residue distance metric; one of \code{"all"}, \code{"ca"},
#'   \code{"backbone"}, or \code{"sidechain"} (default \code{"all"}).
#' @param drop_incomplete_residue Logical; drop residues missing backbone atoms
#'   (default \code{TRUE}).
#' @param rsa_method Method for RSA normalization; one of \code{"rose"},
#'   \code{"miller"}, \code{"theoretical_tien"}, or \code{"empirical_tien"}.
#' @param dist_cutoff Numeric; maximum distance (Å) for defining patch neighbors
#'   (default 15).
#' @param rsa_cutoff Numeric; minimum RSA for defining seed residues (default 0.1).
#' @param sasa_cutoff Optional numeric; minimum SASA for defining seed residues.
#' @param only_exposed_in_patch Logical; if \code{TRUE}, restrict patch members to
#'   exposed residues only (default \code{TRUE}).
#' @param use_rsa_sasa Logical operator for combining RSA and SASA cutoffs; one of
#'   \code{"and"} or \code{"or"} (default \code{"and"}).
#' @param max_patch Optional integer; maximum number of neighbors per patch.
#' @param interface_dist_cutoff Numeric; maximum distance (Å) for defining interface
#'   residue contacts (default 5).
#' @param verbose Integer; verbosity level (0 = silent, >0 prints progress).
#' @param detail_level Integer controlling return content:
#'   \itemize{
#'     \item \code{0}: minimal output
#'     \item \code{1}: include chain identifiers
#'     \item \code{2}: include residue distance matrix
#'   }
#' @param force_file_type Optional; override automatic detection of PDB vs CIF input.
#' @param in_wrapper Logical; internal flag indicating whether the function is
#'   being called from a higher-level wrapper. When \code{TRUE}, certain input
#'   checks or messages may be suppressed.
#' @param patch_mode Placeholder argument; not implemented.
#'
#' @return A list with class \code{"evo3D_pdb_info"} containing:
#' \itemize{
#'   \item \code{pdb}: PDB object (or \code{NULL} when cached externally).
#'   \item \code{chain}: Chain identifiers if \code{detail_level > 0}.
#'   \item \code{seq_set}: Amino acid sequences per chain.
#'   \item \code{residue_dist}: Residue–residue distance matrix if
#'     \code{detail_level > 1}.
#'   \item \code{residue_df}: Data frame of residues with RSA/SASA values, exposure
#'     status, and patch membership.
#' }
#'
#' @export

pdb_to_patch = function(pdb, chain = NA, interface_chain = NA, occlusion_chain = NA,
                        distance_method = 'centroid',
                        drop_incomplete_residue = TRUE, rsa_method = 'rose',
                        dist_cutoff = 15, rsa_cutoff = 0.1,
                        sasa_cutoff = NA, only_exposed_in_patch = TRUE,
                        use_rsa_sasa = 'and',
                        max_patch = NA, interface_dist_cutoff = 5,
                        verbose = 1, detail_level = 1,
                        force_file_type = NULL, in_wrapper = TRUE,
                        patch_mode = 'placeholder'){

  # need to validate more inputs #

  # split chain and occlusion chain (they dont need grouped)
  chain = .split_chain(chain)
  occlusion_chain = .split_chain(occlusion_chain)

  # step 0: validate pdb and other arguments ----
  if(verbose > 0){
    cat('\tpdb_to_patch: Standardizing PDB input\n')
  }

  if(!in_wrapper){
    pdb = .standardize_pdb_input(pdb, force_file_type)
    # check that max_patch or dist_cutoff is set #
    if(is.na(max_patch) && is.na(dist_cutoff)){
      stop('max_patch and dist_cutoff cannot both be NA')
    }
  }

  distance_method = match.arg(distance_method, c('all', 'ca', 'centroid'))
  rsa_method = match.arg(rsa_method, c('rose', 'miller', 'theoretical_tien', 'empirical_tien'))
  use_rsa_sasa = match.arg(use_rsa_sasa, c('and', 'or'))

  # step 1: retrieve sequences for chains of interest ----
  if(verbose > 0){
    cat('\tpdb_to_patch: Extracting sequences for chains\n')
  }

  seq_set = .get_pdb_sequence(pdb, chain = chain, in_module = TRUE)

  # step 2: calculate residue-wise distance matrix ----
  if(verbose > 0){
    cat('\tpdb_to_patch: Calculating residue-wise distance matrix\n')
  }

  residue_dist = .calculate_residue_distance(pdb, chain = chain,
                                            distance_method = distance_method,
                                            in_module = TRUE)

  # step 3: calculate residue-wsie accessibility ----
  if(verbose > 0){
    cat('\tpdb_to_patch: Calculating residue-wise accessibility\n')
  }

  chain_set = c(chain, occlusion_chain)
  chain_set = chain_set[!is.na(chain_set)]

  # probably need a slot of occlusion chain, so its not in resulting residue df #
  residue_df = .calculate_accessibility(pdb, chain = chain_set,
                                       drop_incomplete = drop_incomplete_residue,
                                       method = rsa_method,
                                       in_module = T)

  # drop any occlusion set chains from residue_df #
  occlusion_chain = occlusion_chain[!is.na(occlusion_chain)]
  residue_df = residue_df[!residue_df$orig_chain %in% occlusion_chain,]

  # step 3.5: classifying exposure ----
  residue_df = .is_exposed(residue_df = residue_df,
                           rsa_cutoff = rsa_cutoff,
                           sasa_cutoff = sasa_cutoff,
                           use_rsa_sasa = use_rsa_sasa)


  # step 4: identify surface patches (expands residue_df) ----
  if(verbose > 0){
    cat('\tpdb_to_patch: Identifying patches\n')
  }

  residue_df = .identify_patches(residue_dist,
                                 residue_df,
                                 only_exposed_in_patch = only_exposed_in_patch,
                                 dist_cutoff = dist_cutoff,
                                 max_patch = max_patch)

  # step 5: capture interface patches #
  if(!(length(interface_chain) == 1 && is.na(interface_chain))){

    if(verbose > 0){
      cat('\tpdb_to_patch: Identifying interface patches\n')
    }

    # apply identify interface to all sets of interface chains #
    for(i in 1:length(interface_chain)){
      int_ch = unlist(strsplit(interface_chain[i], split = ''))
      interface_patches = .identify_interface(pdb, chain = chain, interface_chain = int_ch, dist_cutoff = interface_dist_cutoff)

      # add to residue df #
      residue_df[nrow(residue_df)+1,] = NA
      residue_df$residue_id[nrow(residue_df)] = interface_patches$name
      residue_df$patch[nrow(residue_df)] = interface_patches$interf
    }

  }

  # detail level will control return #
  # 0 - just pdb obj, seq_set, and residue_df
  # 1 - pdb obj, seq_set, residue_df, and chain
  # 2 - all

  # return list object
  results = list(
    pdb = pdb,
    chain = if (detail_level > 0) chain else NULL,
    seq_set = seq_set,
    residue_dist = if (detail_level > 1) residue_dist else NULL,
    residue_df = residue_df
  )

  class(results) = 'evo3D_pdb_info'

  return(results)

}

