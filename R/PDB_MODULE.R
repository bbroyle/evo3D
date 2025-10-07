# --------------------------------------------------------------- #
# PDB MODULE
#
# utilities for working with PDB and mmCIF data, returning standardized
# pdb_info objects ready for downstream evo3D analysis
#
# main exported wrapper: pdb_to_patch()
# internal .functions() handle standardization, sequence extraction,
# residue distances, solvent accessibility, exposure, and interface contacts
#
# bonus utilities: download matching PDBs, visualize chain layouts
#
# NOTES:
# - internal helpers keep dot prefix, but some may be exported for modular use
# - hydrogens and hetatoms are removed by default
# - rsa and sasa both supported; sasa recommended where possible
#
# contact: bbroyle@purdue.edu
# --------------------------------------------------------------- #

# helper functions ----
# .find_matching_structures() * intended for an explore_evo3d() mode * ----

#' find matching pdb structures
#'
#' blast a protein sequence against the rcsb pdb to find candidate structures.
#' returns top hits filtered by sequence identity, with option to download pdbs
#' and generate a hit summary plot.
#' useful when no structure is supplied and you want evo3d to suggest starting models.
#'
#' @param pep protein sequence as a single character string
#' @param identity_cutoff minimum sequence identity percentage (default 80)
#' @param max_hits maximum number of structures to return (default 5)
#' @param generate_plot logical, whether to return a ggplot summary (default TRUE)
#' @param download_pdbs logical, whether to download matching pdbs (default TRUE)
#' @param output_dir directory to save downloaded structures (default "retrieved_pdbs")
#'
#' @return a list with elements:
#' \describe{
#'   \item{top_hits}{data frame of top pdb hits passing filters}
#'   \item{all_hits}{data frame of all blast hits with status annotations}
#'   \item{hit_plot}{ggplot object if \code{generate_plot=TRUE}, else NULL}
#'   \item{pdb_table}{table of downloaded pdb files if \code{download_pdbs=TRUE}, else NULL}
#' }
#'
#' @details
#' this is a wrapper for \code{bio3d::blast.pdb()} and \code{bio3d::get.pdb()}.
#' results are filtered on sequence identity, then the top hits are kept.
#' external database queries can be rate limited, so avoid running in tight parallel loops.
#'
#' @seealso \code{\link{download_structures}}, \code{bio3d::blast.pdb}, \code{bio3d::get.pdb}
#'
#' @keywords experimental
#' @keywords internal

.find_matching_structures = function(pep, identity_cutoff = 80, max_hits = 5, generate_plot = T,
                                    download_pdbs = T, output_dir = 'retrieved_pdbs'){
  # If no PDB is provided to program -- this function will blast RCSB pdb for matching structures #
  # We want to return multiple structures for more robust window and solvent accessibility calculations #
  # essentially a wrapper for bio3d::blast.pdb() & bio3d::get.pdb() with helpful plot and table output

  stopifnot(is.character(pep), length(pep) == 1)
  if (!is.numeric(identity_cutoff) || identity_cutoff < 0 || identity_cutoff > 100) {
    stop('`identity_cutoff` must be numeric between 0 and 100')
  }
  if (!is.numeric(max_hits) || max_hits < 1) {
    stop('`max_hits` must be a positive integer')
  }

  message(
    "\n⚠️ This function queries external databases; avoid rapid parallel calls to prevent rate-limiting.\n"
  )

  # blast pep (expecting single protein sequence) against PDB database
  bl = bio3d::blast.pdb(pep)
  hits = bl$hit.tbl
  hits = hits[order(hits$mlog.evalue, decreasing = T),]

  # filter on identity cutoff and then grab top hits #
  filtered_hits = hits[hits$identity >= identity_cutoff,]
  top_hits = filtered_hits[seq_len(min(max_hits, nrow(filtered_hits))),]

  # add status to hits
  hits$pdb_status = 'raw_hit'
  hits$pdb_status[hits$subjectids %in% filtered_hits$subjectids] = 'filtered_hit'
  hits$pdb_status[hits$subjectids %in% top_hits$subjectids] = 'top_hit'

  # generate plot if desired
  if(generate_plot){
    hit_plot = ggplot2::ggplot(hits, ggplot2::aes(alignmentlength, identity, color = pdb_status))+
      ggplot2::geom_jitter(height = 0.1, width = 0.5)+
      ggrepel::geom_text_repel(ggplot2::aes(label = subjectids), show.legend = F)+
      ggplot2::geom_hline(yintercept = identity_cutoff, linetype = 'dashed')+
      ggplot2::theme_bw()+
      ggplot2::scale_color_manual(values = c('raw_hit' = '#d95f02',
                                             'filtered_hit' = '#7570b3',
                                             'top_hit' = '#1b9e77'))+
      ggplot2::labs(title = 'PDB Blast Results',
                    subtitle = paste0('max_hits = ',
                                      max_hits, ', identity_cutoff = ',
                                      identity_cutoff, '%'),
                    x = 'Reported Alignment Length',
                    y = 'Identity (%)')
  } else {
    hit_plot = NULL
  }

  # ---------------------------------------------------------------------------- #
  # check if any hits passed filters -- if not stop here and return data at hand #
  if(nrow(top_hits) == 0){
    message('\n🛑 No hits found >= identity cutoff: ', identity_cutoff, '% 🛑',
            '\nInspect `$all_hits` and `$hit_plot` for more information\n',
            'If you want to continue manually provide `all_hits` to:\n',
            'download_structures(output$all_hits)\n')

    return(list(
      top_hits = top_hits,
      all_hits = hits,
      hit_plot = hit_plot,
      pdb_table = NULL)
    )
  }

  # ---------------------------------------------------------------------------- #
  # if download is desired, call download_structures() #
  # return hit_tables and hit_plot #
  if(download_pdbs){
    hit_table = suppressMessages(download_structures(top_hits, output_dir))

    return(list(
      top_hits = top_hits,
      all_hits = hits,
      hit_plot = hit_plot,
      pdb_table = hit_table)
    )
  } else {
    # return top_hits, all_hits, and hit_plot #
    return(list(
      top_hits = top_hits,
      all_hits = hits,
      hit_plot = hit_plot,
      pdb_table = NULL)
    )
  }

}


#' Download PDB Structures
#'
#' Downloads structures from a table of BLAST hits or PDB IDs.
#'
#' @param hit_table Data frame of BLAST hits or character vector of PDB IDs.
#' @param output_dir Output directory to save downloaded PDB files.
#'
#' @return A data frame with PDB IDs, chain IDs, and local file paths.
#' @export
download_structures = function(hit_table, output_dir = 'retrieved_pdbs'){
  # hit_table is one of the two table outputs of find_matching_structures() #
  # hit_list$all_hits or hit_list$top_hits #
  # or a user provided table with at least subjectids column #
  # or can be a character vector of subjectids #

  # subjectids are four letter pdb codes -- anything longer is trimmed in bio3d::get.pdb() #
  # user may provide unique subjectids that are duplicated in the first 4 characters #
  # thats fine I will add downloaded file path to both #

  # !!! output_dir cannot be empty !!! # -- need to check #
  if(output_dir == ''){
    output_dir = '.'
  } else if (grepl('/$', output_dir)){
    output_dir = substr(output_dir, 1, nchar(output_dir)-1)
  }

  message("\n⚠️ This function makes live queries to external databases ⚠️\n",
          "Please avoid running in parallel or sending requests too quickly\n",
          "You may be rate-limited or blocked\n")

  # check if hit_table is a character vector #
  if(is.character(hit_table)){
    hit_table = data.frame(subjectids = hit_table)
  }

  # check if subjectids column is available #
  if(!'subjectids' %in% colnames(hit_table)){
    message('\n🛑 hit_table must have a subjectids column or be a character vector 🛑\n')
    return(NULL)
  }

  # make output directory (also handled by bio3d::get.pdb())
  if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = T)
  }

  # output is variable depending on content of output directory and success of download #
  # if all successful (or previously downloaded) output is unnamed vector of file paths #
  # if some failed "1" and newly download get "0", while previous download get path #
  file_paths = suppressWarnings(bio3d::get.pdb(hit_table$subjectids, path = output_dir))

  # if all download then unnamed vector -- if some missing then named vector output #
  if(is.null(names(file_paths))){
    names(file_paths) = file_paths
  }

  # use names as paths -- drop paths that failed
  file_paths = file_paths[file_paths != '1']

  # add file_paths to hit_table (force ordering)
  ins = toupper(substr(hit_table$subjectids, 1, 4))
  colord = paste0(output_dir, '/', ins, '.pdb')

  hit_table$file_paths = names(file_paths[colord])

  # return a table of subject ids and file paths
  results = data.frame(
    pdb_id = gsub('_.+', '', hit_table$subjectids),
    chain_id = gsub('^[^_]+|_', '', hit_table$subjectids),
    file_path = hit_table$file_paths,
    stringsAsFactors = F
  )

  return(results)

}


# .plot_chain_map () ----

#' plot mds map of pdb chains
#'
#' computes a 2d mds projection of Cα coordinates from a pdb structure
#' and plots chain positions using base R graphics.
#' each chain is colored separately and labeled at its median coordinate.
#'
#' @param pdb pdb input: file path, bio3d pdb object, or standardized pdb from \code{.standardize_pdb_input()}
#' @param chain optional chain identifier to subset (default NA = all chains)
#' @param in_module logical, internal use flag to skip extra validation (default FALSE)
#'
#' @return no return value. produces a base R plot
#' @export

.plot_chain_map = function(pdb, chain = NA, in_module = FALSE) {

  # local helper to compute coords
  mds_pdb = function(pdb, chain, in_module) {
    pdb = .standardize_pdb_input(pdb)
    ca = pdb$atom[pdb$atom$elety == "CA", ]

    if (!is.na(chain)) {
      ca = ca[ca$chain == chain, ]
    }
    ca$insert = ifelse(is.na(ca$insert), "", ca$insert)

    mds_coords = cmdscale(dist(ca[, c("x", "y", "z")]), k = 2)

    data.frame(
      residue_id = paste0(ca$resno, "_", ca$chain, "_", ca$insert),
      chain = as.character(ca$chain),
      x = mds_coords[, 1],
      y = mds_coords[, 2]
    )
  }

  # compute coords
  plot_df = mds_pdb(pdb, chain, in_module = in_module)

  # median coords for labels
  label_df = aggregate(cbind(x, y) ~ chain, data = plot_df, FUN = median)

  # make sure chain is character
  chain_ids = unique(plot_df$chain)

  if (length(chain_ids) == 0) {
    stop("No chains found for plotting")
  }

  # if only one chain, still assign a color
  cols = if (length(chain_ids) == 1) "black" else rainbow(length(chain_ids))
  col_map = setNames(cols, chain_ids)

  # plot points
  plot(
    plot_df$x, plot_df$y,
    col = col_map[plot_df$chain],
    pch = 19, cex = 0.7,
    xlab = "", ylab = "",
    axes = FALSE, main = "MDS projection of Cα coordinates by chain"
  )

  # add labels at medians
  text(label_df$x, label_df$y, labels = label_df$chain, col = "black", font = 2, cex = 1)

  # add legend if more than one chain
  if (length(chain_ids) > 1) {
    legend("topright", legend = chain_ids, col = cols, pch = 19, bty = "n", cex = 1)
  }
}





# PDB UTILS ----

# .standardize_pdb_input() ----

#' standardize pdb input
#'
#' reads a pdb or cif file (or an existing bio3d object) and returns a cleaned pdb
#' object for downstream evo3d analysis. fixes missing insert codes, replaces
#' illegal '+' characters, and appends a unique residue identifier column.
#'
#' @param pdb pdb input: file path, or object from \code{bio3d::read.pdb()} or \code{bio3d::read.cif()}
#' @param force_file_type optional character; override auto-detection.
#'   one of \code{"cif"} or \code{"pdb"}. default NULL = auto-detect by file extension.
#'
#' @details
#' This function performs several preprocessing steps:
#'
#' - If the PDB contains '+' in residue numbers, chain IDs, or insert codes,
#'   these are replaced with a safe alternate character (one of \code{!,$,\%,\&,~,@}).
#'   If no safe replacement is available, the function stops with an error.
#' - Insert codes with NA are replaced by \code{""} (empty string).
#' - Residue IDs are constructed as \code{"resno_chain_insert"}
#'   and stored in \code{atom$residue_id}.
#' @export

.standardize_pdb_input = function(pdb, force_file_type = NULL){

  # expects single entry #
  input_class = class(pdb)[1]

  # At the end of if/else -- either pdb is read in or returns error #
  if(input_class == 'character'){
    if(!file.exists(pdb)){
      stop('PDB file does not exist: ', pdb)
    }

    if(!is.null(force_file_type)){
      if(force_file_type == 'cif'){
        pdb = bio3d::read.cif(pdb)
      } else if(force_file_type == 'pdb'){
        pdb = bio3d::read.pdb(pdb)
      } else {
        stop('Invalid force_file_type. Use "cif" or "pdb".')
      }
    } else {
      if(grepl('cif$', pdb)){
        pdb = bio3d::read.cif(pdb)
      } else {
        pdb = bio3d::read.pdb(pdb)
      }
    }
  } else if(!input_class == 'pdb'){
    # print i dont know what you have #
    stop('pdb provided is not a recognized format. Please provide a file path, an object from bio3d::read.pdb(), or an object from bio3d::read.cif().')
  }

  # fix insert NA to ''
  pdb$atom$insert[is.na(pdb$atom$insert)] = ''

  # "+" is used internally to combine residues into patches -- make sure "+" doesnt exist to start #
  in_resno = any(grepl('\\+', pdb$atom$resno))
  in_chain = any(grepl('\\+', pdb$atom$chain))
  in_insert = any(grepl('\\+', pdb$atom$insert))

  if(any(in_resno, in_chain, in_insert)){
    # replacements for "+"
    replacements = c("!","$","%","&","~")

    # used characters in set #
    all_chars = paste(pdb$atom$resno, pdb$atom$chain, pdb$atom$insert, collapse = '')

    # which replacement can we use
    repl_char = NULL
    for (cand in replacements) {
      if (!grepl(cand, all_chars, fixed = TRUE)) {
        repl_char = cand
        break
      }
    }

    # need to stop here and tell user I couldnt process this PDB it has "+" and couldn't use
    # any of the replacement characters (all in use)

    if(is.null(repl_char)){
      stop("PDB contains '+' but all replacement characters (!,$,%,& ,~, @) are already in use. Cannot sanitize input safely.")
    }

    # replace those "+"
    pdb$atom$resno = gsub('+', repl_char, pdb$atom$resno, fixed = TRUE)
    pdb$atom$chain = gsub('+', repl_char, pdb$atom$chain, fixed = TRUE)
    pdb$atom$insert = gsub('+', repl_char, pdb$atom$insert, fixed = TRUE)

  }

  # add residue id information #
  pdb$atom$residue_id = paste0(pdb$atom$resno, '_', pdb$atom$chain, '_', pdb$atom$insert)

  # pass along pdb #
  return(pdb)

}


# .get_pdb_sequence() ----

#' Extract PDB Sequences
#'
#' Retrieves amino acid sequences for one or more chains from a PDB structure.
#'
#' Called internally by PDB alignment and chain-detection functions.
#'
#' @param pdb A PDB object from \code{bio3d}, or file path passed to
#'   \code{.standardize_pdb_input()} if \code{in_module = FALSE}.
#' @param chain Optional character vector of chain IDs. If \code{NA}, all
#'   protein chains are extracted.
#' @param in_module Logical. If \code{TRUE}, assumes \code{pdb} has already been
#'   standardized (default \code{FALSE}).
#'
#' @return A named character vector of amino acid sequences, one per chain.
#' @keywords internal

.get_pdb_sequence = function(pdb, chain = NA, in_module = FALSE){

  #if running in module (.standardize_pdb_input() is already run) #
  if(!in_module){
    pdb = .standardize_pdb_input(pdb)
  }

  if(length(chain) == 1 && is.na(chain)) {
    chain = unique(pdb$atom$chain)
  }

  # get pdb sequences with names as chains
  aa_seq = sapply(chain, function(x){
    paste0(
      bio3d::pdbseq(bio3d::trim.pdb(pdb, chain = x, 'protein')),
      collapse = '')
  }, USE.NAMES = T)

  return(aa_seq)
}


# .calculate_residue_distance() ----

#' Residue-wise Distance Matrix
#'
#' Computes a pairwise Euclidean distance matrix between residues of a PDB structure.
#' Distances are based on atom coordinates, grouped by residue, with options to
#' restrict to backbone, sidechain, or Cα atoms.
#'
#' Called internally by patch-definition and exposure functions.
#'
#' @param pdb A PDB object from \code{bio3d}, or file path passed to
#'   \code{.standardize_pdb_input()} if \code{in_module = FALSE}.
#' @param chain Optional character vector of chain IDs to retain. Default NA
#'   keeps all chains.
#' @param distance_method Character string, one of \code{"all"}, \code{"ca"},
#'   \code{"backbone"}, or \code{"sidechain"}. Controls which atoms are used for
#'   distance calculation. Default \code{"all"}.
#' @param in_module Logical. If \code{TRUE}, assumes \code{pdb} has already been
#'   standardized (default \code{FALSE}).
#'
#' @return A symmetric square numeric matrix of residue–residue distances,
#'   with residue IDs as row and column names.
#' @keywords internal

.calculate_residue_distance = function(pdb, chain = NA, distance_method = 'all', in_module = FALSE){

  # distance_method can be 'backbone', 'sidechain', 'ca', 'all'
  # hydrogens always removed

  #if running in module (.standardize_pdb_input() is already run) #
  if(!in_module){
    pdb = .standardize_pdb_input(pdb)
  }

  # trim chains
  if(!(length(chain) == 1 && is.na(chain))) {
    pdb = bio3d::trim.pdb(pdb, chain = chain)
  }

  # pdb should be protein only and no H (what about glycan distance?)
  pdb = bio3d::trim.pdb(pdb, 'protein')
  pdb = bio3d::trim.pdb(pdb, 'noh')

  # METHOD IS TOO MUCH #
  # check method, maybe we are only doing backbone or sidechain
  if(distance_method == 'backbone'){
    pdb = bio3d::trim.pdb(pdb, 'backbone')

  } else if(distance_method == 'sidechain'){
    pdb = bio3d::trim.pdb(pdb, 'sidechain')

  } else if(distance_method == 'ca'){
    pdb = bio3d::trim.pdb(pdb, 'calpha')
  }

  # create residue / atom identifiers
  pdb$atom$insert[is.na(pdb$atom$insert)] = ''
  pdb$atom$residue_id = paste0(pdb$atom$resno, '_', pdb$atom$chain, '_', pdb$atom$insert)

  # get atom distance matrix
  res_dist = bio3d::dm.xyz(pdb$xyz, grpby = pdb$atom$residue_id)

  # set row and column names
  colnames(res_dist) = rownames(res_dist) = unique(pdb$atom$residue_id)

  # set diag to 0 and make symmetrical
  res_dist = pmin(res_dist, t(res_dist), na.rm = TRUE)
  diag(res_dist) = 0

  # return residue wise distance matrix
  return(res_dist)

}

# .calculate_accessibility() ----

#' Calculate Solvent Accessibility
#'
#' calculates per-residue solvent accessible surface area (SASA) using
#' a DSSP-style algorithm reimplemented in R and C++.
#'
#' Hydrogens are removed, residues are first renumbered sequentially to ensure
#' contiguous indexing during calculation, then mapped back to their original
#' PDB residue numbers, chain IDs, and insert codes. Accessibility is reported
#' as both absolute SASA (Å²) and normalized relative solvent accessibility (RSA).
#'
#' Called internally by surface-exposure and patch-definition functions.
#'
#' @param pdb A PDB object from \code{bio3d}, or file path passed to
#'   \code{.standardize_pdb_input()} if \code{in_module = FALSE}.
#' @param chain Optional character vector of chain IDs to retain.
#'   Default NA keeps all chains.
#' @param method Character string specifying normalization scheme,
#'   one of \code{"rose"}, \code{"miller"}, \code{"theoretical_tien"},
#'   or \code{"empirical_tien"} (default \code{"rose"}).
#' @param drop_incomplete Logical. If \code{TRUE} (default), drop residues
#'   missing backbone atoms (N, CA, C, O), matching DSSP behavior.
#' @param in_module Logical. If \code{TRUE}, assumes \code{pdb} has already been
#'   standardized (default \code{FALSE}).
#'
#' @return A data frame with one row per residue, containing:
#' \itemize{
#'   \item \code{residue_index} – sequential index used internally
#'   \item \code{aa} – one-letter amino acid code
#'   \item \code{orig_resno}, \code{orig_chain}, \code{orig_insert} – original PDB identifiers
#'   \item \code{sasa} – absolute solvent accessible surface area (Å²)
#'   \item \code{rsa} – relative solvent accessibility (0–1, capped at 1)
#'   \item \code{residue_id} – concatenated unique ID (\code{orig_resno_chain_insert})
#' }
#' @export

.calculate_accessibility = function(pdb, chain = NA, method = 'rose', drop_incomplete = TRUE, in_module = FALSE){
  # return residue wise solvent accessibility

  #if running in module (.standardize_pdb_input() is already run) #
  if(!in_module){
    pdb = .standardize_pdb_input(pdb)
  }

  # check method is valid #
  if(!method %in% c('rose', 'miller', 'theoretical_tien', 'empirical_tien')){
    message('method not recognized - setting to default: rose')
  }

  # trim chains
  if(!(length(chain) == 1 && is.na(chain))) {
    pdb = bio3d::trim.pdb(pdb, chain = chain)
  }

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

  # Drop hydrogen atoms and non-ATOM records (sometimes data is but into b factor that is too large -- messing up elesy - can we use elety)
  #atoms = atoms[atoms$elesy != 'H' & atoms$type == 'ATOM', ]
  # dropping elety
  atoms = atoms[atoms$elety != 'H' & atoms$type == 'ATOM', ]

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
  #Rcpp::sourceCpp('src/dssp_sasa_rewrite.cpp')  # should move outside
  accessibility = calculateDSSPAccessibility(atom_df, residue_df)

  # Add results to the residue dataframe (round to 2 decimals) ** doesnt save mem **
  residue_df$sasa = round(accessibility, 2)

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

#' Classify Exposed Residues
#'
#' Marks residues as exposed based on solvent accessibility thresholds.
#' Supports filtering by relative solvent accessibility (RSA), absolute
#' solvent accessible surface area (SASA), or both. When both are provided,
#' the combination can be controlled with \code{use_rsa_sasa}.
#'
#' Called internally after \code{.calculate_accessibility()} to add an
#' \code{exposed} column to the residue data frame.
#'
#' @param residue_df Data frame of residues, typically output from
#'   \code{.calculate_accessibility()}, containing \code{rsa} and/or \code{sasa}.
#' @param rsa_cutoff Numeric. RSA threshold for exposure classification
#'   (default \code{NA}, ignored).
#' @param sasa_cutoff Numeric. SASA threshold for exposure classification
#'   (default \code{NA}, ignored).
#' @param use_rsa_sasa Character, either \code{"and"} or \code{"or"}.
#'   Controls how RSA and SASA criteria are combined if both are provided.
#'   Default \code{"and"}; automatically reset to \code{"and"} if either
#'   cutoff is missing.
#'
#' @return The input \code{residue_df} with an added logical column
#'   \code{exposed}, marking residues that meet the exposure criteria.
#' @keywords internal

.is_exposed = function(residue_df, rsa_cutoff = NA,
                       sasa_cutoff = NA, use_rsa_sasa = 'and') {

  # how to combine (and or else)
  if(any(is.na(c(rsa_cutoff, sasa_cutoff)))){
    use_rsa_sasa = 'and'
  }

  # build exposure under rsa (if available)
  if (!is.na(rsa_cutoff)) {
    pass_rsa = residue_df$rsa >= rsa_cutoff
  } else {
    pass_rsa = TRUE
  }

  # build exposure under sasa (if available)
  if (!is.na(sasa_cutoff)) {
    pass_sasa = residue_df$sasa >= sasa_cutoff
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

#' identify surface patches
#'
#' defines residue patches by growing neighborhoods around exposed seed residues
#' using rsa, sasa, and spatial distance criteria. seeds are always taken from
#' residues marked \code{exposed}; if cutoffs are set to 0, all residues qualify
#' as seeds.
#'
#' @param dist_mat residue–residue distance matrix
#' @param residue_df data frame of residues with rsa, sasa, and exposure info
#' @param dist_cutoff numeric. maximum distance in angstroms for neighbors (default 15)
#' @param max_patch numeric. optional maximum number of neighbors per patch (default NA = unlimited)
#' @param only_exposed_in_patch logical. if TRUE (default), restricts patch members to exposed residues only
#'
#' @return input \code{residue_df} with new columns:
#' \itemize{
#'   \item \code{patch} – string of neighbor residue ids joined by '+'
#'   \item \code{patch_len} – number of residues in patch
#'   \item \code{max_dist} – maximum pairwise distance among patch members
#' }
#' @keywords internal

.identify_patches = function(dist_mat, residue_df,
                             dist_cutoff = 15, max_patch = NA,
                             only_exposed_in_patch = TRUE){

  # in case dist mat is too large #
  res_ids = intersect(rownames(dist_mat), residue_df$residue_id)
  dist_mat = dist_mat[res_ids, res_ids]

  # first is to get seed residues #
  # only exposed residues are seeds #
  seeds = residue_df[residue_df$exposed,]

  # If only exposed can be in patch -- shrink to seeds #
  if(only_exposed_in_patch){
    in_set = seeds$residue_id
    dist_mat = dist_mat[in_set, in_set]
  }

  # get neighbors (could be dist_cutoff, or could be max_patch, or both)
  seeds$patch = NA
  seeds$patch_len = NA
  seeds$max_dist = NA

  for(i in 1:nrow(seeds)){
    center = seeds$residue_id[i]

    # get neighbors (sorted by distance) #
    neighbors = sort(dist_mat[center,])

    # filtered by dist_cut #
    if(!is.na(dist_cutoff)){
      neighbors = neighbors[neighbors <= dist_cutoff]
    }

    # filtered by max_patch #
    if(!is.na(max_patch)){
      neighbors = neighbors[1:min(c(length(neighbors), max_patch))]
    }

    # add to seeds


    if (length(neighbors) > 0) {
      seeds$patch[i] = paste0(names(neighbors), collapse = '+')
      seeds$patch_len[i] = length(neighbors)
      seeds$max_dist[i] = round(max(neighbors),2)
    } else {
      patch_len = 0
    }
  }

  # add back to residue_df
  idx = match(seeds$residue_id, residue_df$residue_id)
  residue_df$patch[idx] = seeds$patch
  residue_df$patch_len[idx] = seeds$patch_len
  residue_df$max_dist[idx] = seeds$max_dist



  return(residue_df)
}

# .identify_interfaces() ----

#' identify interface contacts
#'
#' finds residue contacts between two chain sets in a pdb structure using atom distances.
#'
#' @param pdb bio3d pdb object
#' @param chain chain id(s) of interest
#' @param interface_chain chain id(s) to test for contacts with
#' @param dist_cutoff maximum angstrom distance for contact (default 5)
#'
#' @return list with:
#' \itemize{
#'   \item \code{name} – interface name string of the form "interface_[chain]_[interface_chain(s)]"
#'   \item \code{interf} – concatenated residue ids (e.g. "35_A_+42_A_") of contacting residues
#' }
#' @keywords internal

.identify_interface = function(pdb, chain = NULL, interface_chain = NULL, dist_cutoff = 5){

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

#' extract surface patches from pdb
#'
#' high-level wrapper to compute sequences, residue distances,
#' solvent accessibility, exposure, and surface patches from a pdb structure.
#'
#' @param pdb pdb input: file path, bio3d object, or standardized pdb
#' @param chain chain id(s) to analyze (default NA = all chains)
#' @param interface_chain optional chain id(s) to define interface contacts
#' @param occlusion_chain optional chain id(s) used to occlude surface during accessibility calculation
#' @param distance_method residue distance metric: one of "all", "ca", "backbone", "sidechain" (default "all")
#' @param drop_incomplete_residue logical. drop residues missing backbone atoms (default TRUE)
#' @param rsa_method method for rsa normalization ("rose", "miller", "theoretical_tien", "empirical_tien")
#' @param dist_cutoff numeric. maximum distance in angstroms for neighbors in patches (default 15)
#' @param rsa_cutoff minimum rsa for defining seed residues (default 0.1)
#' @param sasa_cutoff optional minimum sasa for seed residues
#' @param only_exposed_in_patch logical. if TRUE, restricts patch members to exposed residues only (default TRUE)
#' @param use_rsa_sasa logical operator for combining rsa and sasa cutoffs: "and" or "or" (default "and")
#' @param max_patch optional maximum number of neighbors per patch
#' @param interface_dist_cutoff maximum distance in angstroms for interface residue contacts (default 5)
#' @param verbose integer. verbosity level (0 silent, >0 prints progress)
#' @param detail_level controls return content: 0 = minimal, 1 = include chain, 2 = include distance matrix
#' @param force_file_type optional. override auto-detection of pdb vs cif input
#' @param patch_mode placeholder argument, not implemented
#'
#' @return list with:
#' \itemize{
#'   \item \code{pdb} – NULL (object cached elsewhere)
#'   \item \code{chain} – chain ids if \code{detail_level > 0}
#'   \item \code{seq_set} – amino acid sequences per chain
#'   \item \code{residue_dist} – residue–residue distance matrix if \code{detail_level > 1}
#'   \item \code{residue_df} – data frame with residues, rsa/sasa values, exposure, and patch membership
#' }
#' @export

pdb_to_patch = function(pdb, chain = NA, interface_chain = NA, occlusion_chain = NA,
                        distance_method = 'all',
                        drop_incomplete_residue = TRUE, rsa_method = 'rose',
                        dist_cutoff = 15, rsa_cutoff = 0.1,
                        sasa_cutoff = NA, only_exposed_in_patch = TRUE,
                        use_rsa_sasa = 'and',
                        max_patch = NA, interface_dist_cutoff = 5, verbose = 1, detail_level = 1,
                        force_file_type = NULL, patch_mode = 'justaholder'){

  # need to validate more inputs #

  # split chain and occlusion chain (they dont need grouped)
  if(!(length(chain) == 1 && is.na(chain))) {
    chain = unlist(strsplit(chain, split = ''))
  }

  if(!(length(occlusion_chain) == 1 && is.na(occlusion_chain))) {
    occlusion_chain = unlist(strsplit(occlusion_chain, split = ''))
  }

  # step 0: validate pdb and return pdb object ----
  if(verbose > 0){
    cat('\tpdb_to_patch: Standardizing PDB input\n')
  }

  pdb = .standardize_pdb_input(pdb, force_file_type)
  # actually dont need because pdbs cached in wrapper #
  # add in_wrapper flag #

  # step 1: retrieve sequences for chains of interest ----
  if(verbose > 0){
    cat('\tpdb_to_patch: Extracting sequences for chains\n')
  }

  seq_set = .get_pdb_sequence(pdb, chain = chain)

  # step 2: calculate residue-wise distance matrix ----
  if(verbose > 0){
    cat('\tpdb_to_patch: Calculating residue-wise distance matrix\n')
  }

  residue_dist = .calculate_residue_distance(pdb, chain = chain,
                                            distance_method = distance_method,
                                            in_module = T)

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

  # residue_dist could have residues not present in residue_df ** if drop_incomplete_residue #
  # just keep in mind #

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
                                residue_df, only_exposed_in_patch = only_exposed_in_patch,
                                dist_cutoff = dist_cutoff,
                                max_patch = max_patch)

  # step 5: capture interface patches #
  # --- 9.30.25 -- may need debug dont know yet
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
  return(list(
    pdb = NULL,  # 8/20/25 - turned null (already stored elsewhere)
    chain = if (detail_level > 0) chain else NULL,
    seq_set = seq_set,
    residue_dist = if (detail_level > 1) residue_dist else NULL,
    residue_df = residue_df
  ))

}

