# ------------------------------------------ #
# ANALYSIS_UTILS.R
# utilities for msa statistics and writing outputs
# Brad Broyles
# ------------------------------------------ #


# write_patch_fastas_slow() ----

#' write patch fasta files
#'
#' writes msa subsets to fasta files, one per patch
#'
#' @param msa_subsets list of msa matrices (as from \code{extract_msa_subsets()})
#' @param output_dir directory to write fasta files (default "patch_fastas")
#'
#' @keywords internal

write_patch_fastas_slow = function(msa_subsets, output_dir = 'patch_fastas') {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  file_paths = file.path(output_dir, paste0(names(msa_subsets), ".fa"))
  for (i in seq_along(msa_subsets)) {
    fa_mat = msa_subsets[[i]]
    seqs = do.call(paste0, as.data.frame(fa_mat, stringsAsFactors = FALSE))
    lines = paste0(">", row.names(fa_mat), "\n", seqs)
    writeLines(lines, con = file_paths[i])
  }
}

# write_patch_fastas() ----

#' write patch-level msa subsets to fasta files
#'
#' writes nucleotide msa subsets (one per structural patch) to fasta files.
#' files are written to a temporary directory, bundled into a compressed tarball,
#' then extracted into the output directory. this ensures safe creation and
#' cleanup of intermediate files.
#'
#' @param msa_subsets named list of character matrices (from \code{aln_msa_to_pdb()}),
#'   each matrix representing a patch
#' @param output_dir directory to extract fasta files into (default "patch_fastas").
#'   created if it does not exist
#'
#' @return invisibly returns \code{NULL}. fasta files are written to disk
#' @export

write_patch_fastas = function(msa_subsets, output_dir = "patch_fastas") {
  temp_dir = tempfile("patch_fastas_")
  dir.create(temp_dir)

  for (i in seq_along(msa_subsets)) {
    fa_mat = msa_subsets[[i]]
    seqs = do.call(paste0, as.data.frame(fa_mat, stringsAsFactors = FALSE))
    lines = paste0(">", row.names(fa_mat), "\n", seqs)
    writeLines(lines, file.path(temp_dir, paste0(names(msa_subsets)[i], ".fa")))
  }

  # Set wd to temp_dir for relative tar
  old_wd = getwd()
  setwd(temp_dir)
  on.exit(setwd(old_wd), add = TRUE)  # protect session even on failure

  tar_path = tempfile(fileext = ".tar.gz")
  utils::tar(tarfile = tar_path, files = list.files(), compression = "gzip")

  setwd(old_wd)  # <-- move back BEFORE untar
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  utils::untar(tarfile = tar_path, exdir = output_dir)

  unlink(temp_dir, recursive = TRUE)
  unlink(tar_path)
}


# write_stat_to_pdb() ----

#' Write selection statistics to PDB B-factors and occupancy
#'
#' Embed selection statistics (e.g., Tajima's D, nucleotide diversity) into the
#' \code{b} (temperature factor) and occupancy (\code{q}) fields of a PDB file for
#' visualization in PyMOL, Chimera, or other structure viewers. In PyMOL, these
#' fields can be visualized with commands such as \code{spectrum b} and
#' \code{spectrum q}.
#'
#' If two statistics are provided, the first is written to \code{b} and the second
#' is written to \code{q}.
#'
#' @param evo3d_results An \code{evo3D_results} object returned by \code{run_evo3d()}.
#'   Must include \code{evo3d_df} and \code{pdb_info_sets}.
#' @param pdb_id Character; PDB identifier to annotate (e.g., \code{"pdb1"}).
#' @param stat_name Character vector of length 1 or 2 naming column(s) in
#'   \code{evo3d_results$evo3d_df}.
#' @param outfile File path for the modified PDB output. If the file already exists,
#'   a timestamp is appended to avoid overwriting.
#' @param mapped_chains_only Logical; if \code{TRUE} (default), retain only chains
#'   mapped in the analysis for the specified PDB.
#' @param adjust_NA_stats Numeric; value written where statistics are missing
#'   (default \code{-99}).
#'
#' @return Invisibly returns the output file path. Writes a PDB file to \code{outfile}
#'   with \code{b} and \code{q} replaced by statistic values.
#'
#' @export

write_stat_to_pdb = function(evo3d_results, pdb_id = 'pdb1', stat_name = 'tajima',
                             outfile = 'test.pdb', mapped_chains_only = TRUE,
                             adjust_NA_stats = -99){

  # validate inputs #
  if (!inherits(evo3d_results, "evo3D_results")) {
    stop("Expecting evo3D_results object")
  }

  if (!is.character(pdb_id) || length(pdb_id) != 1 || is.na(pdb_id)) {
    stop("pdb_id must be a single character id (e.g., 'pdb1')")
  }

  if (!is.character(stat_name) || length(stat_name) < 1 || length(stat_name) > 2) {
    stop("stat_name must be a character vector of length 1 or 2")
  }

  if (!is.numeric(adjust_NA_stats) || length(adjust_NA_stats) != 1 || is.na(adjust_NA_stats)) {
    stop("adjust_NA_stats must be a single numeric value")
  }

  # --- avoid overwrite ---
  if (file.exists(outfile)) {
    base = tools::file_path_sans_ext(outfile)
    ext  = tools::file_ext(outfile)
    tag  = format(Sys.time(), "%Y%m%d%H%M%S")
    outfile2 = if (nzchar(ext)) paste0(base, "_", tag, ".", ext) else paste0(base, "_", tag)
    message(outfile, " already exists, saving to: ", outfile2)
    outfile = outfile2
  }

  # --- check for pdb ---
  if (is.null(evo3d_results$pdb_info_sets[[pdb_id]]$pdb)) {
    stop("No pdb found for pdb_id = ", pdb_id)
  }

  pdb = evo3d_results$pdb_info_sets[[pdb_id]]$pdb
  patch_df = evo3d_results$evo3d_df
  patch_df = patch_df[patch_df$pdb == pdb_id,]

  # --- #
  if(isTRUE(mapped_chains_only)){
    grid = evo3d_results$call_info$run_grid
    chain = unique(grid$chain[which(grid$pdb == pdb_id)])
    chain = chain[!is.na(chain)]
    pdb = bio3d::trim.pdb(pdb, chain = chain)
  }

  # use residue id (resno_chain_ins) to match to pdb #
  pdb$atom$b = patch_df[match(pdb$atom$residue_id, patch_df$residue_id), stat_name[1]]
  pdb$atom$b = round(pdb$atom$b, 2)
  pdb$atom$b = ifelse(is.na(pdb$atom$b), adjust_NA_stats, pdb$atom$b)

  # if length stat_name > 1 we can add the second to q factor
  if(length(stat_name) > 1){
    pdb$atom$o = patch_df[match(pdb$atom$residue_id, patch_df$residue_id), stat_name[2]]
    pdb$atom$o = round(pdb$atom$o, 2)
    pdb$atom$o = ifelse(is.na(pdb$atom$o), adjust_NA_stats, pdb$atom$o)
  }

  # write out #
  bio3d::write.pdb(pdb = pdb, file = outfile)

  invisible(outfile)

}



# run_pegas_three() ----

#' Calculate diversity and neutrality statistics for MSA subsets
#'
#' Compute selection statistics—including nucleotide diversity (pi), Tajima's D,
#' and haplotype diversity—on nucleotide multiple sequence alignments. Accepts a
#' single alignment or a named list of patch-level alignments.
#'
#' Typically used downstream of \code{aln_msa_to_pdb()} to quantify variation within
#' structure-informed windows.
#'
#' @param msa single nucleotide alignment (matrix or DNAbin) or a named list
#'   of alignments, one per patch
#' @param residue_df optional data frame with column \code{msa_subset_id};
#'   if provided, statistics are merged into this data frame
#' @param stat character vector of statistics to compute.
#'   options: \code{"pi"}, \code{"tajima"}, \code{"hap"} (default all)
#'
#' @return data frame with one row per alignment window (\code{msa_subset_id}),
#'   including requested columns:
#'   \itemize{
#'     \item \code{pi} – nucleotide diversity
#'     \item \code{tajima} – tajima's d
#'     \item \code{tajima_pnormal} – normal approximation p-value
#'     \item \code{tajima_pbeta} – beta distribution p-value
#'     \item \code{hap} – haplotype diversity
#'   }
#' @keywords internal

run_pegas_three = function(msa, residue_df = NULL, stat = c('pi', 'tajima', 'hap')) {

  # Convert single MSA to list for consistent handling
  if (!is.list(msa)) {
    msa = list(msa)
  }

  # Now everything is a list - one code path!
  seqs = ape::as.DNAbin.list(msa)

  # Handle naming consistently
  if (is.null(names(msa))) {
    names(msa) = paste0('msa_', 1:length(msa))
  }

  # if residue_df is provided we can add results to residue_df // else make a dataframe #
  if (is.null(residue_df)) {
    residue_df = data.frame(msa_subset_id = names(msa), stringsAsFactors = F)
  }

  # check for pi #
  if('pi' %in% stat) {
    pi = lapply(seqs, pegas::nuc.div)
    names(pi) = names(msa)

    residue_df$pi = NA
    residue_df$pi[match(names(pi), residue_df$msa_subset_id)] = unlist(pi)
  }

  # check for tajima #
  if('tajima' %in% stat) {
    taj = lapply(seqs, pegas::tajima.test)
    names(taj) = names(msa)

    residue_df$tajima = NA
    residue_df$tajima_pnormal = NA
    residue_df$tajima_pbeta = NA

    residue_df$tajima[match(names(taj), residue_df$msa_subset_id)] = unlist(lapply(taj, function(x) x[1]))
    residue_df$tajima_pnormal[match(names(taj), residue_df$msa_subset_id)] = unlist(lapply(taj, function(x) x[2]))
    residue_df$tajima_pbeta[match(names(taj), residue_df$msa_subset_id)] = unlist(lapply(taj, function(x) x[3]))
  }

  # check for hap #
  if('hap' %in% stat) {
    hap = lapply(seqs, pegas::hap.div) # gives warnings if gaps present (okay?) -- also makes more hap than needed if ambiguity
    names(hap) = names(msa)

    residue_df$hap = NA

    residue_df$hap[match(names(hap), residue_df$msa_subset_id)] = unlist(hap)
  }

  return(residue_df)
}

# calculate_site_entropy() ----

#' Calculate per-site (codon-position) Shannon entropy
#'
#' Computes Shannon entropy at each aligned position for one or more MSAs and
#' merges results into \code{residue_df} by \code{codon_id}. If \code{seqtype} is
#' \code{"nucleotide"}, sequences are translated to amino acids before entropy is
#' computed.
#'
#' This is typically called downstream of MSA–structure mapping, where
#' \code{residue_df} contains \code{msa}, \code{codon_id}, and other per-residue
#' annotations.
#'
#' @param msa_info_sets Named list of MSA info objects containing \code{msa_mat}.
#' @param residue_df Data frame of residue/codon annotations. Must include
#'   \code{msa} (MSA id) and \code{codon_id} (e.g., \code{"12_msa1"}).
#' @param valid_aa_only Logical; if \code{TRUE} (default), compute entropy using
#'   only standard amino acid characters.
#' @param seqtype Character; \code{"nucleotide"} or \code{"peptide"}.
#' @param gencode Integer genetic code used by \code{seqinr::translate()} when
#'   \code{seqtype = "nucleotide"} (default 1).
#'
#' @return \code{residue_df} with a new column \code{site_entropy}, containing
#'   per-position Shannon entropy values (base 2). Rows without a matching
#'   \code{codon_id} remain \code{NA}.
#'
#'@keywords internal

calculate_site_entropy = function(msa_info_sets, residue_df, valid_aa_only = TRUE,
                                  seqtype = 'nucleotide', gencode = 1){

  aa_vector = strsplit('AVILMWYFSTNQCGPRHKDE', '')[[1]]

  ids = unique(residue_df$msa)
  ids = ids[!is.na(ids)]

  # it is extneded data in which we need to cycle through msas #
  residue_df$site_entropy = NA

  for(id in ids){
    msa = msa_info_sets[[id]]$msa_mat

    # if necessary translate #
    if(seqtype == 'nucleotide'){
      aa_set = t(apply(msa, 1, seqinr::translate, numcode = gencode))
    } else {
      aa_set = msa
    }

    x = apply(aa_set, 2, table)

    # Calculate Shannon entropy for each position
    entropy = sapply(x, function(pos_table) {
      if(valid_aa_only) {
        valid_counts = pos_table[names(pos_table) %in% aa_vector]
        valid_counts = valid_counts[valid_counts > 0]
      } else {
        valid_counts = pos_table[pos_table > 0]
      }

      if(length(valid_counts) == 0) return(NA)

      freqs = valid_counts / sum(valid_counts)
      -sum(freqs * log2(freqs))
    })

    # add to residue_df -- only where there are codons #

    # better way -- use codon_id #
    res = data.frame(codon_id = paste0(1:length(entropy), '_', id),
                     site_entropy = entropy
    )


    # -------- #
    # THIS PART OVERWRITING IF MORE THAN ONE MSA #
    # 11.24.25 #
    # -------- #

    if(FALSE){
      # build lookups
      poly_lookup  =  setNames(res$polymorphic, res$codon_id)
      entropy_lookup =  setNames(res$site_entropy, res$codon_id)

      # assign by matching codon_id
      residue_df$polymorphic = poly_lookup[residue_df$codon_id]
      residue_df$site_entropy =  entropy_lookup[residue_df$codon_id]
    }

    for(i in 1:nrow(res)){
      ro = which(residue_df$codon_id == res$codon_id[i])
      residue_df$site_entropy[ro] = res$site_entropy[i]
    }
  }

  return(residue_df)
}

# summarize_stat_at_patch() ----
#' Summarize a site-level statistic over patch membership
#'
#' Compute a per-row patch summary of a site-level statistic (e.g., entropy,
#' diversity) using a patch-membership column that encodes member keys as a
#' \code{"+"}-separated string. For each row, patch members are matched to
#' \code{key} and the requested summary statistic is computed across
#' \code{stat_col}.
#'
#' @param residue_df Data frame containing a site-level statistic and patch
#'   membership information.
#' @param stat_col Character; name of the numeric column to summarize
#'   (default \code{"site_entropy"}).
#' @param patch Character; name of the column containing patch membership as a
#'   \code{"+"}-separated string of keys (default \code{"codon_patch"}).
#' @param key Character; name of the column used to match patch members to rows
#'   (default \code{"codon_id"}).
#' @param method Character; summary method. One of \code{"mean"}, \code{"min"},
#'   \code{"max"}, or \code{"median"} (default \code{"mean"}).
#' @param na.rm Logical; if \code{TRUE} (default), ignore \code{NA} values when
#'   summarizing.
#'
#' @return \code{residue_df} with a new numeric column named
#'   \code{paste0(method, "_", stat_col)} containing patch-level summaries.
#'
#' @export

summarize_stat_at_patch = function(residue_df, stat_col = 'site_entropy',
                                   patch = 'codon_patch', key = 'codon_id',
                                   method = 'mean', na.rm = TRUE){

  method = match.arg(method, c('mean', 'min', 'max', 'median'))

  # see if stat_col exists in df if not exit #
  need = c(stat_col, patch, key)
  miss = setdiff(need, colnames(residue_df))
  if (length(miss)) {
    message("Skipping: missing column(s): ", paste(miss, collapse = ", "))
    return(residue_df)
  }


  new_col = paste0(method, '_', stat_col)
  residue_df[[new_col]] = NA_real_

  # for each patch - grab the rows (how to avoid duplicates grab first ro of each) #
  for(i in seq_len(nrow(residue_df))){
    if(is.na(residue_df[i,patch])) next
    if(nchar(residue_df[i,patch]) == 0) next

    p = residue_df[i,patch]
    p = strsplit(p, '\\+')[[1]]

    ro = match(p, residue_df[[key]])

    ro2 = ro[!is.na(ro)]

    # if length of ro2 != ro then some positions in the patch have no row in dataframe #
    if(length(ro2) != length(ro)){
      cat('some patch positions not found in dataframe')
    }

    vals = residue_df[ro2,stat_col]

    if(method == 'mean'){
      x = mean(vals, na.rm = na.rm)
    } else if (method == 'min'){
      x = min(vals, na.rm = na.rm)
    } else if (method == 'max'){
      x = max(vals, na.rm = na.rm)
    } else {
      x = stats::median(vals, na.rm = na.rm)
    }

    residue_df[[new_col]][i] = x

  }

  return(residue_df)

}


# block entropy ----

#' Compute block entropy
#'
#' Calculates Shannon entropy across haplotypes in an MSA subset. Optionally
#' translates nucleotide codons to amino acids prior to entropy calculation.
#'
#' @param x Character matrix or vector of aligned sequences.
#' @param translate Logical; if \code{TRUE} (default), translate nucleotide codons
#'   to amino acids before computing entropy.
#' @param gencode Integer genetic code for \code{seqinr::translate()} (default 1).
#' @param valid_aa_only Logical; if \code{TRUE} (default), exclude haplotypes
#'   containing non-standard amino acid characters.
#'
#' @return Numeric Shannon entropy value (bits, log base 2).
#' @keywords internal

block_entropy = function(x, translate = TRUE, gencode = 1, valid_aa_only = TRUE){
  # take a suite of haplotypes (msa_subsets) #
  # compute a table of frequencies #
  # compute entropy #

  # translate sequence #
  if(translate){
    x = apply(x, 1, function(seq) {
      seq = seqinr::translate(seq, numcode = gencode)
      paste(seq, collapse = '')
    })
  }

  # remove haplotypes with gaps and non AA characters #
  aavec = '[^AVILMWYFSTNQCGPRHKDE]'
  if(valid_aa_only){
    x = x[!grepl(aavec, x)]
  }

  # compute shannon entropy #
  tab = table(x)
  p = tab / sum(tab)
  entropy = -sum(p * log2(p))

  return(entropy)
}

# calculate_patch_stats() ----

#' Calculate patch-level selection statistics from MSA subsets
#'
#' Module-level wrapper that computes requested statistics on structure-informed
#' MSA subsets and merges results into \code{residue_df}. Statistics are toggled
#' via \code{stat_controls} and include per-site entropy, patch-summarized entropy,
#' block (haplotype) entropy, and nucleotide diversity/neutrality statistics.
#'
#' Site entropy and average patch entropy require \code{msa_info_sets} (output of
#' \code{msa_to_ref()}) to access full-length MSAs and construct \code{codon_id}
#' lookups. Nucleotide statistics (\code{pi}, Tajima's D, haplotype diversity) are
#' computed only when \code{seqtype = "nucleotide"}.
#'
#' @param msa_info_sets Optional; named list returned by \code{msa_to_ref()}.
#'   Required when \code{stat_controls$calc_site_entropy} or
#'   \code{stat_controls$calc_avg_patch_entropy} is \code{TRUE}.
#' @param final_msa_subsets Named list of MSA subsets (one per window/patch),
#'   typically produced by the MSA–structure mapping pipeline. Names should match
#'   \code{residue_df$msa_subset_id}.
#' @param residue_df Data frame of residue/window annotations. Must include
#'   \code{msa_subset_id}. Additional columns (e.g., \code{msa}, \code{codon_id},
#'   and patch membership columns) are required depending on enabled statistics.
#' @param stat_controls Named list of logical flags controlling which statistics
#'   are computed. Common entries include:
#'   \itemize{
#'     \item \code{calc_site_entropy}
#'     \item \code{calc_avg_patch_entropy}
#'     \item \code{calc_block_entropy}
#'     \item \code{calc_pi}
#'     \item \code{calc_tajima}
#'     \item \code{calc_hap}
#'     \item \code{valid_aa_only} (used for entropy calculations)
#'   }
#' @param seqtype Character; \code{"nucleotide"} or \code{"peptide"}.
#' @param gencode Integer genetic code used for translation-based calculations
#'   when \code{seqtype = "nucleotide"} (default 1).
#'
#' @return \code{residue_df} with additional statistic columns appended, depending
#'   on \code{stat_controls}. Potential columns include \code{site_entropy},
#'   \code{mean_site_entropy} (or other summary name from the patch summarizer),
#'   \code{block_entropy}, \code{pi}, \code{tajima}, and \code{hap}.
#'
#' @export

calculate_patch_stats = function(msa_info_sets, final_msa_subsets, residue_df,
                                 stat_controls = list(calc_pi = TRUE),
                                 seqtype = 'nucleotide', gencode = 1){

  # if doing site entropy or avg_patch entropy need both #
  if(isTRUE(stat_controls$calc_site_entropy) || isTRUE(stat_controls$calc_avg_patch_entropy)){

    if(is.null(msa_info_sets)){
      stop('Need to set msa_info_sets to output of msa_to_ref() to calcualte site entropy of avg_patch_entropy')
    }

    residue_df = calculate_site_entropy(
      msa_info_sets,
      residue_df,
      valid_aa_only = stat_controls$valid_aa_only,
      seqtype = seqtype,
      gencode = gencode
    )
  }

  if(isTRUE(stat_controls$calc_avg_patch_entropy)){
    # needs site_entropy
    residue_df = summarize_stat_at_patch(residue_df = residue_df,
                                         stat_col = 'site_entropy')
  }

  if(isTRUE(stat_controls$calc_block_entropy)){
    ents = unlist(lapply(final_msa_subsets, function(x) {
      # compute block entropy #
      translate = (seqtype == 'nucleotide')
      block_entropy(x, translate = translate, gencode = gencode)
    }))

    residue_df$block_entropy = NA
    for(i in 1:nrow(residue_df)){
      # id #
      id = residue_df$msa_subset_id[i]
      if(is.na(id)) next
      residue_df$block_entropy[i] = ents[which(names(ents) == id)]
    }
  }

  stat = c()
  if (isTRUE(stat_controls$calc_pi)) {
    stat = c(stat, 'pi')
  }

  if (isTRUE(stat_controls$calc_tajima)) {
    stat = c(stat, 'tajima')
  }

  if (isTRUE(stat_controls$calc_hap)) {
    stat = c(stat, 'hap')
  }

  # also slow but works
  if(length(stat) > 0){
    if(seqtype == 'nucleotide'){
      residue_df = run_pegas_three(msa = final_msa_subsets,
                                   residue_df = residue_df,
                                   stat = stat)
    } else {
      message('MSA subsets are protein. skipping nucleotide stats --', stat, '--\n')
    }
  }

  return(residue_df)


}

# generate_null_model ----

#' Generate a null model of patch-level MSA subsets
#'
#' Generate a null distribution of fixed-length patch windows by sampling codon
#' identifiers from the observed results. Sampling can optionally be weighted to
#' match the empirical codon frequency observed across all patches. The resulting
#' null patches are converted into MSA subsets for downstream statistic
#' calculation.
#'
#' This function is typically used to construct null distributions for comparison
#' against structure-informed patch statistics.
#'
#' @param evo3d_results An \code{evo3D_results} object containing patch definitions,
#'   MSA information, and run metadata.
#' @param n Integer; number of null patches to generate (default 10,000).
#' @param len Integer; number of codons per patch (default 15).
#' @param seed Integer random seed for reproducibility.
#' @param match_codon_frequency Logical; if \code{TRUE} (default), sample codons
#'   with probability proportional to their frequency in observed patches.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{null_df}: Data frame containing generated null patch definitions
#'     (\code{codon_patch}, \code{msa_subset_id}, \code{msa}).
#'   \item \code{msa_subsets}: Named list of extracted MSA subsets corresponding
#'     to the null patches.
#' }
#'
#' @export
generate_null_model = function(evo3d_results, n = 10000,
                               len = 15, seed = 1219, match_codon_frequency = TRUE){
  # generate n haplotypes according to codon usage of results #
  # designed to operate on fixed length #

  # grab all patches from results #
  patches = evo3d_results$evo3d_df$codon_patch
  patches = patches[!is.na(patches)]
  patches = paste0(patches, collapse = '+')
  patches = strsplit(patches, '\\+')

  x = table(patches)
  ids = names(x)

  if(match_codon_frequency){
    weight = as.numeric(x)
  } else {
    weight = rep(1, length(x))
  }

  # sample
  set.seed(seed)

  samples = replicate(
    n,
    sample(ids, size = len, replace = FALSE, prob = weight),
    simplify = FALSE
  )

  df = data.frame(
    codon_patch = sapply(samples, function(x) paste(x, collapse = "+")),
    stringsAsFactors = FALSE
  )

  df$msa_subset_id = 1:nrow(df)

  # add msa id holder -- .extract uses msa id tied to codon in codon patch - msa column just place hoder -- #
  msas = unique(evo3d_results$evo3d_df$msa)
  df$msa[1:length(msas)] = msas
  df$msa[(length(msas)+1):nrow(df)] = msas[1]

  msa_subsets = .extract_msa_subsets(evo3d_results$msa_info_sets,
                                     df,
                                     use_sample_names = evo3d_results$call_info$aln_controls$use_sample_names,
                                     genetic_code = evo3d_results$call_info$msa_controls$genetic_code)

  return(list(
    null_df = df,
    msa_subsets = msa_subsets)
  )
}

# filter-overlaps ----

#' Filter overlapping patch windows
#'
#' Reduce a set of patch definitions by removing patches that overlap too strongly
#' with previously retained patches. Overlap is quantified as the proportion of
#' shared codon identifiers relative to the size of the candidate patch.
#'
#' Patches are evaluated sequentially in the order provided. A patch is retained
#' only if its overlap with all previously retained patches is less than or equal
#' to the specified threshold.
#'
#' @param df Data frame containing a \code{codon_patch} column, where each entry is
#'   a \code{"+"}-separated string of codon identifiers.
#' @param overlap Numeric; maximum allowed fractional overlap with an existing
#'   patch (default \code{1/3}).
#'
#' @return A subset of \code{df} containing only retained (non-overlapping) patches,
#'   in their original order.
#'
#' @export
filter_overlaps = function(df, overlap = 1/3) {

  # keeping a reduced subset of original patches #

  keep = logical(nrow(df))
  ref_sets = list()

  for (i in seq_len(nrow(df))) {
    codon_set = strsplit(df$codon_patch[i], "\\+")[[1]]

    # overlap test against every patch already kept
    clash = vapply(ref_sets, function(s)
      (length(intersect(s, codon_set)) / length(codon_set)) > overlap,
      logical(1))

    if (!any(clash)) {
      keep[i] = TRUE
      ref_sets[[length(ref_sets) + 1]] = codon_set
    }
  }

  # return reduced set #
  return(df[keep, , drop = FALSE])
}
