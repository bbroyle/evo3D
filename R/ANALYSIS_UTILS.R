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

  file_paths <- file.path(output_dir, paste0(names(msa_subsets), ".fa"))
  for (i in seq_along(msa_subsets)) {
    fa_mat <- msa_subsets[[i]]
    seqs <- do.call(paste0, as.data.frame(fa_mat, stringsAsFactors = FALSE))
    lines <- paste0(">", row.names(fa_mat), "\n", seqs)
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

write_patch_fastas <- function(msa_subsets, output_dir = "patch_fastas") {
  temp_dir <- tempfile("patch_fastas_")
  dir.create(temp_dir)

  for (i in seq_along(msa_subsets)) {
    fa_mat <- msa_subsets[[i]]
    seqs <- do.call(paste0, as.data.frame(fa_mat, stringsAsFactors = FALSE))
    lines <- paste0(">", row.names(fa_mat), "\n", seqs)
    writeLines(lines, file.path(temp_dir, paste0(names(msa_subsets)[i], ".fa")))
  }

  # Set wd to temp_dir for relative tar
  old_wd <- getwd()
  setwd(temp_dir)
  on.exit(setwd(old_wd), add = TRUE)  # protect session even on failure

  tar_path <- tempfile(fileext = ".tar.gz")
  utils::tar(tarfile = tar_path, files = list.files(), compression = "gzip")

  setwd(old_wd)  # <-- move back BEFORE untar
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  utils::untar(tarfile = tar_path, exdir = output_dir)

  unlink(temp_dir, recursive = TRUE)
  unlink(tar_path)
}

# block entropy ----

#' compute block entropy
#'
#' calculates shannon entropy across haplotypes in an msa subset.
#' optionally translates nucleotide codons to amino acids first.
#'
#' @param x character matrix or vector of aligned sequences
#' @param translate logical. if TRUE (default), translate nucleotide codons
#'   to amino acids before computing entropy
#'
#' @return numeric shannon entropy value (bits, log base 2)
#' @export

block_entropy = function(x, translate = TRUE){
    # take a suite of haplotypes (msa_subsets) #
    # compute a table of frequencies #
    # compute entropy #

    # translate sequence #
    if(translate){
      x = apply(x, 1, function(seq) {
        seq = seqinr::translate(seq)
        paste(seq, collapse = '')
      })
    }

    # remove haplotypes with gaps and non AA characters #
    # 9.22.25 need to adjust if not translating #
    aavec = '[^AVILMWYFSTNQCGPRHKDE]'
    x = x[!grepl(aavec, x)]

    # compute shannon entropy #
    tab = table(x)
    p = tab / sum(tab)
    entropy = -sum(p * log2(p)) # 9.22.25 changed to log2 from log #

    return(entropy)
  }

# .write_patch_pymol() ----

#' write pymol commands for patches
#'
#' generates pymol command strings to visualize structural patches and centroids.
#' each patch gets a unique selection (patch_1, patch_2, ...) and a random color.
#' centroid residues are grouped and colored black.
#'
#' @param patches data frame with at least:
#'   \itemize{
#'     \item \code{patch} – residue groups formatted "resno_chain"
#'     \item \code{residue_id} – ids of centroid residues
#'   }
#'
#' @return named list of character vectors:
#'   \itemize{
#'     \item \code{patches} – pymol \code{select} commands for each patch
#'     \item \code{colors} – pymol \code{color} commands to assign colors
#'     \item \code{centroid} – commands to select and color centroids black
#'   }
#' @export

write_patch_pymol = function(patches){
  # expecting a data.frame with patch column #
  # patches will be converted to pymol selection commands
  # patch numbered to order in patches data.frame

  # for each row generate pymol command for selecting patch
  cmd_list = c()
  for(i in 1:nrow(patches)){
    patch = patches$patch[i]
    patch = data.frame(id = unlist(strsplit(patch, '\\+')))
    patch$resno = gsub('_.*', '', patch$id)
    patch$chain = gsub('^[^_]+_', '', patch$id)
    patch$ins = gsub('^[^_]+_', '', patch$chain)
    patch$chain = gsub('_.*', '', patch$chain)

    # each resi set should be paired with chain and patch
    patch_cmd = c()
    for(ch in unique(patch$chain)){
      resi = paste0(patch$resno[patch$chain == ch], patch$ins[patch$chain == ch], collapse = '+')
      cmd = paste0('resi ', resi, ' and chain ', ch)
      patch_cmd[length(patch_cmd)+1] = cmd
    }

    # combine patch_cmd into final output
    cmd = paste0('select patch_', i, ', ', paste0('(', patch_cmd, ')', collapse = ' or '))

    cmd_list[i] = cmd
  }

  # copy and paste this into pymol patches #
  #cat(cmd_list, sep = '\n')

  # get color per patch #
  colors <- grDevices::hcl.colors(nrow(patches), palette = 'Dynamic')
  colors = gsub('#', '', colors)
  colors = colors[sample(1:length(colors))]
  color_cmd = paste0('color 0x', colors, ', patch_', 1:nrow(patches))

  # copy and paste this into pymol for colors #
  #cat(color_cmd, sep = '\n')

  # last generate color black for all centroids #
  centroid = data.frame(id = patches$residue_id)
  centroid$resno = gsub('_.*', '', centroid$id)
  centroid$chain = gsub('^[^_]+_', '', centroid$id)
  centroid$ins = gsub('^[^_]+_', '', centroid$chain)
  centroid$chain = gsub('_.*', '', centroid$chain)

  # could exist on different chains #
  centroid_cmd = list()
  for(ch in unique(centroid$chain)){
    resi = paste0(centroid$resno[centroid$chain == ch], centroid$ins[centroid$chain == ch], collapse = '+')
    cmd = paste0('resi ', resi, ' and chain ', ch)
    centroid_cmd[length(centroid_cmd)+1] = cmd
  }

  # combine patch_cmd into final output
  centroid_cmd = paste0('select centroid, ', paste0('(', centroid_cmd, ')', collapse = ' or '))
  centroid_cmd[2] = 'color black, centroid'

  # copy and paste this into pymol patches #
  #cat(centroid_cmd, sep = '\n')

  # return these cmds for pymol #
  return(list(patches = cmd_list,
              colors = color_cmd,
              centroid = centroid_cmd))
}


# .write_stat_to_bfactor() ----

#' write selection stats to pdb b-factors
#'
#' embeds codon-aligned selection statistics (e.g. tajima's d, nucleotide diversity)
#' into the b-factor field of a pdb file for visualization in pymol, chimera, or
#' other structure viewers.
#'
#' @param evo3d_results list returned by \code{run_evo3d()}, must include
#'   \code{evo3d_df} and \code{pdb_info_sets}
#' @param pdb_id numeric index of pdb to annotate (default 1)
#' @param stat_name character. name of statistic column in \code{evo3d_df}
#'   (e.g. "tajima", "pi")
#' @param outfile file path for modified pdb output
#' @param mapped_chains_only logical. if TRUE (default), only chains mapped in
#'   the analysis are retained
#' @param scale_up_pi logical. if TRUE and \code{stat_name == "pi"}, rescales
#'   very small values upward for visibility in pymol
#' @param adjust_NA_stats numeric value written to residues with missing stats
#'   (default -99)
#'
#' @return no R object returned. writes a pdb file to \code{outfile}
#'   with b-factors replaced by statistic values
#' @export

write_stat_to_bfactor = function(evo3d_results, pdb_id = 1, stat_name = 'tajima', outfile = 'test.pdb',
                                  mapped_chains_only = TRUE, scale_up_pi = FALSE, adjust_NA_stats = -99){

  # check if pdb column has pdb1 tags or not #
  multi_run = any(grepl('pdb1', evo3d_results$evo3d_df))

  # set up pdb column to grab #
  if(!multi_run){
    pdb_col = 'residue_id'
    pdb_id = 1
  } else {
    pdb_col = paste0('pdb', pdb_id, '_residue_id')
  }

  # grab pdb #
  pdb_name = paste0('pdb', pdb_id)
  pdb = evo3d_results$pdb_info_sets[[pdb_name]]$pdb

  # --- #
  patch_df = evo3d_results$evo3d_df

  # --- #
  if(mapped_chains_only){
    grid = evo3d_results$call_info$run_grid
    chain = grid$chain[which(grid$pdb == pdb_name)]
    pdb = bio3d::trim.pdb(pdb, chain = chain)
  }

  # if stat is NA convert to 0 #
  #patch_df[is.na(patch_df[stat_name]), stat_name] = 0

  # set pdb insert to '' if NA
  pdb$atom$ins[is.na(pdb$atom$ins)] = ''
  pdb$atom$residue_id = paste(pdb$atom$resno, pdb$atom$chain, pdb$atom$ins, sep = '_')

  # get the stat values
  stat = patch_df[match(pdb$atom$residue_id, patch_df[[pdb_col]]), stat_name]

  # if pi is too low for pymol scale up -- print message #
  if(scale_up_pi && stat_name == 'pi') {
    # what is non zero non NA min #
    # dictates scaling factor #
    min_pi = min(stat[!is.na(stat) & stat > 0])

    if(min_pi < 0.01){
      scale_factor = -floor(log10(min_pi))
      scale_factor = scale_factor - 2 # just want to move into 0.01
      stat = stat * (10^scale_factor)
      message(paste0('Scaling up pi values by 10^', scale_factor, ' for visualization.'))
    }
  }

  # use residue id (resno_chain_ins) to match to pdb #
  pdb$atom$b = patch_df[match(pdb$atom$residue_id, patch_df[[pdb_col]]), stat_name]
  pdb$atom$b = round(pdb$atom$b, 2)
  pdb$atom$b = ifelse(is.na(pdb$atom$b), adjust_NA_stats, pdb$atom$b)



  # drop residue ID, and write
  pdb$atom = pdb$atom[,!names(pdb$atom) %in% c('residue_id')]
  bio3d::write.pdb(pdb = pdb, b = pdb$atom$b, file = outfile)

}

# run_pegas_three() ----

#' calculate diversity and neutrality stats for msa subsets
#'
#' computes nucleotide diversity (pi), tajima's d, and haplotype diversity
#' on nucleotide alignments. accepts a single alignment or a list of
#' patch-level alignments.
#'
#' typically used downstream of \code{aln_msa_to_pdb()} to quantify
#' variation within structure-informed windows.
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
#' @export

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
    #names(seqs) = names(msa)  # sync the names
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

# calculate_polymorphic_residue() ----

#' flag polymorphic codons and compute site entropy
#'
#' translates codon-aligned nucleotide msas to amino acids and flags codons
#' with amino acid polymorphism. also computes shannon entropy per codon
#' site from amino acid frequencies.
#'
#' works with single-msa or multi-msa evo3d results.
#'
#' @param msa_info_sets named list of msa info objects (each with \code{msa_mat})
#' @param residue_df data frame of patch residues (from \code{aln_msa_to_pdb()}),
#'   must include \code{codon_id}, \code{codon}, and optionally \code{msa}
#' @param valid_aa_only logical. if TRUE (default), only standard amino acids
#'   (AVILMWYFSTNQCGPRHKDE) are considered when testing polymorphism
#'
#' @return input \code{residue_df} with two added columns:
#'   \itemize{
#'     \item \code{polymorphic} – 1 if codon site has >1 amino acid, 0 otherwise
#'     \item \code{site_entropy} – shannon entropy (bits, log base 2) of aa distribution
#'   }
#' @export

calculate_polymorphic_residue = function(msa_info_sets, residue_df, valid_aa_only = TRUE){

  aa_vector = strsplit('AVILMWYFSTNQCGPRHKDE', '')[[1]]

  ids = unique(residue_df$msa)
  ids = ids[!is.na(ids)]

  # it is extneded data in which we need to cycle through msas #
  residue_df$polymorphic = NA
  residue_df$site_entropy = NA

  for(id in ids){
    msa = msa_info_sets[[id]]$msa_mat

    # handle residue_df null later #
    aa_set = t(apply(msa, 1, seqinr::translate))

    x = apply(aa_set, 2, table)

    # Check polymorphism for each codon position
    polymorphic = sapply(x, function(pos_table) {
      if(valid_aa_only) {
        valid_counts = pos_table[names(pos_table) %in% aa_vector]
        length(valid_counts[valid_counts > 0]) > 1
      } else {
        length(pos_table[pos_table > 0]) > 1
      }
    })

    # Calculate Shannon entropy for each position
    entropy = sapply(x, function(pos_table) {
      if(valid_aa_only) {
        valid_counts = pos_table[names(pos_table) %in% aa_vector]
        valid_counts = valid_counts[valid_counts > 0]
      } else {
        valid_counts = pos_table[pos_table > 0]
      }

      if(length(valid_counts) == 0) return(0)

      freqs = valid_counts / sum(valid_counts)
      -sum(freqs * log2(freqs))
    })

    # add to residue_df -- only where there are codons #

    # better way -- use codon_id #
    res = data.frame(codon_id = paste0(1:length(polymorphic), '_', id),
                     polymorphic = as.integer(polymorphic),
                     site_entropy = entropy
                     )


    # build lookups
    poly_lookup   <- setNames(res$polymorphic, res$codon_id)
    entropy_lookup <- setNames(res$site_entropy, res$codon_id)

    # assign by matching codon_id
    residue_df$polymorphic <- poly_lookup[residue_df$codon_id]
    residue_df$site_entropy <- entropy_lookup[residue_df$codon_id]

  }
}

# calculate_patch_entropy() ----

#' calculate mean amino acid entropy per patch
#'
#' translates codon-aligned nucleotide msas to amino acids and computes mean
#' shannon entropy across codon positions in each patch. this summarizes
#' amino acid diversity within structure-defined windows.
#'
#' @param msa single codon-aligned msa (matrix) or a named list of msas (one per patch)
#' @param residue_df optional data frame with column \code{msa_subset_id};
#'   if provided, results are merged into this data frame
#' @param valid_aa_only logical. if TRUE (default), restrict entropy to standard
#'   amino acids (ignores stop or X)
#'
#' @return data frame with one row per patch (\code{msa_subset_id}),
#'   including column \code{patch_entropy}
#' @export

calculate_patch_entropy = function(msa, residue_df = NULL, valid_aa_only = TRUE){
  # Convert single MSA to list for consistent handling
  if (!is.list(msa)) {
    msa = list(msa)
  }

  # Handle naming consistently
  if (is.null(names(msa))) {
    names(msa) = paste0('msa_', 1:length(msa))
    #names(seqs) = names(msa)  # sync the names
  }

  # if residue_df is provided we can add results to residue_df // else make a dataframe #
  if (is.null(residue_df)) {
    residue_df = data.frame(msa_subset_id = names(msa), stringsAsFactors = F)
  }

  seq_set = lapply(msa, function(x){
    t(apply(x, 1, seqinr::translate))
  })

  aa_vector = strsplit('AVILMWYFSTNQCGPRHKDE', '')[[1]]

  # Calculate Shannon entropy for each MSA subset (averaged over columns)
  entropy = lapply(seq_set, function(aa_matrix) {
    if(ncol(aa_matrix) == 0) return(0)

    # Calculate entropy for each column (position)
    col_entropies = apply(aa_matrix, 2, function(col) {
      if(valid_aa_only) {
        valid_aas = col[!is.na(col) & col %in% aa_vector]
      } else {
        valid_aas = col[!is.na(col)]
      }

      if(length(valid_aas) == 0) return(0)

      aa_counts = table(valid_aas)
      aa_freqs = aa_counts / sum(aa_counts)

      # Shannon entropy
      -sum(aa_freqs * log2(aa_freqs))
    })

    # Return average entropy across positions in patch
    mean(col_entropies, na.rm = TRUE)
  })


  # Add entropy results to residue_df
  residue_df$patch_entropy = NA
  residue_df$patch_entropy[match(names(entropy), residue_df$msa_subset_id)] = unlist(entropy)

  return(residue_df)

}

