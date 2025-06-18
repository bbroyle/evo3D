# ------------------------------------------ #
# ANALYSIS_UTILS.R 
# utilities for msa statistics and writing outputs
# Brad Broyles 
# ------------------------------------------ #


# .write_patch_fastas() ----

#' Write Patch-Level MSA Subsets to FASTA Files
#'
#' Writes nucleotide MSA subsets for each structural patch to separate FASTA files.
#' Useful for downstream tools or manual inspection.
#'
#' @param msa_subsets A named list of character matrices (typically from \code{aln_msa_to_pdb()}), each representing one patch.
#' @param output_dir Output directory where FASTA files will be saved. Created if it doesn't exist.
#'
#' @return Invisibly returns \code{NULL}. Files are written to disk.
#' @export
write_patch_fastas = function(msa_subsets, output_dir = 'patch_fastas'){
  # nuc_patches is from map_patches_to_nucleotides()
  # fasta dir for saving files
  # save.file to save files or if F return list of fastas
  
  # create dir if it doesn't exist
  if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = T)
  }
  
  for(i in 1:length(msa_subsets)){
    fa = msa_subsets[[i]]
    fa = apply(fa, 1, paste, collapse = '')
    fname = paste0(output_dir, '/', names(msa_subsets)[[i]], '.fa')
    
    # fix fasta seqs for writing
    fa = paste0('>', names(fa), '\n', fa, '\n')
    cat(fa, sep = '', file = fname)
  }
  
}


# .write_patch_pymol() ----

#' Generate PyMOL Selection and Coloring Commands for Structural Patches
#'
#' Converts structural residue patches into PyMOL command strings for selection, coloring, and centroid annotation.
#' Each patch is assigned a unique selection (e.g., \code{patch_1}) and a randomly shuffled color.
#'
#' These commands can be copied directly into a PyMOL session to visualize surface patches and their centroids.
#'
#' @param patches A data frame containing at least \code{patch} (residue groups, formatted as "residue_chain") and \code{residue_id} columns.
#'
#' @return A named list of character vectors:
#' \describe{
#'   \item{patches}{PyMOL \code{select} statements for each structural patch.}
#'   \item{colors}{PyMOL \code{color} commands to visually differentiate patches.}
#'   \item{centroid}{A \code{select} and \code{color} command to highlight centroid residues in black.}
#' }
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

#' Write Selection Statistics to PDB B-Factors (Codon-Based Mapping)
#'
#' Maps codon-aligned selection statistics (e.g., Tajima's D, nucleotide diversity) to B-factors in a PDB structure file.
#' This enables visualization of selection results in molecular viewers such as PyMOL or Chimera.
#'
#' This function operates on output from \code{run_evo3d()} and uses codon-patch-to-structure mappings to embed
#' a selected statistic into the B-factor field of a target PDB.
#'
#' @param evo3d_results Output from \code{run_evo3d()}, including \code{evo3d_df} and structure info.
#' @param pdb_id Numeric index for the target PDB (e.g., 1 for \code{pdb1}). Defaults to 1.
#' @param stat_name Name of the statistic column in \code{evo3d_df} to embed (e.g., \code{"tajima"}, \code{"pi"}).
#' @param outfile Output path for the modified PDB file.
#' @param mapped_chains_only Logical; if \code{TRUE}, only atoms from chains used in the analysis will be retained.
#' @param scale_up_pi Logical; if \code{TRUE} and \code{stat_name == "pi"}, small values will be scaled for PyMOL visibility.
#' @param adjust_NA_stats Numeric value to assign to residues with missing statistic values (default: -10).
#'
#' @return No R return value. A PDB file is written to \code{outfile} with modified B-factor values.
#' @export
write_stat_to_bfactor = function(evo3d_results, pdb_id = 1, stat_name = 'tajima', outfile = 'test.pdb', 
                                  mapped_chains_only = TRUE, scale_up_pi = FALSE, adjust_NA_stats = -10){
  
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
#' Calculate Diversity and Neutrality Statistics for MSA Subsets
#'
#' Computes per-patch summary statistics for nucleotide alignments, including nucleotide diversity (\code{pi}),
#' Tajima's D, and haplotype diversity. Accepts either a single MSA or a list of patch-level MSAs.
#'
#' This function is typically used downstream of \code{aln_msa_to_pdb()} to quantify evolutionary variation
#' within each structure-informed alignment window.
#'
#' @param msa Either a single nucleotide alignment (matrix or DNAbin) or a named list of alignments (one per patch).
#' @param residue_df Optional data frame containing at minimum a column \code{msa_subset_id}.
#'                   If provided, statistics will be merged into this data frame.
#' @param stat Character vector indicating which statistics to compute. Options: \code{"pi"}, \code{"tajima"}, \code{"hap"}.
#'
#' @return A data frame with one row per alignment window (matching \code{msa_subset_id}),
#' including columns for each requested statistic.
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

#' Flag Polymorphic Codons and Compute Site Entropy
#'
#' Translates codon-aligned nucleotide MSAs into amino acids and identifies codon sites with
#' evidence of amino acid polymorphism. Also calculates per-site Shannon entropy based on amino acid usage.
#'
#' Works with both single-MSA and extended multi-MSA evo3D results (via \code{msa_info_sets}).
#'
#' @param msa_info_sets A named list of MSA info objects (e.g., \code{msa1}, \code{msa2}, ...) each with \code{msa_mat}.
#' @param residue_df A data frame of patch residues, typically from \code{aln_msa_to_pdb()}; must include \code{codon} and optionally \code{msa}.
#' @param valid_aa_only Logical; if \code{TRUE}, only standard amino acids are considered polymorphic (excludes X, *, etc.).
#'
#' @return An updated \code{residue_df} with added columns:
#' \describe{
#'   \item{polymorphic}{Binary indicator: 1 if codon is polymorphic (multiple amino acids observed), 0 otherwise.}
#'   \item{site_entropy}{Shannon entropy of amino acid distribution at the codon site.}
#' }
#' @export
calculate_polymorphic_residue = function(msa_info_sets, residue_df, valid_aa_only = TRUE){

  aa_vector = strsplit('AVILMWYFSTNQCGPRHKDE', '')[[1]]
  
  # is it extended? #
  if(!'msa' %in% colnames(residue_df)) {
    # only one info set #
    msa = msa_info_sets$msa1$msa_mat
    
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
    codon_ro = which(!is.na(residue_df$codon))
    
    residue_df$polymorphic = NA
    residue_df$polymorphic[codon_ro] = as.integer(polymorphic)
    
    residue_df$site_entropy = NA
    residue_df$site_entropy[codon_ro] = entropy
  } else {
    # it is extneded data in which we need to cycle through msas #
    residue_df$polymorphic = NA
    residue_df$site_entropy = NA
    
    msa_ids = unique(residue_df$msa)
    msa_ids = msa_ids[!is.na(msa_ids)]
    for(id in msa_ids){
      msa_name = paste0('msa', id)
      msa = msa_info_sets[[msa_name]]$msa_mat
      
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
      codon_ro = which(!is.na(residue_df$codon))
      msa_ro = which(residue_df$msa == id)
      codon_ro = intersect(codon_ro, msa_ro)
      codon_ro = sort(codon_ro)
      
      residue_df$polymorphic[codon_ro] = as.integer(polymorphic)
      
      residue_df$site_entropy[codon_ro] = entropy
      
    }
  }
  
  return(
    residue_df
  )
}

# calculate_patch_entropy() ----

#' Calculate Mean Amino Acid Entropy per Patch
#'
#' Translates codon-aligned nucleotide MSAs to amino acids and computes mean Shannon entropy
#' across all codon positions within each patch. This summarizes the overall amino acid diversity
#' of each structural window.
#'
#' @param msa Either a single MSA (matrix) or a named list of nucleotide MSAs (one per patch).
#'            Each should be codon-aligned (i.e., divisible by 3 in width).
#' @param residue_df Optional data frame with \code{msa_subset_id}; if provided, entropy is merged in.
#' @param valid_aa_only Logical; if \code{TRUE}, only standard amino acids (no stop/X) are included in entropy calculation.
#'
#' @return A data frame with one row per patch (matching \code{msa_subset_id}), including a new column \code{patch_entropy}.
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

