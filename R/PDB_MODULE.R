# --------------------------------------------------------------- #
# PDB MODULE goal is to return pdb_info list
# 1. pdb object
# 2. 
# 
# DEV NOTES:
# 1. should mds proj be saved so I dont have to remake them 
# 2. 
# email me: bbroyle@purdue.edu
# --------------------------------------------------------------- #



# DOWNLOAD PDB UTILS ----
#' Find Matching PDB Structures
#'
#' Uses BLAST to find and optionally download PDB structures matching a query peptide sequence.
#'
#' @param pep A protein sequence (character string).
#' @param identity_cutoff Minimum identity percentage (default 80).
#' @param max_hits Maximum number of structures to return (default 5).
#' @param generate_plot Logical, whether to return a ggplot summary.
#' @param download_pdbs Logical, whether to download matching PDBs.
#' @param output_dir Directory to save downloaded structures.
#'
#' @return A list with BLAST results, plots, and optionally downloaded structure paths.
#' @export
find_matching_structures = function(pep, identity_cutoff = 80, max_hits = 5, generate_plot = T,
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


# EXPLORE PDB UTILS ----
.mds_pdb = function(pdb, chain = NA, in_module = F){
  
  # check if pdb is a file path or object
  pdb = .standardize_pdb_input(pdb)
  
  # maybe just ca 
  ca = pdb$atom[pdb$atom$elety == 'CA', ]
  
  # if chain is provided, filter by chain
  if(!is.na(chain)){
    ca = ca[ca$chain == chain, ]
  }
  
  ca$insert = ifelse(is.na(ca$insert), "", ca$insert)
  
  # Use atom-level distances
  mds_coords <- cmdscale(dist(ca[, c("x", "y", "z")]), k = 2)
  
  # Add chain info for coloring
  mds_df <- data.frame(
    residue_id = paste0(ca$resno, '_', ca$chain, '_', ca$insert),
    chain = ca$chain,
    x = mds_coords[, 1],
    y = mds_coords[, 2]
  )
  
  # return
  return(mds_df)
}

.plot_chain_map = function(pdb, chain = NA, in_module = F){
  
  # get mds coords
  plot_df = .mds_pdb(pdb, chain, in_module = in_module)
  
  # drawing contact:: 
  #geom_segment(data = edge_df, aes(x = x1, y = y1, xend = x2, yend = y2), color = "gray", alpha = 0.5)
  
  # get median x and y per chain
  label_df <- aggregate(cbind(x, y) ~ chain, data = plot_df, FUN = median)
  
  p1 <- ggplot2::ggplot(plot_df, ggplot2::aes(x, y, color = chain)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_label(data = label_df,
                        ggplot2::aes(label = chain, fill = chain),
                        color = "black", show.legend = FALSE) +
    ggplot2::theme_void() +
    ggplot2::ggtitle("MDS Projection of C\u03b1 Coordinates by Chain") +
    ggplot2::theme(plot.margin = ggplot2::margin(10, 10, 10, 10, unit = "mm"))
  
  return(p1)
  
}

.auto_detect_chain = function(pep, pdb, k = 4, in_module = F){
  
  # Changed to coverage instead of jaccard
  kmer_coverage <- function(pdb_seq, msa_seq) {
    # seq is too short for kmer - just return 0
    if (nchar(pdb_seq) < k || nchar(msa_seq) < k) return(0)
    
    pdb_kmers <- substring(pdb_seq, 1:(nchar(pdb_seq) - k + 1), k:(nchar(pdb_seq)))
    msa_kmers <- substring(msa_seq, 1:(nchar(msa_seq) - k + 1), k:(nchar(msa_seq)))
    
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

# PDB UTILS ----
#' Standardize PDB Input
#'
#' Reads a PDB or CIF file and returns a trimmed structure for specified chains.
#'
#' @param pdb A PDB object, or path to pdb or mmCIF file.
#' @param chain Character vector of chain IDs, or \code{"all"}. # REMOVED 5/4
#'
#' @return A trimmed \code{bio3d} PDB object containing only selected chains.
#' @export
#' 
#6/6 -- i dont think we will handle chain here #
.standardize_pdb_input = function(pdb){
  
  # expects single entry #
  input_class = class(pdb)[1]
  
  if(input_class == 'character'){
    # !! add a check if file exists !! #
    # file.exists(msa)
    if(grepl('cif$', pdb)){
      pdb = bio3d::read.cif(pdb)
    } else {
      pdb = bio3d::read.pdb(pdb)
    }
  } else if(!input_class == 'pdb'){
    stop('NOT ONE OF THE TWO PDB OPTIONS')
  }
  
  # add residue id information #
  pdb$atom$insert[is.na(pdb$atom$insert)] = ''
  pdb$atom$residue_id = paste0(pdb$atom$resno, '_', pdb$atom$chain, '_', pdb$atom$insert)
  
  # pass along pdb #
  return(pdb)

}



#' Extract Sequence from PDB
#'
#' Retrieves the amino acid sequence from a PDB file or object.
#'
#' @param pdb A PDB object from \code{bio3d}.
#' @param chain Optional chain ID(s).
#' @param pdb_path Path to a PDB file.
#'
#' @return A named character vector of sequences, one per chain.
#' @export
.get_pdb_sequence = function(pdb, chain = NA, in_module = F){
  
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
      bio3d::pdbseq(bio3d::trim.pdb(pdb, chain = x)),
      collapse = '')
  }, USE.NAMES = T)

  return(aa_seq)
}



#' Residue-wise Distance Matrix
#'
#' Computes a pairwise distance matrix between residues based on 3D coordinates.
#'
#' @param pdb PDB object or NULL.
#' @param chain Chain ID(s).
#' @param pdb_path Path to a PDB file.
#' @param distance_method One of \code{'all'}, \code{'ca'}, \code{'backbone'}, or \code{'sidechain'}.
#'
#' @return A square numeric matrix of distances.
#' @export
.calculate_residue_distance = function(pdb, chain = NA, distance_method = 'all', in_module = F){

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
  res_dist <- pmin(res_dist, t(res_dist), na.rm = TRUE)
  diag(res_dist) <- 0

  # return residue wise distance matrix
  return(res_dist)

}



#' Calculate Solvent Accessibility
#'
#' Estimates residue-wise solvent accessibility using a DSSP rewrite.
#'
#' @param pdb A bio3d PDB object.
#' @param chain Chain ID(s).
#' @param pdb_path Path to the PDB file.
#' @param method One of \code{'rose'}, \code{'miller'}, \code{'theoretical_tien'}, or \code{'empirical_tien'}.
#' @param drop_incomplete Logical, drop residues missing backbone atoms.
#'
#' @return A data frame with residue indices, exposure values, and metadata.
#' @export
.calculate_accessibility = function(pdb, chain = NULL, method = 'rose', drop_incomplete = T, in_module = F){
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

  # Drop hydrogen atoms and non-ATOM records
  atoms = atoms[atoms$elesy != 'H' & atoms$type == 'ATOM', ]

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

#' Identify Surface Patches
#'
#' Defines residue patches around surface-exposed centroids based on RSA, SASA, and spatial distance.
#'
#' @param dist_mat A residue-residue distance matrix.
#' @param accessibility_df Data frame of residue accessibility values (e.g., from \code{calculate_accessibility}).
#' @param rsa_cutoff Minimum RSA for centroid selection (default = 0.1).
#' @param sasa_cutoff Optional minimum SASA for centroid selection.
#' @param dist_cutoff Maximum distance (in Å) for neighbors.
#' @param max_patch Optional maximum number of neighbors per patch.
#' @param only_exposed_in_patch If TRUE, restricts patch members to surface residues only.
#'
#' @return Updated \code{accessibility_df} with a \code{patch} column listing neighbors.
#' @export
.identify_patches = function(dist_mat, accessibility_df,
                            rsa_cutoff = 0.1, sasa_cutoff = NULL,
                            dist_cutoff = 15, max_patch = NULL,
                            only_exposed_in_patch = F){

  # returns patch around each residue #
  # dist_mat is residue_dist from calculate_distance_matrix
  # number of ways to define a window #
  # patch centroids defined by
  #     rsa_cutoff, or sasa_cutoff, or both NULL all residue
  # neighbors defined by
  #    dist_cutoff, or min_patch, or max_patch, or NULL all residues
  #    further filtered by only_exposed_in_patch

  # !! THERE ARE OTHER METHODS TO ADD -- COULD SEPERATE THESE PATCHERS
  # INTO DIFFERENT FUNCTIONS. OR Supply centroid_method = rsa, sasa, none
  # acc_cutoff. Window method = dist, min_patch, max_patch, none #


  # expect dist_mat and accessibility_df to have same amount of residues
  # could be different if residues dropped in calculate_accessibility()
  # dont really need this res_id checking ***
  res_ids = intersect(rownames(dist_mat), accessibility_df$residue_id)

  # get centroids
  if(!is.null(rsa_cutoff)){
    centroids = accessibility_df[accessibility_df$rsa >= rsa_cutoff &
                                   accessibility_df$residue_id %in% res_ids,]
  } else if(!is.null(sasa_cutoff)){
    centroids = accessibility_df[accessibility_df$sasa >= sasa_cutoff &
                                   accessibility_df$residue_id %in% res_ids,]
  } else {
    centroids = accessibility_df[accessibility_df$residue_id %in% res_ids,]
  }

  # shrink dist mat to only exposed residues if desired
  if(only_exposed_in_patch){
    in_set = centroids$residue_id
    dist_mat = dist_mat[in_set, in_set]
  }

  # get neighbors (in dist_cut and also check max_patch)
  centroids$patch = NA

  for(i in 1:nrow(centroids)){
    center = centroids$residue_id[i]

    # get neighbors
    neighbors = dist_mat[center,
                         names(which(dist_mat[center,] <= dist_cutoff))]

    # check max_patch
    if(!is.null(max_patch) && length(neighbors) > max_patch){
      # sort neighbors and take the first max_patch
      neighbors = sort(neighbors)
      neighbors = neighbors[1:max_patch]
    }

    # add to centroids
    centroids$patch[i] = paste0(names(neighbors), collapse = '+')
  }

  # add back to accessibility_df
  accessibility_df$patch = NA

  # join back
  accessibility_df$patch[match(centroids$residue_id, accessibility_df$residue_id)] = centroids$patch

  #accessibility_df$patch_method = paste0('-rsa ', rsa_cutoff, ' -sasa ', sasa_cutoff,
  #                                       ' -dist ', dist_cutoff, ' -max ', max_patch,
  #                                       ' -only_exposed ', only_exposed_in_patch)
  #

  return(accessibility_df)
}


#' Identify Interface Contacts
#'
#' Finds epitope and paratope residues from a PDB structure using chain-based atom distances.
#'
#' @param pdb A \code{bio3d} PDB object.
#' @param ag_chain Antigen chain ID.
#' @param h_chain Heavy chain ID.
#' @param l_chain Light chain ID.
#' @param dist_cutoff Maximum Å distance for contact (default = 5).
#'
#' @return A list with:
#' \item{epitope}{Concatenated residue IDs (e.g., "35_A_ 42_A_") of antigen residues near antibody.}
#' \item{paratope_h}{Heavy chain contact residues.}
#' \item{paratope_l}{Light chain contact residues.}
#' \item{contacts}{Matrix of all atom-level contacts.}
#' @export
.identify_both_interface = function(pdb, chain = NULL, interface_chain = NULL, dist_cutoff = 5){
  
  # remove H and HETATM
  pdb = bio3d::trim.pdb(pdb, 'protein')
  pdb = bio3d::trim.pdb(pdb, 'noh')
  
  # Grab chains of interest #
  ag_c = ag_chain
  h_c = h_chain
  l_c = l_chain
  
  pdb = bio3d::trim.pdb(pdb, chain = c(ag_c, h_c, l_c))
  
  pdb$atom$extra = paste(pdb$atom$resno, pdb$atom$chain, pdb$atom$insert, sep = '_')
  
  dist_mat = bio3d::dm.xyz(pdb$xyz, grpby = pdb$atom[,'extra'])
  
  colnames(dist_mat) = rownames(dist_mat) = unique(pdb$atom$extra)
  
  # grab close contacts
  close = which(dist_mat <= dist_cutoff, arr.ind = T)
  # replace pos with names
  close[, 1] = rownames(dist_mat)[close[, 1]]
  close[, 2] = colnames(dist_mat)[as.numeric(close[, 2])]
  
  # filter out heavy and light from column 1
  # filter out epi from column 2
  # ** probably break if l_c or h_c missing
  close = close[-grep(paste0('_', h_c, '_'), close[, 1]),]
  close = close[-grep(paste0('_', l_c, '_'), close[, 1]),]
  close = close[-grep(paste0('_', ag_c, '_'), close[, 2]),]
  rownames(close) = NULL
  
  # grab epi and paratope
  epi = paste0(unique(close[,1]), collapse = '+')
  para_h = paste0(close[grepl(paste0('_', h_c, '_'), close[, 2]), 2], collapse = '+')
  para_l = paste0(close[grepl(paste0('_', l_c, '_'), close[, 2]), 2], collapse = '+')
  
  return(list(
    epitope = epi, 
    paratope_h = para_h, 
    paratope_l = para_l,
    contacts = close)
  )
}

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


#' Extract Surface Patches from PDB
#'
#' High-level wrapper for computing sequences, distances, solvent accessibility, and patches from a PDB structure.
#'
#' @param pdb A \code{bio3d} PDB object.
#' @param chain Chain ID(s) or \code{"all"}.
#' @param pdb_path Path to a PDB or CIF file.
#' @param distance_method One of \code{"all"}, \code{"ca"}, \code{"backbone"}, \code{"sidechain"}.
#' @param drop_incomplete_residue Logical. Drop residues with incomplete backbone atoms.
#' @param rsa.method Method for RSA normalization (e.g., \code{"rose"}).
#' @param patch.dist.cutoff Max Å distance between residues in a patch.
#' @param patch.rsa.cutoff RSA cutoff for defining patch centroids.
#' @param patch.sasa.cutoff Optional SASA cutoff for patch centroids.
#' @param patch.only.exposed Logical, restrict patches to exposed residues.
#' @param max.patch Optional max number of neighbors per patch.
#'
#' @return A list with \code{pdb}, \code{seq_set}, \code{residue_dist}, and \code{residue_df} (with patches).
#' @export
pdb_to_patch = function(pdb, chain = NA, interface_chain = NA, occlusion_chain = NA,
                        distance_method = 'all',
                        drop_incomplete_residue = T, rsa.method = 'rose',
                        patch.dist.cutoff = 15, patch.rsa.cutoff = 0.1,
                        patch.sasa.cutoff = NA, patch.only.exposed = T,
                        max.patch = NA, interface.dist.cutoff = 5){
  
  # check if interface or occlusion is entered as 'all'
  #if(!is.null(interface_chain) && interface_chain == 'all'){
  #  interface_chain = NULL
  #}
  
  # split chain and occlusion chain (they dont need grouped)
  if(!(length(chain) == 1 && is.na(chain))) {
    chain <- unlist(strsplit(chain, split = ''))
  }
  
  if(!(length(occlusion_chain) == 1 && is.na(occlusion_chain))) {
    occlusion_chain <- unlist(strsplit(occlusion_chain, split = ''))
  }
  
  # step 0: validate pdb and return pdb object
  pdb = .standardize_pdb_input(pdb)

  # step 1: retrieve sequences for chains of interest
  seq_set = .get_pdb_sequence(pdb, chain = chain)

  # step 2: calculate residue-wise distance matrix
  residue_dist = .calculate_residue_distance(pdb, chain = chain,
                                            distance_method = distance_method,
                                            in_module = T)
  
  # step 3: calculate residue-wsie accessibility
  chain_set = c(chain, occlusion_chain)
  chain_set = chain_set[!is.na(chain_set)]
  residue_df = .calculate_accessibility(pdb, chain = chain_set,
                                       drop_incomplete = drop_incomplete_residue,
                                       method = rsa.method,
                                       in_module = T)

  # step 4: identify surface patches (expands residue_df)
  #print('Step 4: Identifying surface patches')
  residue_df = .identify_patches(residue_dist,
                                residue_df, only_exposed_in_patch = patch.only.exposed,
                                dist_cutoff = patch.dist.cutoff,
                                rsa_cutoff = patch.rsa.cutoff,
                                sasa_cutoff = patch.sasa.cutoff)
  
  # step 5: capture interface patches #
  if(!(length(interface_chain) == 1 && is.na(interface_chain))){
    # apply identify interface to all sets of interface chains #
    for(i in 1:length(interface_chain)){
      int_ch = unlist(strsplit(interface_chain[i], split = ''))
      interface_patches = .identify_interface(pdb, chain = chain, interface_chain = int_ch, dist_cutoff = interface.dist.cutoff)
      
      # add to residue df #
      residue_df[nrow(residue_df)+1,] = NA
      residue_df$residue_id[nrow(residue_df)] = interface_patches$name
      residue_df$patch[nrow(residue_df)] = interface_patches$interf
    }
    
  }

  # return list object
  return(list(
    pdb = pdb,
    chain = chain,
    seq_set = seq_set,
    residue_dist = residue_dist,
    residue_df = residue_df
  ))

}
