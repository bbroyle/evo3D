# ------------------------------------------ #
# EXTENSION_UTILS.R
# utilities for incorporating multiple msa and pdbs
# also special logic for homomultimers so each codon recieves one patch
# Brad Broyles
# ------------------------------------------ #

# extend_pdb_homomultimer() ----

#' Extend Patch Definitions Across Homomultimer Chains
#'
#' For homomultimeric proteins, identical codons may occur on multiple chains. This function merges
#' residue patches that map to the same codon number across chains, ensuring each codon position is
#' associated with a single patch. It updates MSA subsets accordingly and consolidates metadata.
#'
#' @param aln_info_sets A named list of \code{aln_df} outputs (e.g., from \code{aln_msa_to_pdb()}), one per chain or structure.
#' @param msa_info Output of \code{WRAPPER_msa_to_ref()}, containing the full codon-aligned nucleotide MSA matrix.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{msa_subsets}}{A named list of updated MSA patch subsets (post-merging).}
#'   \item{\code{aln_df}}{A combined residue-level metadata table after patch unification.}
#' }
#' @export
extend_pdb_homomultimer = function(aln_info_sets, msa_info){

  # just doing union patches #
  # step 0 -- set up patch joining methods (just union for now) # ----
  union_codon_patches <- function(patch1, patch2) {
    # split by '+', get unique numbers, sort, and rejoin
    if(is.na(patch1)) patch1 = ''
    if(is.na(patch2)) patch2 = ''
    nums1 = as.numeric(strsplit(patch1, "\\+")[[1]])
    nums2 = as.numeric(strsplit(patch2, "\\+")[[1]])
    union_nums = sort(unique(c(nums1, nums2)))
    paste(union_nums, collapse = "+")
  }

  # expand codon patches ----

  # finding residue_ids that span chains #
  # rbind all these aln_info_sets together #
  set_names = names(aln_info_sets)  # get all set names dynamically
  dfs = lapply(set_names, function(x) aln_info_sets[[x]]$aln_df)
  combined_df = do.call(rbind, dfs)

  # keep track of same codon different chain #
  # for each codon -- find its residue_ids #
  codons = combined_df[!is.na(combined_df$codon), 'codon']

  resi_list = list()
  for(codon in unique(codons)){
    # get all residue_ids for this codon #
    resi_ids = combined_df$residue_id[combined_df$codon == codon]
    # remove NA #
    resi_ids = resi_ids[!is.na(resi_ids)]

    # if more than one, save it #
    if(length(resi_ids) > 1){
      codon_name = paste('codon_', codon, sep = '')
      resi_list[[codon_name]] = resi_ids
    }
  }

  # which residue_id have multiple entries #
  #id_count = table(combined_df$residue_id)
  #multi_entry_ids = names(id_count[id_count > 1])
  # drop '-'
  #multi_entry_ids = multi_entry_ids[multi_entry_ids != '-']

  # for each of these loop through and union the patches #
  # save new patch, save new msa, remove extra rows #
  for(codon in names(resi_list)){

    # pull all residues at this codon #
    resi = resi_list[[codon]]

    # update residue_id with 2+ residues #
    # nothing to update but residue_id and pdb_aa #
    real_codon = gsub('codon_', '', codon)
    ro = which(combined_df$codon == real_codon)

    new_residue_id = paste(resi, collapse = '+')
    new_pdb_aa = combined_df$pdb_aa[ro]
    new_pdb_aa = new_pdb_aa[!is.na(new_pdb_aa)]
    new_pdb_aa = paste(new_pdb_aa, collapse = '+')

    # filter out '-' there wont be a patch #
    resi = resi[resi != '-']

    # grab rows #
    if(length(resi) == 0){

      combined_df$residue_id[ro[1]] = new_residue_id
      combined_df$pdb_aa[ro[1]] = new_pdb_aa

      # remove extra rows #
      combined_df = combined_df[-ro[-1],]
    } else {
      # actually update patch and msa subset #
      ro = which(combined_df$residue_id %in% resi)
      ro = unique(c(ro, which(combined_df$codon == real_codon)))

      # get all patches for this residue_id #
      patches = combined_df$codon_patch[ro]

      # union them #
      combined_patch = patches[1]
      for(i in 2:length(patches)){
        combined_patch = union_codon_patches(combined_patch, patches[i])
      }

      # patch built -- update row with codon as msa_subset_id #
      msa_subset_id = combined_df$msa_subset_id[ro]
      update_pos = which(grepl('^codon', msa_subset_id))[1] # if its 0 -- then its interface just do first one

      if(is.na(update_pos)) update_pos = 1 # case of interface #
      combined_df$codon_patch[ro[update_pos]] = combined_patch

      # lets grab new subset -- should do outside of loop (later) #
      new_msa = .extract_msa_subsets(msa_info$msa_mat, combined_df[ro[update_pos], , drop = FALSE])

      # replace in aln_info_sets1
      aln_info_sets[[1]]$msa_subsets[[msa_subset_id[update_pos]]] = new_msa[[1]]

      # update residue_id and pdb_aa #
      combined_df$residue_id[ro[update_pos]] = new_residue_id
      combined_df$pdb_aa[ro[update_pos]] = new_pdb_aa

      # remove the rest of rows for this msa_subset_id #
      combined_df = combined_df[-ro[-update_pos], ]

      # remove reside_X entries from aln_set_infos[[1]]$msa_subsets #
      resi_names = paste('residue_', resi, sep = '')
      aln_info_sets[[1]]$msa_subsets[resi_names] = NULL

    }

  }

  # save any extra interface msa_subsets #
  set1_names = names(aln_info_sets[[1]]$msa_subsets)

  set_len = length(aln_info_sets)
  for(i in 2:set_len){
    # get names of msa_subsets in this set #
    seti_names = names(aln_info_sets[[i]]$msa_subsets)

    # just looking for interfaces
    seti_names = seti_names[grepl('^interface', seti_names)]

    # see if interface doesnt exist in set1 #
    seti_names = seti_names[!seti_names %in% set1_names]

    # if there are any move over to aln_info_sets[[1]]$msa_subsets #
    if(length(seti_names) > 0){

      # create new entries
      aln_info_sets[[1]]$msa_subsets[seti_names] = aln_info_sets[[i]]$msa_subsets[seti_names]

    }

  }

  # combined_df is alignment -- # still need to add residues ...
  # aln_info_set[[1]] is msa_subset #

  return(list(
    msa_subsets = aln_info_sets[[1]]$msa_subsets,
    aln_df = combined_df
  ))

}


# extend_pdb() ----
#' Merge Codon Patch Definitions Across Structures
#'
#' This function combines patch metadata from two aligned PDB structures (or chain combinations)
#' that share a reference MSA. For codons present in both inputs, it performs a union of the
#' patch windows and updates the corresponding MSA subsets. It supports recursive application
#' to successively integrate additional models.
#'
#' @param aln_info_set1 A named list from \code{aln_msa_to_pdb()} or previous \code{extend_pdb()} output.
#' @param aln_info_set2 A new alignment info set (output of \code{aln_msa_to_pdb()} for a new structure or chain).
#' @param msa_info Output of \code{WRAPPER_msa_to_ref()}, containing the full codon-aligned MSA matrix.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{msa_subsets}}{A named list of updated MSA patch subsets, merged across structures.}
#'   \item{\code{aln_df}}{A unified residue-level metadata table including all merged structures.}
#' }
#'
#' @details
#' This function allows progressive patch merging across structural models by tagging each new set
#' of PDB columns (e.g., \code{pdb2_residue_id}, \code{pdb3_residue_id}, etc.). Interface and orphan
#' rows are preserved. Patch subsets from \code{aln_info_set2} are added to the master list unless
#' already present.
#'
#' @export
extend_pdb = function(aln_info_set1, aln_info_set2, msa_info){

  # take in two aln_info sets that share an underlying msa #
  # -- look for codon_patches built on shared codon #
  # -- constructu union of codon_patches and repull msa subsets #
  # -- save orphan (touches this msa, but codon is from other chain) #
  # -- save interface #
  # this function can be called recursivly,
  # but aln_info_set2 is always new (from aln_msa_to_pdb)

  # step 0 -- set up patch joining methods (just union for now) # ----
  union_codon_patches <- function(patch1, patch2) {
    # split by '+', get unique numbers, sort, and rejoin
    nums1 = as.numeric(strsplit(patch1, "\\+")[[1]])
    nums2 = as.numeric(strsplit(patch2, "\\+")[[1]])
    union_nums = sort(unique(c(nums1, nums2)))
    paste(union_nums, collapse = "+")
  }

  # step 1 -- see if aln_info_set1 is previous output from extend_pdb() -- track pdb number # ----
  extended_aln_df = aln_info_set1$aln_df
  additional_aln_df = aln_info_set2$aln_df

  cols = colnames(extended_aln_df)
  if('residue_id' %in% cols){
    # this is output from aln_msa_to_pdb() #
    next_pdb = 2

    # update column names #
    colnames(extended_aln_df)[which(cols == 'residue_id')] = 'pdb1_residue_id'
    colnames(extended_aln_df)[which(cols == 'pdb_aa')] = 'pdb1_aa'

    # since this is first extend_pdb() run lets add pdb1 tag to residue and interface codon subsets #
    msa_names = names(aln_info_set1$msa_subsets)
    pos = grep('^residue|^interface', msa_names)
    if(length(pos) > 0) {
      names(aln_info_set1$msa_subsets)[pos] = paste0(msa_names[pos], '_pdb1')
    }

    # NEW ADDITIONS #
    # UPDATE msa_subset_id for interfaces/residues to match the renamed msa_subsets
    interface_rows = grep('^interface_', extended_aln_df$msa_subset_id)
    residue_rows = grep('^residue_', extended_aln_df$msa_subset_id)

    if(length(interface_rows) > 0) {
      extended_aln_df$msa_subset_id[interface_rows] = paste0(extended_aln_df$msa_subset_id[interface_rows], '_pdb1')
    }
    if(length(residue_rows) > 0) {
      extended_aln_df$msa_subset_id[residue_rows] = paste0(extended_aln_df$msa_subset_id[residue_rows], '_pdb1')
    }

  } else {
    # this is output from previous extend_pdb() -- what count are we on #
    pdb_cols = cols[grep("^pdb[0-9]+_res", colnames(extended_aln_df))]
    next_pdb = length(pdb_cols) + 1

    # no need to add _pdbx tag -- should already be present from previous runs #
  }

  # lets update additional_aln_df
  cols = colnames(additional_aln_df)
  colnames(additional_aln_df)[which(cols == 'residue_id')] = paste0('pdb', next_pdb, '_residue_id')
  colnames(additional_aln_df)[which(cols == 'pdb_aa')] = paste0('pdb', next_pdb, '_aa')

  # add pdbx tag to msa_subsets
  msa_names = names(aln_info_set2$msa_subsets)
  pos = grep('^residue|^interface', msa_names)
  if(length(pos) > 0) {
    names(aln_info_set2$msa_subsets)[pos] = paste0(msa_names[pos], '_pdb', next_pdb)
  }

  # UPDATE msa_subset_id for interfaces/residues in additional_aln_df (add to previous section)
  interface_rows = grep('^interface_', additional_aln_df$msa_subset_id)
  residue_rows = grep('^residue_', additional_aln_df$msa_subset_id)

  if(length(interface_rows) > 0) {
    additional_aln_df$msa_subset_id[interface_rows] = paste0(additional_aln_df$msa_subset_id[interface_rows], '_pdb', next_pdb)
  }
  if(length(residue_rows) > 0) {
    additional_aln_df$msa_subset_id[residue_rows] = paste0(additional_aln_df$msa_subset_id[residue_rows], '_pdb', next_pdb)
  }

  # step 2 -- hold on to extra data (orphans and interfaces) to add after merge ----

  # store extra info that will not be merged #
  extra1 = extended_aln_df[is.na(extended_aln_df$codon),]
  extra2 = additional_aln_df[is.na(additional_aln_df$codon),]

  # only process extras if we actually have some
  if(nrow(extra1) > 0 || nrow(extra2) > 0) {
    # get all unique columns from both dataframes
    all_cols = unique(c(names(extra1), names(extra2)))

    # add missing columns to extra1 (will be NA)
    if(nrow(extra1) > 0) {
      missing_in_extra1 = setdiff(all_cols, names(extra1))
      extra1[missing_in_extra1] = NA
    }

    # add missing columns to extra2 (will be NA)
    if(nrow(extra2) > 0) {
      missing_in_extra2 = setdiff(all_cols, names(extra2))
      extra2[missing_in_extra2] = NA
    }

    # now rbind them
    extra_combined = rbind(extra1, extra2)
  } else {
    extra_combined = NULL
  }

  # remove these extras just from extended_aln_df #
  extended_aln_df = extended_aln_df[!is.na(extended_aln_df$codon),]
  additional_aln_df = additional_aln_df[!is.na(additional_aln_df$codon),]

  # step 3 -- build union pdb patches # ----

  # easy enough just find codon positions that match and have data (codon_patch in both) #
  # 6/14/25 -- why not check extended_aln_df for codon_patch?
  codons1 = extended_aln_df$codon[which(!is.na(extended_aln_df$codon))]
  codons2 = additional_aln_df$codon[which(!is.na(additional_aln_df$codon) & !is.na(additional_aln_df$codon_patch))]
  shared_codons = intersect(codons1, codons2)

  # so we need to extend patches for these shared codons -- take union of two codon_patches #
  updated_codons = c()

  for(i in 1:length(shared_codons)){
    # get the patches for this codon in both dataframes
    patch1 = extended_aln_df$codon_patch[which(extended_aln_df$codon == shared_codons[i])]
    patch2 = additional_aln_df$codon_patch[which(additional_aln_df$codon == shared_codons[i])]

    # handle different cases
    if(is.na(patch1)) {
      # patch1 is NA, patch2 has data (guaranteed by codons2 filter) - use patch2
      new_patch = patch2
      updated_codons = c(updated_codons, shared_codons[i])
    } else {
      # both have data - take union
      new_patch = union_codon_patches(patch1, patch2)
      if(new_patch != patch1) {
        updated_codons = c(updated_codons, shared_codons[i])
      } else {
        next  # no change needed
      }
    }

    # update the extended_aln_df with the new patch
    extended_aln_df$codon_patch[which(extended_aln_df$codon == shared_codons[i])] = new_patch
  }

  # merge new pdb and pdb aa columns into extended_aln_df #
  # merge on codon, keeping all rows from both # what will happen if pdb introduce gap in msa 6/14/25?
  extended_aln_df = merge(extended_aln_df, additional_aln_df[, !names(additional_aln_df) %in% c("codon_patch", "ref_aa", 'msa_subset_id')], by = "codon", all = TRUE)
  extended_aln_df = extended_aln_df[order(as.numeric(extended_aln_df$codon)),]

  # step 4 -- rebuild msa_subsets ----
  if(length(updated_codons) > 0){
    ro = which(extended_aln_df$codon %in% updated_codons)
    msa_subsets = .extract_msa_subsets(msa_info$msa_mat, extended_aln_df[ro,])

    # replace or write to aln_info_set1$msa_subsets #
    aln_info_set1$msa_subsets[names(msa_subsets)] = msa_subsets
  }

  # propogate any bonus sets from aln_info_set2 to aln_info_set1 #
  set1 = names(aln_info_set1$msa_subsets)
  set2 = names(aln_info_set2$msa_subsets)
  bonus_sets = set2[!set2 %in% set1]

  aln_info_set1$msa_subsets[bonus_sets] = aln_info_set2$msa_subsets[bonus_sets]

  # step 5 -- add extra back to extended_aln_df ----
  if(!is.null(extra_combined)) {
    extended_aln_df = rbind(extended_aln_df, extra_combined)
  }

  # lets reearrange columns before return #
  col_order = c('codon', 'ref_aa', 'msa_subset_id',
                grep('pdb[0-9]+_residue', colnames(extended_aln_df), value = TRUE),
                grep('pdb[0-9]+_aa', colnames(extended_aln_df), value = TRUE),
                'codon_patch'
  )

  col_order = c(col_order, setdiff(colnames(extended_aln_df), col_order))
  extended_aln_df = extended_aln_df[, col_order]

  # return
  return(list(msa_subsets = aln_info_set1$msa_subsets,
              aln_df = extended_aln_df))
}

# debug #
#aln_info_set1 = extended_result
#aln_info_set2 = working_aln_sets[[3]]
#msa_info = msa_info_sets[[msa_id]]

# extend_msa() ----

#' Merge Patch MSAs Across Structural Chains or Regions
#'
#' This function extends residue- and patch-level alignment data across independently aligned chains
#' or MSA regions that share a common codon reference. It merges patch metadata and corresponding
#' MSA subsets across input alignment sets, resolving codon patches, rescuing orphaned mappings, and
#' concatenating per-codon MSA slices across chains or models.
#'
#' @param aln_info_set1 A named list from \code{aln_msa_to_pdb()} or previous \code{extend_msa()} output.
#' @param aln_info_set2 A second alignment info set from another chain or region (also from \code{aln_msa_to_pdb()} or \code{extend_pdb()}).
#' @param msa_info_set A named list of full MSA matrices from \code{WRAPPER_msa_to_ref()}, containing \code{msa1}, \code{msa2}, etc.
#' @param use_sample_names Logical; if \code{TRUE}, sequences will be matched by FASTA header during MSA merging. If \code{FALSE}, rows are merged by order.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{msa_subsets}}{A combined list of MSA patch subsets, merged across MSAs and structural chains.}
#'   \item{\code{aln_df}}{A residue-level metadata table that includes merged patch and codon information across inputs.}
#' }
#'
#' @details
#' This function enables chain-aware patch merging by propagating patch definitions across residue matches
#' and codon positions shared across independently aligned structures. Codon-level MSA subsets are updated
#' to reflect the union of observed residues across chains. Interface rows are preserved, and orphans
#' (residues with MSA mapping but no associated patch) are rescued using cross-reference with the alternate chain.
#'
#' Repeated calls to \code{extend_msa()} allow recursive merging of additional chains or MSA blocks.
#'
#' @seealso \code{\link{extend_pdb}}, \code{\link{aln_msa_to_pdb}}, \code{\link{WRAPPER_msa_to_ref}}
#' @export
extend_msa = function(aln_info_set1, aln_info_set2, msa_info_set, use_sample_names = TRUE) {

  # take in two aln_info sets that share an underlying pdb structure #
  # -- look for shared residue_ids across pdb columns #
  # -- if msa_info provided: rescue orphans (codon_patch exists, codon is NA) #
  # -- concatenate msa_subsets for matching residue positions #
  # -- preserve interfaces and handle multi-pdb column structure #
  # This function can be used recursively and can be used on outputs from aln_msa_to_pdb() or extend_pdb() #
  # both those dataframes have one codon column -- and no need to track msa #
  # first thing we will do here is add msa column so we can append other msas #

  # -- step 1 ~ check if this is first extend_msa() call or recursive and adding msa tags # ----
  extended_aln_df = aln_info_set1$aln_df
  additional_aln_df = aln_info_set2$aln_df

  # Check if this is first extend_msa() call or recursive
  cols = colnames(extended_aln_df)
  if(!'msa' %in% cols){
    # this is output from aln_msa_to_pdb() or extend_pdb() #
    next_msa = 2

    # add msa column to track source #
    extended_aln_df$msa = 1

    # add msa1 tag to codon msa subsets #
    msa_names = names(aln_info_set1$msa_subsets)
    pos = grep('^codon', msa_names)
    if(length(pos) > 0) {
      names(aln_info_set1$msa_subsets)[pos] = paste0('msa1_', msa_names[pos])
    }

    # add msa1 tag to msa_subset_id in extended_aln_df #
    codon_rows = grep('^codon_', extended_aln_df$msa_subset_id)
    if(length(codon_rows) > 0) {
      extended_aln_df$msa_subset_id[codon_rows] = paste0('msa1_', extended_aln_df$msa_subset_id[codon_rows])
    }

  } else {
    # this is output from previous extend_msa() -- what count are we on #
    max_msa = max(extended_aln_df$msa, na.rm = TRUE)
    next_msa = max_msa + 1

  }

  # lets update additional_aln_df ( always expected to be from extend_pdb() or from aln_msa_to_pdb() )
  additional_aln_df$msa = next_msa

  # add msa tag to msa_subsets
  msa_names = names(aln_info_set2$msa_subsets)
  pos = grep('^codon', msa_names)
  if(length(pos) > 0) {
    names(aln_info_set2$msa_subsets)[pos] = paste0('msa', next_msa, '_', msa_names[pos])
  }

  # add msa tag to msa_subset_id in additional_aln_df #
  codon_rows = grep('^codon_', additional_aln_df$msa_subset_id)
  if(length(codon_rows) > 0) {
    additional_aln_df$msa_subset_id[codon_rows] = paste0('msa', next_msa, '_', additional_aln_df$msa_subset_id[codon_rows])
  }

  # all patches in these two data frames are from $msa
  extended_aln_df$patch_msa = extended_aln_df$msa
  additional_aln_df$patch_msa = additional_aln_df$msa

  # store interfaces (which are not tied to codons) we will merge on common pdbX_residue id # ----
  cols1 = colnames(extended_aln_df)[grep("residue_id$", colnames(extended_aln_df))]
  cols2 = colnames(additional_aln_df)[grep("residue_id$", colnames(additional_aln_df))]
  int1 = which(rowSums(sapply(cols1, function(col) grepl("^interface", extended_aln_df[[col]]))) > 0)
  int2 = which(rowSums(sapply(cols2, function(col) grepl("^interface", additional_aln_df[[col]]))) > 0)

  interfaces = rbind(extended_aln_df[int1, , drop = FALSE],
                     additional_aln_df[int2, , drop = FALSE])

  if(length(int1) > 0) {
    extended_aln_df = extended_aln_df[-c(int1),]
  }

  if(length(int2) > 0) {
    additional_aln_df = additional_aln_df[-c(int2),]
  }


  # step 2 -- extend unmapped residues (orphans) across pdbs ----

  # store those codon and interface rows #
  ready1 = which(!is.na(extended_aln_df$codon))
  ready2 = which(!is.na(additional_aln_df$codon))

  # Check for orphans
  has_orphans1 = (nrow(extended_aln_df) > length(ready1))
  has_orphans2 = (nrow(additional_aln_df) > length(ready2))

  if(has_orphans1){

    # WE WILL MERGE ORPHANS ACROSS PEBS USING THEIR CODON POSITION #
    # THESE ADDITIONS NEED TO BE STORED TO LATER ADD TO extended_aln_df #
    # SOME ORPHANS MAP TO DIFFERENT MSA SET -- in this case they should be saved for later #

    # grab orphaned data
    orphan_df = extended_aln_df[-ready1, ]
    orphan_df$track_id = 1:nrow(orphan_df)  # add track_id to keep track of orphans
    used_ids = c()

    # look for residue id matches in proper pdb column -- save patch to that row #
    codon_map = additional_aln_df[!is.na(additional_aln_df$codon), c("codon", "msa", cols2, 'msa_subset_id', 'patch_msa')]
    codon_map$patch_msa = NA
    filled_codon_map1 = codon_map[0,]

    # for each column in cols 1 grab its codon_patch and see if we can find its residue in that column in codon map
    msa_save1 = list()
    for(msa_num in unique(orphan_df$msa)) {

      # filter orphan_df for current msa_id
      orphan_sub = orphan_df[orphan_df$msa == msa_num, ]

      # store codon_map for current msa_num #
      hold_codon_map = codon_map

      # builds codon_map for current msa_num -- across pdbs #
      # aka residue touched msa1 and now we want its msa2 codon #
      for (col in cols1) {
        if (!col %in% cols2) next # skip if col not in additional_aln_df #

        # find residue_id in codon_map #
        # i dont think != '-' is needed here? how can a gap have a patch
        df = orphan_sub[!is.na(orphan_sub[[col]]) & orphan_sub[[col]] != '-', c(col, 'codon_patch', 'track_id')]
        if (nrow(df) == 0){
          next
        }

        # keep track of which will actually merge in #
        will_merge = which(df[[col]] %in% hold_codon_map[[col]])
        used_ids = c(used_ids, df$track_id[will_merge])
        # drop track id
        df$track_id = NULL

        colnames(df)[2] = paste0(gsub('[_]*residue_id', '', col), '_codon_patch')

        hold_codon_map = merge(hold_codon_map, df, by = col, all.x = TRUE)
      }

      # Keep rows where at least one codon_patch column has data
      patches = grep('codon_patch', colnames(hold_codon_map), value = TRUE)
      hold_codon_map = hold_codon_map[rowSums(sapply(patches, function(col) !is.na(hold_codon_map[[col]]))) > 0,]

      # now for each row -- do union patch # ##NOTE IF EVER ADDING INTERSECT PATCH TO extend_pdb() -- we would want to do it here too #
      hold_codon_map$codon_patch = NA
      for(i in 1:nrow(hold_codon_map)) {
        combine_patch = hold_codon_map[i, patches]
        combine_patch = combine_patch[!is.na(combine_patch)]
        combine_patch = paste(combine_patch, collapse = '+')
        combine_patch = unique(unlist(strsplit(combine_patch, '\\+')))

        # remove previous tags # - or do they need stored?
        combine_patch = gsub('_[0-9]+$', '', combine_patch)

        combine_patch = sort(as.numeric(combine_patch))

        hold_codon_map$codon_patch[i] = paste0(combine_patch, collapse = '+')
      }

      # now we have final patches for all the orphans of df1
      save_sub = .extract_msa_subsets(msa_info_set[[paste0('msa', msa_num)]]$msa_mat, hold_codon_map)
      names(save_sub) = paste0(names(save_sub), '_pulling_msa', msa_num)

      # save msas (grows across msas in extend_aln_df)
      msa_save1 = c(msa_save1, save_sub)

      # add hold_codon_map to filled_codon_map #
      hold_codon_map$patch_msa = msa_num
      filled_codon_map1 = rbind(filled_codon_map1, hold_codon_map[c('codon', 'msa', 'codon_patch', 'patch_msa', 'msa_subset_id')])

    }

    # which orphans were used and which are waiting for different msa? #
    still_orphan1 = orphan_df[!orphan_df$track_id %in% used_ids,]
    still_orphan1$track_id = NULL

    # may want to keep all orphan data -- and signal if it was mapped or not #

  }

  if(has_orphans2){

    orphan_df = additional_aln_df[-ready2, ]
    orphan_df$track_id = 1:nrow(orphan_df)  # add track_id to keep track of orphans
    used_ids = c()

    # this is msa extension but we need to pair residue_id from df1 to their df2 residue and codon #
    codon_map = extended_aln_df[!is.na(extended_aln_df$codon), c("codon", "msa", cols1, 'msa_subset_id', 'patch_msa')]  # msa 2 but we will use to find resi from msa1 df
    codon_map$patch_msa = next_msa
    #filled_codon_map2 = codon_map[0,] -- since all patch_msa source are msa_next - we can collect at the end #

    # for each column in cols 1 grab its codon_patch and see if we can find its residue in that column in codon map
    for(col in cols2){
      if(!col %in% cols1) next

      # find residue_id in codon_map #
      df = orphan_df[!is.na(orphan_df[[col]]) & orphan_df[[col]] != '-', c(col, 'codon_patch', 'track_id')]
      if(nrow(df) == 0) next

      # keep track of which will actually merge in #
      will_merge = which(df[[col]] %in% codon_map[[col]])
      used_ids = c(used_ids, df$track_id[will_merge])
      # drop track id
      df$track_id = NULL

      colnames(df)[2] = paste0(gsub('[_]*residue_id', '', col), '_codon_patch')

      codon_map = merge(codon_map, df, by = col, all.x = TRUE)
    }

    # now we have codon_map1 with codon_patch and residue_id from df1 #

    # Keep rows where at least one codon_patch column has data
    patches = grep('codon_patch', colnames(codon_map), value = TRUE)
    codon_map = codon_map[rowSums(sapply(patches, function(col) !is.na(codon_map[[col]]))) > 0,]

    # now for each row -- do union patch # ##NOTE IF EVER ADDING INTERSECT PATCH TO extend_pdb() -- we would want to do it here too #
    codon_map$codon_patch = NA
    for(i in 1:nrow(codon_map)) {
      combine_patch = codon_map[i, patches]
      combine_patch = combine_patch[!is.na(combine_patch)]
      combine_patch = paste(combine_patch, collapse = '+')
      combine_patch = unique(unlist(strsplit(combine_patch, '\\+')))
      combine_patch = sort(as.numeric(combine_patch))

      combine_patch = gsub('_[0-9]+$', '', combine_patch)

      codon_map$codon_patch[i] = paste0(combine_patch, collapse = '+')
    }

    # now go by msa_num in codon_map ~ to grab msa subsets #
    msa_save2 = list()
    for(msa_num in unique(codon_map$msa)){

      codon_sub = codon_map[codon_map$msa == msa_num, ]

      # now we have final patches for all the orphans of df1
      # save as msaX_codon_X_msa1 and tack on too where?
      save_sub = .extract_msa_subsets(msa_info_set[[paste0('msa', next_msa)]]$msa_mat, codon_sub)
      names(save_sub) = paste0(names(save_sub), '_pulling_msa', next_msa)

      msa_save2 = c(msa_save2, save_sub)
    }

    # which orphans were used and which are waiting for different msa? #
    still_orphan2 = orphan_df[!orphan_df$track_id %in% used_ids,]
    still_orphan2$track_id = NULL

    # need codon_map for additional patch info #
    filled_codon_map2 = codon_map[, c("codon", "msa", "codon_patch", "patch_msa", 'msa_subset_id')]
  }

  # step3 --- ready to rbind data frames and concatonate msa_subsets ----

  # at this point all data is ready to merge #
  # 1. extended_aln_df[ready1,] -- having codon X msa Y and its patch
  #    -- msa subsets are in aln_info_set1$msa_subsets (msaY_codonX)
  # 2. additional_aln_df[ready2,] -- having codox X msa Z and its patch
  #    -- msa subsets are in aln_info_set2$msa_subsets (msaZ_codonX)
  # 3. filled_codon_map1 -- having pdb residues that touch msa Y but come from msa Z
  #    -- msa subsets are in msa_save1
  # 4. filled_codon_map2 -- having pdb residues that touch msa Z but come from msa Y (Y can be multiple msa)
  #    -- msa subsets are in msa_save2
  # 5. still_orphan1 -- pdb residues that touch Y but have not found their home msa
  # 6. still_orphan2 -- pdb residues that touch Z but have not found their home msa

  # so 1 and 4 will go together (focal residue is codonX msaY)
  # and 2 and 3 will go together (focal residue is codonX msaZ)

  # update patch for these top 4 datasets -- their patch to msa pull is done #
  add_msa_tag = function(patch, tag){
    if(is.na(patch)) return(NA)
    if(grepl('_', patch)) return(patch)  # already has msa tag
    patch = strsplit(patch, '\\+')[[1]]
    patch = paste0(patch, '_', tag)
    patch = paste(patch, collapse = '+')
    return(patch)
  }

  clean1 = extended_aln_df[ready1,]
  clean1$codon_patch = mapply(add_msa_tag, clean1$codon_patch, clean1$patch_msa, USE.NAMES = FALSE)

  clean2 = additional_aln_df[ready2,]
  clean2$codon_patch = mapply(add_msa_tag, clean2$codon_patch, clean2$patch_msa, USE.NAMES = FALSE)

  filled_codon_map1$codon_patch = mapply(add_msa_tag, filled_codon_map1$codon_patch, filled_codon_map1$patch_msa, USE.NAMES = FALSE)
  filled_codon_map2$codon_patch = mapply(add_msa_tag, filled_codon_map2$codon_patch, filled_codon_map2$patch_msa, USE.NAMES = FALSE)

  interfaces$codon_patch = mapply(add_msa_tag, interfaces$codon_patch, interfaces$patch_msa, USE.NAMES = FALSE)

  # clean1 and clean2 have full codon to msa maps #
  # add codon_patch from filled2 to clean1
  for(i in 1:nrow(filled_codon_map2)) {
    codon = filled_codon_map2$codon[i]
    msa = filled_codon_map2$msa[i]
    patch = filled_codon_map2$codon_patch[i]
    patch_msa = filled_codon_map2$patch_msa[i]

    # find codon in clean1 and join this data #
    idx = which(clean1$codon == codon & clean1$msa == msa)
    if(length(idx) > 0) {
      clean1$codon_patch[idx] = paste(clean1$codon_patch[idx], patch, sep = '+')
    }

    # pdbs can also be concatonated #
    msa1 = aln_info_set1$msa_subsets[[paste0('msa', msa, '_codon_', codon)]]
    msa2 = msa_save2[[paste0('msa', msa, '_codon_', codon, '_pulling_msa', patch_msa)]]

    # if use_names_true -- then we concate on common name #
    # other wise just concatenate the two #
    if(use_sample_names) {
      names1 = rownames(msa1)
      names2 = rownames(msa2)
      reorder_msa2 = match(names1, names2)

      # see if any data is dropped #
      valid_match = !is.na(reorder_msa2)
      names1_valid = names1[valid_match]
      reorder_msa2_valid = reorder_msa2[valid_match]

      # print (later store) -- if any data is dropped #
      if(length(valid_match) != length(names1)){
        print('Warning: some data is dropped when merging msa subsets:')
        print('use_sample_names = TRUE ... so msa\'s are matched by fasta headers')
        print('this ensures valid msa subsetting for sequence based statistics')
        print('it could be turned off for single site stats like entropy -- but')
        print('you risk pulling gene information from entirely different genomic samples')
        print(paste0('dropped ', length(names1) - length(names1_valid), ' sequences'))
        print(paste0('for msa subset built on ', paste0('msa', msa, '_codon_', codon) ))
      }

      msa1 = msa1[names1_valid,]
      msa2 = msa2[reorder_msa2_valid,]

      # now we can merge the two #
      merged_msa = cbind(msa1, msa2)

      # overwrite the msa subset in aln_info_set1 #
      aln_info_set1$msa_subsets[[paste0('msa', msa, '_codon_', codon)]] = merged_msa
    } else {
      # just combine all available rows #
      min_ro = min(nrow(msa1), nrow(msa2))
      merged_msa = cbind(msa1[1:min_ro,], msa2[1:min_ro,])
      aln_info_set1$msa_subsets[[paste0('msa', msa, '_codon_', codon)]] = merged_msa
    }

  }

  # add codon_patch from filled1 to clean2 #should be seq_len()
  for(i in 1:nrow(filled_codon_map1)) {
    codon = filled_codon_map1$codon[i]
    msa = filled_codon_map1$msa[i]
    patch = filled_codon_map1$codon_patch[i]
    patch_msa = filled_codon_map1$patch_msa[i]

    # find codon in clean1 and join this data #
    idx = which(clean2$codon == codon & clean2$msa == msa)
    if(length(idx) > 0) {
      clean2$codon_patch[idx] = paste(clean2$codon_patch[idx], patch, sep = '+')
    }

    # pdbs can also be concatonated #
    msa1 = aln_info_set2$msa_subsets[[paste0('msa', msa, '_codon_', codon)]]
    msa2 = msa_save1[[paste0('msa', msa, '_codon_', codon, '_pulling_msa', patch_msa)]]

    # if use_names_true -- then we concate on common name #
    # other wise just concatenate the two #
    if(use_sample_names) {
      names1 = rownames(msa1)
      names2 = rownames(msa2)
      reorder_msa2 = match(names1, names2)

      # see if any data is dropped #
      valid_match = !is.na(reorder_msa2)
      names1_valid = names1[valid_match]
      reorder_msa2_valid = reorder_msa2[valid_match]

      # print (later store) -- if any data is dropped #
      if(length(valid_match) != length(names1)){
        print('Warning: some data is dropped when merging msa subsets:')
        print('use_sample_names = TRUE ... so msa\'s are matched by fasta headers')
        print('this ensures valid msa subsetting for sequence based statistics')
        print('it could be turned off for single site stats like entropy -- but')
        print('you risk pulling gene information from entirely different genomic samples')
        print(paste0('dropped ', length(names1) - length(names1_valid), ' sequences'))
        print(paste0('for msa subset built on ', paste0('msa', msa, '_codon_', codon) ))
      }

      msa1 = msa1[names1_valid,]
      msa2 = msa2[reorder_msa2_valid,]

      # now we can merge the two #
      merged_msa = cbind(msa1, msa2)

      # overwrite the msa subset in aln_info_set1 #
      aln_info_set2$msa_subsets[[paste0('msa', msa, '_codon_', codon)]] = merged_msa
    } else {
      # just combine all available rows #
      min_ro = min(nrow(msa1), nrow(msa2))
      merged_msa = cbind(msa1[1:min_ro,], msa2[1:min_ro,])
      aln_info_set2$msa_subsets[[paste0('msa', msa, '_codon_', codon)]] = merged_msa
    }

  }

  # step 4: merge interfaces -- check names too (needs work) ----
  unique_interface_ids = unique(interfaces$msa_subset_id)

  # Handle common interfaces between the two msa_subset lists
  names1 = names(aln_info_set1$msa_subsets)
  names2 = names(aln_info_set2$msa_subsets)

  # Find interfaces that exist in both
  common_interfaces = intersect(names1[grep('^interface_', names1)],
                                names2[grep('^interface_', names2)])

  # Combine MSA subsets for common interfaces
  for(interface_name in common_interfaces) {
    msa1 = aln_info_set1$msa_subsets[[interface_name]]
    msa2 = aln_info_set2$msa_subsets[[interface_name]]

    if(use_sample_names) {
      names1 = rownames(msa1)
      names2 = rownames(msa2)
      reorder_msa2 = match(names1, names2)

      valid_match = !is.na(reorder_msa2)
      names1_valid = names1[valid_match]
      reorder_msa2_valid = reorder_msa2[valid_match]

      if(length(valid_match) != length(names1)){
        print(paste0('Warning: dropped ', length(names1) - length(names1_valid), ' sequences for interface ', interface_name))
      }

      msa1 = msa1[names1_valid,]
      msa2 = msa2[reorder_msa2_valid,]

      combined_msa = cbind(msa1, msa2)
    } else {
      min_rows = min(nrow(msa1), nrow(msa2))
      combined_msa = cbind(msa1[1:min_rows,], msa2[1:min_rows,])
    }

    # Store combined version
    aln_info_set1$msa_subsets[[interface_name]] = combined_msa
  }

  # Collapse duplicate interfaces by msa_subset_id
  if(nrow(interfaces) > 0) {

    # Group by msa_subset_id and collapse
    unique_interface_ids = unique(interfaces$msa_subset_id)
    collapsed_interfaces = data.frame()

    for(interface_id in unique_interface_ids) {
      interface_rows = interfaces[interfaces$msa_subset_id == interface_id, ]

      if(nrow(interface_rows) == 1) {
        # Single row, just add it
        collapsed_row = interface_rows
      } else {
        # Multiple rows, combine patches
        collapsed_row = interface_rows[1, ]  # take first row as template

        # Combine all codon_patches
        all_patches = interface_rows$codon_patch
        all_patches = all_patches[!is.na(all_patches)]
        if(length(all_patches)>0){
          combined_patch = paste(all_patches, collapse = '+')
        } else {
          combined_patch = NA
        }
          collapsed_row$codon_patch = combined_patch
      }

      # Set msa and patch_msa to NA for interfaces
      collapsed_row$msa = NA
      collapsed_row$patch_msa = NA

      collapsed_interfaces = rbind(collapsed_interfaces, collapsed_row)
    }

    # Replace interfaces with collapsed version
    interfaces = collapsed_interfaces
  }

  # not all columns are present so use set_diff #
  all_cols = unique(c(colnames(clean1),
                      colnames(clean2),
                      colnames(interfaces),
                      colnames(still_orphan1),
                      colnames(still_orphan2)))

  # Add missing columns to each dataframe
  # if interfaces, still orphan1 or still orphan2

  if(nrow(still_orphan1) > 0) {
    still_orphan1[setdiff(all_cols, colnames(still_orphan1))] = NA
  } else{
    still_orphan1 = data.frame(matrix(NA, nrow = 0, ncol = length(all_cols)))
    colnames(still_orphan1) = all_cols
  }

  if(nrow(still_orphan2) > 0) {
    still_orphan2[setdiff(all_cols, colnames(still_orphan2))] = NA
  } else{
    still_orphan2 = data.frame(matrix(NA, nrow = 0, ncol = length(all_cols)))
    colnames(still_orphan2) = all_cols
  }

  if(nrow(interfaces) > 0) {
    interfaces[setdiff(all_cols, colnames(interfaces))] = NA
  } else{
    interfaces = data.frame(matrix(NA, nrow = 0, ncol = length(all_cols)))
    colnames(interfaces) = all_cols
  }

  if(nrow(clean1) > 0) {
    clean1[setdiff(all_cols, colnames(clean1))] = NA
  } else{
    clean1 = data.frame(matrix(NA, nrow = 0, ncol = length(all_cols)))
    colnames(clean1) = all_cols
  }

  if(nrow(clean2) > 0) {
    clean2[setdiff(all_cols, colnames(clean2))] = NA
  } else{
    clean2 = data.frame(matrix(NA, nrow = 0, ncol = length(all_cols)))
    colnames(clean2) = all_cols
  }

  #clean1[setdiff(all_cols, colnames(clean1))] = NA
  #clean2[setdiff(all_cols, colnames(clean2))] = NA
  #interfaces[setdiff(all_cols, colnames(interfaces))] = NA
  #still_orphan1[setdiff(all_cols, colnames(still_orphan1))] = NA
  #still_orphan2[setdiff(all_cols, colnames(still_orphan2))] = NA

  # rbind all dataframes together #
  full_set = rbind(clean1, clean2, interfaces, still_orphan1, still_orphan2)

  full_set = full_set[order(as.numeric(full_set$msa), as.numeric(full_set$codon)),]

  # reorder columns - include msa_subset_id #
  pdb_residue_cols = grep('pdb[0-9]+_residue', colnames(full_set), value = TRUE)
  pdb_aa_cols = grep('pdb[0-9]+_aa', colnames(full_set), value = TRUE)

  col_order = c('msa', 'codon', 'ref_aa', 'msa_subset_id',
                pdb_residue_cols, pdb_aa_cols,
                'codon_patch', 'patch_msa')

  col_order = c(col_order, setdiff(colnames(full_set), col_order))
  full_set = full_set[, col_order]

  # copy over msa sets #
  full_msa = aln_info_set1$msa_subsets
  full_msa = c(full_msa, aln_info_set2$msa_subsets)

  full_msa = full_msa[!grepl('^residue', names(full_msa))]  # remove residue msa's #

  # any interfaces (NA missed) #

  #

  return(
    list(
      msa_subsets = full_msa,
      aln_df = full_set
    )
  )

}

# indetify_epitopes() ----

#' Identify Antibody–Antigen Contacts
#'
#' Identifies epitope and paratope residues from a PDB structure based on inter-chain atom distances.
#'
#' @param pdb A \code{bio3d} PDB object.
#' @param ag_chain Chain ID of the antigen.
#' @param h_chain Chain ID of the antibody heavy chain.
#' @param l_chain Chain ID of the antibody light chain.
#' @param dist_cutoff Maximum distance (in Å) to define a contact. Default is 5.
#'
#' @return A list with:
#' \item{epitope}{Residue IDs (e.g., \code{"35_A_+42_A_"}) on the antigen within contact range.}
#' \item{paratope_h}{Heavy chain residues contacting the antigen.}
#' \item{paratope_l}{Light chain residues contacting the antigen.}
#' \item{contacts}{Matrix of all atom-level contacts (row = antigen residue, col = antibody residue).}
#'
#' @export
identify_epitopes = function(pdb, ag_chain = NULL, h_chain = NULL, l_chain = NULL, dist_cutoff = 5){

  # remove H and HETATM
  pdb = bio3d::trim.pdb(pdb, 'protein')
  pdb = bio3d::trim.pdb(pdb, 'noh')

  # set insert to '' if NA
  pdb$atom$insert[is.na(pdb$atom$insert)] = ''

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
