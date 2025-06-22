# ------------------------------------------ #
# WRAPPER_MODULE.R
# utilities for setting up run_grid, evo3d defaults, and wrapper script
# Brad Broyles
# ------------------------------------------ #

# set default controls ----
{
  .evo3d_defaults = list(
    default_msa_controls = list(
      ref_method = 1, # ** ref_method can be numeric, least_gap, or consensus -- see MSA_MODULE
      force_seqtype = NULL, # ** can be protein or nucleic acid -- NULL means to autodetect
      genetic_code = 1 # ** seqinr can use different genebank codes
    ),

    default_pdb_controls = list(
      distance_method = 'all',  # ** all atom - can be side chain, backbone, ca
      drop_incomplete_residue = T, # ** for dssp style sasa -- drop incomplete residues
      rsa.method = 'rose',         # ** sasa to rsa
      patch.dist.cutoff = 15,      # ** defualt angstrom
      patch.rsa.cutoff = 0.1,      # for surface def
      patch.sasa.cutoff = NULL,    # if using sasa instead of rsa cuts
      patch.only.exposed = T,      # seed is exposed but patch can be also buried
      max.patch = NULL,            # max aa in patch
      interface.dist.cutoff = 5    # interface dist cut --
    ),

    default_aln_controls = list(),

    default_stat_controls = list(
      calc_pi = T,
      calc_tajima = T,
      calc_hap = T,
      calc_polymorphic = T,
      calc_patch_entropy = F,
      calc_site_entropy = F, # needs flag
      valid_aa_only = F
    )
  )
}

# .setup_controls ----

#' Validate and Merge Control Lists
#'
#' Internal helper to merge user-provided control lists with defaults, warning on unrecognized keys.
#'
#' @param user_controls A named list of user-supplied parameters.
#' @param default_controls A named list of default parameters.
#' @param module_name Optional string for messaging context.
#'
#' @return A list combining defaults and valid user entries.
#' @keywords internal
.setup_controls = function(user_controls, module_name = "") {

  # switch default based on module name ----
  default_controls = switch(module_name,
                            'msa' = .evo3d_defaults$default_msa_controls,
                            'pdb' = .evo3d_defaults$default_pdb_controls,
                            'aln' = .evo3d_defaults$default_aln_controls,
                            'stat' = .evo3d_defaults$default_stat_controls
                            )

  # check keys sent by user ----
  unknown_keys = setdiff(names(user_controls), names(default_controls))

  if (length(unknown_keys) > 0) {
    message(sprintf("[%s] Unrecognized keys: %s", module_name, paste(unknown_keys, collapse = ", ")))
    user_controls = user_controls[!names(user_controls) %in% unknown_keys]
  }

  # update and return list ----
  return(
    modifyList(default_controls, user_controls)
  )

  # WOULD BE BETTER TO LOAD DEFAULTS ONCE (outside function) -- but for now is okay #
}

# .show_evo3d_defaults ----

#' Show Default evo3D Controls
#'
#' Print default evo3D module control settings (MSA, PDB, alignment, statistics).
#'
#' @param module_name Optional string: one of "msa", "pdb", "aln", or "stat". If `NULL`, shows all.
#' @return None. Prints output to console.
#' @export
show_evo3d_defaults = function(module_name = NULL){
  # module options #
  modules = c("msa", "pdb", "aln", "stat")

  # for nice print formatting #
  nice_print = function(name) {
    header = paste0("[ ", toupper(name), " controls ]")
    cat("\n", paste(rep("-", nchar(header)), collapse = ""), "\n")
    cat(header, "\n")
    cat(paste(rep("-", nchar(header)), collapse = ""), "\n\n")
    str(.evo3d_defaults[[paste0("default_", name, "_controls")]])
  }

  if (!is.null(module_name)) {
    module_name = tolower(module_name)
    if (!module_name %in% modules) {
      stop(sprintf("Invalid module_name: '%s'. Choose from: %s",
                   module_name, paste(modules, collapse = ", ")))
    }
    nice_print(module_name)
  } else {
    for (m in modules) nice_print(m)
  }

  invisible(NULL)  # to avoid printing the function definition
}

# .setup_chain_mapping ----

# handle run grid set ups
#' Setup Chain Mapping for PDB–MSA Alignment
#'
#' Standardizes chain input for main chains, interface chains, or occlusion chains across multiple PDBs and MSAs.
#'
#' This function ensures each PDB has an explicit chain mapping, even when user input is partial or symbolic (e.g. \code{"auto"}, \code{"all"}).
#'
#' @param chain Chain input; can be a character string, vector, or list (e.g. \code{"auto"}, \code{"all"}, \code{c("A", "B")}, or a list of vectors).
#' @param chain_type One of \code{"main"}, \code{"interface"}, or \code{"occlusion"}, determining the mapping logic.
#' @param pdb_count Total number of structural models.
#' @param msa_count Number of MSAs per PDB (e.g. for homomultimers).
#'
#' @return A named list of length \code{pdb_count}, where each element is a chain mapping vector or character.
#' @keywords internal
.setup_chain_mapping = function(chain, chain_type = 'main', pdb_count = 1, msa_count = 1) {
  # Input validation
  if (!chain_type %in% c('main', 'interface', 'occlusion')) {
    stop("chain_type must be one of: 'main', 'interface', 'occlusion'")
  }

  if (pdb_count < 1 || msa_count < 1) {
    stop("pdb_count and msa_count must be positive integers")
  }

  # Helper function to extend list to required length
  extend_list = function(lst, target_length, fill_value = NA) {
    if (length(lst) < target_length) {
      extra = target_length - length(lst)
      return(c(lst, rep(list(fill_value), extra)))
    }
    return(lst)
  }

  # Helper function to extend vector to required length
  extend_vector = function(vec, target_length, fill_value = "auto") {
    if (length(vec) < target_length) {
      extra = target_length - length(vec)
      return(c(vec, rep(fill_value, extra)))
    }
    return(vec)
  }

  # 1. Handle main chain type (most complex)
  if (chain_type == 'main') {

    if (is.list(chain)) {
      # Case 1: List input - extend to pdb_count and fill each entry to msa_count
      chain = extend_list(chain, pdb_count, "auto")

      out_chain = lapply(chain, function(x) {
        if (is.character(x) && length(x) == 1 && x == "auto") {
          return(rep("auto", msa_count))
        } else if (is.vector(x)) {
          return(extend_vector(x, msa_count, "auto"))
        } else {
          return(rep("auto", msa_count))  # fallback for unexpected input
        }
      })

    } else if (is.character(chain) && length(chain) == 1 && chain == "auto") {
      # Case 2: Global auto
      out_chain = replicate(pdb_count, rep("auto", msa_count), simplify = FALSE)

    } else if (is.vector(chain)) {
      # Case 3: Vector input (single PDB) - extend to msa_count and wrap in list
      extended_chain = extend_vector(chain, msa_count, "auto")
      out_chain = list(extended_chain)

      # If we need more PDBs, replicate with "auto"
      if (pdb_count > 1) {
        additional_pdbs = replicate(pdb_count - 1, rep("auto", msa_count), simplify = FALSE)
        out_chain = c(out_chain, additional_pdbs)
      }

    } else {
      # Fallback: treat as auto
      out_chain = replicate(pdb_count, rep("auto", msa_count), simplify = FALSE)
    }
  }

  # 2. Handle interface chain type (simplest)
  else if (chain_type == 'interface') {

    if (is.list(chain)) {
      # Case 1: List input - extend to pdb_count with NA
      out_chain = extend_list(chain, pdb_count, NA)

    } else if (is.character(chain) && length(chain) == 1 && chain == "all") {
      # Case 2: Global "all"
      out_chain = rep(list("all"), pdb_count)

    } else if (is.vector(chain)) {
      # Case 3: Vector input (single PDB)
      if (pdb_count == 1) {
        out_chain = list(chain)
      } else {
        # Multiple PDBs: use chain for first, NA for rest
        out_chain = c(list(chain), rep(list(NA), pdb_count - 1))
      }

    } else {
      # Fallback
      out_chain = rep(list(NA), pdb_count)
    }
  }

  # 3. Handle occlusion chain type
  else if (chain_type == 'occlusion') {

    if (is.character(chain) && length(chain) == 1 && chain == "all") {
      # Global "all"
      out_chain = rep(list("all"), pdb_count)

    } else if (is.list(chain)) {
      # List input - extend to pdb_count and process each entry
      chain = extend_list(chain, pdb_count, NA)

      out_chain = lapply(chain, function(x) {
        if (is.null(x) || (length(x) == 1 && is.na(x))) {
          return(NA)
        } else if (is.character(x) && length(x) == 1 && x == "all") {
          return("all")
        } else if (is.character(x)) {
          # Split multi-character entries unless it's "all"
          if ("all" %in% x) {
            return("all")
          } else {
            return(unlist(strsplit(x, split = "")))
          }
        } else {
          return(x)
        }
      })

    } else if (is.vector(chain)) {
      # Vector input (single PDB)
      if ("all" %in% chain) {
        processed_chain = "all"
      } else if (is.character(chain)) {
        processed_chain = unlist(strsplit(chain, split = ""))
      } else {
        processed_chain = chain
      }

      if (pdb_count == 1) {
        out_chain = list(processed_chain)
      } else {
        # Multiple PDBs: use processed chain for first, NA for rest
        out_chain = c(list(processed_chain), rep(list(NA), pdb_count - 1))
      }

    } else {
      # Fallback
      out_chain = rep(list(NA), pdb_count)
    }
  }

  # add pdb names to outchain
  names(out_chain) = paste0('pdb', seq_len(pdb_count))
  return(out_chain)
}

# .setup_multi_run_info ----

#' Setup evo3D Multi-Run Grid
#'
#' Expands combinations of MSA and PDB inputs into a run grid for multi-model, multi-chain evo3D analyses.
#' Automatically handles chain assignment and resolves homomultimer expansion where applicable.
#'
#' @param msa A single MSA (file path, matrix, or fasta object) or a list of such inputs.
#' @param pdb A single PDB (bio3d object) or a list of PDBs.
#' @param chain Main chain(s) of interest. Can be a character, vector, or nested list (e.g. \code{"A"}, \code{c("A","B")}, or \code{list(c("A","B"), "C")}).
#' @param interface_chain Chains to use for interface calculations. Can be \code{"all"}, a character vector, or list.
#' @param occlusion_chain Chains to use for occlusion in RSA/SASA. Same input rules as \code{interface_chain}.
#'
#' @return A list with standardized inputs and a \code{run_grid} data frame detailing each MSA–PDB–chain combination.
#' @keywords internal
.setup_multi_run_info = function(msa, pdb, chain, interface_chain, occlusion_chain){

  # goal is to setup run grid and return msa/pdb/chain/interface_chain/occlusion_chain objects #

  # NEEDS TO BE A CHECK THAT FILE PATHS ARE IN LIST AND NOT VECTOR #
  # -- OR MOVE VECTOR TO LIST -- #

  # 1. unpack msa information
  if(class(msa) == 'list'){
    msa_count = length(msa)
  } else {
    # assumes 1 entry #
    msa_count = 1
    msa = list(msa)
  }

  names(msa) = paste0('msa', seq_len(msa_count))

  # 2. unpack pdb information
  if(class(pdb)[1] == 'list'){
    pdb_count = length(pdb)
  } else {
    # assumes 1 entry #
    pdb_count = 1
    pdb = list(pdb)
  }

  names(pdb) = paste0('pdb', seq_len(pdb_count))

  # 3. set up msa x pdb run grid
  run_grid = expand.grid(msa = names(msa), pdb = names(pdb), stringsAsFactors = FALSE)

  run_grid$chain = NA
  #run_grid$interface_chain = NA ~ keep out of run grid? they only map to pdb
  #run_grid$occlusion_chain = NA ~ keep out of run grid? they only map to pdb

  # 4. format chains per run grid line #
  chain = .setup_chain_mapping(chain, chain_type = 'main', msa_count = msa_count, pdb_count = pdb_count)
  interface_chain = .setup_chain_mapping(interface_chain, chain_type = 'interface', msa_count = msa_count, pdb_count = pdb_count)
  occlusion_chain = .setup_chain_mapping(occlusion_chain, chain_type = 'occlusion', msa_count = msa_count, pdb_count = pdb_count)

  # unpack to run grid #
  run_grid$chain = unlist(chain)
  #run_grid$interface_chain = I(interface_chain)
  #run_grid$occlusion_chain = I(occlusion_chain)

  # now if nchar(chain) > 1 -- this is homodimer/trimer complex -- need to split into three rows? #
  # or handle later so i only do one msa to pdb mapping? #
  multichain <- which(nchar(run_grid$chain) > 1 & run_grid$chain != 'auto')
  if(length(multichain) > 0){
    expanded_rows <- list()
    for(i in multichain) {
      chains <- strsplit(run_grid$chain[i], "")[[1]]
      for(chain in chains) {
        new_row <- run_grid[i, ]
        new_row$chain <- chain
        expanded_rows <- append(expanded_rows, list(new_row))
      }
    }

    # Remove original multichain rows and add expanded ones
    run_grid <- rbind(run_grid[-multichain, ], do.call(rbind, expanded_rows))
  }

  # return
  return(list(
    msa = msa,
    pdb = pdb,
    interface_chain = interface_chain,
    occlusion_chain = occlusion_chain,
    run_grid = run_grid
  ))

}


# run_evo3D ----

#' Run evo3D Workflow
#'
#' Full evo3D wrapper to align MSA(s) and PDB(s), generate 3D-defined patches, codon mappings, and selection statistics.
#' Handles homomultimers and multi-model structures automatically.
#'
#' @param msa A matrix, file path, or list of MSAs (character matrices, fasta objects, or file paths).
#' @param pdb A \code{bio3d} PDB object or list of such objects.
#' @param chain Chain ID(s) to analyze. Can be "auto", a character, vector, or nested list.
#' @param interface_chain Chain(s) to include in interface-based patching (optional).
#' @param occlusion_chain Chain(s) to include in occlusion masking for RSA (optional).
#' @param run_selection Logical; whether to calculate nucleotide diversity, Tajima’s D, etc. (default: \code{TRUE}).
#' @param auto_chain_threshold Similarity cutoff used to match MSA peptide to PDB chains (default: \code{0.2}).
#' @param write_patch_fastas Logical; whether to save individual patch-level MSA FASTA files.
#' @param write_evo3d_df Logical; whether to save the final evo3D dataframe to CSV.
#' @param output_dir Directory path for writing outputs if above flags are \code{TRUE}.
#' @param msa_controls List of control parameters for MSA preprocessing (e.g. \code{ref_method}, \code{force_seqtype}).
#' @param pdb_controls List of control parameters for patch definition (e.g. \code{patch.dist.cutoff}, \code{rsa.method}).
#' @param aln_controls List of control parameters for MSA–structure alignment (if applicable).
#' @param stat_controls List of control parameters for patch-level statistic calculations (e.g. \code{calc_pi}, \code{calc_tajima}).
#'
#' @return A list containing:
#' \item{evo3d_df}{Final data frame with codon-level mappings and statistics.}
#' \item{final_msa_subsets}{Patch-level MSAs used in downstream statistics.}
#' \item{msa_info_sets}{Reference and peptide information for each input MSA.}
#' \item{pdb_info_sets}{Structure and patch data for each input PDB.}
#' \item{aln_info_sets}{Codon-level alignment mappings from MSA to PDB.}
#' \item{call_info}{Cached input metadata and control parameters.}
#' @export
run_evo3d = function(msa, pdb, chain = 'auto', interface_chain = NA, occlusion_chain = NA,
                     run_selection = TRUE, auto_chain_threshold = 0.2,
                     write_patch_fastas = FALSE, write_evo3d_df = FALSE, output_dir = NULL,
                     msa_controls = list(), pdb_controls = list(), aln_controls = list(), stat_controls = list()){

  #0 SETUP MODULE BEHAVIOR AND BUILD RUN INFO FROM USER INPUTS ----

  cat('STEP 0: Setting up run information and controls...\n\n')

  msa_controls = .setup_controls(msa_controls, 'msa')
  pdb_controls = .setup_controls(pdb_controls, 'pdb')
  aln_controls = .setup_controls(aln_controls, 'aln')
  stat_controls = .setup_controls(stat_controls, 'stat')

  # SETUP RUN INFO #
  run_info = .setup_multi_run_info(msa, pdb, chain, interface_chain, occlusion_chain)

  # extract run_grid and chain info #
  run_grid = run_info$run_grid
  interface_chain = run_info$interface_chain
  occlusion_chain = run_info$occlusion_chain

  #1 MODULE 1 msa_to_ref() ----

  cat('STEP 1: Converting MSAs to reference peptide sequences...\n\n')

  msa_info_sets = list()
  for(msa_name in names(run_info$msa)) {
    call_args = list(msa = run_info$msa[[msa_name]])
    msa_info_sets[[msa_name]] = do.call(msa_to_ref, c(call_args, msa_controls))
  }

  #1.5 CACHING PDBS AND FILLING IN 'auto' CHAINS ----

  #cat('STEP 1.5: Caching PDBs and resolving auto chains...\n')

  pdb_cache = list()
  chain_mappings = list()

  # keep track of similarity for quick inspection of auto chains latter #
  run_grid$auto_similarity = NA

  # Load all PDBs once and resolve chain mappings (really only need to map chains for autos)
  for(pdb_name in names(run_info$pdb)) {
    pdb_cache[[pdb_name]] = .standardize_pdb_input(pdb = run_info$pdb[[pdb_name]])

    # For each MSA that maps to this PDB, detect chains
    chain_mappings[[pdb_name]] = list()

    for(msa_name in names(run_info$msa)) {
      chain_mappings[[pdb_name]][[msa_name]] = .auto_detect_chain(
        pep = msa_info_sets[[msa_name]]$pep,
        pdb = pdb_cache[[pdb_name]],
        in_module = TRUE
      )
    }
  }

  # UPDATE RUN_GRID WITH RESOLVED CHAINS
  for(i in 1:nrow(run_grid)) {
    if(is.na(run_grid$chain[i])) next
    if(run_grid$chain[i] == "auto") {
      pdb_name = run_grid$pdb[i]
      msa_name = run_grid$msa[i]

      # Get best chain match
      best_chain = names(chain_mappings[[pdb_name]][[msa_name]])[1]
      run_grid$chain[i] = best_chain

      # add similarity scores
      sim = chain_mappings[[pdb_name]][[msa_name]][1]
      run_grid$auto_similarity[i] = sim
    }
  }

  # need to filter auto threshold -- off for now #

  #2 MODULE 2 pdb_to_patch() ----

  # COULD BE GRACEFUL ERROR -- one of your pdbs has no MSA mapping in any chain? #

  cat('STEP 2: Converting PDBs to patches...\n\n')

  # GATHER CHAIN, INTERFACE, AND OCCLUSION FOR EACH PDB #
  pdb_info_sets = list()
  for (pdb_name in unique(run_grid$pdb)) {
    chain_set = run_grid$chain[run_grid$pdb == pdb_name]
    chain_set = unique(chain_set)
    call_args = list(
      pdb = pdb_cache[[pdb_name]],
      chain = chain_set,
      interface_chain = interface_chain[[pdb_name]],
      occlusion_chain = occlusion_chain[[pdb_name]]
    )
    pdb_info_sets[[pdb_name]] = do.call(pdb_to_patch, c(call_args, pdb_controls))

  }

  # CLEAR CACHE and save call to reduce memory usage ----

  # capture call info #
  call_info = list(
    msa = lapply(run_info$msa, function(x) if(is.character(x) && length(x) == 1) x else NA),
    pdb = lapply(run_info$pdb, function(x) if(is.character(x) && length(x) == 1) x else NA),
    run_grid = run_grid,
    interface_chain = interface_chain,
    occlusion_chain = occlusion_chain,
    msa_controls = msa_controls,
    pdb_controls = pdb_controls,
    aln_controls = aln_controls,
    stat_controls = stat_controls
  )

  # remove run_info and pdb_cache
  rm(run_info, pdb_cache)
  invisible(gc())

  #3 MODULE 3 aln_msa_to_pdb ----

  cat('STEP 3: Aligning MSAs to PDBs...\n\n')

  # if chain is NA (no pdb to msa mapping -- add empty set) #
  # this strategy keeps pdb and msa numbering from extend_pdb()
  # and extend_msa() -- exact #

  aln_info_sets = list()
  for(i in seq_len(nrow(run_grid))) {

    # IF CHAIN IS NA -- CREATE EMPTY ALN_INFO_SET #
    # THIS ENSURES PDB1 in evo3d_df is pdb1 in run_info #
    if(is.na(run_grid$chain[i])){

      # we can make fake aln_info_set # just pull codon and msa info #
      # if pdb mapped to no chains -- pdb_to_patch() would have failed #

      df = data.frame(
        residue_id = NA,
        codon_patch = NA,
        codon = NA,
        pdb_aa = NA,
        ref_aa = strsplit(msa_info_sets[[run_grid$msa[i]]]$pep, "")[[1]]
      )

      df$codon = 1:nrow(df)

      empty_set = list(
        msa_subsets = NA,
        aln_df = df,
        aln_coverage = NA
      )

      aln_name <- paste(run_grid$msa[i], run_grid$pdb[i], run_grid$chain[i], sep="_")
      aln_info_sets[[aln_name]] = empty_set
      next
    }

    call_args = list(
      msa_info = msa_info_sets[[run_grid$msa[i]]],
      pdb_info = pdb_info_sets[[run_grid$pdb[i]]],
      chain = run_grid$chain[i]
    )

    aln_name <- paste(run_grid$msa[i], run_grid$pdb[i], run_grid$chain[i], sep="_")
    aln_info_sets[[aln_name]] = do.call(aln_msa_to_pdb, c(call_args, aln_controls))
  }

  #3.5 (msa and pdb patch extensions) ----

  # Work on a copy so original is preserved for inspection
  working_aln_sets = aln_info_sets
  working_run_grid = run_grid

  # Check for homomultimers to extend ---
  msa_pdb_counts = table(paste(working_run_grid$msa, working_run_grid$pdb, sep="_"))
  needs_homomultimer_extension = names(msa_pdb_counts)[msa_pdb_counts > 1]

  # Extend homomultimers
  for(multimer in needs_homomultimer_extension){

    cat('STEP 3.5: Extending homomultimers...\n')

    # can do all at once -- then update run grid so pdb extend doesnt run on these #
    rows_for_multimer = which(paste(working_run_grid$msa, working_run_grid$pdb, sep="_") == multimer)
    msa_id = working_run_grid$msa[rows_for_multimer[1]]

    extended_result = extend_pdb_homomultimer(working_aln_sets[rows_for_multimer],
                                              msa_info_sets[[msa_id]])

    # collapse these rows so extend_pdb doesnt try to run (paste chains back to one) #
    chains = unique(working_run_grid$chain[rows_for_multimer])
    chains = paste0(chains, collapse = "")
    working_run_grid$chain[rows_for_multimer] = chains
    working_run_grid = working_run_grid[-rows_for_multimer[-1], ]

    # Replace first entry, remove others
    working_aln_sets[[rows_for_multimer[1]]] = extended_result
    working_aln_sets[rows_for_multimer[-1]] = NULL
  }

  # Check for PDB extension needs (same MSA --> multiple PDbs)
  msa_counts = table(working_run_grid$msa)
  needs_pdb_extension = names(msa_counts)[msa_counts > 1]

  # First: Handle patch extensions (same MSA --> multiple PDBs)
  for(msa_id in needs_pdb_extension) {

    cat('STEP 3.5: Extending complimentary PDB info...\n')

    rows_for_msa = which(working_run_grid$msa == msa_id)

    # Start with the first PDB's result
    extended_result = working_aln_sets[[rows_for_msa[1]]]

    # Iteratively extend with each additional PDB
    for(i in 2:length(rows_for_msa)) {
      extended_result = extend_pdb(extended_result,
                                      working_aln_sets[[rows_for_msa[i]]],
                                      msa_info_sets[[msa_id]])
    }

    # Replace first entry, remove others
    working_aln_sets[[rows_for_msa[1]]] = extended_result
    working_aln_sets[rows_for_msa[-1]] = NULL

    # Update working_run_grid
    working_run_grid$pdb[rows_for_msa[1]] = paste(working_run_grid$pdb[rows_for_msa], collapse = "+")
    working_run_grid$chain[rows_for_msa[1]] = paste(working_run_grid$chain[rows_for_msa], collapse = "+")
    working_run_grid = working_run_grid[-rows_for_msa[-1], ]
  }

  # Second: Handle MSA extensions
  pdb_counts = table(working_run_grid$pdb)
  needs_msa_extension = names(pdb_counts)[pdb_counts > 1]

  for(pdb_id in needs_msa_extension) {

    cat('STEP 3.5: Extending multi-chain info...\n')

    rows_for_pdb = which(working_run_grid$pdb == pdb_id)

    # Start with first, extend with rest
    extended_result = working_aln_sets[[rows_for_pdb[1]]]
    first_msa = working_run_grid$msa[rows_for_pdb[1]]

    msa_set = c(first_msa)
    for(i in 2:length(rows_for_pdb)) {
      current_msa = working_run_grid$msa[rows_for_pdb[i]]
      msa_set = c(msa_set, current_msa)
      extended_result = extend_msa(extended_result, working_aln_sets[[rows_for_pdb[i]]],
                                   msa_info_sets[c(msa_set)],
                                   use_sample_names = use_sample_names
                                   )
    }

    # Replace first, remove others
    working_aln_sets[[rows_for_pdb[1]]] = extended_result
    working_aln_sets[rows_for_pdb[-1]] = NULL
    working_run_grid = working_run_grid[-rows_for_pdb[-1], ]
  }

  # Final result is working_aln_sets[[1]] (should be only one left)
  final_result = working_aln_sets[[1]]

  #4 MODULE 4 calculate_patch_stats ----

  cat('STEP 4: Calculating patch statistics...\n\n')

  # build evo3d -- might not have selection data
  evo3d_df = final_result$aln_df

  if (run_selection) {

    # see which stats are on #
    valid_aa_only = stat_controls$valid_aa_only

    # add polymorphic sites -- needs wrapped across msa subsets #
    if(stat_controls$calc_polymorphic){
      evo3d_df = calculate_polymorphic_residue(msa_info_sets, evo3d_df)
    }

    if(stat_controls$calc_patch_entropy){
      evo3d_df = calculate_patch_entropy(
        msa = final_result$msa_subsets,
        residue_df = evo3d_df,
        valid_aa_only = valid_aa_only
      )
    }

    stat = c()
    if (stat_controls$calc_pi) {
      stat = c(stat, 'pi')
    }

    if (stat_controls$calc_tajima) {
      stat = c(stat, 'tajima')
    }

    if (stat_controls$calc_hap) {
      stat = c(stat, 'hap')
    }

    if(length(stat) > 0){
    evo3d_df = run_pegas_three(msa = final_result$msa_subsets,
                               residue_df = evo3d_df,
                               stat = stat)
    }

  }


  #5 saving/writing -- what to return ----

  cat('STEP 5: Saving results...\n\n')

  if (write_patch_fastas){

    # if output_dir is not provided, use current #
    if (is.null(output_dir)) {
      output_dir = '.'
    }

    fasta_dir = file.path(output_dir, "patch_fastas")

    write_patch_fastas(final_result$msa_subsets, output_dir = fasta_dir)
  }

  # reorder the columns of evo3d_df #
  # msa (if available), codon, msa_subset_id, ref_aa, pdbX_aa, pdbY_aa, ..., pdbX_residue_id, pdbY_residue_id, codon_patch, everything else #
  codon_info =  intersect(c("msa","codon","msa_subset_id","ref_aa"), names(evo3d_df))
  aa_cols = grep("^pdb.*_aa$",          names(evo3d_df), value = TRUE)
  id_cols = grep(".*residue_id$",  names(evo3d_df), value = TRUE)
  patch_col = "codon_patch"
  other = setdiff(names(evo3d_df), c(codon_info, aa_cols, id_cols, patch_col))
  col_order = c(codon_info, aa_cols, id_cols, patch_col, other)
  evo3d_df = evo3d_df[, col_order, drop = FALSE]

  if (write_evo3d_df) {

    # if output_dir is not provided, use current #
    if (is.null(output_dir)) {
      output_dir = '.'
    }

    csv_path = file.path(output_dir, '/evo3d_df.csv')

    write.csv(evo3d_df, file = csv_path, row.names = FALSE, quote = FALSE)

  }

  return(list(
    evo3d_df = evo3d_df,
    final_msa_subsets = final_result$msa_subsets,
    msa_info_sets = msa_info_sets,
    pdb_info_sets = pdb_info_sets,
    aln_info_sets = aln_info_sets,
    call_info = call_info
    )
  )

}
