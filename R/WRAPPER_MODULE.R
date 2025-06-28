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
      patch.only.exposed = TRUE,      # seed is exposed but patch can be also buried
      max.patch = NULL,            # max aa in patch
      interface.dist.cutoff = 5    # interface dist cut --
    ),

    default_aln_controls = list(
      use_sample_names = TRUE,     # actually used in in extend_msa() - stored here
      auto_chain_threshold = 0.4,   # also used before aln module but conceptually fits
      kmer_size = 4                # kmer size for auto chain mapping
    ),

    default_stat_controls = list(
      calc_pi = FALSE,
      calc_tajima = FALSE,
      calc_hap = TRUE,
      calc_polymorphic = TRUE,
      calc_patch_entropy = FALSE,
      calc_site_entropy = FALSE, # needs flag
      valid_aa_only = FALSE # removes non-standard amino acids from stats (entropy and polymorphic)
    ),

    default_output_controls = list(
      output_dir = NULL, # ** if NULL, then we dont save output files
      write_msa_subsets = TRUE, # ** write msa subsets to output_dir/msa_subsets
      write_evo3d_df = TRUE, # ** write evo3d_df to output dir
      write_call_info = TRUE, # ** write call info to output dir
      write_module_intermediates = FALSE, # ** write intermediate files to output dir (as .rds)
      prefix = '' # prefix for files and folders
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
                            'stat' = .evo3d_defaults$default_stat_controls,
                            'output' = .evo3d_defaults$default_output_controls,
                            )

  # check keys sent by user ----
  unknown_keys = setdiff(names(user_controls), names(default_controls))

  # I THINK WE SHOULD STOP AT UNKNOWN KEYS -- THAT WAY WE DONT WASTE TIME RUNNING UNWANTED PARAMETERS #
  # 6/26/25 -- still need to stop #

  if (length(unknown_keys) > 0) {
    message(sprintf("[%s] Unrecognized keys: %s", module_name, paste(unknown_keys, collapse = ", ")))
    user_controls = user_controls[!names(user_controls) %in% unknown_keys]
  }

  # update and return list ----
  return(
    modifyList(default_controls, user_controls)
  )

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
  modules = c("msa", "pdb", "aln", "stat", "output")

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



# .safe_save ----
#' Safely Resolve a Non-Clashing Output Path
#'
#' Avoids overwriting files/directories by appending numeric suffixes.
#' If no safe name is found, uses a fallback with a consistent timestamp.
#'
#' @param path Character. Desired path.
#' @param is_dir Logical. Is it a directory?
#' @param max_tries Integer. How many numbered attempts before fallback to POSIXct timestamp
#' @param tag Character. Optional tag to append to the path.
#'
#' @return A list with:
#'   \item{path}{Safe file path}
#'   \item{systime}{Used timestamp (in fallback or NULL if not used)}
#'
#' @keywords internal
#' @keywords internal
.safe_save <- function(path, is_dir = FALSE) {
  base <- tools::file_path_sans_ext(path)
  ext  <- tools::file_ext(path)
  ext  <- if (nzchar(ext)) paste0(".", ext) else ""

  # If safe, return original
  if (!file.exists(path)) return(list(path = path, tag = NULL))
  if (is_dir && length(list.files(path)) == 0) return(list(path = path, tag = NULL))
  if (!is_dir && file.info(path)$size == 0) return(list(path = path, tag = NULL))

  # Use systime + pid fallback
  systime <- Sys.time()
  pid <- Sys.getpid()
  tag <- paste0(format(systime, "%Y%m%d%H%M%S"), "_pid", pid)
  tagged_path <- if (is_dir) paste0(base, "_", tag) else paste0(base, "_", tag, ext)

  list(path = tagged_path, tag = tag)
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
#' @param compute_stats Logical; whether to calculate nucleotide diversity, Tajima’s D, etc. (default: \code{TRUE}).
#' @param detail_level Controls the level of detail in the returned R object.
#'
#' \describe{
#'   \item{0}{Return only `evo3d_df` and `pdb_info_sets`. Minimal memory usage. Suitable for post hoc visualization and writing to PDB files, but not for restart or stat recomputation.}
#'   \item{1}{Adds `msa_info_sets`, enabling restart with user-provided alignments. Does not include full alignment-to-structure mappings or subsets.}
#'   \item{2}{Adds full `msa_subsets` and MSA-to-PDB patch mappings, allowing for patch-level stat recomputation or post hoc reruns on individual structures.}
#' }
#'
#' @param verbose Controls the verbosity of printed output.
#'
#' \describe{
#'   \item{0}{Silent mode. Only errors or critical messages are printed.}
#'   \item{1}{Basic progress updates showing which module is running.}
#'   \item{2}{Detailed output, including all submodular steps.}
#'  }
#'
#' @param restart_run A previous evo3d run object used to restart the pipeline (e.g., after alignment failure or to reuse mapped PDB metadata).
#' @param user_aln A user-provided MSA aligned to the reference sequence. Can be a single alignment or a list (for multi-run pipelines). Used in conjunction with \code{restart_run}.
#' @param msa_controls List of control parameters for MSA preprocessing (e.g. \code{ref_method}, \code{force_seqtype}).
#' @param pdb_controls List of control parameters for patch definition (e.g. \code{patch.dist.cutoff}, \code{rsa.method}).
#' @param aln_controls List of control parameters for MSA–structure alignment (if applicable).
#' @param stat_controls List of control parameters for patch-level statistic calculations (e.g. \code{calc_pi}, \code{calc_tajima}).
#' @param output_controls List of control parameters for output formatting (e.g. \code{output_format}, \code{output_dir}).
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
                     compute_stats = TRUE, detail_level = 1, verbose = 1,
                     restart_run = NULL, user_aln = NULL,
                     msa_controls = list(), pdb_controls = list(), aln_controls = list(),
                     stat_controls = list(), output_controls = list()){

  # COULD ADD 1 more step -- validating input types (quick - is detail and verbose here ...) #
  # 6/26/25 need to add #

  #0 SETUP MODULE BEHAVIOR AND BUILD RUN INFO FROM USER INPUTS ----
  if(verbose > 0){
    cat('STEP 0: Setting up run information and controls...\n')
    if(verbose > 1){
      cat('\tUpdating controls with defaults and user inputs\n')
    }
  }

  # setup custom or just default parameters #
  msa_controls = .setup_controls(msa_controls, 'msa')
  pdb_controls = .setup_controls(pdb_controls, 'pdb')
  aln_controls = .setup_controls(aln_controls, 'aln')
  stat_controls = .setup_controls(stat_controls, 'stat')
  output_controls = .setup_controls(output_controls, 'output')

  # logically these group with aln_controls but are not called in aln_msa_to_pdb() #
  use_sample_names = aln_controls$use_sample_names
  auto_chain_threshold = aln_controls$auto_chain_threshold
  kmer_size = aln_controls$kmer_size

  aln_controls$use_sample_names = NULL
  aln_controls$auto_chain_threshold = NULL
  aln_controls$kmer_size = NULL

  # setup run info #
  if(verbose > 1){
    cat('\tBuilding run grid\n')
  }

  run_info = .setup_multi_run_info(msa, pdb, chain, interface_chain, occlusion_chain)

  # extract run_grid and chain info #
  run_grid = run_info$run_grid
  interface_chain = run_info$interface_chain
  occlusion_chain = run_info$occlusion_chain

  #1 MODULE 1 msa_to_ref() ----
  if(verbose > 0){
    cat('STEP 1: Converting MSAs to reference peptide sequences...\n')
  }

  msa_info_sets = list()
  for(msa_name in names(run_info$msa)) {
    call_args = list(msa = run_info$msa[[msa_name]],
                     verbose = verbose - 1)
    msa_info_sets[[msa_name]] = do.call(msa_to_ref, c(call_args, msa_controls))
  }

  #1.5 CACHING PDBS AND FILLING IN 'auto' CHAINS ----
  auto_ros = which(run_grid$chain == 'auto')
  run_grid$kmer_match = NA
  pdb_cache = list()

  if (verbose > 0) {
    if (length(auto_ros) > 0) {
      cat("STEP 1.5: Caching PDBs and resolving auto chains...\n")
    } else {
      cat("STEP 1.5: Caching PDBs (no auto chains to resolve)...\n")
    }
  }

  # cache pdbs #
  for(pdb_name in names(run_info$pdb)) {
    pdb_cache[[pdb_name]] = .standardize_pdb_input(pdb = run_info$pdb[[pdb_name]])
  }

  # add auto detection if needed #
  if(F){
  for (i in auto_ros) {
    # get pdb and msa names
    pdb_name = run_grid$pdb[i]
    msa_name = run_grid$msa[i]

    # Run auto detect #
    chain_mappings = .auto_detect_chain(
      pep = msa_info_sets[[msa_name]]$pep,
      pdb = pdb_cache[[pdb_name]],
      in_module = TRUE,
      k = kmer_size
    )

    # check for passing auto_chain_threshold #
    if(!any(chain_mappings >= auto_chain_threshold)) {
      if (verbose > 0) {
        message(sprintf("In auto chain mapping: No chains in '%s' passed the k-mer threshold for '%s' (%.2f).\n",
                    pdb_name, msa_name, auto_chain_threshold))
      }
      run_grid$chain[i] = NA
      next
    }

    # get best chain based on kmer similarity #
    run_grid$chain[i] = names(chain_mappings)[1]

    # add similarity scores #
    run_grid$kmer_match[i] = round(chain_mappings[1], 2)
  }
  }

  # auto chain for all #
  for (i in seq_len(nrow(run_grid))) {
    # get pdb and msa names
    pdb_name = run_grid$pdb[i]
    msa_name = run_grid$msa[i]

    # Run auto detect #
    chain_mappings = .auto_detect_chain(
      pep = msa_info_sets[[msa_name]]$pep,
      pdb = pdb_cache[[pdb_name]],
      in_module = TRUE,
      k = kmer_size
    )

    # if in auto -- try to pass autochain threshold #
    if(i %in% auto_ros){
      # check for passing auto_chain_threshold #
      if(!any(chain_mappings >= auto_chain_threshold)) {
        if (verbose > 0) {
          message(sprintf("In auto chain mapping: No chains in '%s' passed the k-mer threshold for '%s' (%.2f).\n",
                          pdb_name, msa_name, auto_chain_threshold))
        }
        run_grid$chain[i] = NA
        next
      }

      # get best chain based on kmer similarity #
      run_grid$chain[i] = names(chain_mappings)[1]

      # add similarity scores #
      run_grid$kmer_match[i] = round(chain_mappings[1], 2)
    }

    # if not see if specified chain is valid #
    if(run_grid$chain[i] %in% names(chain_mappings)){
      # check for passing auto_chain_threshold #
      if(!any(chain_mappings >= auto_chain_threshold)) {
        if (verbose > 0) {
          message(sprintf("In specified chain mapping: Chain '%s' in '%s' did not pass the k-mer threshold for '%s' (%.2f).\n",
                          run_grid$chain[i], pdb_name, msa_name, auto_chain_threshold))
        }
        run_grid$chain[i] = NA
        next
      }

      # add similarity scores #
      run_grid$kmer_match[i] = round(chain_mappings[run_grid$chain[i]], 2)
    }
  }

  if (verbose > 1) {
    cat("\tRun grid for rest of analysis:\n\n")
    out <- capture.output(print(run_grid))
    cat(paste0("\t", out), sep = "\n")
    cat('\n')
  }

  # !!! STOP RUN IF A PDB OR MSA DOESNT HAVE ANY MAPPINGS? !!! ----

  # save call_info so it can be exported #
  call_info = list(
    msa = lapply(run_info$msa, function(x) if(is.character(x) && length(x) == 1) x else NA),
    pdb = lapply(run_info$pdb, function(x) if(is.character(x) && length(x) == 1) x else NA),
    run_grid = run_grid,
    interface_chain = interface_chain,
    occlusion_chain = occlusion_chain,
    msa_controls = msa_controls,
    pdb_controls = pdb_controls,
    aln_controls = c(aln_controls,
                     use_sample_names = use_sample_names,
                     auto_chain_threshold = auto_chain_threshold,
                     kmer_size = kmer_size), # pack back in previosly removed aln_controls #
    stat_controls = stat_controls,
    output_controls = output_controls
  )

  # stop if no mappings (nonsense data) #
  msa_valid = tapply(!is.na(run_grid$chain), run_grid$msa, any)
  bad_msas = names(msa_valid)[!msa_valid]

  pdb_valid = tapply(!is.na(run_grid$chain), run_grid$pdb, any)
  bad_pdbs = names(pdb_valid)[!pdb_valid]

  if (length(bad_msas) > 0 || length(bad_pdbs) > 0) {
    stop_msg <- "!!! STOPPING EVO3D RUN !!!\nThe following MSAs or PDBs had no valid mappings:\n"
    if (length(bad_msas) > 0) {
      stop_msg <- paste0(stop_msg, sprintf("MSAs: %s\n", paste(bad_msas, collapse = ", ")))
    }
    if (length(bad_pdbs) > 0) {
      stop_msg <- paste0(stop_msg, sprintf("PDBs: %s\n", paste(bad_pdbs, collapse = ", ")))
    }
    stop_msg <- paste0(stop_msg, "Please check inputs and call_info\n")
    stop_msg <- paste0(stop_msg, "Not advised but set auto_chain_threshold = 0 to skip this check")
    message(paste0(stop_msg, "\n"))

    return(list(
      call_info = call_info
    ))
  }

  #2 MODULE 2 pdb_to_patch() ----

  if(verbose > 0) {
    cat('STEP 2: Converting PDBs to patches...\n')
  }

  # GATHER CHAIN, INTERFACE, AND OCCLUSION FOR EACH PDB #
  pdb_info_sets = list()
  for (pdb_name in unique(run_grid$pdb)) {
    # gather all relavant chains for this PDB #
    chain_set = run_grid$chain[run_grid$pdb == pdb_name]

    # remove NA and get unique chains #
    chain_set = chain_set[!is.na(chain_set)]
    chain_set = unique(chain_set)

    # if there are no chains (this pdb doesnt map to any MSA) -- already cleared out #
    #if(length(chain_set) == 0) {
    #  if(verbose > 0) {
    #    message(sprintf("No chains found for PDB '%s'. Skipping...\n", pdb_name))
    #  }
    #  next
    #}

    # gather interface and occlusion chains #
    call_args = list(
      pdb = pdb_cache[[pdb_name]],
      chain = chain_set,
      interface_chain = interface_chain[[pdb_name]],
      occlusion_chain = occlusion_chain[[pdb_name]],
      verbose = verbose - 1, # reduce by one level
      detail_level = detail_level
    )

    # rub pdb_to_patch() with the gathered chains #
    pdb_info_sets[[pdb_name]] = do.call(pdb_to_patch, c(call_args, pdb_controls))
  }

  # CLEAR CACHE reduce memory usage ----

  # remove run_info and pdb_cache
  rm(run_info, pdb_cache, call_args, chain_set, pdb_name, bad_pdbs, bad_msas, msa_valid, pdb_valid,
     chain_mappings, auto_chain_threshold, i, auto_ros, msa_name, interface_chain, occlusion_chain,
     msa_controls, pdb_controls)

  invisible(gc())

  #3 MODULE 3 aln_msa_to_pdb ----
  if (verbose > 0) {
    cat('STEP 3: Aligning MSAs to PDBs...\n')
  }

  aln_info_sets = list()
  for(i in seq_len(nrow(run_grid))) {

    # IF CHAIN IS NA -- CREATE EMPTY ALN_INFO_SET #
    # THIS ENSURES PDB1 in evo3d_df is pdb1 in run_info #
    if(is.na(run_grid$chain[i])){

      # print skipping one row if verbose > 1
      if(verbose > 1) {
        cat(sprintf("\tSkipping run grid row %d: No MSA mapping for PDB\n", i))
      }

      # create empty alignment info set for this row #
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
      chain = run_grid$chain[i],
      verbose = verbose - 1
    )

    aln_name <- paste(run_grid$msa[i], run_grid$pdb[i], run_grid$chain[i], sep="_")
    aln_info_sets[[aln_name]] = do.call(aln_msa_to_pdb, c(call_args, aln_controls))
  }

  suppressWarnings(rm(aln_name, call_args, i))

  #3.5 (msa and pdb patch extensions) ----

  # Work on a copy so original is preserved for inspection
  working_aln_sets = aln_info_sets
  working_run_grid = run_grid

  # Check for homomultimers to extend #
  msa_pdb_counts = table(paste(working_run_grid$msa, working_run_grid$pdb, sep="_"))
  needs_homomultimer_extension = names(msa_pdb_counts)[msa_pdb_counts > 1]

  # Extend homomultimers
  for(multimer in needs_homomultimer_extension){
    if (verbose > 0) {
      cat('STEP 3.5: Extending homomultimers...\n')
    }

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

  suppressWarnings(rm(msa_pdb_counts, needs_homomultimer_extension, multimer, rows_for_multimer,
                      msa_id, extended_result, chains))

  # Handle patch extensions across PDBs
  msa_counts = table(working_run_grid$msa)
  needs_pdb_extension = names(msa_counts)[msa_counts > 1]

  for(msa_id in needs_pdb_extension) {

    # 6/26/25 !! really need to propogate occlusions #
    # think e1e2 complex for 1 pdb, but just e2 for the other #
    # that e1 occlusion should be kept #
    if (verbose > 0) {
      cat('STEP 3.5: Extending complimentary PDB info...\n')
    }

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

  suppressWarnings(rm(msa_counts, needs_pdb_extension, msa_id, rows_for_msa, extended_result, i))

  # Handle patch extensions across MSAs
  pdb_counts = table(working_run_grid$pdb)
  needs_msa_extension = names(pdb_counts)[pdb_counts > 1]

  for(pdb_id in needs_msa_extension) {
    if (verbose > 0) {
      cat('STEP 3.5: Extending multi-chain info...\n')
    }

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

  suppressWarnings(rm(pdb_counts, needs_msa_extension, pdb_id, rows_for_pdb, first_msa,
                      extended_result, i, msa_set, current_msa, working_aln_sets))

  # CLEAN UP INTERMEDIATES TO DESIRED LEVEL ----

  # rest of detail_level clean ups after rights and computing stats #
  if(detail_level < 3){
    # clear out intermediate msa_subsets #
    aln_info_sets = lapply(aln_info_sets, function(x) {
      x$msa_subsets = NULL
      x
    })
  }

  invisible(gc())

  #4 MODULE 4 calculate_patch_stats ----
  evo3d_df = final_result$aln_df

  if (!compute_stats){
    if (verbose >= 1) {
      cat('STEP 4: Skipping stats calculation...\n')
    }
  } else {

    if (verbose >= 1) {
      cat('STEP 4: Calculating patch stats...\n')
    }

    # add polymorphic sites -- needs wrapped across msa subsets #
    if(stat_controls$calc_polymorphic){
      evo3d_df = calculate_polymorphic_residue(
        msa_info_sets,
        evo3d_df,
        valid_aa_only = stat_controls$valid_aa_only
        )
    }

    if(stat_controls$calc_patch_entropy){
      evo3d_df = calculate_patch_entropy(
        msa = final_result$msa_subsets,
        residue_df = evo3d_df,
        valid_aa_only = stat_controls$valid_aa_only
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

  # CLEAN UP INTERMEDIATES TO DESIRED LEVEL ----
  if(detail_level < 1){
    # just enough to plot #
    msa_info_sets = NULL
    aln_info_sets = NULL
  }

  #5 saving to disk ----

  # reorder the columns of evo3d_df #
  # msa (if available), codon, msa_subset_id, ref_aa, pdbX_aa, pdbY_aa, ..., pdbX_residue_id, pdbY_residue_id, codon_patch, everything else #
  codon_info =  intersect(c("msa","codon","msa_subset_id","ref_aa"), names(evo3d_df))
  aa_cols = grep("^pdb.*_aa$", names(evo3d_df), value = TRUE)
  id_cols = grep(".*residue_id$", names(evo3d_df), value = TRUE)
  patch_col = "codon_patch"
  other = setdiff(names(evo3d_df), c(codon_info, aa_cols, id_cols, patch_col))
  col_order = c(codon_info, aa_cols, id_cols, patch_col, other)
  evo3d_df = evo3d_df[, col_order, drop = FALSE]

  # if output_dir is not null then write results to disk #
  if (!is.null(output_controls$output_dir)) {
    output_dir = output_controls$output_dir

    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }

    if (verbose >= 1) {
      cat(sprintf('STEP 5: Writing results to %s/ ...\n', output_dir))
    }

    prefix = output_controls$prefix
    prefix = ifelse(prefix == "", prefix, paste0(prefix, "_"))

    # check if any of the files i will write already exist #
    # if so i want to write all the following with the same tag #
    fasta_dir = file.path(output_dir, paste0(prefix, "msa_subsets"))
    csv_path = file.path(output_dir, paste0(prefix, 'evo3d_df.csv'))
    call_path = file.path(output_dir, paste0(prefix, 'call_info.json'))
    intermediates_path = file.path(output_dir, paste0(prefix, "evo3d_intermediates.rds"))

    # check for safe save versions #
    fasta_dir_safe = .safe_save(fasta_dir, is_dir = TRUE)
    csv_path_safe = .safe_save(csv_path, is_dir = FALSE)
    call_path_safe = .safe_save(call_path, is_dir = FALSE)
    intermediates_path_safe = .safe_save(intermediates_path, is_dir = FALSE)

    tags = c(fasta_dir_safe$tag, csv_path_safe$tag, call_path_safe$tag, intermediates_path_safe$tag)

    if (!is.null(tags)) {
      # just take first tag #
      tag = paste0('_', tags[1])

      if (verbose >= 1) {
        message(
          sprintf(
            'WARNING: Some output files already exist, using tag "%s" to avoid overwriting.\n',
            tag
          )
        )
      }

    } else {
      tag = NULL
    }

    if (output_controls$write_msa_subsets) {
      fasta_dir = file.path(paste0(fasta_dir, tag))

      # fast writes with tar to temp_dir, but could possibly fail on HPC #
      tryCatch({
        write_patch_fastas(final_result$msa_subsets, output_dir = fasta_dir)
      }, error = function(e) {
        warning("Fast write failed — falling back to slow file_writes mode.")
        write_patch_fastas_slow(final_result$msa_subsets, output_dir = fasta_dir)
      })


    }

    if (output_controls$write_evo3d_df) {

      csv_path = file.path(paste0(tools::file_path_sans_ext(csv_path), tag, '.csv'))

      write.csv(evo3d_df,
                file = csv_path,
                row.names = FALSE,
                quote = FALSE)

    }

    if (output_controls$write_call_info) {
      call_path = file.path(paste0(tools::file_path_sans_ext(call_path), tag, '.json'))


      jsonlite::write_json(
        call_info,
        path = call_path,
        auto_unbox = TRUE,
        pretty = TRUE
      )

    }

    if (output_controls$write_module_intermediates){

      intermediates_path = file.path(paste0(tools::file_path_sans_ext(intermediates_path), tag, '.rds'))

      saveRDS(list(msa_info_sets, pdb_info_sets, aln_info_sets))
    }

  } else {
    if (verbose >= 1) {
      cat('STEP 5: No output directory specified in `output_controls$output_dir`, skipping writing results.\n')
    }
  }

  # CLEAN UP INTERMEDIATES TO DESIRED LEVEL ----
  if(detail_level < 2){
    # not returning msa subsets #
    final_result$msa_subsets = NULL
    invisible(gc())
  }

  if(detail_level < 0){
    # no return at all assumes data of intereset is written to disk #
    return(invisible())
  }

  #6 return results ----
  if (verbose > 0) {
    cat('---- RUN COMPLETE ----\n')
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
