# ------------------------------------------ #
# WRAPPER_MODULE.R
# utilities for setting up run_grid, evo3d defaults, and wrapper script
# Brad Broyles
# ------------------------------------------ #

# set default controls ----
{
  .evo3d_defaults = list(

    default_msa_controls = list(
      ref_method = 'consensus', # ** ref_method can be numeric, least_gap, or consensus -- see MSA_MODULE
      force_seq_type = NULL, # ** can be protein or nucleic acid -- NULL means to autodetect
      detect_sequence_threshold = 0.8,
      detect_sequence_len = 100,
      genetic_code = 1 # ** seqinr can use different genebank codes
    ),

    default_pdb_controls = list(
      distance_method = 'ca',  # ** ca atom - can be side chain, backbone, ca, all, centroid
      drop_incomplete_residue = TRUE, # ** for dssp style sasa -- drop incomplete residues
      rsa_method = 'rose',         # ** sasa to rsa
      dist_cutoff = 15,      # **  angstrom distance
      rsa_cutoff = 0.1,      # for surface def
      sasa_cutoff = NA,    # if using sasa
      use_rsa_sasa = 'or',         # if both provided how to use
      only_exposed_in_patch = TRUE,      # seed is exposed but patch can be also buried
      max_patch = NA,            # max aa in patch
      interface_dist_cutoff = 5,    # interface dist cut --
      force_file_type = NULL, # can be pdb or cif
      patch_mode = 'codon'  # will codons be duplicated and count as residues or should codons define windows -- used in rebuild_patches()
    ),

    default_aln_controls = list(
      use_sample_names = TRUE,     # actually used in in extend_msa() - stored here
      auto_chain_threshold = 0.2,   # also used before aln module but conceptually fits
      kmer_size = 4,              # kmer size for auto chain mapping
      #next_chain_tolerance = 0.8, # used with autochain to keep good chain matches
      merge_type = 'exposure_distance',    # change this to union / by_distance / by_exposure_distance
      merge_exposure = 0.5       # how often does a residue need exposure to count as exposed in codon collapse
    ),

    default_stat_controls = list(
      calc_pi = FALSE,
      calc_tajima = FALSE,
      calc_hap = FALSE,
      calc_polymorphic = TRUE,
      calc_patch_entropy = FALSE,
      calc_site_entropy = FALSE, # needs flag
      valid_aa_only = TRUE # removes non-standard amino acids from stats (entropy and polymorphic)
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

#' validate and merge control lists
#'
#' internal helper to merge user-supplied control parameters with defaults.
#' unknown keys stop with an error so unwanted arguments do not pass through.
#'
#' @param user_controls named list of user-supplied parameters
#' @param module_name character. module context used to select defaults
#'   (one of "msa","pdb","aln","stat","output")
#'
#' @return named list of controls, with defaults updated by valid user entries
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
    stop(sprintf("[%s] Unrecognized keys: %s", module_name, paste(unknown_keys, collapse = ", ")))
    user_controls = user_controls[!names(user_controls) %in% unknown_keys]
  }

  # update and return list ----
  return(
    modifyList(default_controls, user_controls)
  )

}

# .show_evo3d_defaults ----

#' show default evo3D controls
#'
#' prints the default evo3D control settings for one module
#' or all modules if none specified
#'
#' @param module_name character. one of "msa","pdb","aln","stat","output".
#'   if NULL (default) prints all
#'
#' @return none. prints to console
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

#' setup chain mapping for pdb–msa alignment
#'
#' standardizes chain input across pdbs and msas for main, interface,
#' or occlusion chains. ensures every pdb has an explicit mapping even
#' when user input is partial or symbolic (e.g. "auto", "all").
#'
#' @param chain character, vector, or list. chain input, e.g. "auto",
#'   "all", c("A","B"), or a list of vectors
#' @param chain_type character. one of "main","interface","occlusion"
#' @param pdb_count integer. total number of pdbs (default 1)
#' @param msa_count integer. number of msas per pdb (default 1)
#'
#' @return named list of length \code{pdb_count}, each entry a chain
#'   mapping vector or character
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

#' setup evo3D multi-run grid
#'
#' expands combinations of msa and pdb inputs into a run grid for multi-model,
#' multi-chain evo3D analyses. handles main chain assignment and resolves
#' homomultimer expansion where a chain string contains multiple ids.
#'
#' @param msa single msa (file path, matrix, or fasta object) or list of such
#' @param pdb single pdb (bio3d object) or list of pdbs
#' @param chain main chains of interest. character, vector, or list
#'   (e.g. "A", c("A","B"), list(c("A","B"),"C"))
#' @param interface_chain chains used for interface calculations. can be "all",
#'   character vector, or list
#' @param occlusion_chain chains used for occlusion in rsa/sasa. same rules as
#'   interface_chain
#'
#' @return list with standardized msa, pdb, interface_chain, occlusion_chain,
#'   and a run_grid data frame of msa–pdb–chain combinations
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
  multichain = which(nchar(run_grid$chain) > 1 & run_grid$chain != 'auto')
  if(length(multichain) > 0){
    expanded_rows = list()
    for(i in multichain) {
      chains = strsplit(run_grid$chain[i], "")[[1]]
      for(chain in chains) {
        new_row = run_grid[i, ]
        new_row$chain = chain
        expanded_rows = append(expanded_rows, list(new_row))
      }
    }

    # Remove original multichain rows and add expanded ones
    run_grid = rbind(run_grid[-multichain, ], do.call(rbind, expanded_rows))
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

#' safely resolve non-clashing output path
#'
#' avoids overwriting by appending a numeric or timestamp tag.
#' returns a safe path for saving a file or directory
#'
#' @param path character. desired path
#' @param is_dir logical. if TRUE, treat as directory (default FALSE)
#'
#' @return list with
#'   \item{path}{safe path string}
#'   \item{tag}{timestamp tag if fallback used, NULL otherwise}
#' @keywords internal
#
.safe_save = function(path, is_dir = FALSE) {
  base = tools::file_path_sans_ext(path)
  ext  = tools::file_ext(path)
  ext  = if (nzchar(ext)) paste0(".", ext) else ""

  # If safe, return original
  if (!file.exists(path)) return(list(path = path, tag = NULL))
  if (is_dir && length(list.files(path)) == 0) return(list(path = path, tag = NULL))
  if (!is_dir && file.info(path)$size == 0) return(list(path = path, tag = NULL))

  # Use systime + pid fallback
  systime = Sys.time()
  pid = Sys.getpid()
  tag = paste0(format(systime, "%Y%m%d%H%M%S"), "_pid", pid)
  tagged_path = if (is_dir) paste0(base, "_", tag) else paste0(base, "_", tag, ext)

  list(path = tagged_path, tag = tag)
}

# new run_evo3d() fixed size patches----

#' run evo3D end-to-end
#'
#' High-level pipeline that links sequence and structure to build codon-aligned,
#' structure-aware windows and (optionally) compute diversity/neutrality stats.
#' Internally calls the MSA, PDB, alignment, codon-collapse, and stats modules.
#'
#' @param msa MSA input(s): a file path, matrix, or \code{bio3d::read.fasta()} object;
#'   or a list of any mix of these (named automatically as \code{msa1}, \code{msa2}, …).
#' @param pdb PDB/mmCIF input(s): a \code{bio3d} PDB object or a file path; or a list of these
#'   (named \code{pdb1}, \code{pdb2}, …). Files are read via \code{.standardize_pdb_input()}.
#' @param chain Main chain(s) to analyze. Accepts:
#'   \itemize{
#'     \item \code{"auto"} (default) – auto-detect per MSA–PDB via k-mer coverage.
#'     \item single chain ID (\code{"A"}) or concatenated homomultimer string (\code{"AB"}) – expanded per PDB.
#'     \item list form for per-PDB/per-MSA control (see \code{.setup_chain_mapping()}).
#'   }
#' @param interface_chain Chains used to define interface patches. Can be \code{NA},
#'   \code{"all"}, a character vector, or list per PDB.
#' @param occlusion_chain Chains used to occlude solvent accessibility (RSA/SASA). Same forms as \code{interface_chain}.
#' @param analysis_mode One of \code{"codon"} (default). Controls downstream collapse mode.
#' @param calculate_stats Logical. If \code{TRUE}, computes per-window stats (pi, Tajima’s D, haplotype div., etc.).
#' @param detail_level Integer verbosity of returned intermediates:
#'   \itemize{
#'     \item \code{<0}: invisible, no return (write-to-disk only).
#'     \item \code{0}: minimal (skip large intermediates).
#'     \item \code{1} (default): keep \code{pdb} chains; drop some heavy objects.
#'     \item \code{2}: keep final MSA subsets.
#'     \item \code{3}: keep most intermediates.
#'   }
#' @param verbose Integer console messaging level (0–2+).
#' @param restart_run Reserved for future restart support (not implemented).
#' @param user_aln Optional pre-aligned input placeholder (not implemented; alignment is performed internally).
#'
#' @param msa_controls Named list of MSA-module controls. Typical fields:
#'   \code{ref_method} (\code{"consensus"}|\code{"most_complete"}|row index),
#'   \code{force_seq_type} (\code{"protein"}|\code{"nucleotide"}|NULL),
#'   \code{detect_sequence_threshold} (0–1), \code{detect_sequence_len},
#'   \code{reading_frame}, \code{reading_sens}, \code{genetic_code}.
#' @param pdb_controls Named list of PDB-module controls. Typical fields:
#'   \code{distance_method} (\code{"all"}|\code{"ca"}|\code{"backbone"}|\code{"sidechain"}),
#'   \code{drop_incomplete_residue}, \code{rsa_method} (\code{"rose"}|\code{"miller"}|\code{"theoretical_tien"}|\code{"empirical_tien"}),
#'   \code{dist_cutoff} (Å) or \code{max_patch} (integers; at least one must be non-\code{NA}),
#'   \code{rsa_cutoff}, \code{sasa_cutoff}, \code{only_exposed_in_patch}, \code{use_rsa_sasa} (\code{"and"}|\code{"or"}),
#'   \code{interface_dist_cutoff}, \code{force_file_type} (\code{"pdb"}|\code{"cif"}|NULL), \code{patch_mode} (\code{"codon"}|\code{"residue"}).
#' @param aln_controls Named list of alignment/merge controls. Typical fields:
#'   \code{merge_type} (\code{"distance"}|\code{"exposure_distance"}),
#'   \code{merge_exposure} (0–1 threshold for exposed consensus),
#'   \code{use_sample_names} (logical), \code{auto_chain_threshold} (k-mer coverage),
#'   \code{kmer_size} (integer).
#' @param stat_controls Named list of statistics controls. Typical fields:
#'   \code{calc_pi}, \code{calc_tajima}, \code{calc_hap}, \code{calc_patch_entropy},
#'   \code{calc_polymorphic}, \code{valid_aa_only} (exclude X/* etc. when TRUE).
#' @param output_controls Named list of output controls. Typical fields:
#'   \code{output_dir} (path or \code{NULL}), \code{prefix} (string),
#'   \code{write_msa_subsets}, \code{write_evo3d_df}, \code{write_call_info},
#'   \code{write_module_intermediates}.
#'
#' @details
#' \strong{Flow:}
#' \enumerate{
#'   \item MSA \textrightarrow{} reference peptide (\code{msa_to_ref()}).
#'   \item PDB \textrightarrow{} distances, accessibility, patches (\code{pdb_to_patch()}).
#'   \item Align ref–PDB, map residues \textrightarrow{} codons (\code{aln_msa_to_pdb()}).
#'   \item Merge multi-PDB/multimer contexts to codon windows (\code{collapse_to_codon()}).
#'   \item Extract codon-aligned nucleotide windows (\code{.extract_msa_subsets()}).
#'   \item (Optional) Compute stats (\code{run_pegas_three()}, \code{calculate_patch_entropy()},
#'         \code{calculate_polymorphic_residue()}).
#' }
#'
#' \strong{Auto chain detection:} \code{chain = "auto"} resolves per MSA–PDB using k-mer coverage
#' (see \code{aln_controls$kmer_size}, \code{aln_controls$auto_chain_threshold}). Passing chains are added per-PDB.
#'
#' \strong{3D windows:} Provide either \code{pdb_controls$dist_cutoff} (variable-size) or
#' \code{pdb_controls$max_patch} (fixed-size). Both may be set; fixed-size enforces quotas.
#'
#' \strong{Exposure logic:} RSA and SASA filters define seeds and (optionally) members.
#' \code{use_rsa_sasa = "and"} by default. SASA is recommended; RSA is provided for convenience.
#'
#' \strong{Not implemented (placeholders):} \code{restart_run}, \code{user_aln}.
#'
#' @return A list:
#' \describe{
#'   \item{evo3d_df}{Per-codon (and interface) table with \code{msa}, \code{codon}, \code{codon_id},
#'     reference AA, per-PDB AA/residue IDs, merged \code{codon_patch}, lengths, exposure, distances,
#'     and (if requested) stats (\code{pi}, \code{tajima}, \code{hap}, \code{patch_entropy}, polymorphism flags).}
#'   \item{final_msa_subsets}{Named list of nucleotide MSA windows per patch (may be \code{NULL} if \code{detail_level < 2}).}
#'   \item{msa_info_sets}{Per-MSA references and matrices.}
#'   \item{pdb_info_sets}{Per-PDB sequences, residue distances, accessibility, and patches.}
#'   \item{aln_info_sets}{Per-PDB alignment coverage and mapped tables.}
#'   \item{call_info}{Inputs, resolved run grid, and all control settings (for reproducibility).}
#' }
#'
#' @note If \code{output_controls$output_dir} is set, results can be written:
#'   MSA subsets (FASTA per patch), \code{evo3d_df.csv}, \code{call_info.json}, and intermediates \code{.rds}.
#'   Existing paths are made safe via \code{.safe_save()} to avoid overwrites.
#'
#' @export

run_evo3d = function(msa, pdb, chain = 'auto', interface_chain = NA, occlusion_chain = NA,    # how to handle pdb info #
                     analysis_mode = 'codon', calculate_stats = TRUE, detail_level = 1, verbose = 1, # how to handle analysis
                     restart_run = NULL, user_aln = NULL, # these need work -- how to restart a run?? # -- ... what does user need to give?
                     msa_controls = list(), pdb_controls = list(), aln_controls = list(),
                     stat_controls = list(), output_controls = list()){ # 5 control lists -- msa, pdb (patches), aln_controls -- how to merge, autochain thresh, stats, outputs


  # COULD ADD 1 more step -- validating input types (quick - is detail and verbose here ...) #
  # 6/26/25 need to add #

  #STEP 0 -- SETUP MODULE BEHAVIOR AND BUILD RUN INFO FROM USER INPUTS ----
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

  #STEP 0.5 -- quick stops to confine param space ----
  if(is.na(pdb_controls$dist_cutoff) && is.na(pdb_controls$max_patch)){
    message('STOPPING!! In pdb_controls: Must provide either a distance cutoff (dist_cutoff) or fixed size (max_patch) for 3D window defintion')
    return(NULL)
  }


  #STEP 1 MODULE 1 msa_to_ref() ----
  if(verbose > 0){
    cat('STEP 1: Converting MSAs to reference peptide sequences...\n')
  }

  # read MSA and convert to reference sequence for aligning to PDB #
  msa_info_sets = list()
  for(msa_name in names(run_info$msa)) {

    call_args = list(
      msa = run_info$msa[[msa_name]],
      verbose = verbose - 1)

    msa_info_sets[[msa_name]] = do.call(msa_to_ref, c(call_args, msa_controls))
  }

  #STEP 2 -- Caching PDBS, auto-detecting chains  ----

  # find matching chains -- for 'auto' rows #
  auto_ros = which(run_grid$chain == 'auto')
  run_grid$kmer_match = NA

  if (verbose > 0) {
    if (length(auto_ros) > 0) {
      cat("STEP 1.5: Caching PDBs and resolving auto chains...\n")
    } else {
      cat("STEP 1.5: Caching PDBs (no auto chains to resolve)...\n")
    }
  }

  # cache pdbs #
  pdb_cache = list()
  for(pdb_name in names(run_info$pdb)) {
    pdb_cache[[pdb_name]] = .standardize_pdb_input(pdb = run_info$pdb[[pdb_name]])
  }

  # all msa-pdb pairs get some kmer threshold #
  run_grid_ro = nrow(run_grid)
  for (i in seq_len(run_grid_ro)) {
    pdb_name = run_grid$pdb[i]
    msa_name = run_grid$msa[i]

    # Run auto detect
    chain_mappings = .auto_detect_chain(
      pep = msa_info_sets[[msa_name]]$pep,
      pdb = pdb_cache[[pdb_name]],
      in_module = TRUE,
      k = kmer_size
    )

    # if in auto -- try to pass autochain threshold
    if (i %in% auto_ros) {
      passing = chain_mappings[chain_mappings > auto_chain_threshold]

      if (length(passing) == 0) {
        if (verbose > 0) {
          message(sprintf(
            "In auto chain mapping: No chains in '%s' passed the k-mer threshold (%.2f) for '%s'.\n",
            pdb_name, auto_chain_threshold, msa_name
          ))
        }
        # mark as NA to show this msa pdb pairing had no matches
        run_grid$chain[i] = NA
        next
      }

      # add all passing chains as *new rows* at the bottom
      new_rows = data.frame(
        msa = msa_name,
        pdb = pdb_name,
        chain = names(passing),
        kmer_match = round(as.numeric(passing), 3),
        stringsAsFactors = FALSE
      )

      run_grid = rbind(run_grid, new_rows)
    }
    if (!(i %in% auto_ros)) {
      chain_choice = run_grid$chain[i]

      if (chain_choice %in% names(chain_mappings)) {
        # always record score, even if it's below threshold
        run_grid$kmer_match[i] = round(chain_mappings[chain_choice], 3)
      } else {
        stop(sprintf(
          "Chain '%s' was specified for PDB '%s' but not found.
       Please check your input or specify a valid chain.",
          chain_choice, run_grid$pdb[i]
        ))
      }
    }

  }

  # drop any remaining auto rows, but keep NA chains that replaced 'auto' #
  run_grid = run_grid[is.na(run_grid$chain) | run_grid$chain != 'auto', , drop = FALSE]

  # possible one pdb chain passed auto_threshold against two msa's #
  run_grid = run_grid[order(run_grid$pdb, run_grid$msa, run_grid$chain), , drop = FALSE]
  run_grid = run_grid[!duplicated(paste(run_grid$pdb, run_grid$chain)), , drop = FALSE]

  # reorder on msa and pdb
  run_grid = run_grid[order(run_grid$msa, run_grid$pdb), , drop = FALSE]
  row.names(run_grid) = NULL

  if (verbose > 1) {
    cat("\tRun grid for rest of analysis:\n\n")
    out = capture.output(print(run_grid))
    cat(paste0("\t", out), sep = "\n")
    cat('\n')
  }

  #STEP 3 -- VALIDATE input set !!! STOP RUN IF A PDB OR MSA DOESNT HAVE ANY MAPPINGS? !!! ----

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
                     kmer_size = kmer_size), # pack back in previously removed aln_controls #
    stat_controls = stat_controls,
    output_controls = output_controls
  )

  # stop if no mappings (nonsense data) #
  msa_valid = tapply(!is.na(run_grid$chain), run_grid$msa, any)
  bad_msas = names(msa_valid)[!msa_valid]

  pdb_valid = tapply(!is.na(run_grid$chain), run_grid$pdb, any)
  bad_pdbs = names(pdb_valid)[!pdb_valid]

  if (length(bad_msas) > 0 || length(bad_pdbs) > 0) {
    stop_msg = "!!! STOPPING EVO3D RUN !!!\nThe following MSAs or PDBs had no valid mappings:\n"
    if (length(bad_msas) > 0) {
      stop_msg = paste0(stop_msg, sprintf("MSAs: %s\n", paste(bad_msas, collapse = ", ")))
    }
    if (length(bad_pdbs) > 0) {
      stop_msg = paste0(stop_msg, sprintf("PDBs: %s\n", paste(bad_pdbs, collapse = ", ")))
    }
    stop_msg = paste0(stop_msg, "Please check inputs and call_info\n")
    stop_msg = paste0(stop_msg, "Not advised but set auto_chain_threshold = 0 to skip this check")
    message(paste0(stop_msg, "\n"))

    return(list(
      call_info = call_info
    ))
  }

  #STEP 4 MODULE 2 pdb_to_patch() ----

  if(verbose > 0) {
    cat('STEP 2: Converting PDBs to patches...\n')
  }

  # GATHER CHAIN, INTERFACE, AND OCCLUSION FOR EACH PDB #
  pdb_info_sets = list()
  for (pdb_name in unique(run_grid$pdb)) {
    # gather all relevant chains for this PDB #
    chain_set = run_grid$chain[run_grid$pdb == pdb_name]

    # remove NA and get unique chains #
    chain_set = chain_set[!is.na(chain_set)]
    chain_set = unique(chain_set)

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

  # CLEAR CACHE reduce memory usage ---- #

  # remove run_info and pdb_cache #
  rm(run_info, pdb_cache, call_args, chain_set, pdb_name, bad_pdbs, bad_msas, msa_valid, pdb_valid,
     chain_mappings, auto_chain_threshold, i, auto_ros, msa_name, interface_chain, occlusion_chain)

  invisible(gc())

  #STEP 5 MODULE 3 aln_msa_to_pdb ----
  if (verbose > 0) {
    cat('STEP 3: Aligning MSAs to PDBs...\n')
  }

  # now for each pdb -- resolve its residue windows into a codon window #
  aln_info_sets = list()
  pdb_names = unique(run_grid$pdb)
  for(i in seq_len(length(pdb_names))) {

    # gather all the msa's that this pdb touched
    pdb_name = pdb_names[i]
    sub_grid = run_grid[run_grid$pdb == pdb_name,]

    msa_set = unique(sub_grid$msa)

    call_args = list(
      msa_info = msa_info_sets[msa_set],
      pdb_info = pdb_info_sets[pdb_name][[1]],
      chain = NA,
      subset_msa = FALSE,
      verbose = verbose - 1,
      run_grid = sub_grid,
      drop_gap_map = TRUE,
      fix_gap_map_and_dedup_codons = TRUE,
      patch_mode = pdb_controls$patch_mode,
      max_patch = pdb_controls$max_patch,
      dist_cutoff = pdb_controls$dist_cutoff,
      only_exposed_in_patch = pdb_controls$only_exposed_in_patch
    )

    # aln_controls needs reworked - should have codon vs residue for patch mode
    aln_info_sets[[pdb_name]] = do.call(aln_msa_to_pdb, c(call_args, aln_controls))
  }

  # make residue df by rbinding all these aln_df's together #
  # thus every residue in analysis has a patch and codon window #

  residue_df = do.call(rbind, lapply(aln_info_sets, function(x) x$aln_df))
  # do i want to remove aln_df in aln_info_sets or leave for later inspection #
  # leave because individual aln_df maybe collapsed to codon level #

  # remove unneeded columns from residue_df #
  residue_df$patch = NULL
  residue_df$patch_len = NULL


  #STEP 6 MODULE 4 collapse_to_codon()  ----
  if(analysis_mode == 'codon'){
    # collapse_to_codon()
    # merge on distance or distance+exposure #
    # can build variable or fixed-size patches #
    residue_df = collapse_to_codon(residue_df = residue_df,
                                    merge_type = aln_controls$merge_type,
                                    merge_exposure = aln_controls$merge_exposure,
                                    max_patch = pdb_controls$max_patch,
                                    only_exposed_in_patch = pdb_controls$only_exposed_in_patch,
                                    patch_mode = pdb_controls$patch_mode,
                                    dist_cutoff = pdb_controls$dist_cutoff,
                                    merge_interface_surface = TRUE, # need to place this
                                    pdb_info_sets = pdb_info_sets
                                    )
  }

  #STEP 7 -- generate_msa_subsets (final_msa_subsets) # ----

  # generate set
  final_msa_subsets = .extract_msa_subsets(msa_info_sets, residue_df)


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

  #STEP 8 MODULE 5 calculate_patch_stats ----
  evo3d_df = residue_df

  if (!calculate_stats){
    if (verbose >= 1) {
      cat('STEP 4: Skipping stats calculation...\n')
    }
  } else {

    if (verbose >= 1) {
      cat('STEP 4: Calculating patch stats...\n')
    }

    # add polymorphic sites -- needs wrapped across msa subsets #
    # broken?
    if(stat_controls$calc_polymorphic){
      evo3d_df = calculate_polymorphic_residue(
        msa_info_sets,
        evo3d_df,
        valid_aa_only = stat_controls$valid_aa_only
      )
    }

    # very slow
    if(stat_controls$calc_patch_entropy){
      evo3d_df = calculate_patch_entropy(
        msa = final_msa_subsets,
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

    # also slow but works
    if(length(stat) > 0){
      evo3d_df = run_pegas_three(msa = final_msa_subsets,
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

  #5 saving to disk (work on saving) ----

  # reorder the columns of evo3d_df #
  if(F){
  # msa (if available), codon, msa_subset_id, ref_aa, pdbX_aa, pdbY_aa, ..., pdbX_residue_id, pdbY_residue_id, codon_patch, everything else #
  codon_info =  intersect(c("msa","codon","msa_subset_id","ref_aa"), names(evo3d_df))
  aa_cols = grep("^pdb.*_aa$", names(evo3d_df), value = TRUE)
  id_cols = grep(".*residue_id$", names(evo3d_df), value = TRUE)
  patch_col = "codon_patch"
  other = setdiff(names(evo3d_df), c(codon_info, aa_cols, id_cols, patch_col))
  col_order = c(codon_info, aa_cols, id_cols, patch_col, other)
  #evo3d_df = evo3d_df[, col_order, drop = FALSE]  ** drop 8/25/25 **
  }

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
    final_msa_subsets = NULL
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
    final_msa_subsets = final_msa_subsets,
    msa_info_sets = msa_info_sets,
    pdb_info_sets = pdb_info_sets,
    aln_info_sets = aln_info_sets,
    call_info = call_info
  )
  )

}


