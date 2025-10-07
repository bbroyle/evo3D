# --------------------------------------------------------------- #
# ALN MODULE — align msa to pdb and derive codon-level windows
# internal .functions() are wrapped by exported aln_msa_to_pdb()
#
# NOTES:
# 1. alignment uses msa::msa(…, method = "ClustalOmega", type = "protein")
#    default sub matrix BLOSUM65; order = "input"
# 2. mismatches are computed on aa-to-aa positions only (ignores '-' and 'X')
# 3. .map_aln_to_positions() overlays indices:
#    ref codon positions and pdb residue_ids are added as extra columns
# 4. patches are converted to codon ids via .map_patches_to_codons()
#    codon_id format: "<codon>_<msaLabel>"; patch strings are '+'-delimited
# 5. gap-map handling:
#    - if a residue aligns to '-', it is considered gap-mapped
#    - drop_gap_map = TRUE removes gap-centered patches and prunes gap members
# 6. cleanup/rebuild step:
#    .fix_gap_map_and_dedup_codons() deduplicates codons and, when max_patch is set,
#    rebuilds patches from nearest valid residues using pdb_info$residue_dist
#    (optional dist_cutoff and only_exposed_in_patch respected)
# 7. patch modes:
#    - "codon": size/uniqueness tracked at codon level
#    - "residue": size built on residues, codons assigned afterward
# 8. interfaces:
#    residue_ids beginning with "interface" are held out during rebuild
#    and reattached unchanged at the end
# 9. multi-MSA support:
#    patches may span multiple MSAs (e.g., "1_msa1+2_msa1+1_msa2");
#    subsets are extracted per MSA and cbind()-merged
# 10. outputs:
#    - coverage: per-run alignment metadata (no aln matrix)
#    - aln_df: residue ↔ codon mapping with patch/codon strings and metadata
#    - msa_subsets: nucleotide windows per patch (samples × nucleotides)
#
# email: bbroyle@purdue.edu
# --------------------------------------------------------------- #

# .auto_detect_chain() ----

#' auto detect matching pdb chain
#'
#' compare a peptide sequence to each chain in a pdb structure
#' and return similarity scores based on k-mer coverage.
#' higher scores indicate closer matches between the query sequence and pdb chain.
#'
#' @param pep amino acid sequence as a single character string
#' @param pdb pdb structure input. can be a file path to a pdb, a pdb object
#'   from \code{bio3d::read.pdb()}, or a standardized pdb from \code{standardize_pdb_input()}.
#'   cannot be a pdb id alone.
#' @param k integer k-mer size used for coverage calculation (default 4)
#' @param in_module logical, for internal use. if TRUE, skips pdb validation
#'   (default FALSE)
#'
#' @return a named numeric vector of k-mer coverage values for each pdb chain,
#' sorted in descending order. names correspond to pdb chain identifiers.
#'
#' @details
#' coverage is defined as the fraction of pdb k-mers also present in the peptide sequence.
#'
#' @seealso \code{\link{find_matching_structures}}, \code{\link{standardize_pdb_input}},
#'   \code{bio3d::read.pdb}
#'
#' @examples
#' \dontrun{
#' pep = "MKTFFVAGVILLLVATATGVHS"
#' pdb = standardize_pdb_input("my_structure.pdb")
#' auto_detect_chain(pep, pdb)
#' }
#'
#' @export


.auto_detect_chain = function(pep, pdb, k = 4, in_module = F){

  # Changed to coverage instead of jaccard
  kmer_coverage = function(pdb_seq, msa_seq) {
    # seq is too short for kmer - just return 0
    if (nchar(pdb_seq) < k || nchar(msa_seq) < k) return(0)

    pdb_kmers = substring(pdb_seq, 1:(nchar(pdb_seq) - k + 1), k:(nchar(pdb_seq)))
    msa_kmers = substring(msa_seq, 1:(nchar(msa_seq) - k + 1), k:(nchar(msa_seq)))

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


# .calculate_coverage() ----

#' Calculate Aligned Coverage Ranges
#'
#' Identifies contiguous non-gap regions in each sequence of an alignment matrix
#' and formats them as start:end ranges. Optionally summarizes contiguous runs
#' of mismatch positions across the alignment.
#'
#' Called internally by \code{.align_sequences()}.
#'
#' @param aln_mat A character matrix representing a sequence alignment
#'   (rows = sequences, columns = positions).
#' @param mismatch Optional numeric vector of mismatch positions to summarize.
#'   If not supplied or empty, the returned \code{mismatch} element is \code{NA}.
#'
#' @return A named list. Each element corresponds to one sequence, containing a
#'   character vector of ranges (e.g., \code{"5:25"}). Includes an additional
#'   element \code{mismatch} if mismatches were provided.
#' @keywords internal

.calculate_coverage = function(aln_mat, mismatch) {
  covered_regions = list()

  # convert alignments into ranges
  for (i in 1:nrow(aln_mat)) {
    seq_name = rownames(aln_mat)[i]
    # ignore gap positions
    positions = which(aln_mat[i, ] != "-")

    if (length(positions) > 0) {
      # find breaks in sequence
      breaks = c(0, which(diff(positions) > 1), length(positions))
      ranges = character()

      for (j in 1:(length(breaks) - 1)) {
        start_idx = breaks[j] + 1
        end_idx = breaks[j + 1]
        start_pos = positions[start_idx]
        end_pos = positions[end_idx]
        ranges = c(ranges, paste0(start_pos, ":", end_pos))
      }

      covered_regions[[seq_name]] = ranges
    } else {
      covered_regions[[seq_name]] = character(0)
    }
  }

  # convert mismatches into ranges
  if (length(mismatch) > 0) {
    # find breaks in sequence
    breaks = c(0, which(diff(mismatch) > 1), length(mismatch))
    ranges = character()

    for (j in 1:(length(breaks) - 1)) {
      start_idx = breaks[j] + 1
      end_idx = breaks[j + 1]
      start_pos = mismatch[start_idx]
      end_pos = mismatch[end_idx]
      ranges = c(ranges, paste0(start_pos, ":", end_pos))
    }

    covered_regions$mismatch = ranges
  } else {
    covered_regions$mismatch = NA
  }


  return(covered_regions)
}

# .align_sequences() ----

#' Align Reference and Structure Sequences
#'
#' Aligns a reference amino acid sequence against a structure-derived sequence
#' using ClustalOmega (via the \pkg{msa} package, default BLOSUM65). Computes
#' positional mismatches and summarizes aligned coverage ranges.
#'
#' Intended for internal use in mapping reference sequences to PDB-derived
#' sequences prior to downstream patch analyses.
#'
#' @param sequences A named character vector of length two, containing the
#'   reference sequence first and the structure-derived sequence second.
#'   Example: \code{c(ref = "MKT...", pdb = "MK-...")}.
#' @param user_supplied_alignment Reserved for future use. Currently ignored.
#'
#' @return A list with components:
#' \describe{
#'   \item{aln_mat}{A two-column character matrix with aligned reference and
#'     structure sequences.}
#'   \item{coverage}{Named list of coverage ranges (see
#'     \code{.calculate_coverage}). Includes mismatch ranges if applicable.}
#' }
#' @seealso \code{.calculate_coverage}
#' @keywords internal

.align_sequences = function(sequences, user_supplied_alignment = NA){

  # user_supplied_alignment not supported yet -- soon #

  # need to track original gap (from msa)
  # then introduced gaps (from alignment)
  # also track matching positions
  if (!requireNamespace("msa", quietly = TRUE)) {
    stop(
      paste(
        "The 'msa' package is required for alignment.",
        "",
        "To install, run:",
        "  if (!requireNamespace(\"BiocManager\", quietly = TRUE))",
        "    install.packages(\"BiocManager\")",
        "  BiocManager::install(\"msa\")",
        sep = "\n"
      ),
      call. = FALSE
    )
  }

  # use msa::msa() ~ with baked in clustal omega (defualt to GONNET sub matrix)
  # could make invisible()
  aln = suppressMessages(msa::msa(sequences, method = 'ClustalOmega', type = 'protein', order = 'input', substitutionMatrix = 'BLOSUM65'))

  aln_set = as.character(aln@unmasked)
  aln_chars = strsplit(aln_set, '')

  # which position are identical
  mismatch = which(aln_chars[[1]] != aln_chars[[2]])

  # only keep mismatch if they are aa to aa (not - or X)
  mismatch = mismatch[!aln_chars[[1]][mismatch] %in% c('-', 'X')]
  mismatch = mismatch[!aln_chars[[2]][mismatch] %in% c('-')]

  # unpack alignment into matrix
  aln_mat = do.call(cbind, aln_chars)
  colnames(aln_mat) = c('ref_aa', 'pdb_aa')

  # calculate coverage and plot #
  coverage = .calculate_coverage(t(aln_mat), mismatch)

  # return alignment matrix and coverage plot
  return(list(
    aln_mat = aln_mat,
    coverage = coverage
  ))
}

# .map_aln_to_positions() ----

#' Map Alignment Columns to Codon and Residue Indices
#'
#' Annotates an alignment matrix by overlaying reference codon positions and
#' PDB residue indices onto aligned amino acids. This preserves the aligned
#' sequence characters while adding positional context for downstream mapping.
#'
#' Typically used after \code{.align_sequences()} and \code{identify_patches()}
#' to connect reference codons with structure-derived residues.
#'
#' @param aln_mat A character matrix returned by \code{.align_sequences()},
#'   with two rows (reference and structure-derived sequences).
#' @param residue_df Data frame of patch or interface residues, such as that
#'   returned by \code{identify_patches()}, containing at least
#'   \code{residue_id} and \code{orig_chain}.
#' @param chain Optional chain identifier (string or vector). If supplied,
#'   restricts mapping to residues belonging to the given chain(s).
#'
#' @return A modified character matrix with additional columns:
#' \describe{
#'   \item{codon}{Reference codon indices for aligned positions.}
#'   \item{residue_id}{PDB residue indices aligned to the structure sequence.}
#' }
#' @keywords internal

.map_aln_to_positions = function(aln_mat, residue_df, chain = NA){

  # aln_mat is from .align_sequences() aln_mat$aln_mat
  # residue_df is from pdb_to_patch()

  # filter for chain of interest #
  if(!any(is.na(chain))){
    residue_df = residue_df[residue_df$orig_chain %in% chain,]
  }

  # replace ref aa residues with codon position #
  # assuming first is ref sequence #

  # lets keep aa info in tack -- otherwise we just merge later #
  aln_mat = cbind(aln_mat, '-')

  ref_pos = which(aln_mat[,1]!='-')
  aln_mat[ref_pos,5] = 1:length(ref_pos)

  # for each set in residue_df update that row in aln_mat
  pos = residue_df$residue_id

  aln_mat = cbind(aln_mat, '-')
  # fill position in pdb seq with position #
  aa_pos = which(aln_mat[,2]!='-')
  aln_mat[aa_pos,6] = pos

  # rename
  colnames(aln_mat)[5:6] = c('codon', 'residue_id')

  return(aln_mat)
}

# .map_patches_to_codons() ----

#' Convert Residue Patches to Codon-Based Coordinates
#'
#' Maps structural patch definitions from residue space into codon-aligned
#' coordinates, using a MSA-PDB alignment table from
#' \code{.map_aln_to_positions()}. This produces codon-based patch identifiers
#' that can be used for downstream nucleotide-level extraction.
#'
#' Each residue inherits codon patch membership, with options to remove
#' gap-mapped residues. When gaps are dropped, both gap-seeded patches and
#' membership of gap residues in other patches are removed. If a residue
#' distance matrix is supplied, \code{max_dist} is recalculated after
#' gap-removal.
#'
#' @param aln_table Alignment table from \code{.map_aln_to_positions()},
#'   containing at least columns \code{codon}, \code{msa}, \code{pdb},
#'   \code{residue_id}, \code{ref_aa}, and \code{pdb_aa}.
#' @param residue_df Data frame of patch-level residues (e.g. from
#'   \code{pdb_to_patch()}), with columns including \code{residue_id},
#'   \code{patch}, \code{patch_len}, \code{exposed}, and \code{max_dist}.
#' @param drop_gap_map Logical; whether to drop gap-mapped residues. If
#'   \code{TRUE} (default), patches centered on gap residues are removed,
#'   and gap residues are excluded from other patches.
#' @param dist_mat Optional residue–residue distance matrix (from
#'   \code{pdb_to_patch()}). If supplied and \code{drop_gap_map = TRUE},
#'   \code{max_dist} is updated after gap removal.
#'
#' @return A data frame aligning residues with codon-based patch information,
#'   including:
#'   \describe{
#'     \item{codon_id}{Codon identifier string (e.g. \code{"15_msa1"}).}
#'     \item{codon_patch}{“+”-delimited string of codon IDs making up the patch.}
#'     \item{codon_len}{Number of codons in the patch.}
#'     \item{unique_codon}{Number of unique codons represented.}
#'     \item{gap_map_count}{Number of gap-mapped residues removed from the patch.}
#'     \item{patch, patch_len}{Residue-based patch membership (cleaned).}
#'     \item{exposed, max_dist}{Exposure and distance metadata (updated if applicable).}
#'   }
#' @keywords internal

.map_patches_to_codons = function(aln_table, residue_df, drop_gap_map = TRUE, dist_mat = NA){

  # in aln_msa_to_pdb() always drop gap map #
  # in run_evo3d() always drop gap map, but dont touch dist_mat - later function #
  # as manual function drop_gap_map can be turned off #

  # copy over original residue patch, patch len, exposure, and max_dist #
  aln_table$orig_patch = residue_df$patch[match(aln_table$residue_id, residue_df$residue_id)]
  aln_table$orig_patch_len = residue_df$patch_len[match(aln_table$residue_id, residue_df$residue_id)]
  aln_table$exposed = residue_df$exposed[match(aln_table$residue_id, residue_df$residue_id)]
  aln_table$max_dist = residue_df$max_dist[match(aln_table$residue_id, residue_df$residue_id)]

  # Open new columns
  aln_table$codon_patch = NA
  aln_table$codon_len = NA
  aln_table$unique_codon = NA
  aln_table$gap_map_count = NA

  aln_table$patch = NA
  aln_table$patch_len = NA

  # Set up codon_id
  aln_table$codon_id = paste0(aln_table$codon, '_', aln_table$msa)
  aln_table$codon_id[is.na(aln_table$codon)] = NA

  # check for gap_map residues #
  gap_ids = aln_table$residue_id[which(aln_table$codon == '-')]

  # if gap_map is dropped remove gap_map seeds #
  if(drop_gap_map){
    aln_table$orig_patch[aln_table$residue_id %in% gap_ids] = NA
  }

  # for each row - copy patch, *remove gaps*, return codons #
  for(i in seq_len(nrow(aln_table))){
    if(is.na(aln_table$orig_patch[i])) next

    p = strsplit(aln_table$orig_patch[i], '\\+')[[1]]

    # count gap map #
    g = which(p %in% gap_ids)

    if(drop_gap_map){
      p = setdiff(p, gap_ids)
    }

    # find codon positions #
    c = aln_table$codon_id[match(p, aln_table$residue_id)]

    # fill table #
    aln_table$patch[i] = paste0(p, collapse = '+')
    aln_table$patch_len[i] = length(p)

    aln_table$codon_patch[i] = paste0(c, collapse = '+')
    aln_table$codon_len[i] = length(c)
    aln_table$unique_codon[i] = length(unique(c))

    aln_table$gap_map_count[i] = length(g)

    # update max_dist if available #
    if(length(g) > 0 && is.matrix(dist_mat) && drop_gap_map){
      # replace gap map dist #
      id = aln_table$residue_id[i]

      if(grepl('^interface', id)) next

      aln_table$max_dist[i] = max(dist_mat[id, p], na.rm = T)
    }

  }

  # clean up table #

  # drop original patch info #
  aln_table = aln_table[ , !(names(aln_table) %in% c("orig_patch", "orig_patch_len")) ]

  # reorder rest of columns #
  # codon_id, codon, msa, pdb, residue_id, ref_aa, pdb_aa, codon_patch,
  # codon_len, unique_codon, gap_map_count, exposed, max_dist, patch, patch_len #

  # and anything else #
  # your desired order
  wanted_order = c(
    "codon_id", "codon", "msa", "pdb", "residue_id", "ref_aa", "pdb_aa",
    "codon_patch", "codon_len", "unique_codon", "gap_map_count",
    "exposed", "max_dist", "patch", "patch_len"
  )

  # columns in data that match your desired order
  present = intersect(wanted_order, names(aln_table))

  # everything else (not in your desired order)
  extra = setdiff(names(aln_table), wanted_order)

  # final reordering
  aln_table = aln_table[ , c(present, extra)]


  # last clean up if codon_patch is '' set all codon info to NA #
  # possible they are rebuilt - so save them #
  #aln_table$codon_patch[aln_table$codon_patch == ''] = NA
  #aln_table$codon_len[aln_table$codon_patch == ''] = NA
  #aln_table$unique_codon[aln_table$codon_patch == ''] = NA

  return(aln_table)
}

# .extract_msa_subsets_single() ----

#' Extract a Codon-Aligned Subset from a single MSA
#'
#' Internal helper to \code{.extract_msa_subsets()}. Converts codon indices
#' (without MSA labels) into nucleotide positions (3 bases per codon) and
#' subsets a single MSA accordingly.
#'
#' @param msa A nucleotide alignment matrix (samples × sites).
#' @param codon_patches Data frame with \code{msa_subset_id} and
#'   \code{codon_patch} (a “+”-delimited string of codon indices).
#'
#' @return A list of one or more nucleotide MSA subsets, each labeled by
#'   \code{msa_subset_id}.
#' @keywords internal

.extract_msa_subsets_single = function(msa, codon_patches){
  # works across MSA's #
  # drop patches that dont have codon positions #
  patches = codon_patches[!is.na(codon_patches$codon_patch),]

  msa_subset = list()
  for(i in 1:nrow(patches)){
    codons = patches$codon_patch[i]
    codon_pos = unlist(strsplit(codons, '\\+'))

    # convert codon positions to nucleotide positions #
    nuc_pos = c()
    for(codon in codon_pos){
      codon_num = as.numeric(codon)
      # each codon spans 3 nucleotides: codon 1 = nucs 1:3, codon 2 = nucs 4:6, etc.
      nuc_range = ((codon_num - 1) * 3 + 1):(codon_num * 3)
      nuc_pos = c(nuc_pos, nuc_range)
    }

    # subset #
    fasta = msa[, nuc_pos]
    msa_subset[[i]] = fasta

    # use the pre-computed msa_subset_id
    names(msa_subset)[i] = patches$msa_subset_id[i]
  }

  # return list of subsets #
  return(msa_subset)
}

# .extract_msa_subsets() ----

#' Extract Codon-Aligned Nucleotide Subsets Across Multiple MSAs
#'
#' Subsets one or more nucleotide multiple sequence alignments into
#' codon-aligned windows corresponding to structural patches. Codon patches
#' may reference multiple MSAs (e.g., \code{"1_msa1+2_msa1+1_msa2"}), in which
#' case codons are extracted from each MSA separately and then merged.
#'
#' Typically used after \code{.map_patches_to_codons()} to build
#' nucleotide-level MSAs for patch-based diversity and selection analyses.
#'
#' @param msa_info A named list of MSA entries. Each entry must contain
#'   \code{msa_mat}, a nucleotide alignment matrix (samples × sites), usually
#'   from \code{ape::as.DNAbin()} or similar.
#' @param codon_patches Data frame with at least \code{codon_patch} and
#'   \code{msa_subset_id} columns, as produced by
#'   \code{.map_patches_to_codons()}. Codon patches are represented as
#'   “+”-delimited strings of codon indices, each suffixed by its MSA label
#'   (e.g., \code{"15_msa1+22_msa1+7_msa2"}).
#' @param use_sample_names Logical; whether to collapse cross-MSA or cross-chain
#'   patches so that resulting subsets share the same FASTA identifiers.
#'   Default is \code{TRUE}.
#'
#' @return A named list of nucleotide alignment subsets, one per patch. Each
#'   element is a merged matrix (samples × nucleotides) containing the codons
#'   drawn from one or more MSAs according to the patch definition.
#' @keywords internal

.extract_msa_subsets = function(msa_info, codon_patches, use_sample_names = TRUE){
  # works across MSA's #
  # schedules extract_msa_subsets_single #

  # still need to incorporate use_sample_names #

  # two options --
  # 1. pull apart patches into msa1, msa2, msaX columns and schedule blocks to extract_msa_subsets #
  # 2. send row by row -- still pulling apart #
  # -- doing option 2

  # add new column for each identifier #
  msa_subsets = list()
  for(i in seq_len(nrow(codon_patches))){
    msa_subsets[[i]] = NULL
    p = codon_patches$codon_patch[i]
    if(is.na(p)) next

    p = strsplit(p, '\\+')[[1]]
    id = codon_patches$msa_subset_id[i]

    # separate codons from msa_ids # (assumes id always present)
    codons = gsub('_.+', '', p)
    msas = gsub('.+_', '', p)

    msa_set = unique(msas)
    msa_list = list()
    for(j in seq_len(length(msa_set))){
      msa_name = msa_set[j]
      codons2 = codons[msas == msa_name]
      msa_list[[j]] = .extract_msa_subsets_single(msa_info[[msa_name]]$msa_mat, data.frame(msa_subset_id = id, codon_patch = paste(codons2, collapse = '+')))[[1]]
    }

    # right here -- there needs to be a checking for sample_names
    msa_subset = do.call(cbind, msa_list)

    msa_subsets[[i]] = msa_subset
    names(msa_subsets)[i] = id

  }

  # clean up set to only the ones with pulls #
  keep = !vapply(msa_subsets, is.null, logical(1))
  msa_subsets = msa_subsets[keep]

  # just return
  return(msa_subsets)
}


# .fix_gap_map_and_dedup_codons ----

#' Fix Gap-Mapped Residues and Deduplicate Codons in Patches
#'
#' Cleans and rebuilds patch definitions by removing gap-mapped residues and
#' ensuring codons are not duplicated within codon-based patches. When a maximum
#' patch size (\code{max_patch}) is specified, patches are rebuilt from the
#' residue distance matrix to satisfy size constraints.
#'
#' Behavior depends on \code{patch_mode}:
#' \itemize{
#'   \item In \code{"codon"} mode, codons within each \code{codon_patch} are
#'   deduplicated. If \code{max_patch} is set, patches smaller than the target
#'   are rebuilt using nearest valid codons until the limit is reached.
#'   \item In \code{"residue"} mode, patches are rebuilt at the residue level,
#'   with codons reassigned from residue membership.
#' }
#'
#' Interface residues (with IDs beginning \code{"interface"}) are held out and
#' reattached after processing.
#'
#' @param residue_df Data frame of patch-level residues, typically from
#'   \code{.map_patches_to_codons()}, with columns \code{residue_id},
#'   \code{patch}, \code{codon_patch}, \code{codon}, \code{pdb_aa},
#'   \code{patch_len}, \code{codon_len}, \code{unique_codon}, and
#'   \code{exposed}.
#' @param patch_mode Either \code{"codon"} or \code{"residue"}, indicating how
#'   patches are rebuilt and deduplicated.
#' @param dist_cutoff Numeric; optional maximum residue–residue distance used
#'   when rebuilding patches. Distances greater than this are excluded.
#' @param max_patch Integer; maximum allowed patch size. If \code{NA}, no patch
#'   rebuilding is performed.
#' @param pdb_info PDB metadata list containing at least \code{residue_dist}, a
#'   residue–residue distance matrix.
#' @param only_exposed_in_patch Logical; whether only exposed residues can be
#'   included when rebuilding patches.
#'
#' @return The input \code{residue_df} with updated patch definitions. Columns
#'   \code{patch}, \code{patch_len}, \code{codon_patch}, \code{codon_len}, and
#'   \code{unique_codon} are rebuilt as needed. Interface residues are preserved
#'   and reattached unchanged.
#' @keywords internal

.fix_gap_map_and_dedup_codons = function(residue_df, patch_mode, dist_cutoff,
                                    max_patch, pdb_info, only_exposed_in_patch) {

  # if patch_mode = codon, dedup whole table #
  if (patch_mode == "codon") {
    # Deduplicate codons within each codon_patch
    for (i in seq_len(nrow(residue_df))) {
      p = residue_df$codon_patch[i]
      if (is.na(p)) next
      p = unique(strsplit(p, "\\+")[[1]])
      residue_df$unique_codon[i] = length(p)
      residue_df$codon_patch[i] = paste(p, collapse = "+")
    }
  }


  # ---- Variable patch scheme ----
  if (is.na(max_patch)) {
    # nothing to do
    return(residue_df)
  }

  # ---- Fixed patch scheme (max_patch is set) ----

  # Step 0: separate interfaces
  ro = grep("^interface", residue_df$residue_id)
  interface_df = NULL
  if (length(ro)) {
    interface_df = residue_df[ro, ]
    residue_df = residue_df[-ro, ]
  }

  # Step 1: rebuild patches depending on mode
  if (patch_mode == "codon") {
    ro = which(!is.na(residue_df$codon_patch) & residue_df$unique_codon < max_patch)

    if(length(ro)){
      # grab valid replacements #
      valid = residue_df[residue_df$codon != "-" &
                            residue_df$pdb_aa != "-" &
                            (if (only_exposed_in_patch) residue_df$exposed else TRUE),
                          c("codon_id", "residue_id"), drop = FALSE]

      # loop through each row and rebuild with valid sets #
      if(nrow(valid)){

        # trim dist_mat
        d = pdb_info$residue_dist
        d = d[, valid$residue_id, drop = FALSE]

        # make mapping
        res2cod = setNames(valid$codon_id, valid$residue_id)

        # loop through ro #
        for (i in ro) {
          resi = residue_df$residue_id[i]
          if (!resi %in% rownames(d)) next

          d2 = sort(d[resi, , drop = TRUE])
          if (!is.na(dist_cutoff)) d2 = d2[d2 <= dist_cutoff]
          if (!length(d2)) next

          cod_seq = res2cod[names(d2)]
          uniq_cod = cod_seq[!duplicated(cod_seq)]
          set = uniq_cod[seq_len(min(max_patch, length(uniq_cod)))]

          residue_df$patch_len[i] = length(set)
          residue_df$patch[i] = paste(names(set), collapse = "+")
          residue_df$codon_patch[i] = paste(set, collapse = "+")
          residue_df$codon_len[i] = length(set)
          residue_df$unique_codon[i]= length(unique(set)) # codon_len and unique_codon should be the same #
        }
      }


      }
    } else if (patch_mode == "residue") {
      ro = which(!is.na(residue_df$codon_patch) & residue_df$codon_len < max_patch)

      if(length(ro)){
        # grab valid replacements #
        valid = residue_df[residue_df$codon != "-" &
                              residue_df$pdb_aa != "-" &
                              (if (only_exposed_in_patch) residue_df$exposed else TRUE),
                            c("codon_id", "residue_id"), drop = FALSE]

        # loop through each row and rebuild with valid sets #
        if(nrow(valid)){

          # trim dist_mat
          d = pdb_info$residue_dist
          d = d[, valid$residue_id, drop = FALSE]

          # make mapping
          res2cod = setNames(valid$codon_id, valid$residue_id)

          # loop through ro #
          for (i in ro) {
            resi = residue_df$residue_id[i]
            if (!resi %in% rownames(d)) next

            d2 = sort(d[resi, , drop = TRUE])
            if (!is.na(dist_cutoff)) d2 = d2[d2 <= dist_cutoff]
            if (!length(d2)) next

            take_k = min(length(d2), max_patch)
            set = names(d2)[seq_len(take_k)]

            residue_df$patch_len[i] = length(set)
            residue_df$patch[i] = paste(set, collapse = "+")
            residue_df$codon_patch[i] = paste(res2cod[set], collapse = "+")
            residue_df$codon_len[i] = length(res2cod[set])
            residue_df$unique_codon[i] = length(unique(res2cod[set]))
          }
        }
      }
    }

  # Step 2: reattach interfaces
  if (!is.null(interface_df)) {
    residue_df = rbind(residue_df, interface_df)
  }

  return(residue_df)
}
# aln_msa_to_pdb() ----

#' Align Structure to Reference MSA and Generate Codon-Level Windows
#'
#' Aligns a structure-derived protein sequence to the reference sequence in one
#' or more nucleotide MSAs, maps structural patches to codon positions, and
#' extracts corresponding nucleotide windows. This function links 3D structure
#' to sequence diversity by chaining together multiple submodules:
#' \enumerate{
#'   \item Sequence alignment of PDB vs. MSA reference (\code{.align_sequences()}).
#'   \item Mapping aligned positions to codon and residue indices
#'         (\code{.map_aln_to_positions()}).
#'   \item Converting structural patches to codon-level definitions
#'         (\code{.map_patches_to_codons()}).
#'   \item Optional patch cleanup: gap-map removal and codon deduplication
#'         (\code{.fix_gap_map_and_dedup_codons()}).
#'   \item Extraction of nucleotide MSA subsets for each patch
#'         (\code{.extract_msa_subsets()}).
#' }
#'
#' This is the core alignment module in the evo3D workflow and is called
#' internally by \code{run_evo3d()}.
#'
#' @param msa_info A named list from \code{msa_to_ref()}, containing
#'   at least \code{msa_mat} (nucleotide MSA matrix) and \code{pep}
#'   (reference peptide sequence).
#' @param pdb_info A named list from \code{WRAPPER_pdb_to_patch()}, including
#'   \code{seq_set} (PDB-derived peptide sequences), \code{residue_df}
#'   (residue/patch annotations), and \code{residue_dist} (residue distance
#'   matrix).
#' @param chain Character string specifying the PDB chain to analyze. Ignored if
#'   \code{run_grid} is provided.
#' @param subset_msa Logical; if \code{TRUE}, nucleotide MSA subsets are
#'   extracted for each patch. Default \code{TRUE}.
#' @param verbose Integer; level of console reporting. Default \code{1}.
#' @param run_grid Optional data frame specifying combinations of MSA, PDB, and
#'   chain to process. If supplied, overrides \code{chain}.
#' @param drop_gap_map Logical; if \code{TRUE}, drops patches centered on
#'   gap-mapped residues. Default \code{TRUE}.
#' @param fix_gap_map_and_dedup_codons Logical; if \code{TRUE}, calls
#'   \code{.fix_gap_map_and_dedup_codons()} to remove gap residues from patch
#'   membership and ensure codon deduplication. Default \code{TRUE}.
#' @param patch_mode Either \code{"codon"} or \code{"residue"}, controlling
#'   how patches are deduplicated and rebuilt.
#' @param max_patch Integer; maximum patch size. If \code{NA}, patches are
#'   variable-length.
#' @param dist_cutoff Numeric; maximum residue–residue distance when expanding
#'   patches. Default \code{NULL} (no cutoff).
#' @param only_exposed_in_patch Logical; if \code{TRUE}, restricts patch
#'   rebuilding to exposed residues only.
#' @param merge_type Character; strategy for merging overlapping patches.
#'   Typically \code{"hold"}. Reserved for internal use.
#' @param merge_exposure Character; strategy for merging exposure information.
#'   Typically \code{"hold"}. Reserved for internal use.
#'
#' @return A list with three components:
#' \describe{
#'   \item{coverage}{List of alignment metadata for each MSA–PDB–chain run.}
#'   \item{aln_df}{Data frame mapping residues to codons, amino acids, and patch
#'         identifiers, after gap/deduplication cleanup.}
#'   \item{msa_subsets}{Named list of nucleotide MSA matrices, one per patch
#'         (if \code{subset_msa = TRUE}).}
#' }
#' @export

aln_msa_to_pdb = function(msa_info, pdb_info, chain = 'auto',
                           subset_msa = TRUE,
                           verbose = 1,
                           run_grid = NULL,
                           drop_gap_map = TRUE, # aln will drop patches centered on gap map residues
                           fix_gap_map_and_dedup_codons = TRUE, # aln will drop gap map residues from patch membership
                           patch_mode, max_patch, dist_cutoff, only_exposed_in_patch,
                           merge_type = 'hold', merge_exposure = 'hold'){

  # msa_info ~ must be list object from msa_to_ref() # ~ can be list of msa_infos
  # pdb_info ~ must be list object from ONE pdb_to_patch() #

  # step 0: prep data ----
  # in module always a run grid -- outside module can be chain specified #

  # If run grid supplied - ignore chain info #
  if(!is.null(run_grid)){
    chain = 'run_grid'
  }

  # at the end lets make sure to have a run_grid format
  if(F){
    # probably dont need chain mapping here 9.30.25 #
    # chain must be specified or 'auto'
    if(chain == 'auto') {
      chain = .auto_detect_chain(msa_info$pep, pdb_info$pdb)
      cat('\tDetected chain:', names(chain)[1], '\n')
      cat('\tAt', round(chain[1], 3) * 100, '% identity', '\n')
      chain = names(chain)[1]  # use the first chain detected
    }
  }

  # step 1: align sequences ----
  # run alignment for each row of run_grid
  if(verbose > 0){
    cat('\taln_msa_to_pdb: aligning msa(s) to pdb ')
  }

  aln_sets = list()
  for(i in 1:nrow(run_grid)){
    pep = msa_info[[run_grid$msa[i]]]$pep
    seq = pdb_info$seq_set[run_grid$chain[i]]

    if(any(is.na(seq))){
      stop(paste0('No sequence found for chain ', chain, '. Check PDB input.'))
    }

    aln_sets[[i]] = .align_sequences(c(pep, seq))

    # add msa and pdb info here (pdb_id and msa_id can pull from wrapper level -- as module it just is pdb1, and msa1-n)
    aln_sets[[i]]$aln_mat = cbind(aln_sets[[i]]$aln_mat, msa = run_grid$msa[i], pdb = run_grid$pdb[i])

    names(aln_sets)[i] = paste0(run_grid$msa[i], '_', run_grid$pdb[i], '_', run_grid$chain[i])
  }

  # could stack same chain now -- or later #

  # step 2: map alignment to positions ----
  if(verbose > 0){
    cat('\taln_msa_to_pdb: mapping alignment to codon and pdb positions\n')
  }

  # grab residue df #
  residue_df = pdb_info$residue_df
  pos_sets = list()
  for(i in 1:nrow(run_grid)){
    pos_sets[[i]] = .map_aln_to_positions(aln_sets[[i]]$aln_mat, residue_df, chain = run_grid$chain[i])
    names(pos_sets)[i] = paste0(run_grid$msa[i], '_', run_grid$pdb[i], '_', run_grid$chain[i])
  }

  # step 2.5: create large rbind() list -- every residue_id has row ----
  aln_table = do.call(rbind, pos_sets)
  aln_table = as.data.frame(aln_table, stringsAsFactors = FALSE)

  # at this point need to recover "interfaces" from residue_df and place them as residue_id #
  # any interfaces or drug binding sites (can be coded as interface) #
  ro = grep('^interface', residue_df$residue_id)
  if(length(ro)>0){
    interface_table = data.frame(
      residue_id = residue_df$residue_id[ro],
      ref_aa = NA,
      pdb_aa = NA,
      msa = NA,
      pdb = run_grid$pdb[1],
      codon = NA
    )

    aln_table = rbind(aln_table, interface_table)
  }


  # step 3: map patches to codons ----
  if(verbose > 0){
    cat('\taln_msa_to_pdb: converting pdb patches to codon\n')
  }

  codon_patches = .map_patches_to_codons(aln_table, residue_df,
                                          drop_gap_map = TRUE)

  # STEP 4: add MSA subset ids ----
  # here tied to residue_id #
  codon_patches$msa_subset_id = paste0(codon_patches$residue_id, '_', codon_patches$pdb)
  codon_patches$msa_subset_id[is.na(codon_patches$codon_patch)] = NA

  # STEP 5: fix gap map and dedup codons ----
  if(fix_gap_map_and_dedup_codons){
    codon_patches = .fix_gap_map_and_dedup_codons(residue_df = codon_patches,
                                             patch_mode = patch_mode,
                                             max_patch = max_patch,
                                             dist_cutoff = dist_cutoff,
                                             only_exposed_in_patch = only_exposed_in_patch,
                                             pdb_info = pdb_info)
  }

  table(codon_patches$codon_len)
  table(codon_patches$unique_codon)

  # step 4: grab subsets of MSA ----
  msa_subsets = NULL
  if(subset_msa){
    msa_subsets = .extract_msa_subsets(msa_info, codon_patches)
  }

  # step 5: gather coverage data for individual aliments and return ----
  aln_sets = lapply(aln_sets, function(x) {
    x$aln_mat = NULL
    x
  })

  # also fixed exposure to FALSE if NA
  #codon_patches$exposed[is.na(codon_patches$exposed)] = FALSE

  # return data #
  return(list(
    coverage = aln_sets,
    aln_df = codon_patches,
    msa_subsets = msa_subsets)
  )

}


