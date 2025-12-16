# --------------------------------------------------------------- #
# MSA MODULE
# Utilities for processing multiple sequence alignment (MSA) data
# and returning standardized msa_info objects.
#
# Internal helper functions (.functions) are wrapped by the
# user-facing msa_to_ref() interface.
#
# NOTES:
# 1. Gaps in peptide reference sequences are represented as 'X'.
#    This distinguishes original alignment gaps from gaps introduced
#    later in the pipeline.
# 2. Consensus and most-complete reference methods consider only
#    valid DNA or amino acid characters.
#    For example, if a DNA column contains 25 'R' and 1 'A',
#    the consensus method returns 'A'.
# 3. Ambiguous characters (e.g., 'R', 'Y', 'B') are not resolved
#    to specific amino acids (e.g., 'GG[RBY]' -> Glycine).
#
# Contact: bbroyle@purdue.edu
# --------------------------------------------------------------- #

# .standardize_msa_input() ----

#' Standardize MSA Input
#'
#' Converts multiple sequence alignments provided in different formats into a
#' standardized character matrix suitable for evo3D analyses. Accepts FASTA
#' files, objects returned by \code{bio3d::read.fasta()}, or pre-loaded matrices.
#'
#' This function ensures that input alignments are in uppercase, contain unique
#' row names, and have the format \code{[samples × positions]} expected by all
#' downstream evo3D modules.
#'
#' @param msa Multiple sequence alignment, provided as:
#'   \itemize{
#'     \item A character string: file path to a FASTA alignment.
#'     \item A \code{bio3d} fasta object (from \code{bio3d::read.fasta()}).
#'     \item A character matrix with sequences as rows and alignment positions
#'           as columns.
#'   }
#'
#' @return A character matrix with sequences as rows and alignment positions as
#'   columns. Guarantees:
#'   \itemize{
#'     \item All letters are converted to uppercase.
#'     \item Row names are present and unique (defaulting to
#'           \code{"seq_1"}, \code{"seq_2"}, ... if not).
#'   }
#'
#' @examples
#' \dontrun{
#' # From FASTA file
#' msa_mat <- .standardize_msa_input("example_alignment.fasta")
#'
#' # From bio3d fasta object
#' fasta_obj <- bio3d::read.fasta("example_alignment.fasta")
#' msa_mat <- standardize_msa_input(fasta_obj)
#'
#' # From pre-loaded matrix (rows ~ samples, aligned positions ~ columns)
#' fasta_obj <- bio3d::read.fasta("example_alignment.fasta")
#' my_matrix <- fasta_obj$ali
#' msa_mat <- .standardize_msa_input(my_matrix)
#' }
#'
#' @export

.standardize_msa_input = function(msa){
  if(inherits(msa, 'fasta')){
    # object from bio3d::read.fasta()
    msa = msa$ali
  } else if(length(msa) == 1L && is.character(msa)){
    # file path #
    if (!file.exists(msa)) {
      stop("MSA file does not exist: ", msa)
    }
    msa = bio3d::read.fasta(msa, rm.dup = FALSE)$ali
  } else if(is.matrix(msa)){
    # nothing to adjust #
  } else {
    stop('MSA provided is not a recognized format.\n Please provide a file path, an object from bio3d::read.fasta(), or a matrix [samples, positions]')
  }

  # force uppercase #
  msa = toupper(msa)

  # force sample (row) names #
  if(is.null(rownames(msa))){
    rownames(msa) = paste0('seq_', seq_len(nrow(msa)))
    message('Row names were not provided in MSA Assigning default names: seq_1, seq_2, ...')
  }

  # check for unique row names #
  if(any(duplicated(rownames(msa)))){
    message('MSA has non-unique row names. Reassigning default names: seq_1, seq_2, ...')
    rownames(msa) = paste0('seq_', seq_len(nrow(msa)))
  }

  return(msa)
}

# .detect_sequence_type() ----

#' Detect Sequence Type
#'
#' determines whether a given sequence is nucleotide or protein based on the proportion of standard nucleotide characters in the first \code{detect_sequence_len} positions.
#' characters '-' and 'X' are dropped because they do not add information about either sequence type
#'
#' @param seq a character string representing a single biological sequence
#' @param detect_sequence_threshold proportion of characters that must be nucleotide-like (A, T, C, G) to classify the sequence as \code{"nucleotide"}. default is 0.8
#' @param detect_sequence_len maximum number of non-gap '-' or 'X' characters (from the start of the sequence) to consider when computing the proportion. default is 100
#'
#' @return a string: either \code{"nucleotide"} or \code{"protein"}
#' @keywords internal

.detect_sequence_type = function(seq, detect_sequence_threshold = 0.8, detect_sequence_len = 100) {

  # remove '-' and 'X' characters from the sequence
  seq = gsub("[-X]", "", seq)

  # if sequence was all gap or X - print message and return 'nucleotide' #
  if (nchar(seq) == 0L) {
    message("Automatic sequence detection failed (all gaps/X); defaulting to nucleotide")
    return("nucleotide")
  }

  # get sequence substring
  seq = substr(seq, 1, min(nchar(seq), detect_sequence_len))

  # Count how many are A/T/C/G
  # FUTURE implementation could accomidate 'U' for RNA
  nuc_like = sum(strsplit(seq, "")[[1]] %in% c("A", "T", "C", "G"))
  prop = nuc_like / nchar(seq)

  # Check if the proportion of nucleotide-like characters is above the threshold
  seq_type = if(prop >= detect_sequence_threshold) "nucleotide" else "protein"

  # return the sequence type
  return(seq_type)
}

# .get_reference_sequence() ----

#' Get Reference Sequence from MSA
#'
#' Extract a reference sequence from a multiple sequence alignment (MSA) using
#' one of three strategies: a specified row index, the most complete sequence,
#' or a consensus across all sequences.
#'
#' Consensus is built from the most frequent valid characters
#' (A/T/C/G for nucleotide; 20 standard amino acids for protein).
#' If no valid character is found in a column, the position becomes 'X'.
#'
#' @param msa A character matrix of sequences (rows) by alignment positions (columns).
#'   Should be standardized using \code{.standardize_msa_input()}.
#' @param ref_method Reference selection method: \code{"most_complete"},
#'   \code{"consensus"}, or a numeric row index.
#' @param force_seq_type Optional sequence type override: \code{"nucleotide"} or
#'   \code{"protein"}. If \code{NULL}, the type is auto-detected.
#' @param detect_sequence_threshold Proportion of A/T/C/G characters required to
#'   classify a sequence as nucleotide during auto-detection (default: 0.8).
#' @param detect_sequence_len Maximum number of non-gap (\code{-}) or \code{X}
#'   positions (from the start of the sequence) considered during auto-detection
#'   (default: 100).
#'
#' @return a list with two elements:
#'   \item{ref}{the reference sequence as a named character string}
#'   \item{seq_type}{sequence type, either \code{"nucleotide"} or \code{"protein"}}
#' @keywords internal

.get_reference_sequence = function(msa, ref_method = 'consensus', force_seq_type = NULL, detect_sequence_threshold = 0.8, detect_sequence_len = 100){
  # grab the reference sequence based on the method provided #

  # three passes #
  # 1. lenient (most_complete or consensus)
  # 2. detect sequence type
  # 3. strict (most_complete or consensus)
  # -- if method is numeric or force_seq_type is not NULL, skip the first two passes #

  # STEP 0 validating data ----
    # validate ref_method
    if (length(ref_method) != 1L) {
      stop('ref_method must be a single value')
    }

    if (is.numeric(ref_method)) {
      if (is.na(ref_method) || ref_method <= 0 || ref_method %% 1 != 0) {
        stop('numeric ref_method must be a positive integer row index')
      }
      if (ref_method > nrow(msa)) {
        stop('numeric ref_method exceeds number of sequences in MSA')
      }
    } else if (!ref_method %in% c('most_complete', 'consensus')) {
      stop('ref_method must be "most_complete", "consensus", or a numeric row index')
    }

  # check that seq_type is valid #
  if (!is.null(force_seq_type)) {
    if (length(force_seq_type) != 1L ||
        !force_seq_type %in% c("nucleotide", "protein")) {
      stop('force_seq_type must be NULL, "nucleotide", or "protein"')
    }
  }

  # STEP 1 & 2 lenient pass + detect sequence type ----
  # if method is most_complete or consensus and force_seq_type is NULL, do lenient first pass #
  if (is.null(force_seq_type) && ref_method %in% c('most_complete', 'consensus')) {

    avoid = c('X', '-')

    # if most_complete #
    if(ref_method == 'most_complete'){
      complete = apply(msa, 1, function(x) sum(!x %in% avoid))
      most_complete = which.max(complete)

      ref = msa[most_complete, ]
      ref = paste0(ref, collapse = '')
      names(ref) = paste0('ref.',rownames(msa)[most_complete])

    } else {
      # if consensus #
      chars = unique(as.vector(msa))
      chars = chars[!chars %in% avoid]

      if (length(chars) == 0L) {
        ref = rep("X", ncol(msa))
      } else {
        counts = vapply(chars, function(ch) colSums(msa == ch, na.rm = TRUE), numeric(ncol(msa)))
        ref = chars[max.col(counts, ties.method = "first")]
        ref[is.na(ref)] = "X"
      }

      ref = paste0(ref, collapse = '')
      names(ref) = 'ref.consensus'
    }

    # detect sequence type #
    seq_type = .detect_sequence_type(ref,
                                     detect_sequence_threshold = detect_sequence_threshold,
                                     detect_sequence_len = detect_sequence_len)
  }

  # step 3 strict pass ----

  # if method is numeric grab row -- detect_seq_type if needed #
  # otherwise construct most complete or consensus
  if(is.numeric(ref_method)){

    # return sequence by number
    ref = paste0(msa[ref_method,], collapse = '')
    names(ref) = paste0('ref.', rownames(msa)[ref_method])

    # if force_seq_type is NULL, detect seq_type #
    if(is.null(force_seq_type)){
      seq_type = .detect_sequence_type(ref,
                                       detect_sequence_threshold = detect_sequence_threshold,
                                       detect_sequence_len = detect_sequence_len)
    } else {
      seq_type = force_seq_type
    }

  } else {

    # strict pass for most_complete or consensus #
    # set up characters that count as complete #

    # if forcing seq_type, use that #
    if(!is.null(force_seq_type)){
      seq_type = force_seq_type
    }

    # use seq_type to get valid chars #
    if(seq_type == 'nucleotide'){
      chars = c('A', 'T', 'C', 'G')
    } else {
      chars = strsplit('AVILMWYFSTNQCGPRHKDE', '')[[1]]
    }

    # if most complete #
    if(ref_method == 'most_complete'){
      complete = apply(msa, 1, function(x) sum(x %in% chars))

      # return sequence with least_gaps
      most_complete = which.max(complete)
      ref = msa[most_complete, ]
      ref = paste0(ref, collapse = '')
      names(ref) = paste0('ref.',rownames(msa)[most_complete])
    } else {
      # method is consensus #
      counts = vapply(chars, function(ch) colSums(msa == ch, na.rm=TRUE), numeric(ncol(msa)))

      # add column X at 0.5 (if full column was gaps or ambiguous bases) -- take X #
      counts = cbind(counts, X = 0.5)
      chars = c(chars, 'X')

      # grab the max count for each column #
      consensus = chars[max.col(counts, ties.method="first")]
      ref = paste0(consensus, collapse = '')
      names(ref) = 'ref.consensus'
    }

  }

  # retun results ----
  return(list(ref = ref, seq_type = seq_type))
}

# .translate_dna_to_protein() ----

#' Translate DNA to amino acids
#'
#' Translate a DNA sequence to an amino acid sequence using \code{seqinr::translate()}.
#' Ambiguous codons (including gaps) translate to \code{"X"}. Internal stop codons are
#' reported via \code{message()}, but translation continues.
#'
#' @param seq A named character string containing a single DNA sequence.
#' @param frame Reading frame (0, 1, or 2). Fixed to 0 in the main pipeline, retained for flexibility.
#' @param sens Translation direction: \code{"F"} (forward) or \code{"R"} (reverse). Fixed to \code{"F"} in the main pipeline.
#' @param genetic_code NCBI genetic code number (default: 1, standard code).
#'
#' @return A named character string containing the translated amino acid sequence.
#' @keywords internal

.translate_dna_to_protein = function(seq, frame = 0, sens = 'F', genetic_code = 1){

  # translate (NNN and --- are treated the same 'X')
  pep = seqinr::translate(strsplit(seq, '')[[1]],
                          frame = frame,
                          sens = sens,
                          numcode = genetic_code,
                          ambiguous = FALSE
                          )

  # check for internal stops (can proceed)
  if (length(pep) > 1L && any(pep[-length(pep)] == "*")) {
    message("Internal stop codon(s) detected in the reference sequence; translation will proceed.\n Consider using a different ref_method or genetic_code.")
  }

  # return as character vector
  pep = paste0(pep, collapse = '')
  names(pep) = names(seq)

  return(pep)
}

# msa_to_ref() ----

#' Extract reference and peptide sequences from an MSA
#'
#' Standardize MSA input, extract a reference sequence, detect or apply sequence
#' type, and translate the reference to a peptide sequence if needed.
#'
#' This can be called directly in modular workflows, or internally by
#' \code{run_evo3d()}, where \code{reading_frame} and \code{reading_sens} are fixed.
#'
#' @param msa MSA input: a FASTA file path, a character matrix, or an object returned
#'   by \code{bio3d::read.fasta()}.
#' @param ref_method Reference selection method: \code{"most_complete"},
#'   \code{"consensus"}, or a numeric row index (default: \code{"consensus"}).
#' @param force_seq_type Optional sequence type override: \code{"protein"},
#'   \code{"nucleotide"}, or \code{NULL} to auto-detect (default).
#' @param verbose Integer; print progress messages if > 0 (default: 0).
#' @param detect_sequence_threshold Proportion of A/T/C/G required to classify as
#'   nucleotide during auto-detection (default: 0.8).
#' @param detect_sequence_len Number of leading characters used for auto-detection
#'   (default: 100).
#' @param reading_frame Reading frame for translation (0, 1, or 2). Fixed to 0 in
#'   \code{run_evo3d()}.
#' @param reading_sens Translation direction: \code{"F"} (forward) or \code{"R"}
#'   (reverse). Fixed to \code{"F"} in \code{run_evo3d()}.
#' @param genetic_code NCBI genetic code number (default: 1, standard code).
#'
#' @return A list with class \code{"evo3D_msa_info"} containing:
#' \itemize{
#'   \item \code{msa_mat}: Standardized alignment matrix.
#'   \item \code{ref}: Reference sequence (DNA or protein).
#'   \item \code{pep}: Peptide sequence (translated if nucleotide input; otherwise
#'     \code{ref} with \code{-} replaced by \code{X}).
#'   \item \code{seq_type}: Detected or forced sequence type.
#' }
#' @export

msa_to_ref = function(msa, ref_method = 'consensus', force_seq_type = NULL, verbose = 0,
                      genetic_code = 1,
                      detect_sequence_threshold = 0.8,
                      detect_sequence_len = 100,
                      reading_frame = 0,
                      reading_sens = 'F'){

  # standardize msa input ----

  if(verbose > 0) {
    cat('\tmsa_to_ref: Standardizing MSA input\n')
  }

  msa = .standardize_msa_input(msa)

  # grab reference sequence ----
  if(verbose > 0) {
    cat('\tmsa_to_ref: Extracting reference sequence\n')
  }

  ref = .get_reference_sequence(msa = msa,
                               ref_method = ref_method,
                               force_seq_type = force_seq_type,
                               detect_sequence_threshold = detect_sequence_threshold,
                               detect_sequence_len = detect_sequence_len
                              )

  # translate to protein sequence if needed ----
  if(verbose > 0) {
    cat('\tmsa_to_ref: Building reference peptide sequence\n')
  }

  if(ref$seq_type == 'nucleotide'){
    pep = .translate_dna_to_protein(ref$ref,
                                    frame = reading_frame,
                                    sens = reading_sens,
                                    genetic_code = genetic_code
                                   )
  } else {
    # replace '-' with 'X' for protein
    pep = ref$ref
    pep = gsub('-', 'X', pep)
  }

  # return results ----
  results = list(msa_mat = msa,
              ref = ref$ref,
              pep = pep,
              seq_type = ref$seq_type)

  class(results) = 'evo3D_msa_info'

  return(results)
}


