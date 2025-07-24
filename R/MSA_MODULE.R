# --------------------------------------------------------------- #
# MSA MODULE -- utilities for taking multiple sequence alignment (MSA) data and returning standardized msa_info #
# internal .functions() are wrapped with external msa_to_ref() module
#
# NOTES:
# 1. gaps in peptide reference are represented as 'X' -- this helps distinguish from new introduced gaps later in the pipeline #
# 2. Consensus and most complete reference methods only count valid DNA and AA characters
#   -- meaning if a DNA column has 25 'R' and 1 'A' the consensus method returns 'A' for that column
# 3. No attempt is made to salvage ambigious characters 'GG[RBY]' -> 'Glycine'
#
# email me: bbroyle@purdue.edu
# --------------------------------------------------------------- #

#' Standardize MSA Input
#'
#' Accepts various forms of MSA input (file path, matrix, or bio3d fasta object) and standardizes it into a character matrix.
#'
#' @param msa A character file path to a FASTA file, a matrix with sequences as rows and alignment positions as columns, or an object returned by \code{bio3d::read.fasta()}.
#'
#' @return A character matrix with sequences as rows and alignment positions as columns. Letter case is standardized to uppercase, and row names are assigned if not present.
#' @keywords internal
.standardize_msa_input = function(msa){
  # Take a variety of msa input types and return standardized matrix #
  # expecting:
  # 1. file path (character)
  # -- or two formats from bio3d::read.fasta() -- #
  # 2. fasta object
  # 3. fasta_object$ali (a matrix)

  input_class = class(msa)[1]

  if(input_class == 'fasta'){
    msa = msa$ali
  } else if(input_class == 'character'){
    # if character, assume it's a file path # 7/2/25 (stops may need to be message() for batch mode)
    if(!file.exists(msa)){
      stop('MSA file does not exist: ', msa)
    }

    msa = bio3d::read.fasta(msa, rm.dup = FALSE)$ali
  } else if(!input_class == 'matrix'){
    # print i dont know what you have #
    stop('msa provided is not a recognized format. Please provide a file path, an object from bio3d::read.fasta(), or a matrix [samples, positions]')
  }

  # force uppercase #
  msa = toupper(msa)

  # force sample (row) names #
  if(is.null(rownames(msa))){
    rownames(msa) = paste0('seq_', seq_len(nrow(msa)))
    message('Row names were not provided in msa. Assigning default names: seq_1, seq_2, ...')
  }

  ## NEED TO CHECK FOR UNIQUE ROWNAMES ##
  if(length(unique(rownames(msa))) != nrow(msa)){
    message('HOLD YOUR HORSES! MSA has non-unique row names. This may cause issues later in the pipeline. Please ensure row names are unique.')
    rownames(msa) = paste0('seq_', seq_len(nrow(msa)))
  }

  return(msa)
}


#' Detect Sequence Type
#'
#' Determines whether a given sequence is nucleotide or protein based on the proportion of standard nucleotide characters in the first \code{max_len} positions.
#' characters '-' and 'X' are dropped because they do not add information about either sequence type
#'
#' @param seq A character string representing a single biological sequence.
#' @param detect_sequence_threshold Proportion of characters that must be nucleotide-like (A, T, C, G) to classify the sequence as \code{"nucleotide"}. Default is 0.9.
#' @param detect_sequence_len Maximum number of non-gap '-' or 'X' characters (from the start of the sequence) to consider when computing the proportion. Default is 100.
#'
#' @return A string: either \code{"nucleotide"} or \code{"protein"}.
#' @keywords internal
.detect_sequence_type <- function(seq, detect_sequence_threshold = 0.8, detect_sequence_len = 100) {
  # check first 100 characters of sequence for nucleotides #
  # if >=90% are ATCG, return nucleotide # ~ similar to muscle approach ~ #

  # remove '-' and 'X' characters from the sequence
  seq = gsub("[-X]", "", seq)

  # get sequence substring
  seq = substr(seq, 1, min(nchar(seq), detect_sequence_len))

  # Count how many are A/T/C/G/N/U- # 7/2/25 (should i really count N and U and -) / maybe count not [VILMWYFSTNQPRHKDE]
  nuc_like = sum(strsplit(seq, "")[[1]] %in% c("A", "T", "C", "G"))
  prop = nuc_like / nchar(seq)

  # Check if the proportion of nucleotide-like characters is above the threshold
  seq_type = if(prop >= detect_sequence_threshold) "nucleotide" else "protein"

  # return the sequence type
  return(seq_type)
}

#' Get Reference Sequence from MSA
#'
#' Extracts a reference sequence from a multiple sequence alignment (MSA) using one of several strategies: a specified row index, the most complete sequence, or the consensus across all sequences.
#' Consensus is built from most frequent ATCG for nucleotide or 20 standardize amino acids for protein data. If none of these characters are found position becomes 'X'
#'
#' @param msa A character matrix of sequences (rows) by alignment positions (columns). Should be standardized using \code{.standardize_msa_input()}.
#' @param ref_method Either a character string (\code{"most_complete"} or \code{"consensus"}) or a numeric value indicating the row number to use as the reference.
#' @param force_seq_type Optional. Force sequence type to be either \code{"nucleotide"} or \code{"protein"}; if \code{NULL}, type is auto-detected.
#' @param ... Additional arguments passed to .detect_sequence_type()
#'
#' @return A list with two elements: \code{ref}, the reference sequence as a named character string, and \code{seq_type}, either \code{"nucleotide"} or \code{"protein"}.
#' @keywords internal
.get_reference_sequence = function(msa, ref_method = 'consensus', force_seq_type = NULL, detect_sequence_threshold = 0.8, detect_sequence_len = 100){
  # grab the reference sequence based on the method provided #

  # three passes #
  # 1. lenient (most_complete or consensus)
  # 2. detect sequence type
  # 3. strict (most_complete or consensus)
  # -- if method is numeric or force_seq_type is not NULL, skip the first two passes #

  # STEP 0 validating data ----
  # check that method is valid #
  if (!ref_method %in% c('most_complete', 'consensus') & !is.numeric(ref_method)){
    stop('ref_method must be one of "most_complete", "consensus", or a numeric value (row number)')
  }

  # check that seq_type is valid #
  if (!is.null(force_seq_type) && !force_seq_type %in% c('nucleotide', 'protein')) {
    stop('force_seq_type must be NULL, "nucleotide", or "protein"')
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

      counts = vapply(chars, function(ch) colSums(msa == ch, na.rm=TRUE), numeric(ncol(msa)))

      # grab the max count for each column #
      ref = chars[max.col(counts, ties.method="first")]
      ref = paste0(ref, collapse = '')
      names(ref) = 'ref.consensus'
    }

    # detect sequence type #
    seq_type = .detect_sequence_type(ref, detect_sequence_threshold = detect_sequence_threshold, detect_sequence_len = detect_sequence_len)
  }

  # step 3 strict pass ----

  # if method is numeric grab row -- detect_seq_type if needed #
  if(is.numeric(ref_method)){
    # check that row exists
    if (ref_method <= 0) stop("Reference method position must be a positive integer")
    if (ref_method > nrow(msa)) stop("Reference method position exceeds number of sequences in MSA")

    # return sequence by number
    ref = paste0(msa[ref_method,], collapse = '')
    names(ref) = paste0('ref.', rownames(msa)[ref_method])

    # if force_seq_type is NULL, detect seq_type #
    if(is.null(force_seq_type)){
      seq_type = .detect_sequence_type(ref)
    } else {
      seq_type = force_seq_type
    }

    return(list(ref = ref, seq_type = seq_type))
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
      ref = paste(consensus, collapse = '')
      names(ref) = 'ref.consensus'
    }

  }

  # retrun results ----
  return(list(ref = ref, seq_type = seq_type))
}

#' Translate DNA to Amino Acids
#'
#' Translates a DNA sequence to its corresponding amino acid sequence.
#' Ambiguous positions are represented as 'X'.
#'
#' @param seq A named character string representing a DNA sequence.
#' @param frame parameter for seqinr::translate, indicating the reading frame (0, 1, or 2).
#' @param sens parameter for seqinr::translate, indicating the sense of translation ('F' for forward, 'R' for reverse).
#' @param numcode parameter for seqinr::translate, indicating the ncbi genetic code to use (default is 1 for standard genetic code).
#'
#' @return A named character string representing the translated amino acid sequence.
#' @keywords internal
.translate_dna_to_protein = function(seq, frame = 0, sens = 'F', numcode = 1){

  # translate (NNN and --- are treated the same 'X')
  pep = seqinr::translate(strsplit(seq, '')[[1]],
                          frame = frame,
                          sens = sens,
                          numcode = numcode,
                          ambiguous = FALSE  # this is default but must be FALSE because gaps are also treated as X #
                          )

  # check for internal stops (can proceed)
  # consider stopping (but could be false stop from consensus building)
  if(any(pep[-length(pep)] == '*')){
    message('Internal stop codon(s) found in reference sequence\ntry different ref seq or try different frame')
  }

  # return as character vector
  pep = paste0(pep, collapse = '')
  names(pep) = names(seq)

  return(pep)
}

#' Extract Reference and Peptide Sequence from MSA
#'
#' Wrapper function that standardizes input, extracts a reference sequence from an MSA, detects or applies a sequence type, and translates the reference to peptide if needed.
#'
#' @param msa A character file path to a FASTA file, a matrix, or an object returned by \code{bio3d::read.fasta()}.
#' @param ref_method Reference extraction method: one of \code{"most_complete"}, \code{"consensus"}, or a numeric row index.
#' @param force_seq_type Optional sequence type: \code{"protein"}, \code{"nucleotide"}, or \code{NULL} to auto-detect.
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item \code{msa_mat}: The standardized alignment matrix.
#'   \item \code{ref}: The reference sequence (DNA or protein).
#'   \item \code{pep}: The translated peptide sequence (if nucleotide input).
#'   \item \code{seq_type}: The detected or specified sequence type.
#' }
#' @export
msa_to_ref = function(msa, ref_method = 'consensus', force_seq_type = NULL, verbose = 0,
                      detect_sequence_threshold = 0.8, detect_sequence_len = 100,
                      reading_frame = 0, reading_sens = 'F', genetic_code = 1){

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
    cat('\tmsa_to_ref: Translating reference sequence to peptide\n')
  }

  if(ref$seq_type == 'nucleotide'){
    pep = .translate_dna_to_protein(ref$ref,
                                    frame = reading_frame,
                                    sens = reading_sens,
                                    numcode = genetic_code
                                   )
  } else {
    # replace '-' with 'X' for protein
    pep = ref$ref
    pep = gsub('-', 'X', pep)
  }

  # return results ----
  return(list(msa_mat = msa,
              ref = ref$ref,
              pep = pep,
              seq_type = ref$seq_type)
         )
}

# TESTING -- remove before push #
if(F){
aavec = strsplit('AVIL-X', '')[[1]]
aasamp = sample(aavec, 100, replace = TRUE)
aamat = matrix(aasamp, ncol = 25, byrow = TRUE)

msa_aa = msa_to_ref(aamat, ref_method = 2, verbose = 1, reading_frame = 0, genetic_code = 1, frame = 0, detect_sequence_threshold = 0.8)

nucvec = strsplit('ATGCATGCX-RU', '')[[1]]
nucsamp = sample(nucvec, 100, replace = TRUE)
nucmat = matrix(nucsamp, ncol = 25, byrow = TRUE)

msa_nuc = msa_to_ref(nucmat, ref_method = 2, verbose = 1, reading_frame = 0, genetic_code = 1, frame = 0, detect_sequence_threshold = 0.8)
msa_nuc
}
#
