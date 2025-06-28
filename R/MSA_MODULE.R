# --------------------------------------------------------------- #
# MSA MODULE goal is to return msa_info list
# 1. msa as matrix
# 2. reference sequence
# 3. peptide reference sequence
# 4. msa sequence type
#
# ** input to msa can be matrix, fasta file, or bio3d::read.fasta() object
# ** sequence type is autodetected based on the first 100 characters
# ** reference sequence can be specified by row number, most complete, or consensus
# ** sequence type can be forced to nucleotide or protein
#
# email me: bbroyle@purdue.edu
# --------------------------------------------------------------- #

#' Standardize MSA Input
#'
#' Accepts various forms of MSA input (file path, matrix, or bio3d fasta object) and standardizes it into a character matrix.
#'
#' @param msa A character file path to a FASTA file, a matrix, or an object returned by \code{bio3d::read.fasta()}.
#'
#' @return A character matrix with sequences as rows and alignment positions as columns.
#' @keywords internal
.standardize_msa_input = function(msa){
  # Take a variety of msa input types and return standardized matrix #
  # expecting file path (character) | matrix | fasta classes #

  input_class = class(msa)[1]

  if(input_class == 'fasta'){
    msa = msa$ali
  } else if(input_class == 'character'){
    # !! add a check if file exists !! #
    # file.exists(msa)
    msa = bio3d::read.fasta(msa)$ali
  } else if(!input_class == 'matrix'){
    # print i dont know what you have #
    stop('NOT ONE OF THREE MSA OPTIONS')
  }

  # force uppercase #
  msa = toupper(msa)

  return(msa)
}


#' Detect Sequence Type
#'
#' Heuristically determines whether a given sequence is nucleotide or protein based on the proportion of standard nucleotide characters in the first \code{max_len} positions.
#'
#' @param seq A character string representing a single biological sequence.
#' @param threshold Proportion of characters that must be nucleotide-like (A, T, C, G, N, U, -) to classify the sequence as \code{"nucleotide"}. Default is 0.9.
#' @param max_len Maximum number of characters (from the start of the sequence) to consider when computing the proportion. Default is 100.
#'
#' @return A string: either \code{"nucleotide"} or \code{"protein"}.
#' @keywords internal
.detect_sequence_type <- function(seq, threshold = 0.9, max_len = 100) {
  # check first 100 characters of sequence for nucleotides #
  # if >90% are A/T/C/G/N/U/-, return nucleotide #

  # get sequence substring
  seq = substr(seq, 1, min(nchar(seq), max_len))

  # Count how many are A/T/C/G/N/U-
  nuc_like = sum(strsplit(seq, "")[[1]] %in% c("A", "T", "C", "G", "N", "U", "-"))
  prop = nuc_like / nchar(seq)

  # Check if the proportion of nucleotide-like characters is above the threshold
  seqtype = if(prop >= threshold) "nucleotide" else "protein"

  # return the sequence type
  return(seqtype)
}


#' Get Reference Sequence from MSA
#'
#' Extracts a reference sequence from a multiple sequence alignment (MSA) using one of several strategies: a specified row index, the most complete sequence, or the consensus across all sequences.
#'
#' @param msa A character matrix of sequences (rows) by alignment positions (columns). Should be standardized using \code{.standardize_msa_input()}.
#' @param method Either a character string (\code{"most_complete"} or \code{"consensus"}) or a numeric value indicating the row number to use as the reference.
#' @param force_seqtype Optional. Force sequence type to be either \code{"nucleotide"} or \code{"protein"}; if \code{NULL}, type is auto-detected.
#'
#' @return A list with two elements: \code{ref}, the reference sequence as a named character string, and \code{msa_type}, either \code{"nucleotide"} or \code{"protein"}.
#' @keywords internal
.get_reference_sequence = function(msa, method = 1, force_seqtype = NULL){
  # grab the reference sequence based on the method provided #

  # assuming .standardize_msa_input() has been called #

  # check that method is valid #
  if (!method %in% c('most_complete', 'consensus') & !is.numeric(method)){
    stop('method must be one of "most_complete", "consensus", or a numeric value (row number)')
  }

  # detect msa type -- depends on the set of characters to use #
  if (is.null(force_seqtype)){
    seq_type = .detect_sequence_type(paste0(msa[1,], collapse = ''))
  } else {

    if(!force_seqtype %in% c('nucleotide', 'protein')){
      stop('force_seqtype must be either "nucleotide" or "protein"')
    }

    seq_type = force_seqtype
  }

  if (is.numeric(method)) {

    # check that row exists
    if (method > nrow(msa)) stop("Reference method position exceeds number of sequences in MSA")

    # return sequence by number
    ref = paste0(msa[method,], collapse = '')
    names(ref) = paste0('ref.', rownames(msa)[method])
    return(list(ref = ref, seq_type = seq_type))

  } else if (method == "most_complete") {

    # format characters that count as complete #
    if(seq_type == 'nucleotide'){
      chars = c('A', 'T', 'C', 'G')
    } else {
      chars = strsplit('AVILMWYFSTNQCGPRHKDE', '')[[1]]
    }

    # count resolved characters in each sequence #
    complete = apply(msa, 1, function(x) sum(x %in% chars))

    # return sequence with least_gaps
    most_complete = which.max(complete)
    ref = paste0(msa[most_complete, ], collapse = '')
    names(ref) = paste0('ref.',rownames(msa)[most_complete])
    return(list(ref = ref, seq_type = seq_type))

  } else if (method == "consensus") {

    # format characters that count as complete #
    if(seq_type == 'nucleotide'){
      chars = c('A', 'T', 'C', 'G', 'U', 'N') # add U and N so they are not replaced with gap
    } else {
      chars = strsplit('AVILMWYFSTNQCGPRHKDE', '')[[1]]
    }

    # count the number of matches for each char in column #
    counts <- vapply(chars,
                     function(ch) colSums(msa == ch, na.rm=TRUE),
                     numeric(ncol(msa)))

    # grab the max count for each column #
    consensus <- chars[max.col(counts, ties.method="first")]

    # if row sum for counts = 0. no valid characters found, put gap
    gap = which(rowSums(counts) == 0)
    if(length(gap) > 0){
      consensus[gap] = '-'
    }

    # paste consensus sequence
    ref = paste(consensus, collapse = '')
    names(ref) = 'ref.consensus'
    return(list(ref = ref, seq_type = seq_type))
  }
}



#' Translate DNA to Amino Acids
#'
#' Translates a DNA sequence to its corresponding amino acid sequence using the standard genetic code.
#' Assumes the input sequence is in-frame and codons are contiguous. Gaps or ambiguous bases may result in 'X' or stop codons.
#'
#' @param seq A named character string representing a DNA sequence.
#'
#' @return A named character string representing the translated amino acid sequence.
#' @keywords internal
.translate_dna_to_protein = function(seq){
  # translate a DNA sequence into an amino acid sequence
  # assumes that the sequence is in the correct frame

  # remove gaps if present (might not be good) ('A-T' goes to X)
  #ref = gsub('-', '', ref)

  # translate (NNN and --- are treated the same 'X')
  pep = seqinr::translate(strsplit(seq, '')[[1]])

  # check for internal stops (can proceed)
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
#' @param force_seqtype Optional sequence type: \code{"protein"}, \code{"nucleotide"}, or \code{NULL} to auto-detect.
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item \code{msa_mat}: The standardized alignment matrix.
#'   \item \code{ref}: The reference sequence (DNA or protein).
#'   \item \code{pep}: The translated peptide sequence (if nucleotide input).
#'   \item \code{seq_type}: The detected or specified sequence type.
#' }
#' @export
msa_to_ref = function(msa, ref_method = 1, force_seqtype = NULL, genetic_code = 1, verbose = 0){

  # standardize msa input #
  if(verbose > 0) {
    cat('\tmsa_to_ref: Standardizing MSA input\n')
  }

  msa = .standardize_msa_input(msa)

  # grab reference sequence #
  if(verbose > 0) {
    cat('\tmsa_to_ref: Extracting reference sequence\n')
  }
  ref = .get_reference_sequence(msa, ref_method, force_seqtype = force_seqtype)

  # grab protein sequence if needed #
  if(verbose > 0) {
    cat('\tmsa_to_ref: Translating reference sequence to peptide\n')
  }

  if(ref$seq_type == 'nucleotide'){
    pep = .translate_dna_to_protein(ref$ref)
  } else {
    # replace '-' with 'X' for protein
    pep = ref$ref
    pep = gsub('-', 'X', pep)
  }

  return(list(msa_mat = msa,
              ref = ref$ref,
              pep = pep,
              seq_type = ref$seq_type)
         )
}
