# --------------------------------------------------------------- #
# MSA MODULE -- utilities for taking multiple sequence alignment (MSA) data and returning standardized msa_info #
# internal .functions() are wrapped with external msa_to_ref() module
#
# NOTES:
# 1. gaps in peptide reference are represented as 'X' -- this helps distinguish from new introduced gaps later in the pipeline #
# 2. Consensus and most complete reference methods only count valid DNA and AA characters
#   -- meaning if a DNA column has 25 'R' and 1 'A' the consensus method returns 'A' for that column
# 3. No attempt is made to salvage ambiguous characters 'GG[RBY]' -> 'Glycine'
#
# email: bbroyle@purdue.edu
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
#' msa_mat <- standardize_msa_input("example_alignment.fasta")
#'
#' # From bio3d fasta object
#' fasta_obj <- bio3d::read.fasta("example_alignment.fasta")
#' msa_mat <- standardize_msa_input(fasta_obj)
#'
#' # From pre-loaded matrix (rows ~ samples, aligned positions ~ columns)
#' fasta_obj <- bio3d::read.fasta("example_alignment.fasta")
#' my_matrix <- fasta_obj$ali
#' msa_mat <- standardize_msa_input(my_matrix)
#' }
#'
#' @export

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

# .get_reference_sequence() ----

#' Get Reference Sequence from MSA
#'
#' extracts a reference sequence from a multiple sequence alignment (MSA) using one of three strategies:
#' a specified row index, the most complete sequence, or the consensus across all sequences.
#' consensus is built from the most frequent valid characters (ATCG for nucleotide; 20 standard amino acids for protein).
#' if no valid character is found in a column, the position becomes 'X'.
#'
#' @param msa a character matrix of sequences (rows) by alignment positions (columns).
#'   should be standardized using \code{.standardize_msa_input()}
#' @param ref_method either \code{"most_complete"}, \code{"consensus"}, or a numeric row index
#' @param force_seq_type optional. force sequence type to \code{"nucleotide"} or \code{"protein"}.
#'   if \code{NULL}, type is auto-detected
#' @param detect_sequence_threshold proportion of ATCG characters required to classify as nucleotide when auto-detecting. default is 0.8
#' @param detect_sequence_len maximum number of non-gap '-' or 'X' positions (from the start of the sequence) to consider for auto-detection. default is 100
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

# .translate_dna_to_protein() ----

#' Translate dna to amino acids
#'
#' translates a dna sequence to its amino acid sequence using \code{seqinr::translate}.
#' ambiguous triplets (including gaps) are translated as 'X'.
#' internal stop codons are reported with a message but translation continues.
#'
#' @param seq a named character string dna sequence
#' @param frame reading frame (0, 1, or 2). always fixed to 0 in full pipeline but retained here for flexibility
#' @param sens translation sense ('F' forward, 'R' reverse). always fixed to 'F' in full pipeline
#' @param numcode ncbi genetic code number to use (default 1, standard code)
#'
#' @return a named character string of the translated amino acid sequence
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

# msa_to_ref() ----

#' extract reference and peptide sequence from msa
#'
#' wrapper to standardize msa input, extract a reference sequence,
#' detect or apply sequence type, and translate the reference to peptide if needed.
#'
#' can be called directly as part of a modular workflow, or internally by
#' \code{run_evo3d()}, where frame and sens are fixed and only \code{genetic_code}
#' may be user-specified.
#'
#' @param msa path to fasta file, a matrix, or an object returned by \code{bio3d::read.fasta()}
#' @param ref_method method to choose reference: one of \code{"most_complete"},
#'   \code{"consensus"}, or numeric row index (default "consensus")
#' @param force_seq_type optional sequence type: \code{"protein"}, \code{"nucleotide"},
#'   or \code{NULL} to auto-detect (default)
#' @param verbose integer, print progress messages if > 0
#' @param detect_sequence_threshold proportion of atcg required to call nucleotide
#'   when auto-detecting (default 0.8)
#' @param detect_sequence_len number of leading characters used for detection (default 100)
#' @param reading_frame reading frame for translation (0,1,2). fixed to 0 in \code{run_evo3d}
#' @param reading_sens translation sense ('F' forward, 'R' reverse). fixed to 'F' in \code{run_evo3d}
#' @param genetic_code ncbi genetic code number (default 1, standard code)
#'
#' @return a list with
#' \itemize{
#'   \item \code{msa_mat}: standardized alignment matrix
#'   \item \code{ref}: reference sequence (dna or protein)
#'   \item \code{pep}: peptide sequence (translated if nucleotide input,
#'     or reference with '-' replaced by 'X' if protein)
#'   \item \code{seq_type}: detected or forced sequence type
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


