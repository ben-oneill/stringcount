#' Joint Distribution of String-Count and State Variable
#'
#' \code{full.array} returns the array of joint log-probabilities for array of string-count values and state variables.  This is the
#' full object that is computed in the stringcountdist functions (density, cumulative distribution, quantile function) in order to obtain
#' the String Count Distribution.  Those other functions allow the user to input the pre-computed array object, which can be computed and
#' stored using this function.
#'
#' @usage \code{full.array()}
#' @param max.size The maximum size argument (a positive integer)
#' @param string A numeric/character vector
#' @param probs A vector of the symbol probabilities (taken over the symbols in the ```alphabet```)
#' @param alphabet A numeric/character vector containing the alphabet for the analysis
#' @param allow.overlap Logical; if ```TRUE``` then string occurrances are counted even if they overlap with previously counted occurrances
#' @return The full array of log-probability values from the joint distribution of the string count and state variable

full.array <- function(max.size, string, probs, alphabet = NULL, allow.overlap = TRUE) {

  #Check inputs max.size and allow.overlap
  if (!is.vector(max.size))                                    { stop('Error: max.size must be a positive integer') }
  if (!is.numeric(max.size))                                   { stop('Error: max.size must be a positive integer') }
  if (length(max.size) != 1)                                   { stop('Error: max.size must be a single positive integer') }
  if (as.integer(max.size) != max.size)                        { stop('Error: max.size must be a positive integer') }
  if (min(max.size) < 1)                                       { stop('Error: max.size must be a positive integer') }
  if (!is.vector(allow.overlap))                               { stop('Error: allow.overlap must be a single logical value') }
  if (!is.logical(allow.overlap))                              { stop('Error: allow.overlap must be a single logical value') }
  if (length(allow.overlap) != 1)                              { stop('Error: allow.overlap must be a single logical value') }

  #Check input probs
  if (!is.vector(probs))                                       { stop('Error: probs must be a probability vector') }
  if (!is.numeric(probs))                                      { stop('Error: probs must be a probability vector') }
  if (length(probs) == 0)                                      { stop('Error: probs must have at least one value') }
  if (min(probs) < 0)                                          { stop('Error: probs must be a probability vector') }
  if (sum(probs) != 1)                                         { stop('Error: probs must be a probability vector') }

  #Check input alphabet
  if (missing(alphabet)) {
    alphabet <- 1:length(probs) }
  if (!is.vector(alphabet))                                    { stop('Error: alphabet must be a vector') }
  K <- length(unique(alphabet))
  if (length(alphabet) != K) {
    warning('Input alphabet contained duplicate elements')
    alphabet <- unique(alphabet) }
  if (length(probs) != K)                                      { stop('Error: probs must have the same length as alphabet') }

  #Set alphabet type
  TYPE <- NULL
  if (is.numeric(alphabet))   { TYPE <- 'numeric' } else {
    if (is.character(alphabet)) { TYPE <- 'character' } }
  if (is.null(TYPE))                                           { stop('Error: alphabet should be a numeric or character vector') }

  #Check input string
  if (!is.vector(string))                                      { stop('Error: string must be a vector') }
  if ((TYPE ==   'numeric')&&(!is.numeric(string)))            { stop('Error: string must be the same type as alphabet') }
  if ((TYPE == 'character')&&(!is.character(string)))          { stop('Error: string must be the same type as alphabet') }
  m <- length(string)
  for (i in 1:m) {
    if (!(string[i] %in% alphabet))                            { stop(paste0('Error: string element ', i, ' is not in the alphabet')) } }

  #Get string matrix information
  MATS <- stringmatrix(string = string, probs = probs, alphabet = alphabet, allow.overlap = allow.overlap)
  MO   <- MATS$max.overlap
  LOGH <- log(MATS$transition)

  #set values
  max.n <- max.size
  max.r <- max(0, 1+floor((max.n-m)/(m-MO)))

  #Generate log-probabilities
  ARRAY <- array(-Inf, dim = c(max.n+1, max.r+1, m+1),
                 dimnames = list(sprintf('n[%s]', 0:max.n), sprintf('r[%s]', 0:max.r), sprintf('h[%s]', 0:m)))
  ARRAY[1,1,1] <- 0
  for (nn in 1:max.n) {
    UR <- max(0, 1+floor((nn-m)/(m-MO)))
    for (rr in 0:UR) {
      for (hh in 0:(m-1)) {
        ARRAY[nn+1, rr+1, hh+1] <- matrixStats::logSumExp(LOGH[, hh+1] + ARRAY[nn, rr+1, ]) }
      if (rr > 0) { ARRAY[nn+1, rr+1, m+1] <- matrixStats::logSumExp(LOGH[, m+1] + ARRAY[nn, rr, ]) } } }

  #Return output
  list(array = ARRAY, max.size = max.size, string = string, probs = probs, alphabet = alphabet, allow.overlap = allow.overlap) }
