#' Cumulative distribution function for the String Count Distribution
#'
#' \code{pstringcountdist} returns the cumulative distribution function for the String Count Distribution
#'
#' The String Count Distribution is the distribution of the string-count for a specified string vector ```string``` in a random
#' vector of IID categorical variables from an alphabet ```alphabet``` with probability vector ```probs```.  The function allows the user
#' to specify the ```alphabet``` for the analysis, with the default alphabet being the natural numbers up to the length of the probability
#' vector.  (Note: The user can give either a numeric vector or a character vector for the ```string``` and ```alphabet```, but the elements
#' in the string must be in the alphabet, and both vectors must be the same type.)
#'
#' @usage \code{pstringcountdist()}
#' @param q The argument value(s) of the string-count quantile in the cumulative distribution function
#' @param size The size argument (either a scalar or a vector with the same length as ```q```)
#' @param string A numeric/character vector
#' @param probs A vector of the symbol probabilities (taken over the symbols in the ```alphabet```)
#' @param alphabet A numeric/character vector containing the alphabet for the analysis
#' @param allow.overlap Logical; if ```TRUE``` then string occurrances are counted even if they overlap with previously counted occurrances
#' @param lower.tail Logical; if ```TRUE``` the cumulative probabilities are ```P[X ??? q]```; if ```FALSE``` they are ```P[X > q]```
#' @param log.p Logical; if ```TRUE``` the function returns the cumulative log-probability; if ```FALSE``` the function returns the cumulative probability
#' @param full.array Optional input for the full array of joint log-probabilities for the string-count and state variable (from \code{full.array})
#' @return The probability or log-probability values from the cumulative distribution function

pstringcountdist <- function(q, size, string, probs, alphabet = NULL, allow.overlap = TRUE, lower.tail = TRUE, log.p = FALSE, full.array = NULL) {

  #Check inputs q, size, allow.overlap, lower.tail and log
  if (!is.vector(q))                                           { stop('Error: q must be a numeric vector') }
  if (!is.numeric(q))                                          { stop('Error: q must be a numeric vector') }
  RR <- length(q)
  if (RR == 0)                                                 { stop('Error: q must contain at least one value') }
  if (!is.vector(size))                                        { stop('Error: size must be an integer vector') }
  if (!is.numeric(size))                                       { stop('Error: size must be an integer vector') }
  if (as.integer(size) != size)                                { stop('Error: size must be an integer vector') }
  if (min(size) < 1)                                           { stop('Error: size must contain only positive integers') }
  if (length(size) == 1) { size <- rep(size, RR) }
  if (length(size) != RR)                                      { stop('Error: size should either be a scalar or a vector with the same length as x') }
  if (!is.vector(allow.overlap))                               { stop('Error: allow.overlap must be a single logical value') }
  if (!is.logical(allow.overlap))                              { stop('Error: allow.overlap must be a single logical value') }
  if (length(allow.overlap) != 1)                              { stop('Error: allow.overlap must be a single logical value') }
  if (!is.vector(lower.tail))                                  { stop('Error: lower.tail must be a single logical value') }
  if (!is.logical(lower.tail))                                 { stop('Error: lower.tail must be a single logical value') }
  if (length(lower.tail) != 1)                                 { stop('Error: lower.tail must be a single logical value') }
  if (!is.vector(log.p))                                       { stop('Error: log.p must be a single logical value') }
  if (!is.logical(log.p))                                      { stop('Error: log.p must be a single logical value') }
  if (length(log.p) != 1)                                      { stop('Error: log.p must be a single logical value') }

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

  #Check full array (if provided)
  if (!missing(full.array)) {
    if (max(size) > full.array$max.size)                       { stop('Error: The full.array you have provided is not large enough to cover your size') }
    if (!identical(string, full.array$string))                 { stop('Error: The full.array you have provided does not match your string') }
    if (!identical(probs, full.array$probs))                   { stop('Error: The full.array you have provided does not match your probs') }
    if (!identical(alphabet, full.array$alphabet))             { stop('Error: The full.array you have provided does not match your alphabet') }
    if (!identical(allow.overlap, full.array$allow.overlap))   { stop('Error: The full.array you have provided does not match your allow.overlap') } }

  #Get string matrix information
  MATS <- stringmatrix(string = string, probs = probs, alphabet = alphabet, allow.overlap = allow.overlap)
  MO   <- MATS$max.overlap

  #set initial values and log-probability matrix
  max.n <- max(size)
  max.r <- max(0, 1+floor((max.n-m)/(m-MO)))
  LOGPROBS <- matrix(-Inf, nrow = max.n+1, ncol = max.r+1)
  LOGPROBS[1,1] <- 0

  #Compute log-probabilities from string-count distribution (compute full array if not provided)
  if (!missing(full.array)) {

    #Use provided ARRAY to compute LOGPROBS
    ARRAY <- full.array$array
    for (nn in 1:max.n) {
      UR <- max(0, 1+floor((nn-m)/(m-MO)))
      for (rr in 0:(max.r)) {
        LOGPROBS[nn+1, rr+1] <- matrixStats::logSumExp(ARRAY[nn+1, rr+1, ]) } }

  } else {

    #Compute ARRAY and LOGPROBS together
    LOGH  <- log(MATS$transition.probs)
    max.h <- nrow(LOGH)-1
    hh.start <- MATS$state.start
    hh.count <- MATS$state.count
    ARRAY <- array(-Inf, dim = c(max.n+1, max.r+1, max.h+1),
                   dimnames = list(sprintf('n[%s]', 0:max.n), sprintf('r[%s]', 0:max.r), sprintf('h[%s]', 0:m)))
    ARRAY[1, 1, hh.start+1] <- 0
    for (nn in 1:max.n) {
      UR <- max(0, 1+floor((nn-m)/(m-MO)))
      for (rr in 0:UR)    {
        for (hh in 0:max.h) {
          if (hh %in% hh.count) { if (rr > 0) { ARRAY[nn+1, rr+1, hh+1] <- matrixStats::logSumExp(LOGH[, hh+1] + ARRAY[nn, rr, ])   } } else {
            ARRAY[nn+1, rr+1, hh+1] <- matrixStats::logSumExp(LOGH[, hh+1] + ARRAY[nn, rr+1, ]) } }
        LOGPROBS[nn+1, rr+1] <- matrixStats::logSumExp(ARRAY[nn+1, rr+1, ]) } } }

  #Extract cumulative log-probabilities for argument values
  OUT <- rep(-Inf, RR)
  for (i in 1:RR) { if (q[i] >= 0) { OUT[i] <- matrixStats::logSumExp(LOGPROBS[size[i]+1, 1:(min(floor(q[i]), max.r)+1)]) } }

  #Adjust log-probabilities if computing upper tail
  if (!lower.tail) { OUT <- VGAM::log1mexp(-OUT) }

  #Return output
  if (log.p) { OUT } else { exp(OUT) } }
