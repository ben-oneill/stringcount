#' Density function for the String Count Distribution
#'
#' \code{dstringcountdist} returns the density function for the String Count Distribution
#'
#' The String Count Distribution is the distribution of the string-count for a specified string vector ```string``` in a random
#' vector of IID categorical variables from an alphabet ```alphabet``` with probability vector ```probs```.  The function allows the user
#' to specify the ```alphabet``` for the analysis, with the default alphabet being the natural numbers up to the length of the probability
#' vector.  (Note: The user can give either a numeric vector or a character vector for the ```string``` and ```alphabet```, but the elements
#' in the string must be in the alphabet, and both vectors must be the same type.)
#'
#' @usage \code{dstringcountdist()}
#' @param r The argument value(s) of the string-count in the density function
#' @param size The size argument (either a scalar or a vector with the same length as x)
#' @param string A numeric/character vector
#' @param probs A vector of the symbol probabilities (taken over the symbols in the \code{alphabet})
#' @param alphabet A numeric/character vector containing the alphabet for the analysis
#' @param log Logical; if ```TRUE``` the function returns the log-probability; if ```FALSE``` the function returns the probability
#' @param give.state Logical; if ```TRUE``` the output has an attribute that gives the joint log-probabilities for the string-count and state variable
#' @return The probability or log-probability values from the density function

dstringcountdist <- function(r, size, string, probs, alphabet = NULL, log = FALSE, give.state = FALSE) {

  #Check inputs r, size and log
  if (!is.vector(r))                                           { stop('Error: r must be a numeric vector') }
  if (!is.numeric(r))                                          { stop('Error: r must be a numeric vector') }
  RR <- length(r)
  if (RR == 0)                                                 { stop('Error: r must contain at least one value') }
  if (!is.vector(size))                                        { stop('Error: size must be an integer vector') }
  if (!is.numeric(size))                                       { stop('Error: size must be an integer vector') }
  if (as.integer(size) != size)                                { stop('Error: size must be an integer vector') }
  if (min(size) < 1)                                           { stop('Error: size must contain only positive integers') }
  if (length(size) == 1) { size <- rep(size, RR) }
  if (length(size) != RR)                                      { stop('Error: size should either be a scalar or a vector with the same length as x') }
  if (!is.vector(log))                                         { stop('Error: log must be a single logical value') }
  if (!is.logical(log))                                        { stop('Error: log must be a single logical value') }
  if (length(log) != 1)                                        { stop('Error: log must be a single logical value') }
  if (!is.vector(give.state))                                  { stop('Error: give.state must be a single logical value') }
  if (!is.logical(give.state))                                 { stop('Error: give.state must be a single logical value') }
  if (length(give.state) != 1)                                 { stop('Error: give.state must be a single logical value') }

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

  #set values
  max.n <- max(size)
  max.r <- max.n-m+1

  #Get transition log-probability matrix
  H    <- stringmatrix(string = string, probs = probs, alphabet = alphabet)$transition
  LOGH <- log(H)

  #Generate log-probabilities
  ARRAY <- array(-Inf, dim = c(max.n+1, max.r+1, m+1),
                 dimnames = list(sprintf('n[%s]', 0:max.n), sprintf('r[%s]', 0:max.r), sprintf('h[%s]', 0:m)))
  LLL   <- matrix(-Inf, nrow = max.n+1, ncol = max.r+1)
  ARRAY[1,1,1] <- 0
  LLL[1,1]     <- 0
  for (nn in 1:max.n) {
  for (rr in 0:(max.r)) {
  for (hh in 0:(m-1)) {
    ARRAY[nn+1, rr+1, hh+1] <- matrixStats::logSumExp(LOGH[, hh+1] + ARRAY[nn, rr+1, ]) }
    if (rr > 0) { ARRAY[nn+1, rr+1, m+1] <- matrixStats::logSumExp(LOGH[, m+1] + ARRAY[nn, rr, ]) }
    LLL[nn+1, rr+1] <- matrixStats::logSumExp(ARRAY[nn+1, rr+1, ]) } }

  #Extract log-probabilities for argument values
  LOGDENS <- rep(-Inf, RR)
  for (i in 1:RR) { if (r[i] %in% (0:max.r)) { LOGDENS[i] <- LLL[size[i]+1, r[i]+1] } }

  #Give state probabilities if required
  if (give.state) {
    attr(LOGDENS, 'state.logprobs') <- ARRAY[max.n+1, , ] }

  #Return output
  if (log) { LOGDENS } else { exp(LOGDENS) } }
