#' Quantile function for the String Count Distribution
#'
#' \code{qstringcountdist} returns the quantile function for the String Count Distribution
#'
#' The String Count Distribution is the distribution of the string-count for a specified string vector ```string``` in a random
#' vector of IID categorical variables from an alphabet ```alphabet``` with probability vector ```probs```.  The function allows the user
#' to specify the ```alphabet``` for the analysis, with the default alphabet being the natural numbers up to the length of the probability
#' vector.  (Note: The user can give either a numeric vector or a character vector for the ```string``` and ```alphabet```, but the elements
#' in the string must be in the alphabet, and both vectors must be the same type.)
#'
#' @usage \code{qstringcountdist()}
#' @param p The argument cumulative probabilities/log-probabilities in the quantile function
#' @param size The size argument (either a scalar or a vector with the same length as ```p```)
#' @param string A numeric/character vector
#' @param probs A vector of the symbol probabilities (taken over the symbols in the ```alphabet```)
#' @param alphabet A numeric/character vector containing the alphabet for the analysis
#' @param lower.tail Logical; if ```TRUE``` the cumulative probabilities are ```P[X ??? q]```; if ```FALSE``` they are ```P[X > q]```
#' @param log.p Logical; if ```TRUE``` the values in ```p``` are log-probabilities; if ```FALSE``` they are probabilities
#' @return The quantiles from the quantile function

qstringcountdist <- function(p, size, string, probs, alphabet = NULL, lower.tail = TRUE, log.p = FALSE) {

  #Check inputs p, size, lower.tail and log.p
  if (!is.vector(p))                                           { stop('Error: p must be a numeric vector') }
  if (!is.numeric(p))                                          { stop('Error: p must be a numeric vector') }
  RR <- length(p)
  if (RR == 0)                                                 { stop('Error: p must contain at least one value') }
  if (!is.vector(size))                                        { stop('Error: size must be an integer vector') }
  if (!is.numeric(size))                                       { stop('Error: size must be an integer vector') }
  if (as.integer(size) != size)                                { stop('Error: size must be an integer vector') }
  if (min(size) < 1)                                           { stop('Error: size must contain only positive integers') }
  if (length(size) == 1) { size <- rep(size, RR) }
  if (length(size) != RR)                                      { stop('Error: size should either be a scalar or a vector with the same length as x') }
  if (!is.vector(lower.tail))                                  { stop('Error: lower.tail must be a single logical value') }
  if (!is.logical(lower.tail))                                 { stop('Error: lower.tail must be a single logical value') }
  if (length(lower.tail) != 1)                                 { stop('Error: lower.tail must be a single logical value') }
  if (!is.vector(log.p))                                       { stop('Error: log.p must be a single logical value') }
  if (!is.logical(log.p))                                      { stop('Error: log.p must be a single logical value') }
  if (length(log.p) != 1)                                      { stop('Error: log.p must be a single logical value') }
  if (!log.p) {
    if (min(p) < 0)                                            { stop('Error: Elements of p cannot be less than zero') }
    if (max(p) > 1)                                            { stop('Error: Elements of p cannot be greater than one') } }

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
  LOGP  <- if (log.p) { p } else { log(p) }
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

  #Compute quantiles
  QUANTILE <- rep(NA, RR)
  if (lower.tail) {
  for (i in 1:RR) {
    Q <- 0
    LOGPROB <- LLL[size[i]+1, 1]
    while (LOGP[i] > LOGPROB) {
      Q <- Q+1
      LOGPROB <- matrixStats::logSumExp(c(LOGPROB, LLL[size[i]+1, Q+1])) }
    QUANTILE[i] <- Q } }
  if (!lower.tail) {
    for (i in 1:RR) {
      Q <- max.r
      LOGPROB <- 0
      while (LOGP[i] <= LOGPROB) {
        LOGPROB <- matrixStats::logSumExp(c(LOGPROB, LLL[size[i]+1, Q+1]))
        Q <- Q-1 }
      QUANTILE[i] <- Q } }

  #Return output
  QUANTILE }
