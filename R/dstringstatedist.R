#' Density function for the String State Distribution
#'
#' \code{dstringstatedist} returns the density function for the String State Distribution
#'
#' The String State Distribution is the marginal/conditional distribution of the string-state for a specified string vector ```string``` in
#' a random vector of IID categorical variables from an alphabet ```alphabet``` with probability vector ```probs```.  The function allows the
#' user to specify the ```alphabet``` for the analysis, with the default alphabet being the natural numbers up to the length of the probability
#' vector.  (Note: The user can give either a numeric vector or a character vector for the ```string``` and ```alphabet```, but the elements
#' in the string must be in the alphabet, and both vectors must be the same type.)
#'
#' @usage \code{dstringstatedist()}
#' @param x The argument value(s) of the string-state in the density function
#' @param r A vector of conditioning values for the string count (optional)
#' @param size The size argument (either a scalar or a vector with the same length as ```x```)
#' @param string A numeric/character vector
#' @param probs A vector of the symbol probabilities (taken over the symbols in the ```alphabet```)
#' @param alphabet A numeric/character vector containing the alphabet for the analysis
#' @param allow.overlap Logical; if ```TRUE``` then string occurrances are counted even if they overlap with previously counted occurrances
#' @param log Logical; if ```TRUE``` the function returns the log-probability; if ```FALSE``` the function returns the probability
#' @param full.array Optional input for the full array of joint log-probabilities for the string-count and state variable (from \code{full.array})
#' @return The probability or log-probability values from the density function

dstringstatedist <- function(x, r = NULL, size, string, probs, alphabet = NULL, allow.overlap = TRUE, log = FALSE, full.array = NULL) {

  #Check inputs x, size, allow.overlap and log
  if (!is.vector(x))                                           { stop('Error: x must be a numeric vector') }
  if (!is.numeric(x))                                          { stop('Error: x must be a numeric vector') }
  HH <- length(x)
  if (HH == 0)                                                 { stop('Error: x must contain at least one value') }
  if (!is.vector(size))                                        { stop('Error: size must be an integer vector') }
  if (!is.numeric(size))                                       { stop('Error: size must be an integer vector') }
  if (as.integer(size) != size)                                { stop('Error: size must be an integer vector') }
  if (min(size) < 1)                                           { stop('Error: size must contain only positive integers') }
  if (length(size) == 1) { size <- rep(size, HH) }
  if (length(size) != HH)                                      { stop('Error: size should either be a scalar or a vector with the same length as x') }
  if (!is.vector(allow.overlap))                               { stop('Error: allow.overlap must be a single logical value') }
  if (!is.logical(allow.overlap))                              { stop('Error: allow.overlap must be a single logical value') }
  if (length(allow.overlap) != 1)                              { stop('Error: allow.overlap must be a single logical value') }
  if (!is.vector(log))                                         { stop('Error: log must be a single logical value') }
  if (!is.logical(log))                                        { stop('Error: log must be a single logical value') }
  if (length(log) != 1)                                        { stop('Error: log must be a single logical value') }

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

  #Check input r
  if (!missing(r)) {
    if (!is.vector(r))                                         { stop('Error: r must be an integer vector (if supplied)') }
    if (!is.numeric(r))                                        { stop('Error: r must be an integer vector (if supplied)') }
    if (as.integer(r) != r)                                    { stop('Error: r must be an integer vector (if supplied)') }
    if (min(r) < 0)                                            { stop('Error: r must contain only positive integers (if supplied)') }
    if (length(r) == 0)                                        { stop('Error: r must contain at least one value (if supplied)') }
    if (length(r) == 1) { r <- rep(r, HH) }
    if (length(r) != HH)                                       { stop('Error: r should either be a scalar or a vector with the same length as x') } }

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

  #More checks on input r
  if (!missing(r)) {
    for (i in 1:HH) {
      UR <- max(0, 1+floor((size[i]-m)/(m-MO)))
      if (r[i] > UR)                                           { stop(paste0('Error: Element ', i, ' of r is too large to be a feasible string-count value')) } } }
  if (missing(r)) {
    r <- as.numeric(rep(NA, HH)) }

  #set initial values and log-probability matrix
  max.n <- max(size)
  max.r <- max(0, 1+floor((max.n-m)/(m-MO)))
  max.h <- nrow(MATS$transition.table)-1

  #Compute log-probabilities from string-count distribution (compute full array if not provided)
  if (!missing(full.array)) {

    #Use provided ARRAY to compute LOGPROBS
    ARRAY <- full.array$array

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
                                              ARRAY[nn+1, rr+1, hh+1] <- matrixStats::logSumExp(LOGH[, hh+1] + ARRAY[nn, rr+1, ]) } } } } }

  #Extract log-probabilities for argument values
  OUT <- rep(-Inf, HH)
  for (i in 1:HH) {
  if (x[i] %in% (0:max.h)) {
    if (is.na(r[i])) {
      OUT[i] <- matrixStats::logSumExp(ARRAY[size[i]+1, , x[i]+1])
    } else {
      OUT[i] <- ARRAY[size[i]+1, r[i]+1, x[i]+1] - matrixStats::logSumExp(ARRAY[size[i]+1, r[i]+1, ]) } } }

  #Return output
  if (log) { OUT } else { exp(OUT) } }
