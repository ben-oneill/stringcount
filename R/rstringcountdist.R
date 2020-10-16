#' Random generation from the String Count Distribution
#'
#' \code{rstringcountdist} returns the random values from the String Count Distribution
#'
#' The String Count Distribution is the distribution of the string-count for a specified string vector ```string```` in a random
#' vector of IID categorical variables from an alphabet ```alphabet``` with probability vector ```probs```.  This function generates
#' random outputs from this distribution using the underlying pseudo-random number generator (PRNG) in ```R```.  The function allows the
#' user to specify the ```alphabet```` for the analysis, with the default alphabet being the natural numbers up to the length of the
#' probability vector.  (Note: The user can give either a numeric vector or a character vector for the ```string```` and ```alphabet````,
#' but the elements in the string must be in the alphabet, and both vectors must be the same type.)  The output of the function is a
#' vector composed of ```n``` values generated from the distribution.
#'
#' @usage \code{rstringcountdist()}
#' @param n The numer of random variables to generate for the output
#' @param size The size argument (either a scalar or a vector with the same length as x)
#' @param string A numeric/character vector
#' @param probs A vector of the symbol probabilities (taken over the symbols in the \code{alphabet})
#' @param alphabet A numeric/character vector containing the alphabet for the analysis
#' @return The probability or log-probability values from the density function

rstringcountdist <- function(n, size, string, probs, alphabet = NULL) {

  #Check inputs r, size and log
  if (!is.vector(n))                                           { stop('Error: n must be a positive integer') }
  if (!is.numeric(n))                                          { stop('Error: n must be a positive integer') }
  if (length(n) != 1)                                          { stop('Error: n must be a single positive integer') }
  if (as.integer(n) != n)                                      { stop('Error: n must be a positive integer') }
  if (min(n) < 1)                                              { stop('Error: n must be a positive integer') }
  if (!is.vector(size))                                        { stop('Error: size must be a positive integer') }
  if (!is.numeric(size))                                       { stop('Error: size must be a positive integer') }
  if (length(size) != 1)                                       { stop('Error: size must be a single positive integer') }
  if (as.integer(size) != size)                                { stop('Error: size must be a positive integer') }
  if (min(size) < 1)                                           { stop('Error: size must be a positive integer') }

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

  #Match string characters to alphabet
  MATCH <- match(string, alphabet)

  #Generate random string-counts
  COUNT <- rep(NA, n)
  for (i in 1:n) {
    TEXT     <- sample.int(K, size = size, replace = TRUE, prob = probs)
    COUNT[i] <- stringcount(text = TEXT, string = string) }

  #Output the counts
  COUNT }
