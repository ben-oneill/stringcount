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
#' @param start.state.probs A vector of probabilities for the starting state (defaults to point mass on the initial state)
#' @param array Optional input for an array of joint log-probabilities for the string-count and state variable
#' @return The full array of log-probability values from the joint distribution of the string count and state variable

full.array <- function(max.size, string, probs, alphabet = NULL, allow.overlap = TRUE, start.state.probs = NULL, array = NULL) {

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
  if (is.numeric(alphabet))     { TYPE <- 'numeric' } else {
    if (is.character(alphabet)) { TYPE <- 'character' } }
  if (is.null(TYPE))                                           { stop('Error: alphabet should be a numeric or character vector') }

  #Check input string
  if (!is.vector(string))                                      { stop('Error: string must be a vector') }
  if ((TYPE ==   'numeric')&&(!is.numeric(string)))            { stop('Error: string must be the same type as alphabet') }
  if ((TYPE == 'character')&&(!is.character(string)))          { stop('Error: string must be the same type as alphabet') }
  m <- length(string)
  for (i in 1:m) {
    if (!(string[i] %in% alphabet))                            { stop(paste0('Error: string element ', i, ' is not in the alphabet')) } }

  #Check start.state.probs (if not 'stationary')
  if (!missing(start.state.probs)) {
    if (!is.vector(start.state.probs))                         { stop('Error: start.state.probs must be a probability vector (if specified)') }
    if (is.character(start.state.probs)) {
      if (length(start.state.probs) != 1)                      { stop('Error: start.state.probs must be a probability vector (or \'stationary\')') }
      NCHAR <- nchar(start.state.probs[1])
      STATN <- substr('stationary', 1, NCHAR)
      if (start.state.probs != STATN)                          { stop('Error: start.state.probs must be a probability vector (or \'stationary\')') } }
    if (!is.character(start.state.probs)) {
      if (!is.numeric(start.state.probs))                      { stop('Error: start.state.probs must be a probability vector (if specified)') }
      if (length(start.state.probs) != max.h+1)                { stop('Error: start.state.probs must have one entry for each state (if specified)') }
      if (min(start.state.probs) < 0)                          { stop('Error: start.state.probs must be a probability vector (if specified)') }
      if (sum(start.state.probs) != 1)                         { stop('Error: start.state.probs must be a probability vector (if specified)') } } }

  #Check array (if provided)
  if (!missing(array)) {
    if (!identical(string, array$string))                      { stop('Error: The array you have provided does not match your string') }
    if (!identical(probs, array$probs))                        { stop('Error: The array you have provided does not match your probs') }
    if (!identical(alphabet, array$alphabet))                  { stop('Error: The array you have provided does not match your alphabet') }
    if (!identical(allow.overlap, array$allow.overlap))        { stop('Error: The array you have provided does not match your allow.overlap') }
    if (!missing(start.state.probs)) {
      if (!identical(start.state.probs, array$start.state.probs)) { stop('Error: The array you have provided does not match your start.state.probs') } }
    arr.dim   <- dim(array$array)
    arr.max.n <- arr.dim[1]-1
    arr.max.r <- arr.dim[2]-1
    arr.max.h <- arr.dim[3]-1 }

  #Get DFA of language and extraction information
  MATS     <- DFA(string = string, alphabet = alphabet, probs = probs, allow.overlap = allow.overlap)
  LOGH     <- log(MATS$transition.probs)
  MO       <- MATS$max.overlap
  hh.start <- 0
  hh.count <- MATS$state.count

  #Set starting state log-probabilities (HH0)
  #The value m0 is the highest state with a non-zero probability (which affects the permissible range of r)
  if (missing(start.state.probs)) {
    HH0 <- rep(-Inf, nrow(LOGH))
    HH0[hh.start+1] <- 0
  } else {
    if (is.character(start.state.probs)) {
      HH0 <- log(MATS$stationary.probs[1,])
      HH0 <- HH0 - matrixStats::logSumExp(HH0) }
    if (is.numeric(start.state.probs)) {
      HH0 <- log(start.state.probs)
      HH0 <- HH0 - matrixStats::logSumExp(HH0) } }
  m0 <- max(which(HH0 != -Inf))-1

  #set values
  max.n <- max.size
  max.r <- max(0, 1+floor((max.n-m+m0)/(m-MO[1])))
  max.h <- nrow(LOGH)-1
  if (!missing(array)) {
  if (!identical(max.h, arr.max.h))                            { stop('Error: The array you provided does not use the correct number of states for the string') } }

  #Generate array of log-probabilities
  ARRAY <- array(-Inf, dim = c(max.n+1, max.r+1, m+1),
                 dimnames = list(sprintf('n[%s]', 0:max.n), sprintf('r[%s]', 0:max.r), sprintf('h[%s]', 0:m)))

  #Take vaules from array (if provided)
  if (!missing(array)) {

    top.n <- min(max.n, arr.max.n)
    top.r <- min(max.r, arr.max.r)
    ARRAY[1:(top.n+1), 1:(top.r+1), ] <- array$array
    if (max.n > top.n) {
    for (nn in (top.n+1):max.n) {
      UR <- max(0, 1+floor((nn-m+m0)/(m-MO[1])))
      for (rr in 0:UR)    {
        for (hh in 0:max.h) {
          if (hh.count[hh+1, 1] == 1) {
            if (rr > 0) {
              ARRAY[nn+1, rr+1, hh+1] <- matrixStats::logSumExp(LOGH[, hh+1] + ARRAY[nn, rr, ])   } } else {
            ARRAY[nn+1, rr+1, hh+1] <- matrixStats::logSumExp(LOGH[, hh+1] + ARRAY[nn, rr+1, ]) } } }
      ERROR <- matrixStats::logSumExp(ARRAY[nn+1, , ])
      ARRAY[nn+1, , ] <- ARRAY[nn+1, , ] - ERROR } }

  } else {

    ARRAY[1, 1, ] <- HH0
    for (nn in 1:max.n) {
      UR <- max(0, 1+floor((nn-m+m0)/(m-MO[1])))
      for (rr in 0:UR)    {
      for (hh in 0:max.h) {
        if (hh.count[hh+1, 1] == 1) {
          if (rr > 0) {
            ARRAY[nn+1, rr+1, hh+1] <- matrixStats::logSumExp(LOGH[, hh+1] + ARRAY[nn, rr, ])   } } else {
          ARRAY[nn+1, rr+1, hh+1] <- matrixStats::logSumExp(LOGH[, hh+1] + ARRAY[nn, rr+1, ]) } } }
      ERROR <- matrixStats::logSumExp(ARRAY[nn+1, , ])
      ARRAY[nn+1, , ] <- ARRAY[nn+1, , ] - ERROR } }

  #Return output
  list(array = ARRAY, max.size = max.size,
       string = string, probs = probs, alphabet = alphabet,
       allow.overlap = allow.overlap, start.state.probs = start.state.probs) }
