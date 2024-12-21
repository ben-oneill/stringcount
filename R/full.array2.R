#' Joint Distribution of String-Count and String-State Variables
#'
#' \code{full.array} returns the array of joint log-probabilities for a set of string-count and string-state variables.  This is the full
#' object that is computed in the \code{stringcountdist} functions (density, cumulative distribution, quantile function) in order to obtain
#' the probability functions for the String Count Distribution.  Those other functions allow the user to input the pre-computed array object,
#' which can be computed and stored using this function.
#'
#' **Outputs for sublanguages:** By default the funciton will produce the the array for the full language, which means that each string will have one
#' string-count variable in the array.  However, the user may use the \code{sublanguages} input to specify an alternative sublanguage for the output.
#' This input can be specified as a binary matrix with one row for each string in the language and any number of columns, where each column is
#' interpreted as a sublanguage including all strings with a positive binary indicator.  Alternatively, this input can be specified as a list of
#' vectors of string indices, where each vector is interpreted as a sublanguage including all strings with the specified indices.
#'
#' **Note:** By default the function takes the inputs \code{string}, \code{language} and \code{alphabet} in "simple form", where each string is written
#' as a single character string.  This form is fine so long as the symbols in the alphabet are single characters.  If the user wishes to use symbols
#' in the alphabet that are composed of multiple characters, this can be done by using the underlying non-simple form for the language.  To do the
#' latter the user should set \code{simple.form = FALSE} and express all strings/alphabet as vectors with each element being a symbol from the alphabet.
#' Also note that if the \code{alphabet} is unspecified then the function will assume that it is composed of all symbols that appear in the language
#' (and no other symbols).
#'
#' @usage \code{full.array()}
#' @param max.size The maximum size argument (a positive integer)
#' @param string A string expressed as a numeric/character vector
#' @param language A list of strings with each string expressed as a numeric/character vector
#' @param sublanguages A matrix or list specifying a class of sublanguages within the language (see details above)
#' @param alphabet A numeric/character vector containing the alphabet for the analysis
#' @param probs A vector of the symbol probabilities (taken over the symbols in \code{alphabet})
#' @param allow.overlap Logical; if \code{TRUE} then string occurrances are counted even if they overlap with previously counted occurrances
#' @param simple.form Logical; if \code{TRUE} the inputs \code{string}, \code{language} and \code{alphabet} are taken in simple form (see note)
#' @param minimise/minimize Logical; if \code{TRUE} then the DFA is minimised
#' @param start.state.probs A vector of probabilities for the starting state (or input 'stationary' to use the stationary distribution)
#' @param array Optional input for an array of joint log-probabilities for the string-count and state variable
#' @param log Logical; if \code{TRUE} the log-probabilities are returned; if \code{FALSE} the probabilities are returned
#' @return The full array of log-probability/probability values from the joint distribution of the string-count and string-state

full.array <- function(max.size, string = NULL, language = NULL, sublanguages = NULL, alphabet = NULL, probs,
                       allow.overlap = TRUE, simple.form = TRUE, minimise = FALSE, minimize = NULL,
                       start.state.probs = NULL, array = NULL, log = TRUE) {

  #Check input simple.form
  if (!is.vector(simple.form))                                 { stop('Error: Input simple.form should be a single logical value') }
  if (length(simple.form) != 1)                                { stop('Error: Input simple.form should be a single logical value') }
  if (!is.logical(simple.form))                                { stop('Error: Input simple.form should be a single logical value') }

  #Check input max.size
  if (!is.vector(max.size))                                    { stop('Error: max.size must be a positive integer') }
  if (!is.numeric(max.size))                                   { stop('Error: max.size must be a positive integer') }
  if (length(max.size) != 1)                                   { stop('Error: max.size must be a single positive integer') }
  if (as.integer(max.size) != max.size)                        { stop('Error: max.size must be a positive integer') }
  if (min(max.size) < 1)                                       { stop('Error: max.size must be a positive integer') }

  #Check inputs string and language
  if ((missing(string))&(missing(language)))                   { stop('Error: You must input a string or language') }
  if ((!missing(string))&(!missing(language)))                 { stop('Error: Input a string or language, but not both') }
  if (!missing(string))   {
    if (simple.form) { LANGUAGE <- list(sc(string)) } else { LANGUAGE <- list(string) }
    S <- 1
    names(LANGUAGE) <- 'String[1]' }
  if (!missing(language)) {
    if (simple.form) { LANGUAGE <- lapply(language, sc) } else { LANGUAGE <- language }
    S <- length(LANGUAGE)
    names(LANGUAGE) <- sprintf('String[%s]', 1:S) }

  #Check input alphabet
  STRING.SYMBOLS <- sort(unique(unlist(LANGUAGE)))
  if (missing(alphabet)) {
    ALPHABET <- STRING.SYMBOLS } else {
      if (simple.form) { ALPHABET <- sc(alphabet) } else { ALPHABET <- alphabet }
      if (!is.vector(as.vector(ALPHABET)))                     { stop('Error: Input alphabet must be a vector') }
      if (length(ALPHABET) != length(unique(ALPHABET))) {
        warning('Input alphabet contained duplicate elements --- duplicate elements were ignored')
        ALPHABET <- sort(unique(ALPHABET)) }
      for (i in 1:length(STRING.SYMBOLS)) {
        if (!(STRING.SYMBOLS[i] %in% ALPHABET)) {
          stop(paste0('Error: Symbol ', STRING.SYMBOLS[i], ' is in a string, but not in the alphabet')) } } }
  K <- length(ALPHABET)

  #Check input probs
  if (!is.vector(as.vector(probs)))                            { stop('Error: probs must be a probability vector') }
  if (!is.numeric(probs))                                      { stop('Error: probs must be a probability vector') }
  if (length(probs) == 0)                                      { stop('Error: probs must have at least one value') }
  if (min(probs) < 0)                                          { stop('Error: probs must be a probability vector') }
  if (sum(probs) != 1)                                         { stop('Error: probs must be a probability vector') }
  if (length(probs) != K)                                      { stop('Error: probs must have the same length as alphabet') }
  PROBS <- probs

  #Check input allow.overlap
  if (!is.vector(allow.overlap))                               { stop('Error: Input allow.overlap should be a single logical value') }
  if (length(allow.overlap) != 1)                              { stop('Error: Input allow.overlap should be a single logical value') }
  if (!is.logical(allow.overlap))                              { stop('Error: Input allow.overlap should be a single logical value') }
  ALLOW.OVERLAP <- allow.overlap

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

  #Check input minimise
  if ((missing(minimise))&(missing(minimize))) {
    MINIMISE <- minimise }
  if ((!missing(minimise))&(!missing(minimize))) {
    if (!identical(minimise, minimize))                        { stop('Error: Input minimise/minimize, but not both arguments') }
    if (!is.vector(minimise))                                  { stop('Error: Input minimise should be a single logical value') }
    if (length(minimise) != 1)                                 { stop('Error: Input minimise should be a single logical value') }
    if (!is.logical(minimise))                                 { stop('Error: Input minimise should be a single logical value') }
    MINIMISE <- minimise }
  if ((missing(minimise))&(!missing(minimize))) {
    if (!is.vector(minimize))                                  { stop('Error: Input minimize should be a single logical value') }
    if (length(minimize) != 1)                                 { stop('Error: Input minimize should be a single logical value') }
    if (!is.logical(minimize))                                 { stop('Error: Input minimize should be a single logical value') }
    MINIMISE <- minimize }
  if ((!missing(minimise))&(missing(minimize))) {
    if (!is.vector(minimise))                                  { stop('Error: Input minimise should be a single logical value') }
    if (length(minimise) != 1)                                 { stop('Error: Input minimise should be a single logical value') }
    if (!is.logical(minimise))                                 { stop('Error: Input minimise should be a single logical value') }
    MINIMISE <- minimise }

  #Check input sublanguages
  #Convert input to a list if it is not already in this form
  SUB.INCL <- FALSE
  if (missing(sublanguages)) {
    SUBLANGUAGES <- as.list(1:S)
    D <- S
  } else {
    SUB.INCL <- TRUE
    if ((!is.list(sublanguages))&(!is.matrix(sublanguages)))   { stop('Error: Input sublanguages should be a list or binary matrix') }
    if (is.list(sublanguages)) {
      SUBLANGUAGES <- sublanguages
      D <- length(SUBLANGUAGES)
      for (d in 1:D) {
        SUB <- SUBLANGUAGES[[d]]
        if (!is.vector(SUB))                                   {
          stop(paste0('Error: Input sublanguage ', d, ' is not in the proper form')) }
        if (!is.numeric(SUB))                                  {
          stop(paste0('Error: Input sublanguage ', d, ' is not in the proper form')) }
        SUB <- sort(unique(SUB))
        for (i in 1:length(SUB)) {
          if (!(SUB[i] %in% 1:S)) {
            stop(paste0('Error: Input sublanguage ', d, ' refers to string ', SUB[i], ', which is not present in the language')) }
          SUBLANGUAGES[[d]] <- SUB } } }
    if (is.matrix(sublanguages)) {
      if (nrow(sublanguages) != S)                             { stop('Error: Input sublanguages should have one row for each string in the language') }
      if (!all(MM == as.logical(MM)))                          { stop('Error: Input sublanguages should be a binary/logical matrix') }
      D <- ncol(sublanguages)
      SUBLANGUAGES <- vector(mode = 'list', length = D)
      for (d in 1:D) { SUBLANGUAGES[[d]] <- which(sublanguages[, d] == 1) } } }

  #Check array (if provided)
  #FIX THIS FOR MULTIPLE STRINGS
  if (!missing(array)) {
    if (!identical(PROBS, array$probs))                        { stop('Error: The array you have provided does not match your probs') }
    if (!identical(ALPHABET, array$alphabet))                  { stop('Error: The array you have provided does not match your alphabet') }
    if (!identical(LANGUAGE, array$language))                  { stop('Error: The array you have provided does not match your language') }
    if (!identical(SUBLANGUAGES, array$sublanguages))          { stop('Error: The array you have provided does not match your sublanguages') }
    if (!identical(ALLOW.OVERLAP, array$allow.overlap))        { stop('Error: The array you have provided does not match your allow.overlap') }
    if (!identical(MINIMISE, array$minimise))                  { stop('Error: The array you have provided does not match your minimise') }
    if (!missing(start.state.probs)) {
      if (!identical(start.state.probs, array$start.state.probs)) { stop('Error: The array you have provided does not match your start.state.probs') } }
    arr.dim   <- dim(array$array)
    arr.max.n <- arr.dim[1]-1
    arr.max.q <- arr.dim[2]-1
    arr.max.r <- arr.dim[3:length(arr.dim)]-1 }

  #Check input log
  if (!is.vector(log))                                         { stop('Error: Input log should be a single logical value') }
  if (length(log) != 1)                                        { stop('Error: Input log should be a single logical value') }
  if (!is.logical(log))                                        { stop('Error: Input log should be a single logical value') }

  #-----------------------------------------------------------------------------------------------------------
  #Compute the DFA using the function inputs
  #-----------------------------------------------------------------------------------------------------------

  #Get DFA of language and extraction information
  DFA.INFO <- DFA(language = LANGUAGE, sublanguages = SUBLANGUAGES, alphabet = ALPHABET, probs = PROBS,
                  allow.overlap = ALLOW.OVERLAP, simple.form = FALSE, minimise = MINIMISE)
  LANGUAGE <- DFA.INFO$language
  ALPHABET <- DFA.INFO$alphabet
  SUBS     <- DFA.INFO$sublanguages
  MC       <- attributes(SUBS)$minimum.completion
  S        <- length(LANGUAGE)
  D        <- length(SUBS)
  K        <- length(ALPHABET)
  LOGTT    <- log(DFA.INFO$transition.probs)
  FF       <- DFA.INFO$state.count
  MINIMISE <- DFA.INFO$minimise

  #Set starting state log-probabilities (QQ0)
  if (missing(start.state.probs)) {
    QQ0 <- rep(-Inf, nrow(LOGTT))
    QQ0[1] <- 0
  } else {
    if (is.character(start.state.probs)) {
      QQ0 <- log(DFA.INFO$stationary.probs[1,])
      QQ0 <- QQ0 - matrixStats::logSumExp(QQ0) }
    if (is.numeric(start.state.probs)) {
      QQ0 <- log(start.state.probs)
      QQ0 <- QQ0 - matrixStats::logSumExp(QQ0) } }

  #set maximum values
  max.n <- max.size
  max.q <- nrow(LOGTT)-1
  max.r <- rep(0, D)
  for (d in 1:D) {
    max.r[d] <- max(0, 1+floor((max.n-1)/MC[d])) }
  if (!missing(array)) {
    if (!identical(max.q, arr.max.q))                          { stop('Error: The array you provided does not use the correct number of states for the string') }
    if (!identical(max.r, arr.max.r))                          { stop('Error: The array you provided does not use the correct number of states for the string') } }

  #Generate array of log-probabilities
  if (D == 1) {
    COUNT.NAMES <- list(sprintf('r[%s]', 0:max.r[1]))
  } else {
    COUNT.NAMES <- vector(mode = 'list', length = D)
    for (d in 1:D) { COUNT.NAMES[[d]] <- paste0('r', d, sprintf('[%s]', 0:max.r[d])) } }


  #Take values from array (if provided)
  if (!missing(array)) {

    ARRAY <- array(-Inf, dim = c(max.n+1, max.q+1, max.r+1),
                   dimnames = c(list(sprintf('n[%s]', 0:max.n), sprintf('State[%s]', 0:max.q)), COUNT.NAMES))

    top.n <- min(max.n, arr.max.n)
    top.r <- pmin(max.r, arr.max.r)
    ARRAY[1:(top.n+1), , 1:(top.r+1)] <- array$array
    if (max.n > top.n) {
    for (nn in (top.n+1):max.n) {

      #ADD THIS COMPUTATION

      } }

  } else {

    #Set error vector
    ERRORS <- rep(-Inf, max.n)
    names(ERRORS) <- sprintf('n[%s]', 1:max.n)

    #Create the array
    ARRAY <- array(-Inf, dim = c(max.n+1, max.q+1, max.r+1),
                   dimnames = c(list(sprintf('n[%s]', 0:max.n), sprintf('State[%s]', 0:max.q)), COUNT.NAMES))

    #Set array by forward propagation
    IND <- matrix(1, nrow = max.q+1, ncol = D+2)
    IND[,2] <- 0:max.q + 1
    ARRAY[IND] <- QQ0
    if (max.n > 0) {
    for (nn in 0:(max.n-1)) {

      #Set count vector range
      RR <- vector(length = D, mode = 'list')
      names(RR) <- sprintf('r%s', 1:D)
      for (dd in 1:D) { RR[[dd]] <- 0:max(0, 1+floor((nn-1)/MC[dd])) }
      R.VECS <- expand.grid(RR)

      #Propagate forward
      for (qq in 0:max.q) {
        FFF <- unname(FF[qq+1, ])
        for (ii in 0:max.q) {
        for (rr in 1:nrow(R.VECS)) {
          r.vec <- unname(unlist(R.VECS[rr, ]))
          if (all(r.vec+FFF <= max.r)) {
            IND1 <- matrix(c(nn+2, qq+1, r.vec+FFF+1), nrow = 1)
            IND2 <- matrix(c(nn+1, ii+1, r.vec+1), nrow = 1)
            T1   <- ARRAY[IND1]
            T2   <- LOGTT[ii+1, qq+1] + ARRAY[IND2]
            ARRAY[IND1] <- matrixStats::logSumExp(c(T1, T2)) } } } }

      #Ensure each probability vector sums to one
      RR <- vector(length = D+1, mode = 'list')
      names(RR) <- c('q', sprintf('r%s', 1:D))
      RR[[1]] <- 0:max.q
      for (dd in 1:D) { RR[[dd+1]] <- 0:max.r[dd] }
      IND0 <- expand.grid(RR)
      IND  <- as.matrix(cbind(nn+1, IND0))
      ERROR <- matrixStats::logSumExp(ARRAY[IND+1])
      ARRAY[IND+1] <- ARRAY[IND+1] - ERROR
      ERRORS[nn+1] <- ERROR } } }

  #Return output
  if (!log) { ARRAY <- exp(ARRAY) }
  list(array = ARRAY, max.size = max.size,
       language = LANGUAGE, sublanguages = SUBS, alphabet = ALPHABET,
       allow.overlap = ALLOW.OVERLAP, minimise = MINIMISE, probs = PROBS, start.state.probs = start.state.probs,
       errors = ERRORS, log = log) }
