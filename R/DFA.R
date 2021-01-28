#' Deterministic Finite Automaton (DFA) for a string or language
#'
#' \code{DFA} returns the Deterministic Finite Automaton (DFA) for a string or language (a list of strings).
#'
#' Any finite language (list of strings) has a minimal Deterministic Finite Automaton (DFA).  This function takes either a single string input
#' or a finite language and computes the minimal DFA.  This function gives a variation on the DFA where the final-states for each string in the
#' language is returned separately, but these can easily be combined if required.  The output is an object of class \code{'DFA'} which contains the
#' transition table, progress table, and final-state table for the DFA, along with other useful information about the language and alphabet.
#' The states for the DFA are always labelled with consecutive non-negative integers and the starting state is always the zero state.  The user
#' also has the option of giving a probability vector for the symbols in the alphabet; if this is included then the DFA includes this probabilty
#' vector, the resulting stationary probabilities for the states, and the transition-probability matrix for the states.
#'
#' Note: By default the function takes the inputs \code{string}, \code{language} and \code{alphabet} in "simple form", where each string is written
#' as a single character string.  This form is fine so long as the symbols in the alphabet are single characters.  If the user wishes to use symbols
#' in the alphabet that are composed of multiple characters, this can be done by using the underlying non-simple form for the language.  To do the
#' latter the user should set \code{simple.form = FALSE} and express all strings/alphabet as vectors with each element being a symbol from the alphabet.
#' Also note that if the \code{alphabet} is unspecified then the function will assume that it is composed of all symbols that appear in the language
#' (and no other symbols).
#'
#' @usage \code{DFA()}
#' @param string A string expressed as a character vector
#' @param language A list of strings with each string expressed as a character vector
#' @param alphabet A numeric/character vector containing the alphabet for the analysis
#' @param probs A vector of the symbol probabilities (taken over the symbols in the \code{alphabet})
#' @param allow.overlap Logical; if \code{TRUE} then string occurrances are counted even if they overlap with previously counted occurrances
#' @param simple.form Logical; if \code{TRUE} the inputs \code{string}, \code{language} and \code{alphabet} are taken in simple form (see note)
#' @return A list of class 'DFA' containing the deterministic finite automaton (DFA)

DFA <- function(string = NULL, language = NULL, alphabet = NULL, probs = NULL, allow.overlap = TRUE, simple.form = TRUE) {

  #Get input names
  STRING.NAME   <- deparse(substitute(string))
  LANGUAGE.NAME <- deparse(substitute(language))
  ALPHABET.NAME <- deparse(substitute(alphabet))

  #Check input simple.form
  if (!is.vector(simple.form))                                 { stop('Error: Input simple.form should be a single logical value') }
  if (length(simple.form) != 1)                                { stop('Error: Input simple.form should be a single logical value') }
  if (!is.logical(simple.form))                                { stop('Error: Input simple.form should be a single logical value') }

  #Check inputs string and language
  if ((missing(string))&(missing(language)))                   { stop('Error: You must input a string or language') }
  if ((!missing(string))&(!missing(language)))                 { stop('Error: Input a string or language, but not both') }
  if (!missing(string))   {
    LANGUAGE.NAME <- STRING.NAME
    if (simple.form) { LANGUAGE <- list(sc(string)) } else { LANGUAGE <- list(string) }
    S <- 1
    names(LANGUAGE) <- 'String[1]' }
  if (!missing(language)) {
    if (simple.form) { LANGUAGE <- lapply(language, sc) } else { LANGUAGE <- language }
    S <- length(LANGUAGE)
    names(LANGUAGE) <- sprintf('String[%s]', 1:S) }
  attr(LANGUAGE, 'language.name') <- LANGUAGE.NAME

  #Check input alphabet
  STRING.SYMBOLS <- sort(unique(unlist(LANGUAGE)))
  if (missing(alphabet)) {
    ALPHABET.NAME <- NA
    ALPHABET <- STRING.SYMBOLS } else {
      if (simple.form) { ALPHABET <- sc(alphabet) } else { ALPHABET <- alphabet }
      if (!is.vector(as.vector(ALPHABET)))                     { stop('Error: Input alphabet must be a vector') }
      if (length(ALPHABET) != length(unique(ALPHABET))) {
        warning('Input alphabet contained duplicate elements --- duplicate elements were ignored')
        ALPHABET <- unique(ALPHABET) }
      for (i in 1:length(STRING.SYMBOLS)) {
        if (!(STRING.SYMBOLS[i] %in% ALPHABET)) {
          stop(paste0('Error: Symbol ', STRING.SYMBOLS[i], ' is in a string, but not in the alphabet')) } } }
  K <- length(ALPHABET)
  attr(ALPHABET, 'alphabet.name') <- ALPHABET.NAME

  #Check input probs
  if (!missing(probs)) {
    if (!is.vector(as.vector(probs)))                          { stop('Error: probs must be a probability vector') }
    if (!is.numeric(probs))                                    { stop('Error: probs must be a probability vector') }
    if (length(probs) == 0)                                    { stop('Error: probs must have at least one value') }
    if (min(probs) < 0)                                        { stop('Error: probs must be a probability vector') }
    if (sum(probs) != 1)                                       { stop('Error: probs must be a probability vector') }
    if (length(probs) != K)                                    { stop('Error: probs must have the same length as alphabet') } }

  #Check input allow.overlap
  if (!is.vector(allow.overlap))                               { stop('Error: Input allow.overlap should be a single logical value') }
  if (length(allow.overlap) != 1)                              { stop('Error: Input allow.overlap should be a single logical value') }
  if (!is.logical(allow.overlap))                              { stop('Error: Input allow.overlap should be a single logical value') }
  ALLOW.OVERLAP <- allow.overlap

  #Compute individual DFAs, etc., for each string in the language
  DFA <- vector(mode = 'list', length = S)

  for (s in 1:S) {

    #Extract the string and get its length
    STRING <- LANGUAGE[[s]]
    MATCH  <- match(STRING, ALPHABET)
    m      <- length(STRING)

    #Generate the transition matrix
    TT <- matrix(0, nrow = m+1, ncol = K)
    rownames(TT) <- sprintf('State[%s]', 0:m)
    colnames(TT) <- ALPHABET
    for (i in 0:m) {
    for (x in 1:K) {
      t     <- i+1
      BREAK <- FALSE
      while (!BREAK) {
        if (t > 1) { VEC1 <- c(STRING[(i-t+2):i], ALPHABET[x]) } else { VEC1 <- ALPHABET[x] }
        VEC2 <- STRING[1:t]
        if (identical(VEC1, VEC2)) { BREAK <- TRUE } else { t <- t-1 }
        if (t == 0) { BREAK <- TRUE } }
      TT[i+1, x] <- t } }

    #Adjust last row of structure matrix if overlap is disallowed
    if (!ALLOW.OVERLAP) {
      TT[m+1, ] <- 0
      TT[m+1, MATCH[1]] <- 1 }

    #Compute progress matrix
    PP <- matrix(0:m, nrow = m+1, ncol = 1)
    rownames(PP) <- sprintf('State[%s]', 0:m)
    colnames(PP) <- 'String[1]'

    #Compute maximum overlap value
    OVERLAP <- max(TT[m+1,])-1
    names(OVERLAP) <- 'String[1]'

    #Generate final state matrix
    FF           <- matrix(0, nrow = m+1, ncol = 1)
    FF[m+1,]     <- 1
    rownames(FF) <- sprintf('State[%s]', 0:m)
    colnames(FF) <- 'String[1]'

    DFA[[s]] <- list(transition.table = TT, progress = PP, state.count = FF, max.overlap = OVERLAP,
                     language = list('String[1]' = STRING), alphabet = ALPHABET, allow.overlap = ALLOW.OVERLAP)
    class(DFA[[s]]) <- 'DFA' }

  #Combine and minimise DFAs
  DFA.all <- DFA[[1]]
  if (S > 1) {
    for (s in 2:S) { DFA.all <- combine.DFA(DFA.all, DFA[[s]]) } }
  DFA <- DFA.all
  TT  <- DFA$transition.table
  m   <- nrow(TT)-1

  #Generate the transition probability matrix and stationary vectors
  if (!missing(probs)) {

    #Add probabilities
    DFA$probs <- matrix(probs, nrow = 1)
    rownames(DFA$probs) <- 'Probability'
    colnames(DFA$probs) <- ALPHABET

    #Transition probability matrix
    H <- matrix(0, nrow = m+1, ncol = m+1)
    rownames(H) <- sprintf('State[%s]', 0:m)
    colnames(H) <- sprintf('State[%s]', 0:m)
    for (i in 0:m) {
    for (x in 1:K) {
        H[i+1, TT[i+1, x]+1] <- H[i+1, TT[i+1, x]+1] + probs[x] } }

    #Compute the eigendecomposition
    EIGEN <- eigen(t(H))
    VALS  <- EIGEN$values
    VECS  <- EIGEN$vectors

    #Stationary vector(s)
    RES <- zapsmall(Mod(VALS - as.complex(rep(1, m+1))))
    IND <- which(RES == 0)
    SSS <- Re(VECS[, IND])
    STR <- matrix(SSS/sum(SSS), nrow = 1)
    rownames(STR) <- 'Stationary'
    colnames(STR) <- rownames(H)

    #Add transition matrix and stationary vector to DFA
    DFA$stationary.probs <- STR
    DFA$transition.probs <- H }

  #Give output
  DFA }


print.DFA <- function(object) {

  if (class(object) != 'DFA')                    { stop('Error: This function only operates on objects of class \'DFA\'') }

  #Extract information
  LANGUAGE.NAME <- attributes(object$language)$language.name
  LANGUAGE      <- object$language
  STRING.NAMES  <- names(LANGUAGE)
  ALPHABET.NAME <- attributes(object$alphabet)$alphabet.name
  ALPHABET      <- object$alphabet
  ALLOW.OVERLAP <- object$allow.overlap
  OVERLAP       <- object$max.overlap
  PROBS.INCL    <- FALSE
  if (!is.null(object$probs)) {
    PROBS.INCL  <- TRUE
    PROBS       <- object$probs }
  TRANSITIONS   <- object$transition.table
  PROGRESS      <- object$progress
  FINALSTATES   <- object$state.count

  #Set parameters
  S  <- length(LANGUAGE)
  M  <- nrow(TRANSITIONS)
  M0 <- ifelse(PROBS.INCL, M, 6)
  L1 <- paste0(rep('-', floor(9*(M0-4)/2)), collapse = '')
  L2 <- paste0(rep('-', ceiling(9*(M0-4)/2)), collapse = '')

  #Print title and language details
  cat('\n    Deterministic Finite Automaton (DFA) \n \n')
  PP <- ifelse(PROBS.INCL, 'with a specified probability vector', '')
  if (is.na(ALPHABET.NAME)) {
    if (S == 1) { cat('DFA for a single string', LANGUAGE.NAME, 'from an unspecified alphabet', PP, '\n') } else {
                  cat('DFA for a language', LANGUAGE.NAME, 'containing', S, 'strings \n') } }
  if (!is.na(ALPHABET.NAME)) {
    if (S == 1) { cat('DFA for a single string', LANGUAGE.NAME, 'from alphabet', ALPHABET.NAME, PP, '\n') } else {
      cat('DFA for a language', LANGUAGE.NAME, 'containing', S, 'strings from alphabet', ALPHABET.NAME, '\n') } }
  if (!ALLOW.OVERLAP) { cat('(Overlap in strings is not allowed) \n') }
  cat('\n')

  #Print transition table and counting (final) states
  cat(paste0(L1, '--Transitions, progress, and final states---', L2), '\n \n')
  DASHES <- matrix(rep('|', M), nrow = M, dimnames = list(NULL, '|'))
  print(cbind(DASHES, TRANSITIONS, DASHES, PROGRESS, DASHES, FINALSTATES), quote = FALSE)
  cat('\n')

  #Print probability vector, transition probability matrix and stationary state probabilities (if present)
  if (PROBS.INCL) {
    cat(paste0(L1, '------------Symbol probabilities------------', L2), '\n \n')
    print(object$probs)
    cat('\n')
    cat(paste0(L1, '-------Stationary state probabilities-------', L2), '\n \n')
    print(object$stationary.probs)
    cat('\n')
    cat(paste0(L1, '--------Transition probability matrix-------', L2), '\n \n')
    print(object$transition.probs)
    cat('\n') }

  #Print information on alphabet and language
  cat(paste0(L1, '------Alphabet and language information-----', L2), '\n \n')
  SSS <- ifelse(length(ALPHABET) == 1, 'symbol', 'symbols')
  cat(ALPHABET.NAME, 'is an alphabet containing', length(ALPHABET), SSS, '\n', c(ALPHABET), '\n \n')
  for (s in 1:S) {
    SSS <- ifelse(length(LANGUAGE[[s]]) == 1, 'symbol', 'symbols')
    cat(STRING.NAMES[s], 'contains', length(LANGUAGE[[s]]), SSS, 'with maximum overlap', OVERLAP[s], '\n', c(LANGUAGE[[s]]), '\n \n') }
  cat(paste0(L1, '--------------------------------------------', L2), '\n \n') }


