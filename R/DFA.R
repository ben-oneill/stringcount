#' Deterministic Finite Automaton (DFA) for a string or language
#'
#' \code{DFA} returns the Deterministic Finite Automaton (DFA) for a string or language (a list of strings).
#'
#' Any finite language (list of strings) has a minimal Deterministic Finite Automaton (DFA).  This function takes either a single string input
#' or a finite language and computes the minimal DFA.  By default, the DFA is computed for the set of individual strings in the language, but
#' the user may use the \code{sublanguage} input to specify an arbitrary sublanguage, and in this case the DFA is computed for the sublanguage
#' (see details below).  This function actually gives a variation on the DFA where there can be multiple setes of final-states (e.g., one for
#' each string) manifesting in multiple columns of final-state indicators.  The output is an object of class \code{'DFA'} which contains the
#' transition table, progress table, and final-state table for the DFA, along with other useful information about the language and alphabet.
#' The states for the DFA are always labelled with consecutive non-negative integers and the starting state is always the zero state.  The user
#' also has the option of giving a probability vector for the symbols in the alphabet; if this is included then the DFA includes this probabilty
#' vector, the resulting stationary probabilities for the states, and the transition-probability matrix for the states.
#'
#' **Outputs for sublanguages:** By default the funciton will produce the DFA for the full language, which means that each string will have one
#' column in the final-state table.  However, the user may use the \code{sublanguages} input to specify an alternative sublanguage for the output.
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
#' @usage \code{DFA()}
#' @param string A string expressed as a numeric/character vector
#' @param language A list of strings with each string expressed as a numeric/character vector
#' @param sublanguages A matrix or list specifying a class of sublanguages within the language (see details above)
#' @param alphabet A numeric/character vector containing the alphabet for the analysis
#' @param probs A vector of the symbol probabilities (taken over the symbols in the \code{alphabet})
#' @param allow.overlap Logical; if \code{TRUE} then string occurrances are counted even if they overlap with previously counted occurrances
#' @param simple.form Logical; if \code{TRUE} the inputs \code{string}, \code{language} and \code{alphabet} are taken in simple form (see note)
#' @param minimise/minimize Logical; if \code{TRUE} then the DFA is minimised
#' @return A list of class 'DFA' containing the deterministic finite automaton (DFA)

DFA <- function(string = NULL, language = NULL, sublanguages = NULL, alphabet = NULL, probs = NULL,
                allow.overlap = TRUE, simple.form = TRUE, minimise = FALSE, minimize = NULL) {

  #Check input simple.form
  if (!is.vector(simple.form))                                 { stop('Error: Input simple.form should be a single logical value') }
  if (length(simple.form) != 1)                                { stop('Error: Input simple.form should be a single logical value') }
  if (!is.logical(simple.form))                                { stop('Error: Input simple.form should be a single logical value') }

  #Check inputs string and language
  if ((missing(string))&(missing(language)))                   { stop('Error: You must input a string or language') }
  if ((!missing(string))&(!missing(language)))                 { stop('Error: Input a string or language, but not both') }
  if (!missing(string))   {
    if (simple.form) { LANGUAGE <- list(sc(string))     } else { LANGUAGE <- list(string) }
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

  #-----------------------------------------------------------------------------------------------------------
  #Compute the combined DFA
  #This step computes individual DFAs and then combines them, removing impossible states (but not minimising)
  #-----------------------------------------------------------------------------------------------------------

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
    colnames(PP) <- 'P[1]'

    #Generate final state matrix
    FF           <- matrix(0, nrow = m+1, ncol = 1)
    FF[m+1,]     <- 1
    rownames(FF) <- sprintf('State[%s]', 0:m)
    colnames(FF) <- 'F[1]'

    #Compute minimum-completion of the string
    MC <- 0
    STATES <- which(FF[, 1] == 1)-1
    BREAK  <- FALSE
    while (!BREAK) {
      MC <- MC+1
      STATES <- unique(TT[STATES+1, ])
      if (sum(FF[STATES+1, 1]) > 0) { BREAK <- TRUE } }

    #Set language and sublanguage for individual DFA
    LANG1 <- list('String[1]' = STRING)
    SUB1  <- list('Sub[1]' = 1)
    attr(SUB1, 'minimum.completion') <- MC

    #Generate the DFA for the string
    DFA[[s]] <- list(transition.table = TT, progress = PP, state.count = FF,
                     language = LANG1, sublanguages = SUB1, alphabet = ALPHABET,
                     allow.overlap = ALLOW.OVERLAP, minimise = TRUE)
    class(DFA[[s]]) <- 'DFA' }

  #Combine the DFAs and reobtain objects
  DFA.all <- DFA[[1]]
  if (S > 1) {
    for (s in 2:S) { DFA.all <- combine.DFA(DFA.all, DFA[[s]], minimise = FALSE) } }
  DFA <- DFA.all
  TT  <- DFA$transition.table
  PP  <- DFA$progress
  FF  <- DFA$state.count
  m   <- nrow(TT)-1

  #Create output for sublanguages
  if (SUB.INCL) {
    FINAL <- DFA$state.count
    FF    <- matrix(0, nrow = m+1, D)
    COLNAMES <- character(D)
    for (i in 0:m) {
    for (d in 1:D) {
      VALS <- FINAL[i+1, SUBLANGUAGES[[d]]]
      FF[i+1, d] <- 1-prod(1-VALS)
      COLNAMES[d] <- paste0('F[', paste(SUBLANGUAGES[[d]], collapse = ','), ']') } }
    rownames(FF) <- sprintf('State[%s]', 0:m)
    colnames(FF) <- COLNAMES
    DFA$state.count <- FF }
  DFA$sublanguages <- SUBLANGUAGES
  names(DFA$sublanguages) <- sprintf('Sub[%s]', 1:D)

  #Compute the minimum-completion for each sublanguage
  MC <- rep(0, D)
  for (d in 1:D) {
    STATES <- which(FF[, d] == 1)-1
    BREAK  <- FALSE
    while (!BREAK) {
      MC[d] <- MC[d]+1
      STATES <- unique(TT[STATES+1, ])
      if (sum(FF[STATES+1, d]) > 0) { BREAK <- TRUE } } }
  attr(DFA$sublanguages, 'minimum.completion') <- MC

  #-----------------------------------------------------------------------------------------------------------
  #Minimise the DFA
  #This step uses minimisation using the Myhill-Nerode table method
  #-----------------------------------------------------------------------------------------------------------

  if (MINIMISE) {
  if (m > 0)    {

    #Set up Myhill-Nerode table
    #In the completed table, entries marked as TRUE will denote pairs of states that are equivalent
    MN.TABLE <- matrix(nrow = m+1, ncol = m+1)

    #Mark pairs as FALSE if they have different final-state vectors
    for (q1 in 0:m) {
    for (q2 in 0:m) {
    MN.TABLE[q1+1, q2+1] <- identical(FF[q1+1, ], FF[q2+1, ]) } }

    #Recursively mark Myhill-Nerode table to indicate non-equivalent states
    FILL <- TRUE
    while (FILL) {
      FILL <- FALSE
      for (q1 in 0:(m-1))  {
      for (q2 in (q1+1):m) {
        BREAK <- FALSE
        x <- 1
        while ((x <= K)&(!BREAK)) {
          if (MN.TABLE[q1+1, q2+1]) {
            t1 <- TT[q1+1, x]
            t2 <- TT[q2+1, x]
            if (!MN.TABLE[t1+1, t2+1]) {
              MN.TABLE[q1+1, q2+1] <- FALSE
              MN.TABLE[q2+1, q1+1] <- FALSE
              BREAK <- TRUE
              FILL  <- TRUE } }
           x <- x+1 } } } }

    #Create mapping for old and new states
    NEWSTATES <- logical(m+1)
    OLDSTATES <- logical(0)
    USED <- rep(FALSE, m+1)
    qq <- 0
    while (sum(USED) < m+1) {
      q <- min(which(!USED))-1
      EQUIV <- which(MN.TABLE[q+1, ])
      NEWSTATES[EQUIV] <- qq
      OLDSTATES[qq+1]  <- q
      USED[EQUIV]      <- TRUE
      qq <- qq+1 }
    mm <- qq-1

    if (!all(NEWSTATES == 0:m)) {

      #Create minimised transition table
      TT.MIN  <- matrix(0, nrow = mm+1, ncol = K)
      rownames(TT.MIN) <- sprintf('State[%s]', 0:mm)
      colnames(TT.MIN) <- ALPHABET
      FF.MIN  <- matrix(0, nrow = mm+1, ncol = ncol(FF))
      rownames(FF.MIN) <- sprintf('State[%s]', 0:mm)
      colnames(FF.MIN) <- colnames(FF)
      for (qq in 0:mm) {
        TRANS <- rep(0, K)
        FINAL <- rep(0, K)
        STATE <- OLDSTATES[qq+1]
        for (x in 1:K) {
          TRANS[x] <- NEWSTATES[TT[STATE+1, x]+1] }
        TT.MIN[qq+1, ] <- TRANS
        FF.MIN[qq+1, ] <- FF[STATE+1, ] }
      TT <- TT.MIN
      FF <- FF.MIN
      m <- mm

      #Update progress table
      PP.MIN <- PP
      rownames(PP.MIN) <- sprintf('State[%s]', NEWSTATES)
      PP.MIN <- PP.MIN[order(NEWSTATES), , drop = FALSE]
      PP <- PP.MIN

      #Replace tables in DFA object

      DFA$transition.table <- TT
      DFA$progress         <- PP
      DFA$state.count      <- FF } }

    DFA$minimise <- TRUE }

  #-----------------------------------------------------------------------------------------------------------
  #Add alphabet probabilities and transition probabilities (if required)
  #-----------------------------------------------------------------------------------------------------------

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
    for (x in 1:K)  {
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
  LANGUAGE      <- object$language
  STRING.NAMES  <- names(LANGUAGE)
  ALPHABET      <- object$alphabet
  ALLOW.OVERLAP <- object$allow.overlap
  PROBS.INCL    <- FALSE
  if (!is.null(object$probs)) {
    PROBS.INCL  <- TRUE
    PROBS       <- object$probs }
  TRANSITIONS   <- object$transition.table
  PROGRESS      <- object$progress
  FINALSTATES   <- object$state.count
  MINIMISE      <- object$minimise

  #Set parameters
  S  <- length(LANGUAGE)
  M  <- nrow(TRANSITIONS)
  M0 <- ifelse(PROBS.INCL, M, 6)
  L1 <- paste0(rep('-',   floor(9*(M0-4)/2)), collapse = '')
  L2 <- paste0(rep('-', ceiling(9*(M0-4)/2)), collapse = '')

  #Print title and language details
  cat('\n    Deterministic Finite Automaton (DFA) \n \n')
  MM <- ifelse(MINIMISE, 'Minimised ', '')
  PP <- ifelse(PROBS.INCL, 'with a specified probability vector', '')
  if (S == 1) { cat(paste0(MM, 'DFA for a single string'), '\n') } else {
                cat(paste0(MM, 'DFA for a language containing'), S, 'strings \n') }
  if (!ALLOW.OVERLAP) { cat('(String-counts for DFA do not count overlapping occurrences) \n') }
  cat('\n')

  #Print transition table and counting (final) states
  cat(paste0(L1, '--------Transitions and final states--------', L2), '\n \n')
  DASHES <- matrix(rep('|', M), nrow = M, dimnames = list(NULL, '|'))
  print(cbind(DASHES, TRANSITIONS, DASHES, FINALSTATES), quote = FALSE)
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
  cat('The alphabet contains', length(ALPHABET), SSS, '\n', c(ALPHABET), '\n \n')
  for (s in 1:S) {
    SSS <- ifelse(length(LANGUAGE[[s]]) == 1, 'symbol', 'symbols')
    cat(STRING.NAMES[s], 'contains', length(LANGUAGE[[s]]), SSS, '\n', c(LANGUAGE[[s]]), '\n \n') }
  cat(paste0(L1, '--------------------------------------------', L2), '\n \n') }


progress <- function(dfa) {

  if (class(dfa) != 'DFA')                       { stop('Error: This function only operates on objects of class \'DFA\'') }

  #Create output
  OUT <- dfa[c('transition.table', 'progress', 'state.count')]
  class(OUT) <- 'progress.DFA'

  #Return output
  OUT }


print.progress.DFA <- function(object) {

  if (class(object) != 'progress.DFA')           { stop('Error: This function only operates on objects of class \'progress.DFA\'') }

  #Extract information
  TRANSITIONS <- object$transition.table
  PROGRESS    <- object$progress
  FINALSTATES <- object$state.count

  #Produce output tables for transitions and final states
  m  <- nrow(PROGRESS)-1
  m0 <- nrow(TRANSITIONS)-1
  if (m == m0) {
    TT <- TRANSITIONS
    FF <- FINALSTATES
  } else {
    NEWSTATES <- as.integer(gsub(pattern = 'State\\[', '', gsub(pattern = '\\]', replacement = '', x = rownames(PROGRESS))))
    TT <- TRANSITIONS[NEWSTATES+1, , drop = FALSE]
    FF <- FINALSTATES[NEWSTATES+1, , drop = FALSE]
    RR <- rownames(TT)
    for (q in 1:m) {
      if (RR[q+1] == RR[q]) {
        TT[q+1, ] <- NA
        FF[q+1, ] <- NA
        rownames(TT)[q+1] <- ''
        rownames(TT)[q+1] <- '' } } }

  #Produce output matrix
  DASHES <- matrix(rep('|', m+1), nrow = m+1, ncol = 1, dimnames = list(NULL, '|'))
  OUT <- cbind(DASHES, TT, DASHES, PROGRESS, DASHES, FF)
  rownames(OUT) <- rownames(TT)

  #Print the output
  cat('\n--------Transitions, progress and final states-------- \n \n')
  prmatrix(OUT, quote = FALSE, na.print = '')
  cat('\n------------------------------------------------------ \n \n') }


sublanguages <- function(dfa) {

  if (class(dfa) != 'DFA')                       { stop('Error: This function only operates on objects of class \'DFA\'') }

  #Create output
  OUT <- dfa[c('language', 'sublanguages', 'alphabet')]
  class(OUT) <- 'sublanguages.DFA'

  #Return output
  OUT }


print.sublanguages.DFA <- function(object) {

  if (class(object) != 'sublanguages.DFA')       { stop('Error: This function only operates on objects of class \'progress.DFA\'') }

  #Extract information
  LANGUAGE     <- object$language
  SUBLANGUAGES <- object$sublanguages
  ALPHABET     <- object$alphabet
  STRING.NAMES <- names(LANGUAGE)
  SUB.NAMES    <- names(SUBLANGUAGES)
  MC           <- attributes(SUBLANGUAGES)$minimum.completion

  #Set parameters
  S <- length(LANGUAGE)
  D <- length(SUBLANGUAGES)

  #Print information on sublanguages
  cat('--------------------Sublanguage information------------------- \n \n')
  for (d in 1:D) {
    SSS <- ifelse(length(SUBLANGUAGES[[d]]) == 1, 'string', 'strings')
    cat(SUB.NAMES[d], 'contains', length(SUBLANGUAGES[[d]]), SSS, '(minimum completion =', MC[d], 'symbols) \n', c(SUBLANGUAGES[[d]]), '\n \n') }
  cat('-------------------------------------------------------------- \n \n') }


probs <- function(dfa) {

  if (class(dfa) != 'DFA')                       { stop('Error: This function only operates on objects of class \'DFA\'') }

  #Create output
  OUT <- dfa[c('probs', 'transition.probs', 'stationary.probs')]
  class(OUT) <- 'probs.DFA'

  #Return output
  OUT }


print.probs.DFA <- function(object) {

  if (class(object) != 'probs.DFA')           { stop('Error: This function only operates on objects of class \'progress.DFA\'') }

  if (is.null(object$probs)) {

    cat('This DFA does not have specified symbol probabilities\n')

    } else {

    #Print probability vector, transition probability matrix and stationary state probabilities (if present)
    cat('\n-----------------Symbol probabilities----------------- \n \n')
    print(object$probs)
    cat('\n------------Stationary state probabilities------------ \n \n')
    print(object$stationary.probs)
    cat('\n-------------Transition probability matrix------------ \n \n')
    print(object$transition.probs)
    cat('\n------------------------------------------------------ \n \n') } }

