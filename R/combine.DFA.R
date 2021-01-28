#' Combine Deterministic Finite Automata (DFAs)
#'
#' \code{combine.DFA} returns the minimised combined Deterministic Finite Automaton (DFA) for two input DFAs.
#'
#' This function takes two input DFAs and combines these into a single minimised DFA for the union language (i.e., the language for the
#' output is the union of the langauges in the input DFAs).  The input DFAs must be compatible with one another, which requires them to have
#' the same \code{alphabet}, \code{allow.overlap} and \code{probs} (if specified).  Names of these objects need not be the same, and in the
#' event of a difference the names from the first input DFA are used.
#'
#' @usage \code{combine.DFA()}
#' @param dfa1 A deterministic finite automaton (DFA) (see the \code{DFA} function)
#' @param dfa2 dfa1 A deterministic finite automaton (DFA) (see the \code{DFA} function)
#' @return A list of class 'DFA' containing the deterministic finite automaton (DFA)

combine.DFA <- function(dfa1, dfa2) {

  if (class(dfa1) != 'DFA')                                { stop('Error: This function only operates on objects of class \'DFA\'') }
  if (class(dfa2) != 'DFA')                                { stop('Error: This function only operates on objects of class \'DFA\'') }

  #Extract information from DFA1
  LANGUAGE.NAME1 <- attributes(dfa1$language)$language.name
  LANGUAGE1      <- dfa1$language
  STRING.NAMES1  <- names(LANGUAGE1)
  ALPHABET.NAME1 <- attributes(dfa1$alphabet)$alphabet.name
  ALPHABET1      <- dfa1$alphabet
  ALLOW.OVERLAP1 <- dfa1$allow.overlap
  OVERLAP1       <- dfa1$max.overlap
  PROBS.INCL1    <- FALSE
  if (!is.null(dfa1$probs)) {
    PROBS.INCL1  <- TRUE
    PROBS1       <- as.vector(dfa1$probs) }
  TRANSITIONS1   <- dfa1$transition.table
  PROGRESS1      <- dfa1$progress
  FINALSTATES1   <- dfa1$state.count

  #Extract information from DFA2
  LANGUAGE.NAME2 <- attributes(dfa2$language)$language.name
  LANGUAGE2      <- dfa2$language
  STRING.NAMES2  <- names(LANGUAGE2)
  ALPHABET.NAME2 <- attributes(dfa2$alphabet)$alphabet.name
  ALPHABET2      <- dfa2$alphabet
  ALLOW.OVERLAP2 <- dfa2$allow.overlap
  OVERLAP2       <- dfa2$max.overlap
  PROBS.INCL2    <- FALSE
  if (!is.null(dfa2$probs)) {
    PROBS.INCL2  <- TRUE
    PROBS2       <- as.vector(dfa2$probs) }
  TRANSITIONS2   <- dfa2$transition.table
  PROGRESS2      <- dfa2$progress
  FINALSTATES2   <- dfa2$state.count

  #Check inputs to ensure DFAs are compatible
  if (!identical(c(ALPHABET1), c(ALPHABET2)))              { stop('Error: DFAs can only be combined if they use the same alphabet') }
  if (!identical(ALPHABET.NAME1, ALPHABET.NAME2))          { warning('Inconsistent alphabet names --- we will use the first name') }
  if (!identical(ALLOW.OVERLAP1, ALLOW.OVERLAP2))          { stop('Error: DFAs can only be combined if they use the same allow.overlap') }
  if (PROBS.INCL1&PROBS.INCL2) {
  if (!identical(PROBS1, PROBS2))                          { stop('Error: DFAs can only be combined if they use the same alphabet probabilities') } }

  #Remove duplicated strings from the languages
  S1 <- length(LANGUAGE1)
  S2 <- length(LANGUAGE2)
  REMOVE <- rep(FALSE, S2)
  for (s2 in 1:S2) {
    s1 <- 1
    while ((s1 <= S1)&(!REMOVE[s2])) {
      if (identical(LANGUAGE1[[s1]], LANGUAGE2[[s2]])) { REMOVE[s2] <- TRUE }
      s1 <- s1 + 1 } }
  if (sum(REMOVE) == length(REMOVE)) { return(dfa1) }
  if (sum(REMOVE) > 0) {
    S2 <- length(LANGUAGE2)
    LANGUAGE2    <- LANGUAGE2[!REMOVE]
    OVERLAP2     <- OVERLAP2[!REMOVE]
    PROGRESS2    <- PROGRESS2[, !REMOVE]
    FINALSTATES2 <- FINALSTATES2[, !REMOVE] }

  #Create combined objects and set parameters
  m1 <- nrow(TRANSITIONS1)
  m2 <- nrow(TRANSITIONS2)
  SS <- S1 + S2
  ALPHABET <- ALPHABET1
  attr(ALPHABET, 'alphabet.name') <- ALPHABET.NAME1
  K  <- length(ALPHABET)
  LANGUAGE <- c(LANGUAGE1, LANGUAGE2)
  names(LANGUAGE) <- sprintf('String[%s]', 1:SS)
  attr(LANGUAGE, 'language.name') <- LANGUAGE.NAME1
  ALLOW.OVERLAP <- ALLOW.OVERLAP1
  if (!is.null(dfa1$probs)) { PROBS <- PROBS1 } else { if (!is.null(dfa2$probs)) { PROBS <- PROBS2 } }
  OVERLAP <- c(OVERLAP1, OVERLAP2)
  names(OVERLAP)  <- sprintf('String[%s]', 1:SS)
  PROBS.INCL <- (PROBS.INCL1|PROBS.INCL2)
  if (PROBS.INCL) { if (PROBS.INCL1) { PROBS <- PROBS1 } else { PROBS <- PROBS2 } }

  #Generate the combined transition matrix, progress matrix and final-state matrix
  mm <- m1*m2-1
  TT <- matrix(0, nrow = mm+1, ncol = K)
  colnames(TT) <- ALPHABET
  PP <- matrix(0, nrow = mm+1, ncol = SS)
  colnames(PP) <- sprintf('P[%s]', 1:SS)
  FF <- matrix(0, nrow = mm+1, ncol = SS)
  colnames(FF) <- sprintf('F[%s]', 1:SS)

  #Create state-mapping function
  statemap     <- function(i, j) { m2*i + j }

  #Set values in combined states
  for (i in 0:(m1-1)) {
  for (j in 0:(m2-1)) {
    STATE <- statemap(i, j)
    for (x in 1:K)  { TT[STATE+1, x] <- statemap(TRANSITIONS1[i+1, x], TRANSITIONS2[j+1, x]) }
    for (s in 1:S1) {
      PP[STATE+1, s]    <- PROGRESS1[i+1,    s]
      FF[STATE+1, s]    <- FINALSTATES1[i+1, s] }
    for (s in 1:S2) {
      PP[STATE+1, S1+s] <- PROGRESS2[j+1,    s]
      FF[STATE+1, S1+s] <- FINALSTATES2[j+1, s] } } }

  #Remove impossible states via iterative method
  KEEP  <- c(TRUE, rep(FALSE, mm))
  BREAK <- FALSE
  while (!BREAK) {
    KEEP.OLD <- KEEP
    for (i in (0:mm)[KEEP.OLD]) {
    for (x in 1:K) {
      KEEP[TT[i+1, x]+1] <- TRUE } }
    if (identical(KEEP.OLD, KEEP)) { BREAK <- TRUE } }
  mm.new <- sum(KEEP)-1
  TT.new <- TT[KEEP, ]
  PP.new <- PP[KEEP, ]
  FF.new <- FF[KEEP, ]
  rownames(TT.new) <- sprintf('State[%s]', 0:mm.new)
  rownames(PP.new) <- sprintf('State[%s]', 0:mm.new)
  rownames(FF.new) <- sprintf('State[%s]', 0:mm.new)
  INDEX <- (0:mm)[KEEP]
  for (i in 0:mm.new) {
  for (x in 1:K) {
    OLDSTATE <- unname(TT.new[i+1, x])
    NEWSTATE <- which(INDEX == OLDSTATE)-1
    TT.new[i+1,x] <- NEWSTATE } }

  #Set combined DFA
  DFA.COMBINE <- list(transition.table = TT.new, progress = PP.new, state.count = FF.new, max.overlap = OVERLAP,
                      language = LANGUAGE, alphabet = ALPHABET, allow.overlap = ALLOW.OVERLAP)
  class(DFA.COMBINE) <- 'DFA'

  #Generate the transition probability matrix and stationary vectors
  if (PROBS.INCL) {

    #Add probability vector
    DFA.COMBINE$probs <- matrix(PROBS, nrow = 1)
    rownames(DFA.COMBINE$probs) <- 'Probability'
    colnames(DFA.COMBINE$probs) <- ALPHABET

    #Transition probability matrix
    H <- matrix(0, nrow = mm.new+1, ncol = mm.new+1)
    rownames(H) <- sprintf('State[%s]', 0:mm.new)
    colnames(H) <- sprintf('State[%s]', 0:mm.new)
    for (i in 0:mm.new) {
    for (x in 1:K) {
      H[i+1, TT.new[i+1, x]+1] <- H[i+1, TT.new[i+1, x]+1] + PROBS[x] } }

    #Compute the eigendecomposition
    EIGEN <- eigen(t(H))
    VALS  <- EIGEN$values
    VECS  <- EIGEN$vectors

    #Stationary vector(s)
    RES <- zapsmall(Mod(VALS - as.complex(rep(1, mm.new+1))))
    IND <- which(RES == 0)
    SSS <- Re(VECS[, IND])
    STR <- matrix(SSS/sum(SSS), nrow = 1)
    rownames(STR) <- 'Stationary'
    colnames(STR) <- rownames(H)

    #Add transition matrix and stationary vector to DFA
    DFA.COMBINE$stationary.probs <- STR
    DFA.COMBINE$transition.probs <- H }

  #Give the output
  DFA.COMBINE }
