#' Combine Deterministic Finite Automata (DFAs)
#'
#' \code{combine.DFA} returns the minimised combined Deterministic Finite Automaton (DFA) for two input DFAs.
#'
#' This function takes two input DFAs and combines these into a single minimised DFA for the union language (i.e., the language for the
#' output is the union of the langauges in the input DFAs).  The input DFAs must be compatible with one another, which requires them to have
#' the same\code{allow.overlap} and \code{probs} (if specified).
#'
#' @usage \code{combine.DFA()}
#' @param dfa1 A deterministic finite automaton (DFA) (see the \code{DFA} function)
#' @param dfa2 A deterministic finite automaton (DFA) (see the \code{DFA} function)
#' @param minimise/minimize Logical; if \code{TRUE} then the DFA is minimised
#' @return A list of class 'DFA' containing the deterministic finite automaton (DFA)

combine.DFA <- function(dfa1, dfa2, minimise = FALSE, minimize = NULL) {

  if (class(dfa1) != 'DFA')                                { stop('Error: This function only operates on objects of class \'DFA\'') }
  if (class(dfa2) != 'DFA')                                { stop('Error: This function only operates on objects of class \'DFA\'') }

  #Extract information from DFA1
  LANGUAGE1      <- dfa1$language
  STRING.NAMES1  <- names(LANGUAGE1)
  ALPHABET1      <- dfa1$alphabet
  ALLOW.OVERLAP1 <- dfa1$allow.overlap
  MC1            <- attributes(dfa1$sublanguage)$minimum.completion
  PROBS.INCL1    <- FALSE
  if (!is.null(dfa1$probs)) {
    PROBS.INCL1  <- TRUE
    PROBS1       <- as.vector(dfa1$probs) }
  TRANSITIONS1   <- dfa1$transition.table
  PROGRESS1      <- dfa1$progress
  FINALSTATES1   <- dfa1$state.count
  SUBLANGUAGES1  <- dfa1$sublanguages

  #Extract information from DFA2
  LANGUAGE2      <- dfa2$language
  STRING.NAMES2  <- names(LANGUAGE2)
  ALPHABET2      <- dfa2$alphabet
  ALLOW.OVERLAP2 <- dfa2$allow.overlap
  MC2            <- attributes(dfa2$sublanguage)$minimum.completion
  PROBS.INCL2    <- FALSE
  if (!is.null(dfa2$probs)) {
    PROBS.INCL2  <- TRUE
    PROBS2       <- as.vector(dfa2$probs) }
  TRANSITIONS2   <- dfa2$transition.table
  PROGRESS2      <- dfa2$progress
  FINALSTATES2   <- dfa2$state.count
  SUBLANGUAGES2  <- dfa2$sublanguages

  #Check inputs to ensure DFAs are compatible
  if (!identical(ALLOW.OVERLAP1, ALLOW.OVERLAP2))          { stop('Error: DFAs can only be combined if they use the same allow.overlap') }
  if (PROBS.INCL1&PROBS.INCL2) {
    if (!identical(c(ALPHABET1), c(ALPHABET2)))            { stop('Error: If alphabet probabilities are specified, DFAs must use the same alphabet') }
    if (!identical(PROBS1, PROBS2))                        { stop('Error: If alphabet probabilities are specified, DFAs must use the same alphabet probabilities') } }

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

  #-----------------------------------------------------------------------------------------------------------
  #Combine the input DFAs and remove impossible states
  #-----------------------------------------------------------------------------------------------------------

  #Set parameters for combined DFA
  ALPHABET <- sort(unique(ALPHABET1, ALPHABET2))
  K  <- length(ALPHABET)
  ALLOW.OVERLAP <- ALLOW.OVERLAP1
  if (!is.null(dfa1$probs)) { PROBS <- PROBS1 } else { if (!is.null(dfa2$probs)) { PROBS <- PROBS2 } }
  PROBS.INCL <- (PROBS.INCL1|PROBS.INCL2)
  if (PROBS.INCL) { if (PROBS.INCL1) { PROBS <- PROBS1 } else { PROBS <- PROBS2 } }

  #Remove duplicated strings from the languages and merge them
  #This requires changes to the sublanguages to keep track of
  S1 <- length(LANGUAGE1)
  S2 <- length(LANGUAGE2)
  D1 <- length(SUBLANGUAGES1)
  D2 <- length(SUBLANGUAGES2)
  for (d2 in 1:D2) { SUBLANGUAGES2[[d2]] <- SUBLANGUAGES2[[d2]] + S1 }
  REPLACE <- rep(NA, S2)
  for (s2 in 1:S2) {
    s1 <- 1
    while ((s1 <= S1)&(is.na(REPLACE[s2]))) {
      if (identical(LANGUAGE1[[s1]], LANGUAGE2[[s2]])) { REPLACE[s2] <- s1 }
      s1 <- s1 + 1 } }
  if (sum(!is.na(REPLACE)) > 0) {
    for (ss in 1:S2) {
      if (!is.na(REPLACE[ss])) {
      for (d2 in 1:D2) {
        SUB <- SUBLANGUAGES2[[d2]]
        SUB[SUB == ss] <- REPLACE[ss]
        SUBLANGUAGES2[[d2]] <- SUB } } }
    LANGUAGE2 <- LANGUAGE2[is.na(REPLACE), drop = FALSE]
    S2 <- length(LANGUAGE2) }
  SS <- S1 + S2
  LANGUAGE <- c(LANGUAGE1, LANGUAGE2)
  names(LANGUAGE) <- sprintf('String[%s]', 1:SS)

  #Remove duplicated sublanguages
  REMOVE <- rep(FALSE, D2)
  for (d2 in 1:D2) {
    d1 <- 1
    while ((d1 <= D1)&(!REMOVE[d2])) {
      if (identical(SUBLANGUAGES1[[d1]], SUBLANGUAGES2[[d2]])) { REMOVE[d2] <- TRUE }
      d1 <- d1 + 1 } }
  if (sum(REMOVE) > 0) {
    SUBLANGUAGES2 <- SUBLANGUAGES2[!REMOVE, drop = FALSE]
    MC2           <- MC2[!REMOVE]
    FINALSTATES2  <- FINALSTATES2[, !REMOVE, drop = FALSE]
    D2            <- length(SUBLANGUAGES2) }
  DD <- D1 + D2
  SUBLANGUAGES <- c(SUBLANGUAGES1, SUBLANGUAGES2)
  names(SUBLANGUAGES) <- sprintf('Sub[%s]', 1:DD)
  attr(SUBLANGUAGES, 'minimum.completion') <- c(MC1, MC2)

  #Generate the combined transition matrix, progress matrix and final-state matrix
  m1 <- nrow(TRANSITIONS1)-1
  m2 <- nrow(TRANSITIONS2)-1
  mm <- (m1+1)*(m2+1)-1
  TT <- matrix(0, nrow = mm+1, ncol = K)
  colnames(TT) <- ALPHABET
  PP <- matrix(0, nrow = mm+1, ncol = SS)
  colnames(PP) <- sprintf('P[%s]', 1:SS)
  FF <- matrix(0, nrow = mm+1, ncol = DD)
  colnames(FF) <- sprintf('F[%s]', 1:DD)

  #Create state-mapping function
  statemap     <- function(i, j) { (m2+1)*i + j }

  #Set values in combined states
  for (i in 0:m1) {
  for (j in 0:m2) {
    STATE <- statemap(i, j)
    for (x in 1:K)  { TT[STATE+1, x] <- statemap(TRANSITIONS1[i+1, x], TRANSITIONS2[j+1, x]) }
    for (s in 1:S1) { PP[STATE+1, s]    <- PROGRESS1[i+1, s] }
    for (s in 1:S2) { PP[STATE+1, S1+s] <- PROGRESS2[j+1, s] }
    for (d in 1:D1) { FF[STATE+1, d]    <- FINALSTATES1[i+1, d] }
    for (d in 1:D2) { FF[STATE+1, D1+d] <- FINALSTATES2[j+1, d] } } }

  #Remove impossible states via iterative method
  KEEP  <- c(TRUE, rep(FALSE, mm))
  BREAK <- FALSE
  while (!BREAK) {
    KEEP.OLD <- KEEP
    for (i in (0:mm)[KEEP.OLD]) {
    for (x in 1:K) {
      KEEP[TT[i+1, x]+1] <- TRUE } }
    if (identical(KEEP.OLD, KEEP)) { BREAK <- TRUE } }
  m.new  <- sum(KEEP)-1
  TT.new <- TT[KEEP, ]
  PP.new <- PP[KEEP, ]
  FF.new <- FF[KEEP, ]
  rownames(TT.new) <- sprintf('State[%s]', 0:m.new)
  rownames(PP.new) <- sprintf('State[%s]', 0:m.new)
  rownames(FF.new) <- sprintf('State[%s]', 0:m.new)
  INDEX <- (0:mm)[KEEP]
  for (i in 0:m.new) {
  for (x in 1:K) {
    OLDSTATE <- unname(TT.new[i+1, x])
    NEWSTATE <- which(INDEX == OLDSTATE)-1
    TT.new[i+1,x] <- NEWSTATE } }
  m  <- m.new
  TT <- TT.new
  PP <- PP.new
  FF <- FF.new

  #Set combined DFA
  DFA.COMBINE <- list(transition.table = TT, progress = PP, state.count = FF,
                      language = LANGUAGE, sublanguages = SUBLANGUAGES, alphabet = ALPHABET,
                      allow.overlap = ALLOW.OVERLAP, minimise = FALSE)
  class(DFA.COMBINE) <- 'DFA'

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

      DFA.COMBINE$transition.table <- TT
      DFA.COMBINE$progress         <- PP
      DFA.COMBINE$state.count      <- FF } }

    DFA.COMBINE$minimise <- TRUE }

  #-----------------------------------------------------------------------------------------------------------
  #Add alphabet probabilities and transition probabilities (if required)
  #-----------------------------------------------------------------------------------------------------------

  #Generate the transition probability matrix and stationary vectors
  if (PROBS.INCL) {

    #Add probability vector
    DFA.COMBINE$probs <- matrix(PROBS, nrow = 1)
    rownames(DFA.COMBINE$probs) <- 'Probability'
    colnames(DFA.COMBINE$probs) <- ALPHABET

    #Transition probability matrix
    H <- matrix(0, nrow = m+1, ncol = m+1)
    rownames(H) <- sprintf('State[%s]', 0:m)
    colnames(H) <- sprintf('State[%s]', 0:m)
    for (i in 0:m) {
    for (x in 1:K) {
      H[i+1, TT[i+1, x]+1] <- H[i+1, TT[i+1, x]+1] + PROBS[x] } }

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
    DFA.COMBINE$stationary.probs <- STR
    DFA.COMBINE$transition.probs <- H }

  #Give the output
  DFA.COMBINE }

