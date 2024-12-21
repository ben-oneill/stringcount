#' Minimise Deterministic Finite Automata (DFAs)
#'
#' \code{minimise.DFA/minimize.DFA} returns the minimised Deterministic Finite Automaton (DFA).
#'
#' This function takes an input DFA and minimises this by amalgamating any equivalent states.  Minimisation is done using the Myhill-Nerode
#' table method.  The output is a minimised version of the input DFA.
#'
#' @usage \code{minimise.DFA()/minimize.DFA()}
#' @param dfa A deterministic finite automaton (DFA) (see the \code{DFA} function)
#' @return A list of class 'DFA' containing the minimised deterministic finite automaton (DFA)

minimise.DFA <- minimize.DFA <- function(dfa) {

  if (class(dfa) != 'DFA')                                 { stop('Error: This function only operates on objects of class \'DFA\'') }

  #Extract information from DFA
  DFA <- dfa
  ALPHABET <- DFA$alphabet
  TT <- DFA$transition.table
  PP <- DFA$progress
  FF <- DFA$state.count
  m  <- nrow(TT)-1
  K  <- length(ALPHABET)
  PROBS.INCL    <- FALSE
  if (!is.null(DFA$probs)) {
    PROBS.INCL <- TRUE
    PROBS      <- as.vector(DFA$probs) }

  if (m > 0) {

    #-----------------------------------------------------------------------------------------------------------
    #Minimise the DFA
    #This step uses minimisation using the Myhill-Nerode table method
    #-----------------------------------------------------------------------------------------------------------

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
      DFA$state.count      <- FF

      #-----------------------------------------------------------------------------------------------------------
      #Update transition probabilities (if required)
      #-----------------------------------------------------------------------------------------------------------

      #Generate the transition probability matrix and stationary vectors
      if (PROBS.INCL) {

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
        DFA$stationary.probs <- STR
        DFA$transition.probs <- H } } }

  #Update minimisation parameter
  DFA$minimise <- TRUE

  #Give the output
  DFA }

