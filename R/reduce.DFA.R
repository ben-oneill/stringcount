#' Reduce a Deterministic Finite Automaton (DFA)
#'
#' \code{reduce.DFA} returns the reduced Deterministic Finite Automaton (DFA).
#'
#' This function takes an input DFA and reduces it by removing sublanguages.  The user must input a DFA and an input \code{remove} that specifies
#' a vector ofublanguages to remove from the DFA.  The indices in the \code{remove} input must be integers referring to elements of the list
#' of sublanguages in the DFA.  The output is a reduced version of the input DFA where these terms have been removed from the final-state table.
#'
#' @usage \code{reduce.DFA()}
#' @param dfa A deterministic finite automaton (DFA) (see the \code{DFA} function)
#' @return A list of class 'DFA' containing the reduced deterministic finite automaton (DFA)

reduce.DFA <- function(dfa, remove, minimise = FALSE, minimize = NULL) {

  if (class(dfa) != 'DFA')                                     { stop('Error: This function only operates on objects of class \'DFA\'') }

  #Extract information from DFA
  DFA <- dfa
  DFA.REDUCE <- dfa
  ALPHABET <- DFA$alphabet
  SUBLANGUAGES <- DFA$sublanguages
  FF <- DFA$state.count
  K  <- length(ALPHABET)
  D  <- length(SUBLANGUAGES)

  #Check input remove
  if (!is.vector(remove))                                      { stop('Error: Input remove should be a vector') }
  if ((!is.numeric(remove))&(!is.logical(remove)))             { stop('Error: Input remove should be a logical or numeric vector') }
  if (is.numeric(remove)) {
    RR <- length(remove)
    for (r in 1:RR) {
      if (!(remove[r] %in% 1:D)) {
        stop(paste0('Error: Element ', r, ' of input remove is not a valid index for the language')) } }
    REMOVE <- logical(D)
    for (i in 1:RR) { REMOVE[remove[i]] <- TRUE } }
  if (is.logical(remove)) {
    if (length(remove) != D)                                   { stop('Error: Input remove is the wrong length') } }

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

  #Reduce the sublangauges and final-state table
  DD <- sum(!REMOVE)
  SUBLANGUAGES.new <- SUBLANGUAGES[!REMOVE, drop = FALSE]
  names(SUBLANGUAGES.new) <- sprintf('Sub[%s]', 1:DD)
  FF.new <- FF[, !REMOVE, drop = FALSE]
  colnames(FF.new) <- sprintf('F[%s]', 1:DD)

  #Reduce the DFA
  FF <- FF.new
  SUBLANGUAGES <- SUBLANGUAGES.new
  DFA.REDUCE$state.count  <- FF
  DFA.REDUCE$sublanguages <- SUBLANGUAGES

  #-----------------------------------------------------------------------------------------------------------
  #Minimise the DFA
  #This step uses minimisation using the Myhill-Nerode table method
  #-----------------------------------------------------------------------------------------------------------

  if (MINIMISE) {

    #Get parameters
    TT <- DFA$transition.table
    PP <- DFA$progress
    m  <- nrow(TT)-1
    PROBS.INCL    <- FALSE
    if (!is.null(DFA$probs)) {
      PROBS.INCL  <- TRUE
      PROBS       <- as.vector(DFA$probs) }

    if (m > 0)  {

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
        DFA.REDUCE$transition.table <- TT
        DFA.REDUCE$progress         <- PP
        DFA.REDUCE$state.count      <- FF }

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
        DFA.REDUCE$stationary.probs <- STR
        DFA.REDUCE$transition.probs <- H } }

    DFA.REDUCE$minimise <- TRUE }

  #Give the output
  DFA.REDUCE }

