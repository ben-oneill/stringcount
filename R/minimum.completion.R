#' Determines the minimum-completion values for strings/sublanguages in a DFA
#'
#' \code{minimum.completion} returns a matrix of minimum-completion values for all the sublanguages in a DFA
#'
#' This function takes an input DFA produced with the \code{DFA} function and computes a matrix showing the minimum-completion
#' value for each sublanguage in the DFA over each possible starting state.  The minimum-completion value for a sublanguage is
#' the minimum number of additional symbols in the text that are required to complete the sublanguage when it starts in the
#' specified starting state.
#'
#' @usage \code{minimum.completion()}
#' @param dfa A deterministic finite automaton (DFA) produced by the \code{DFA} function
#' @return A matrix of minimum-completion values for all the sublanguages in a DFA over all possible starting states

minimum.completion <- function(dfa) {

  #Check inputs
  if (class(dfa) != 'DFA')                       { stop('Error: This function only operates on objects of class \'DFA\'') }

  #Extract information
  SUBLANGUAGE <- dfa$sublanguages
  TRANSITIONS <- dfa$transition.table
  FINAL       <- dfa$state.count
  D <- length(SUBLANGUAGE)
  K <- ncol(TRANSITIONS)
  m <- nrow(TRANSITIONS)-1

  #Create output
  OUT <- matrix(NA, nrow = m+1, ncol = D)
  rownames(OUT) <- sprintf('StartState[%s]', 0:m)
  colnames(OUT) <- sprintf('MC[%s]', 1:D)

  #Compute minimum completion
  COMP  <- 0
  REACH <- matrix(FALSE, nrow = m+1, ncol = m+1)
  for (i in 0:m) { REACH[i+1, i+1] <- TRUE }
  while (sum(is.na(OUT)) > 0) {

    #Increase completion value
    COMP <- COMP+1

    #Compute reachable states
    OLDREACH <- REACH
    REACH <- matrix(FALSE, nrow = m+1, ncol = m+1)
    for (q in 0:m) {
      STATES <- unname(which(OLDREACH[q+1, ])) - 1
      for (r in STATES) {
        for (x in 1:K) {
          TT <- TRANSITIONS[r+1, x]
          REACH[q+1, TT+1] <- TRUE } } }

    #Check completion for each state/sublanguage
    for (q in 0:m) {
    for (d in 1:D) {
      if (is.na(OUT[q+1, d])) {
        if (sum(REACH[q+1, ]*FINAL[, d] > 0)) { OUT[q+1, d] <- COMP } } } } }

  #Give output
  OUT }

