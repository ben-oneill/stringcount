#' Determines the reachable states in a DFA
#'
#' \code{reachable.states} returns an array of indicators for reachable states in a DFA
#'
#' This function takes an input DFA produced with the \code{DFA} function and a stipulation for the maximum number of symbols of
#' interest.  The function produces a three-dimensional array of indicators showing which states are reachable from other states
#' with a given number of symbols.  For each possible number of symbols from zero up to \code{max.symbol} the array gives a two-
#' dimensional matrix showing the indicators for whether an end-state can be reached from a start-state with exactly that number
#' of symbols.
#'
#' @usage \code{reachable.states()}
#' @param dfa A deterministic finite automaton (DFA) produced by the \code{DFA} function
#' @param max.symbols The maximum number of symbols of interest
#' @return An array of indicators showing which states are reachable from other states with a given number of symbols.

reachable.states <- function(dfa, max.symbols) {

  #Check inputs
  if (class(dfa) != 'DFA')                       { stop('Error: This function only operates on objects of class \'DFA\'') }
  if (!is.numeric(max.symbols))                  { stop('Error: Input max.symbols should be a positive integer') }
  if (length(max.symbols) != 1)                  { stop('Error: Input max.symbols should be a single positive integer') }
  if (as.integer(max.symbols) != max.symbols)    { stop('Error: Input max.symbols should be a positive integer') }
  if (max.symbols < 1)                           { stop('Error: Input max.symbols should be a positive integer') }

  #Extract information
  TRANSITIONS <- dfa$transition.table
  K <- ncol(TRANSITIONS)
  m <- nrow(TRANSITIONS)-1

  #Create output
  OUT <- array(FALSE, dim = c(m+1, m+1, max.symbols+1),
               dimnames = list(sprintf('StartState[%s]', 0:m), sprintf('ReachState[%s]', 0:m), sprintf('n[%s]', 0:max.symbols)))

  #Compute reachable states
  for (i in 0:m) { OUT[i+1, i+1, 1] <- TRUE }
  for (k in 1:max.symbols) {
    for (q in 0:m) {
      STATES <- unname(which(OUT[q+1, , k])) - 1
      for (r in STATES) {
        for (x in 1:K) {
          TT <- TRANSITIONS[r+1, x]
          OUT[q+1, TT+1, k+1] <- TRUE } } } }

  #Return output
  OUT }
