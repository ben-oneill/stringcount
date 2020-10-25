#' Matrices for the state variable in the String Count Distribution
#'
#' \code{stringmatrix} returns the transition table and transition probability matrix for the state variable in the String Count Distribution
#'
#' The hidden state variable used in the String Count Distribution is the number of characters in the string that are presently matched
#' after any given number of characters in the text.  The present function computes the transition table for the state variable and the
#' transition probability matrix for the state variable.  The first of these matrices depends on the string vector ```string``` and the
#' second depends on the string vector and the probability vector ```probs``` for the underlying symbols in the alphabet.  The function
#' allows the user to specify the ```alphabet``` for the analysis, with the default alphabet being the natural numbers up to the length
#' of the probability vector.  (Note: The user can give either a numeric vector or a character vector for the ```string``` and ```alphabet```,
#' but the elements in the string must be in the alphabet, and both vectors must be the same type.)  The state variable is affected by whether
#' or not overlapping occurrances of the string are both included in the count.  If ```allow.overlap``` is set to ```TRUE``` then the string-
#' count includes all occurrances of overlapping strings (i.e., they are both counted even if they share some symbols), and this is reflected
#' in the transition table and transition probability matrix.  Contrarily, if ```allow.overlap``` is set to ```FALSE``` then the string-count
#' excludes any occurrances that share symbols with a previously counted occurrance, and this is again reflected in the transition table
#' and transition probability matrix.  The output also includes the maximum overlap of the string and the starting and counting states for the
#' state variable.  (Note that the transition table, together with the starting and counting states, define the Deterministic Finite Automata
#' (DFA) for the string.)
#'
#' @usage \code{stringmatrix()}
#' @param string A numeric/character vector
#' @param probs A vector of the symbol probabilities (taken over the symbols in the ```alphabet```)
#' @param alphabet A numeric/character vector containing the alphabet for the analysis
#' @param allow.overlap Logical; if ```TRUE``` then string occurrances are counted even if they overlap with previously counted occurrances
#' @return A list containing the maximum overlap, starting and counting states, transition table and transition probability matrix for the state variable

stringmatrix <- function(string, probs, alphabet = NULL, allow.overlap = TRUE) {

  #Check inputs probs and allow.overlap
  if (!is.vector(probs))                                       { stop('Error: probs must be a probability vector') }
  if (!is.numeric(probs))                                      { stop('Error: probs must be a probability vector') }
  if (length(probs) == 0)                                      { stop('Error: probs must have at least one value') }
  if (min(probs) < 0)                                          { stop('Error: probs must be a probability vector') }
  if (sum(probs) != 1)                                         { stop('Error: probs must be a probability vector') }
  if (!is.vector(allow.overlap))                               { stop('Error: allow.overlap must be a single logical value') }
  if (!is.logical(allow.overlap))                              { stop('Error: allow.overlap must be a single logical value') }
  if (length(allow.overlap) != 1)                              { stop('Error: allow.overlap must be a single logical value') }

  #Check input alphabet
  if (is.null(alphabet)) {
    alphabet <- 1:length(probs) }
  if (!is.vector(alphabet))                                    { stop('Error: alphabet must be a vector') }
  K <- length(unique(alphabet))
  if (length(alphabet) != K) {
    warning('Input alphabet contained duplicate elements')
    alphabet <- unique(alphabet) }
  if (length(probs) != K)                                      { stop('Error: probs must have the same length as alphabet') }

  #Set alphabet type
  TYPE <- NULL
  if (is.numeric(alphabet))   { TYPE <- 'numeric' } else {
  if (is.character(alphabet)) { TYPE <- 'character' } }
  if (is.null(TYPE))                                           { stop('Error: alphabet should be a numeric or character vector') }

  #Check input string
  if (!is.vector(string))                                      { stop('Error: string must be a vector') }
  if ((TYPE == 'numeric')&&(!is.numeric(string)))              { stop('Error: string must be the same type as alphabet') }
  if ((TYPE == 'character')&&(!is.character(string)))          { stop('Error: string must be the same type as alphabet') }
  m <- length(string)
  for (i in 1:m) {
    if (!(string[i] %in% alphabet))                            { stop(paste0('Error: string element ', i, ' is not in the alphabet')) } }

  #Match string characters to alphabet
  MATCH <- match(string, alphabet)

  #Generate the structure matrix
  G <- matrix(0, nrow = m+1, ncol = K)
  rownames(G) <- sprintf('State[%s]', 0:m)
  colnames(G) <- alphabet
  for (i in 0:m) {
  for (x in 1:K) {
    g     <- i+1
    BREAK <- FALSE
    while (!BREAK) {
      if (g > 1) { VEC1 <- c(string[(i-g+2):i], alphabet[x]) } else { VEC1 <- alphabet[x] }
      VEC2 <- string[1:g]
      if (identical(VEC1, VEC2)) { BREAK <- TRUE } else { g <- g-1 }
      if (g == 0) { BREAK <- TRUE } }
    G[i+1, x] <- g } }

  #Adjust last row of structure matrix is overlap is disallowed
  if (!allow.overlap) {
    G[m+1, ] <- 0
    G[m+1, MATCH[1]] <- 1 }

  #Set initial and final states
  INITIAL <- 0
  FINAL   <- m

  #Compute overlap value
  OVERLAP <- max(G[m+1,])-1

  #Generate the string advancement matrix
  H <- matrix(0, nrow = m+1, ncol = m+1)
  rownames(H) <- sprintf('State[%s]', 0:m)
  colnames(H) <- sprintf('State[%s]', 0:m)
  for (h in 1:m) { H[h, h+1] <- probs[MATCH[h]] }
  H[1,1] <- 1 - probs[MATCH[1]]
  for (i in 0:m) {
  for (h in 0:i) {
    IND <- (G[i+1, ] == h)
    H[i+1, h+1] <- sum(probs*IND) } }

  #Output the matrix
  list(max.overlap = OVERLAP, state.start = INITIAL, state.count = FINAL, transition.table = G, transition.probs = H) }
