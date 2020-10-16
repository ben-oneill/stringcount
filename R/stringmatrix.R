#' Matrices for the state variable in the String Count Distribution
#'
#' \code{stringmatrix} returns the structure and transition matrices for the state variable in the String Count Distribution
#'
#' The hidden state variable used in the String Count Distribution is the number of characters in the string that are presently matched
#' after any given number of characters in the text.  The present function computes the structure matrix for the state variable and the
#' transition probability matrix for the state variable.  The first of these matrices depends on the string vector ```string``` and the
#' second depends on the string vector and the probability vector ```probs``` for the underlying symbols in the alphabet.  The function
#' allows the user to specify the ```alphabet``` for the analysis, with the default alphabet being the natural numbers up to the length
#' of the probability vector.  (Note: The user can give either a numeric vector or a character vector for the ```string``` and ```alphabet```,
#' but the elements in the string must be in the alphabet, and both vectors must be the same type.)
#'
#' @usage \code{stringmatrix()}
#' @param string A numeric/character vector
#' @param probs A vector of the symbol probabilities (taken over the symbols in the ```alphabet```)
#' @param alphabet A numeric/character vector containing the alphabet for the analysis
#' @return A list containing the structure matrix and transition probability matrix for the state variable

stringmatrix <- function(string, probs, alphabet = NULL) {

  #Check input probs
  if (!is.vector(probs))                                       { stop('Error: probs must be a probability vector') }
  if (!is.numeric(probs))                                      { stop('Error: probs must be a probability vector') }
  if (length(probs) == 0)                                      { stop('Error: probs must have at least one value') }
  if (min(probs) < 0)                                          { stop('Error: probs must be a probability vector') }
  if (sum(probs) != 1)                                         { stop('Error: probs must be a probability vector') }

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
  list(structure = G, transition = H) }
