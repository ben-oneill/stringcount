#' Counts the number of occurrences of a string in text
#'
#' \code{stringcount} returns the number of occurrances of a string in a text
#'
#' This function takes a ```text``` that is either a single character value or a vector of values, and a ```string``` that is either a single
#' character value or a vector of values.  It counts the number of occurrances of the string in the text and returns this count.  By default
#' the function counts overlapping strings as separate occurrances (i.e., they both add to the count), however the user can set ```allow.overlap```
#' to ```FALSE``` to disallow occurrances of the string that overlap with previously counted occurrances (i.e., the latter occurrance does not
#' add to the count).  If the inputs are single character strings then the function checks for the occurrance of the string as a substring of the
#' text.  If the inputs are longer vectors of numeric or character values then the elements of the vectors are treated as individual symbols and
#' the vector gives the string of symbols that is checked in the count.
#'
#' @usage \code{stringcount()}
#' @param text Either a single character string or a numeric/character vector
#' @param string Either a single character string or a numeric/character vector (must have the same type as the ```text```)
#' @param allow.overlap Logical; if ```TRUE``` then string occurrances are counted even if they overlap with previously counted occurrances
#' @return The count value for the number of times the string occurs in the text

stringcount <- function(text, string, allow.overlap = TRUE) {

  #Check input text
  if (!is.vector(text))                                        { stop('Error: text must be a vector') }
  if ((!is.numeric(text))&&(!is.character(text)))              { stop('Error: text must be a numeric or character vector') }
  if (is.numeric(text)) { TYPE <- 'numeric' } else { TYPE <- 'character' }
  n <- length(text)

  #Check input string
  if (!is.vector(string))                                      { stop('Error: string must be a vector') }
  if ((TYPE == 'numeric')&&(!is.numeric(string)))              { stop('Error: string must be the same type as text') }
  if ((TYPE == 'character')&&(!is.character(string)))          { stop('Error: string must be the same type as text') }
  m <- length(string)

  #Check lengths
  if (n == 0)                                                  { stop('Error: test must have at least one symbol') }
  if (m == 0)                                                  { stop('Error: string must have at least one symbol') }

  #Compute string count
  COUNT <- 0
  i <- 0
  while(i <= n-m) {
    if (identical(string, text[(i+1):(i+m)])) {
      COUNT <- COUNT + 1
      if (!allow.overlap) { i <- i+m-1 } }
    i <- i+1 }

  #Output count
  COUNT }
