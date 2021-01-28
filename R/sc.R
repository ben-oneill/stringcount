#' Split/concatenate between a single character string and a vector of its individual symbols
#'
#' \code{sc} can be used to split a single character string into a vector of its individual symbols, or concatenate a character
#' vector with multiple elements into a single character string.  The user may set the logical value \code{split} to specify if they wish to
#' split or concatenate the \code{string} input.  If this value is unspecified then the function will determine whether to split or concatenate
#' by looking at the length of the \code{string} vector; if this vector has a single element it will be split; if it has multiple elements it
#' will be concatenated.  (Note: if the user concatenates a string vector with one or more elements that have multiple characters, they will
#' get a warning to alert the to the fact that concatenation loses the distinction between individual symbols and longer string of symbols.)
#' If the user split the \code{string} then the function will return a character vector showing the individual symbols in the input string; if
#' the user concatenated the \code{string} then the function will return a single character value with the concatenated string.
#'
#' @usage \code{sc()}
#' @param string A character value/vector
#' @param split Logical; if \code{TRUE} the function will split the string; if \code{FALSE} the function will concatenate the string
#' @param warning Logical; if \code{TRUE} the function will give a warning for concatenation of strings with elements that are not single characters
#' @return A character vector giving the split/concatenated string.

sc <- function(string, split = NULL, warning = TRUE) {

  #Check input
  if (!is.character(string))                                   { stop('Error: string must be a character vector') }
  if (length(string) == 0)                                     { stop('Error: string must have at least one element') }
  if (!missing(split)) {
  if (!is.vector(split))                                       { stop('Error: split should be a single logical value (if specified)') }
  if (!is.logical(split))                                      { stop('Error: split should be a single logical value (if specified)') }
  if (length(split) != 1)                                      { stop('Error: split should be a single logical value (if specified)') } }
  if (!missing(warning)) {
  if (!is.vector(warning))                                     { stop('Error: warning should be a single logical value (if specified)') }
  if (!is.logical(warning))                                    { stop('Error: warning should be a single logical value (if specified)') }
  if (length(warning) != 1)                                    { stop('Error: warning should be a single logical value (if specified)') } }

  #Determine transformation type
  if (missing(split)) { SPLIT <- (length(string) == 1) }

  #Split the input
  if (SPLIT) {
    if (length(string) > 1) { warning('You are asking to split a string that is already a vector\n  Only the first element will be used') }
    M <- nchar(string[1], type = 'char')
    OUT <- character(M)
    for (i in 1:M) { OUT[i] <- substr(string[1], start = i, stop = i) } }

  #Concatenate the input
  if (!SPLIT) {
    if (warning) {
    if (!all(nchar(string) == 1)) {
      WARNING <- paste0('At least one element of string was not a single character symbol\n',
                        '  Concatenation loses the distinction between individual symbols and longer strings of symbols\n',
                        '  If you try to reverse the concatenation you will not get back your original vector')
      warning(WARNING) } }
    OUT <- paste(string, collapse = '') }

  #Give the output
  OUT }
