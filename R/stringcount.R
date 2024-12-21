#' Counts the number of occurrences of strings in text
#'
#' \code{stringcount} returns the number of occurrances of a given set of strings in a text
#'
#' This function takes a \code{text} that is either a single character value or a vector of values, and a \code{string} that is either a single
#' character value or a vector of values.  It counts the number of occurrances of the string in the text and returns this count.  By default
#' the function counts overlapping strings as separate occurrances (i.e., they both add to the count), however the user can set \code{allow.overlap}
#' to \code{FALSE} to disallow occurrances of the string that overlap with previously counted occurrances (i.e., the latter occurrance does not
#' add to the count).  If the inputs are single character strings then the function checks for the occurrance of the string as a substring of the
#' text.  If the inputs are longer vectors of numeric or character values then the elements of the vectors are treated as individual symbols and
#' the vector gives the string of symbols that is checked in the count.
#'
#' @usage \code{stringcount()}
#' @param text Either a single character string or a numeric/character vector
#' @param string A string expressed as a numeric/character vector
#' @param language A list of strings with each string expressed as a numeric/character vector
#' @param sublanguages A matrix or list specifying a class of sublanguages within the language (see details above)
#' @param allow.overlap Logical; if \code{TRUE} then string occurrances are counted even if they overlap with previously counted occurrances
#' @param simple.form Logical; if \code{TRUE} the inputs \code{string}, \code{language} and \code{alphabet} are taken in simple form (see note)
#' @return The count value for the number of times the string occurs in the text

stringcount <- function(text, string = NULL, language = NULL, sublanguages = NULL,
                        allow.overlap = TRUE, simple.form = TRUE) {

  #Check input simple.form
  if (!is.vector(simple.form))                                 { stop('Error: Input simple.form should be a single logical value') }
  if (length(simple.form) != 1)                                { stop('Error: Input simple.form should be a single logical value') }
  if (!is.logical(simple.form))                                { stop('Error: Input simple.form should be a single logical value') }

  #Check input text
  if (!is.vector(text))                                        { stop('Error: text must be a vector') }
  if (simple.form) {
    if (length(text) > 1) {
      warning('You have specified that you are using simple form for your inputs, but input text has multiple elements; it does not appear to be in simple form') }
    TEXT <- sc(text) } else { TEXT <- text }
  ALPHABET.TEXT <- sort(unique(TEXT))
  n <- length(TEXT)

  #Check inputs string and language
  if ((missing(string))&(missing(language)))                   { stop('Error: You must input a string or language') }
  if ((!missing(string))&(!missing(language)))                 { stop('Error: Input a string or language, but not both') }
  if (!missing(string))   {
    if (simple.form) { LANGUAGE <- list(sc(string))     } else { LANGUAGE <- list(string) }
    S <- 1
    names(LANGUAGE) <- 'String[1]' }
  if (!missing(language)) {
    if (simple.form) { LANGUAGE <- lapply(language, sc) } else { LANGUAGE <- language }
    S <- length(LANGUAGE)
    names(LANGUAGE) <- sprintf('String[%s]', 1:S) }
  ALPHABET.LANG <- sort(unique(unlist(LANGUAGE)))

  #Generate minimum alphabet for text and language
  ALPHABET <- sort(unique(c(ALPHABET.TEXT, ALPHABET.LANG)))
  K <- length(ALPHABET)

  #Check input sublanguages
  #Convert input to a list if it is not already in this form
  SUB.INCL <- FALSE
  if (missing(sublanguages)) {
    SUBLANGUAGES <- as.list(1:S)
    D <- S
  } else {
    SUB.INCL <- TRUE
    if ((!is.list(sublanguages))&(!is.matrix(sublanguages)))   { stop('Error: Input sublanguages should be a list or binary matrix') }
    if (is.list(sublanguages)) {
      SUBLANGUAGES <- sublanguages
      D <- length(SUBLANGUAGES)
      for (d in 1:D) {
        SUB <- SUBLANGUAGES[[d]]
        if (!is.vector(SUB))                                   {
          stop(paste0('Error: Input sublanguage ', d, ' is not in the proper form')) }
        if (!is.numeric(SUB))                                  {
          stop(paste0('Error: Input sublanguage ', d, ' is not in the proper form')) }
        SUB <- sort(unique(SUB))
        for (i in 1:length(SUB)) {
          if (!(SUB[i] %in% 1:S)) {
            stop(paste0('Error: Input sublanguage ', d, ' refers to string ', SUB[i], ', which is not present in the language')) }
          SUBLANGUAGES[[d]] <- SUB } } }
    if (is.matrix(sublanguages)) {
      if (nrow(sublanguages) != S)                             { stop('Error: Input sublanguages should have one row for each string in the language') }
      if (!all(MM == as.logical(MM)))                          { stop('Error: Input sublanguages should be a binary/logical matrix') }
      D <- ncol(sublanguages)
      SUBLANGUAGES <- vector(mode = 'list', length = D)
      for (d in 1:D) { SUBLANGUAGES[[d]] <- which(sublanguages[, d] == 1) } } }

  #Check input allow.overlap
  if (!is.vector(allow.overlap))                               { stop('Error: Input allow.overlap should be a single logical value') }
  if (length(allow.overlap) != 1)                              { stop('Error: Input allow.overlap should be a single logical value') }
  if (!is.logical(allow.overlap))                              { stop('Error: Input allow.overlap should be a single logical value') }
  ALLOW.OVERLAP <- allow.overlap

  #-----------------------------------------------------------------------------------------------------------
  #Compute the string-counts using the DFA
  #-----------------------------------------------------------------------------------------------------------

  #Compute the DFA for the inputs
  if (SUB.INCL) {
    DFA1 <- DFA(language = LANGUAGE, sublanguages = SUBLANGUAGES, alphabet = ALPHABET,
               allow.overlap = ALLOW.OVERLAP, simple.form = FALSE, minimise = TRUE)
    } else {
    DFA1 <- DFA(language = LANGUAGE, alphabet = ALPHABET,
                 allow.overlap = ALLOW.OVERLAP, simple.form = FALSE, minimise = TRUE) }

  #Match text to alphabet
  MATCH <- match(TEXT, ALPHABET)

  #Compute the string-states and string-counts
  TT <- DFA1$transition.table
  FF <- DFA1$state.count
  STATES <- rep(0, n+1)
  COUNTS <- rep(0, D)
  for (i in 1:n) {
    STATES[i+1] <- TT[STATES[i]+1, MATCH[i]]
    COUNTS <- COUNTS + FF[STATES[i+1]+1, ] }
  COUNTS <- unname(COUNTS)

  #Output counts
  COUNTS }
