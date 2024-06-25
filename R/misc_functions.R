#' Alias for indicator (as.integer) function
#' 
#' @param u Logical input
#' @return BInary (integer) output
In <- as.integer



#' Helper function for debugging; prints timestamps
#' 
#' @param num Number
#' @param msg Message
chk <- function(num, msg="") {
  if (msg=="") {
    str <- paste0("Check ", num, ": ", Sys.time())
  } else {
    str <- paste0("Check ", num, " (", msg, "): ", Sys.time())
  }
  print(str)
}



#' Expit function
#' 
#' @param x Numeric input
#' @return Number
expit <- function(x) {1/(exp(-x)+1)}
