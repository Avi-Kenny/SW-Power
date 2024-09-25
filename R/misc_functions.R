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



# Helper function to calculate power given inputs
calc_power <- function(model, n_sequences, n_clust_per_seq, n_ind_per_cell,
                       effect_size, icc, n_omit, n_wash, n_extra_c=0,
                       n_extra_t=0) {
  
  design <- swDsn(
    clusters=rep(n_clust_per_seq, n_sequences),
    extra.ctl.time = n_extra_c,
    extra.trt.time = n_extra_t
  )
  
  if (model=="IT") {
    H <- NULL
  } else if (model=="ETI") {
    H <- rep(1, n_sequences)
    if (n_omit!=0) {
      for (i in c(1:n_omit)) { H[round(n_sequences-(i-1))] <- 0 }
    }
    if (n_wash!=0) {
      for (i in c(1:n_wash)) { H[i] <- 0 }
    }
  }
  H <- c(H, rep(0, n_extra_t))
  
  power <- suppressWarnings({
    swPwr(
      design = design,
      distn = "gaussian",
      n = n_ind_per_cell,
      mu0 = 0,
      mu1 = effect_size,
      H = H,
      sigma = 1,
      icc = icc,
      # cac = 1,
      alpha = 0.05
    )
  })
  
  return(as.numeric(power))
  
}
