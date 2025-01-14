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
                       effect_size, icc, cac, n_omit, n_wash, n_extra_c=0,
                       n_extra_t=0, staircase=F, n_before_and_after=NA) {
  
  if (staircase) {
    design <- swDsn(
      clusters=rep(n_clust_per_seq, n_sequences),
      extra.ctl.time = n_before_and_after-1,
      extra.trt.time = n_before_and_after-1
    )
    n_row <- nrow(design$swDsn)
    n_col <- ncol(design$swDsn)
    n_matrix <- matrix(NA, nrow=n_row, ncol=n_col)
    for (row in c(1:n_row)) {
      vec <- design$swDsn[row,]
      first_tx <- which.max(vec==1)
      n_vals <- rep(0, n_col)
      ind_1 <- first_tx-n_before_and_after
      ind_2 <- first_tx+n_before_and_after-1
      n_vals[ind_1:ind_2] <- n_ind_per_cell
      n_matrix[row,] <- n_vals
    }
  } else {
    design <- swDsn(
      clusters=rep(n_clust_per_seq, n_sequences),
      extra.ctl.time = n_extra_c,
      extra.trt.time = n_extra_t
    )
    n_matrix <- n_ind_per_cell
  }
  
  if (model=="IT") {
    H <- NULL
  } else if (model=="ETI") {
    if (staircase) {
      n_exp <- sum(design$swDsn.unique.clusters[1,])
      H <- c(rep(1, n_before_and_after), rep(0,n_exp-n_before_and_after))
    } else {
      H <- rep(1, n_sequences)
      if (n_omit!=0) {
        for (i in c(1:n_omit)) { H[round(n_sequences-(i-1))] <- 0 }
      }
      if (n_wash!=0) {
        for (i in c(1:n_wash)) { H[i] <- 0 }
      }
      H <- c(H, rep(0, n_extra_t))
    }
  }
  
  browser() # !!!!!
  power <- suppressWarnings({
    swPwr(
      design = design,
      distn = "gaussian",
      n = n_matrix,
      mu0 = 0,
      mu1 = effect_size,
      H = H,
      sigma = 1,
      icc = icc,
      cac = cac,
      alpha = 0.05
    )
  })
  print(as.numeric(power))
  
  return(as.numeric(power))
  
}
