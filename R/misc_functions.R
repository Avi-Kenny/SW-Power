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
                       n_extra_t=0, staircase=F, n_before_and_after=NA,
                       n_baseline_scale=1) {
  
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
    if (n_baseline_scale==1) {
      n_matrix <- n_ind_per_cell
    } else {
      n_row <- nrow(design$swDsn)
      n_col <- ncol(design$swDsn)
      n_matrix <- matrix(n_ind_per_cell, nrow=n_row, ncol=n_col)
      n_matrix[,1] <- n_ind_per_cell*n_baseline_scale
    }
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
  
  # browser() # !!!!!
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

  return(as.numeric(power))
  
}

# Function to calculate sample size corresponding to desired power
calc_ss2 <- function(power, model, n_sequences, n_ind_per_cell, effect_size,
                     icc, cac, n_omit, n_wash, staircase=F,
                     n_before_and_after=NA) {
  
  power_n <- function(n_clust) {
    calc_power(
      model = model,
      n_sequences = n_sequences,
      n_clust_per_seq = n_clust,
      n_ind_per_cell = n_ind_per_cell,
      effect_size = effect_size,
      icc = icc,
      cac = cac,
      n_omit = n_omit,
      n_wash = n_wash,
      staircase = staircase,
      n_before_and_after = n_before_and_after
    )
  }
  
  n_lo <- 1
  n_up <- 2
  power_lo <- power_n(n_lo)
  power_up <- power_n(n_up)
  if (power_lo>power) { warning("Power sufficient with n_clust==1.") }
  
  # Get initial sample size range that contains 90% power
  while(power_up<power) {
    n_lo <- n_up
    power_lo <- power_up
    n_up <- round(n_up*2)
    power_up <- power_n(n_up)
    if (n_up>10^5) {
      warning("90% power not possible with ", model, " model")
      return(list(n=NA, power=NA))
    }
  }
  
  # Bisect interval until 90% range is found
  while(n_up-n_lo>1) {
    n_mid <- round(mean(c(n_lo,n_up)))
    power_mid <- power_n(n_mid)
    if (power_mid<power) {
      n_lo <- n_mid
      power_lo <- power_mid
    } else {
      n_up <- n_mid
      power_up <- power_mid
    }
  }
  
  if (abs(power-power_lo)<abs(power-power_up)) {
    return(list(n=n_lo, power=power_lo))
  } else {
    return(list(n=n_up, power=power_up))
  }
  
}

# Function to calculate sample size ratio (ETI:IT)
ss_ratio2 <- function(power, n_sequences, n_ind_per_cell, effect_size, icc, cac,
                      n_omit, n_wash, staircase=F, n_before_and_after=NA,
                      quiet=T) {
  
  ss_it <- calc_ss2(
    power = power,
    model = "IT",
    n_sequences = n_sequences,
    n_ind_per_cell = n_ind_per_cell,
    effect_size = effect_size,
    icc = icc,
    cac = cac,
    n_omit = n_omit,
    n_wash = n_wash,
    staircase = staircase,
    n_before_and_after = n_before_and_after
  )
  
  ss_eti <- calc_ss2(
    power = power,
    model = "ETI",
    n_sequences = n_sequences,
    n_ind_per_cell = n_ind_per_cell,
    effect_size = effect_size,
    icc = icc,
    cac = cac,
    n_omit = n_omit,
    n_wash = n_wash,
    staircase = staircase,
    n_before_and_after = n_before_and_after
  )
  
  if (!quiet) {
    cat(paste0("IT sample size: ", ss_it$n, "\n"))
    cat(paste0("ETI sample size: ", ss_eti$n, "\n"))
    cat(paste0("SSR: ", round(ss_eti$n/ss_it$n,1), "\n"))
  }
  
  if (is.na(ss_eti$n) || is.na(ss_it$n)) {
    return(list(ss_eti=NA, ss_it=NA, ssr=NA))
  } else {
    return(list(
      ss_eti = ss_eti$n,
      ss_it = ss_it$n,
      ssr = ss_eti$n/ss_it$n
    ))
  }
  
}

# Function to calculate sample size corresponding to desired power
calc_ss <- function(power, model, n_sequences, n_clust_per_seq, effect_size,
                    icc, cac, n_omit, n_wash, staircase=F,
                    n_before_and_after=NA) {
  
  power_n <- function(n) {
    calc_power(
      model = model,
      n_sequences = n_sequences,
      n_clust_per_seq = n_clust_per_seq,
      n_ind_per_cell = n,
      effect_size = effect_size,
      icc = icc,
      cac = cac,
      n_omit = n_omit,
      n_wash = n_wash,
      staircase = staircase,
      n_before_and_after = n_before_and_after
    )
  }
  
  n_lo <- 1
  n_up <- 2
  power_lo <- power_n(n_lo)
  power_up <- power_n(n_up)
  if (power_lo>power) { warning("Power sufficient with n_ind_per_cell==1.") }
  
  # Get initial sample size range that contains 90% power
  while(power_up<power) {
    n_lo <- n_up
    power_lo <- power_up
    n_up <- round(n_up*2)
    power_up <- power_n(n_up)
    if (n_up>10^5) {
      warning("90% power not possible with ", model, " model")
      return(list(n=NA, power=NA))
    }
  }
  
  # Bisect interval until 90% range is found
  while(n_up-n_lo>1) {
    n_mid <- round(mean(c(n_lo,n_up)))
    power_mid <- power_n(n_mid)
    if (power_mid<power) {
      n_lo <- n_mid
      power_lo <- power_mid
    } else {
      n_up <- n_mid
      power_up <- power_mid
    }
  }
  
  if (abs(power-power_lo)<abs(power-power_up)) {
    return(list(n=n_lo, power=power_lo))
  } else {
    return(list(n=n_up, power=power_up))
  }
  
}

# Function to calculate sample size ratio (ETI:IT)
ss_ratio <- function(power, n_sequences, n_clust_per_seq, effect_size, icc, cac,
                     n_omit, n_wash, staircase=F, n_before_and_after=NA,
                     quiet=T) {
  
  ss_it <- calc_ss(
    power = power,
    model = "IT",
    n_sequences = n_sequences,
    n_clust_per_seq = n_clust_per_seq,
    effect_size = effect_size,
    icc = icc,
    cac = cac,
    n_omit = n_omit,
    n_wash = n_wash,
    staircase = staircase,
    n_before_and_after = n_before_and_after
  )
  
  ss_eti <- calc_ss(
    power = power,
    model = "ETI",
    n_sequences = n_sequences,
    n_clust_per_seq = n_clust_per_seq,
    effect_size = effect_size,
    icc = icc,
    cac = cac,
    n_omit = n_omit,
    n_wash = n_wash,
    staircase = staircase,
    n_before_and_after = n_before_and_after
  )
  
  if (!quiet) {
    cat(paste0("IT sample size: ", ss_it$n, "\n"))
    cat(paste0("ETI sample size: ", ss_eti$n, "\n"))
    cat(paste0("SSR: ", round(ss_eti$n/ss_it$n,1), "\n"))
  }
  
  if (is.na(ss_eti$n) || is.na(ss_it$n)) {
    return(list(ss_eti=NA, ss_it=NA, ssr=NA))
  } else {
    return(list(
      ss_eti = ss_eti$n,
      ss_it = ss_it$n,
      ssr = ss_eti$n/ss_it$n
    ))
  }
  
}
