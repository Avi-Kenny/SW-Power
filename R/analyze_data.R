#' Generate one stepped wedge dataset
#'
#' @param dat A dataset returned by generate_data()
#' @param cal_time One of c("cat", "linear", "NCS"); how the secular
#'     time trend should be modeled
#' @param exp_time One of c("IT", "cat", "NCS"); how the exposure time trend
#'     should be modeled
#' @param re One of c("cluster", "cluster+time"); whether a cluster intercept
#'     should be included ("cluster") or a cluster intercept plus a
#'     cluster-period intercept ("cluster+time")
#' @param estimand One of c("TATE", "PTE-1", "PTE-S"); the target of estimation
#' @param return_curve Boolean; if true, the entire estimated curve is returned
#' @return A list of the form list(est=0, se=0) representing the point estimate
#'     and SE estimate
analyze_data <- function(dat, cal_time, exp_time, re, estimand, return_curve) {
  
  if (!(cal_time %in% c("cat", "linear", "NCS"))) {
    stop("`cal_time` misspecified.")
  }
  if (!(exp_time %in% c("IT", "cat", "NCS"))) {
    stop("`exp_time` misspecified.")
  }
  if (!(re %in% c("cluster", "cluster+time"))) {
    stop("`re` misspecified.")
  }
  
  params <- attr(dat, "params")
  
  # Create exposure time variables (factors or spline basis)
  S <- params$n_sequences
  if (exp_time=="cat") {
    
    for (s in c(1:S)) {
      dat[[paste0("s_", s)]] <- In(dat$s_ij==s)
    }
    
  } else if (exp_time=="NCS") {
    
    knots_exp <- seq(min(dat$s_ij), max(dat$s_ij), length.out=4)
    basis_exp <- splines::ns(
      x = dat$s_ij,
      knots = knots_exp[2:3],
      intercept = TRUE,
      Boundary.knots = knots_exp[c(1,4)]
    )
    # dat$s_1 <- basis_exp[,1]
    # dat$s_2 <- basis_exp[,2]
    # dat$s_3 <- basis_exp[,3]
    # dat$s_4 <- basis_exp[,4]
    
    # !!!!! Modification needed for spline with interept
    dat$s_1 <- basis_exp[,1] * dat$x_ij
    dat$s_2 <- basis_exp[,2] * dat$x_ij
    dat$s_3 <- basis_exp[,3] * dat$x_ij
    dat$s_4 <- basis_exp[,4] * dat$x_ij
    
    rm(basis_exp)
    
  }
  
  # Create calendar time variables (factors or spline basis)
  if (cal_time=="cat") {
    
    dat[["j_fac"]] <- factor(dat$j)
    
  } else if (cal_time=="NCS") {
    
    knots_cal <- seq(min(dat$j), max(dat$j), length.out=5)
    basis_cal <- splines::ns(
      x = dat$j,
      knots = knots_cal[2:4],
      intercept = FALSE,
      Boundary.knots = knots_cal[c(1,5)]
    )
    dat$j_1 <- basis_cal[,1]
    dat$j_2 <- basis_cal[,2]
    dat$j_3 <- basis_cal[,3]
    dat$j_4 <- basis_cal[,4]
    rm(knots_cal,basis_cal)
    
  }

  # Create cluster-time interaction for "random time effect"
  if (re=="cluster+time") {
    dat$ij <- as.integer(factor(paste0(dat$i,"-",dat$j)))
  }
  
  if (exp_time=="IT") {
    
    formula <- "y_ij ~ x_ij + "
    if (cal_time=="cat") { formula <- paste0(formula, "j_fac - 1 + ") }
    if (re=="cluster") {
      formula <- paste0(formula, "(1|i)")
    } else if (re=="cluster+time") {
      formula <- paste0(formula, "(1|i) + (1|ij)")
    }
    
    suppressMessages({ model <- lme4::lmer(formula, data=dat) })
    
    summ <- summary(model)
    est <- summ$coefficients["x_ij", "Estimate"]
    se <- summ$coefficients["x_ij", "Std. Error"]
    
  } else if (exp_time %in% c("cat", "NCS")) {
    
    formula <- "y_ij ~ "
    num_s_terms <- ifelse(exp_time=="cat", params$n_sequences, 4)
    for (s in c(1:num_s_terms)) { formula <- paste0(formula, "s_", s, " + ") }
    if (cal_time=="cat") {
      formula <- paste0(formula, "j_fac - 1 + ")
    } else if (cal_time=="linear") {
      formula <- paste0(formula, "j + ")
    } else if (cal_time=="NCS") {
      formula <- paste0(formula, "j_1 + j_2 + j_3 + j_4 + ")
    }
    if (re=="cluster") {
      formula <- paste0(formula, "(1|i)")
    } else if (re=="cluster+time") {
      formula <- paste0(formula, "(1|i) + (1|ij)")
    }

    suppressMessages({ model <- lme4::lmer(formula, data=dat) })
    
    summ <- summary(model)
    coeff_names <- names(summ$coefficients[,1])
    delta_s_hat <- as.numeric(summary(model)$coefficients[,1])
    sigma_s_hat <- vcov(model)
    indices <- c(1:length(coeff_names))[substr(coeff_names, 1, 2)=="s_"]
    coeff_names <- coeff_names[indices]
    delta_s_hat <- delta_s_hat[indices]
    sigma_s_hat <- sigma_s_hat[indices,indices]
    sigma_s_hat <- as.matrix(sigma_s_hat)
    
    if (exp_time=="NCS") {
      
      B <- splines::ns(
        x = c(1:S),
        knots = knots_exp[2:3],
        intercept = TRUE,
        Boundary.knots = knots_exp[c(1,4)]
      )
      delta_s_hat <- as.numeric(B %*% delta_s_hat)
      sigma_s_hat <- B %*% sigma_s_hat %*% t(B)
      
    }
    
    if (estimand=="TATE") {
      n_omit <- 0 # !!!!!
      H_len <- S - n_omit
      A <- (1/H_len) * matrix(c(rep(1,H_len), rep(0,n_omit)), nrow=1)
    } else if (estimand=="PTE-1") {
      A <- matrix(c(1, rep(0,S-1)), nrow=1)
    } else if (estimand=="PTE-S") {
      A <- matrix(c(rep(0,S-1), 1), nrow=1)
    }
    est <- as.numeric(A %*% delta_s_hat)
    se <- sqrt(as.numeric(A %*% sigma_s_hat %*% t(A)))
    
  }
  
  res <- list(est=est, se=se)
  
  if (return_curve) {
    if (exp_time=="IT") {
      res$curve <- rep(est, S)
    } else {
      res$curve <- delta_s_hat
    }
  }
  
  return(res)
  
}
