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
#' @param n_omit How many time points should be omitted from the end of the TATE
#'     estimator
#' @return A list of the form list(est=0, se=0) representing the point estimate
#'     and SE estimate
analyze_data <- function(dat, cal_time, exp_time, re, n_omit) {
  
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
  
  # Create exposure time variables and calendar time factor
  for (s in c(1:params$n_sequences)) {
    dat[[paste0("s_", s)]] <- In(dat$s_ij==s)
  }
  
  # Create transformed calendar time variables
  if (cal_time=="cat") {
    
    dat[["j_fac"]] <- factor(dat$j)
    
  } else if (cal_time=="NCS") {
    
    knots <- seq(min(dat$j), max(dat$j), length.out=5)
    ns_basis <- splines::ns(
      x = dat$j,
      knots = knots[2:4],
      intercept = F,
      Boundary.knots = knots[c(1,5)]
    )
    dat$j_1 <- ns_basis[,1]
    dat$j_2 <- ns_basis[,2]
    dat$j_3 <- ns_basis[,3]
    dat$j_4 <- ns_basis[,4]
    rm(knots,ns_basis)
    
  }

  # Create cluster-time interaction for "random time effect"
  dat$ij <- as.integer(factor(paste0(dat$i,"-",dat$j)))
  
  if (exp_time=="IT") {
    
    formula <- "y_ij ~ x_ij + "
    if (cal_time=="cat") { formula <- paste0(formula, "j_fac - 1 + ") }
    if (re=="cluster") {
      formula <- paste0(formula, "(1|i)")
    } else if (re=="cluster+time") {
      formula <- paste0(formula, "(1|i) + (1|ij)")
    }
    
    suppressMessages({
      model <- lme4::lmer(formula, data=dat)
    })
    
    summ <- summary(model)
    est <- summ$coefficients["x_ij", "Estimate"]
    se <- summ$coefficients["x_ij", "Std. Error"]
    
  } else if (exp_time=="cat") {
    
    formula <- "y_ij ~ "
    for (s in c(1:params$n_sequences)) {
      formula <- paste0(formula, "s_", s, " + ")
    }
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

    suppressMessages({
      model <- lme4::lmer(formula, data=dat)
    })
    
    summ <- summary(model)
    coeff_names <- names(summ$coefficients[,1])
    delta_s_hat <- as.numeric(summary(model)$coefficients[,1])
    sigma_s_hat <- vcov(model)
    indices <- c(1:length(coeff_names))[substr(coeff_names, 1, 2)=="s_"]
    coeff_names <- coeff_names[indices]
    delta_s_hat <- delta_s_hat[indices]
    sigma_s_hat <- sigma_s_hat[indices,indices]
    sigma_s_hat <- as.matrix(sigma_s_hat)
    
    n_exp_times <- length(delta_s_hat)
    H_len <- n_exp_times - n_omit
    A <- (1/H_len) * matrix(c(rep(1,H_len), rep(0,n_omit)), nrow=1)
    est <- as.numeric(A %*% delta_s_hat)
    se <- sqrt(as.numeric(A %*% sigma_s_hat %*% t(A)))
    
  } else if (exp_time=="NCS") {
    
    
    
    formula <- "y_ij ~ "
    for (s in c(1:4)) {
      formula <- paste0(formula, "s_", s, " + ")
    }
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
    
    suppressMessages({
      model <- lme4::lmer(formula, data=dat)
    })
    
    summ <- summary(model)
    coeff_names <- names(summ$coefficients[,1])
    delta_s_hat <- as.numeric(summary(model)$coefficients[,1])
    sigma_s_hat <- vcov(model)
    indices <- c(1:length(coeff_names))[substr(coeff_names, 1, 2)=="s_"]
    coeff_names <- coeff_names[indices]
    delta_s_hat <- delta_s_hat[indices]
    sigma_s_hat <- sigma_s_hat[indices,indices]
    sigma_s_hat <- as.matrix(sigma_s_hat)
    
    n_exp_times <- length(delta_s_hat)
    H_len <- n_exp_times - n_omit
    A <- (1/H_len) * matrix(c(rep(1,H_len), rep(0,n_omit)), nrow=1)
    est <- as.numeric(A %*% delta_s_hat)
    se <- sqrt(as.numeric(A %*% sigma_s_hat %*% t(A)))
    
  }
  
  return(list(est=est, se=se))
  
}
