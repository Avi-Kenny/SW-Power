####################################.
##### Declare helper functions #####
####################################.

# Setup
cfg2 <- list(
  suppress_title = T,
  d = format(Sys.time(), "%Y-%m-%d")
)

# Function to calculate sample size corresponding to desired power
calc_ss <- function(power, model, n_sequences, n_clust_per_seq, effect_size,
                    icc, n_omit, n_wash) {
  
  power_n <- function(n) {
    calc_power(
      model = model,
      n_sequences = n_sequences,
      n_clust_per_seq = n_clust_per_seq,
      n_ind_per_cell = n,
      effect_size = effect_size,
      icc = icc,
      cac = 1,
      n_omit = n_omit,
      n_wash = n_wash
    )
  }
  
  n_lo <- 1
  n_up <- 2
  power_lo <- power_n(n_lo)
  power_up <- power_n(n_up)
  if (power_lo>power) { stop("Power sufficient with n_ind_per_cell==1.") }
  
  # Get initial sample size range that contains 90% power
  while(power_up<power) {
    n_lo <- n_up
    power_lo <- power_up
    n_up <- round(n_up*2)
    power_up <- power_n(n_up)
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
ss_ratio <- function(power, n_sequences, n_clust_per_seq, effect_size, icc,
                     n_omit, n_wash) {
  
  ss_it <- calc_ss(
    power = power,
    model = "IT",
    n_sequences = n_sequences,
    n_clust_per_seq = n_clust_per_seq,
    effect_size = effect_size,
    icc = icc,
    n_omit = n_omit,
    n_wash = n_wash
  )
  
  ss_eti <- calc_ss(
    power = power,
    model = "ETI",
    n_sequences = n_sequences,
    n_clust_per_seq = n_clust_per_seq,
    effect_size = effect_size,
    icc = icc,
    n_omit = n_omit,
    n_wash = n_wash
  )
  
  return(ss_eti$n/ss_it$n)
  
}

# Plotting function
make_plot <- function(x_lab, which_plot, df) {
  
  if (which_plot=="n_omit") {
    lab_col <- "# time points omitted from end"
    scale_y <- scale_y_continuous(breaks=seq(1,3,0.2), limits=c(1,3))
    theme_ <- theme(legend.position="bottom")
  } else if (which_plot=="n_wash") {
    lab_col <- "# washout time points"
    scale_y <- scale_y_continuous(breaks=seq(1,7,0.5), limits=c(1,7.1))
    theme_ <- theme(legend.position="bottom")
  } else if (which_plot=="basic") {
    lab_col <- ""
    scale_y <- scale_y_continuous(breaks=seq(1,3,0.2), limits=c(1,3))
    theme_ <- theme(legend.position="none")
  } else if (which_plot=="pte") {
    lab_col <- "s, PTE(s)"
    scale_y <- scale_y_continuous(breaks=c(1:12), limits=c(1,12.5))
    theme_ <- theme(legend.position="bottom")
    df %<>% dplyr::mutate(
      n_wo = pmax(2*n_wo,1) # Converts from c(0,1,2,3) to PTE scale c(1,2,4,6)
    )
  }
  
  plot <- ggplot(df, aes(x=x, y=y, color=factor(n_wo))) +
    geom_line() +
    geom_point() +
    scale_color_manual(values=c("#009E73", "#56B4E9", "#CC79A7", "#E69F00")) +
    scale_y +
    labs(x=x_lab, y="ETI/IT sample size ratio", color=lab_col) +
    theme_
  
  return(plot)
  
}



#####################.
##### Plot: SSR #####
#####################.

# for (which_plot in c("basic", "n_omit", "n_wash", "pte")) {
for (which_plot in c("pte")) {
  
  if (which_plot=="basic") {
    ow_vec <- list(c(0,0))
  } else if (which_plot=="n_omit") {
    ow_vec <- list(c(0,0), c(1,0), c(2,0), c(3,0))
  } else if (which_plot=="n_wash") {
    ow_vec <- list(c(0,0), c(0,1), c(0,2), c(0,3))
  } else if (which_plot=="pte") {
    ow_vec <- list(c(5,0), c(4,1), c(2,3), c(0,5))
  }
  
  # Plot component 1: Sequences
  n_sequences <- c(2:10)
  v_ratios_1 <- c()
  for (ow in ow_vec) {
    v_ratios_1 <- c(v_ratios_1, sapply(n_sequences, function(x) {
      if (which_plot=="pte") { ow[1] <- x-(ow[2]+1) } # Hack to get PTE plot to work correctly
      tryCatch(
        expr = {
          return(ss_ratio(
            power = 0.9,
            n_sequences = x,
            n_clust_per_seq = 4,
            effect_size = 0.1,
            icc = 0.1,
            n_omit = ow[1],
            n_wash = ow[2]
          ))
        },
        error = function(e) { return(NA) }
      )
    }))
  }
  p01 <- make_plot(
    x_lab = "# sequences",
    which_plot = which_plot,
    df = data.frame(
      x = rep(n_sequences, length(ow_vec)),
      y = v_ratios_1,
      n_wo = rep(c(1:length(ow_vec))-1, each=length(n_sequences))
    )
  )
  
  # Plot component 2: Effect sizes
  effect_sizes <- seq(0.05,0.4,0.05)
  v_ratios_2 <- c()
  for (ow in ow_vec) {
    v_ratios_2 <- c(v_ratios_2, sapply(effect_sizes, function(x) {
      tryCatch(
        expr = {
          return(ss_ratio(
            power = 0.9,
            n_sequences = 6,
            n_clust_per_seq = 4,
            effect_size = x,
            icc = 0.1,
            n_omit = ow[1],
            n_wash = ow[2]
          ))
        },
        error = function(e) { return(NA) }
      )
    }))
  }
  p02 <- make_plot(
    x_lab = "Effect size",
    which_plot = which_plot,
    df = data.frame(
      x = rep(effect_sizes, length(ow_vec)),
      y = v_ratios_2,
      n_wo = rep(c(1:length(ow_vec))-1, each=length(effect_sizes))
    )
  )
  
  # Plot component 3: ICCs
  iccs <- c(seq(0,0.01,0.002),seq(0.02,0.2,0.01))
  v_ratios_3 <- c()
  for (ow in ow_vec) {
    v_ratios_3 <- c(v_ratios_3, sapply(iccs, function(x) {
      tryCatch(
        expr = {
          return(ss_ratio(
            power = 0.9,
            n_sequences = 6,
            n_clust_per_seq = 4,
            effect_size = 0.1,
            icc = x,
            n_omit = ow[1],
            n_wash = ow[2]
          ))
        },
        error = function(e) { return(NA) }
      )
    }))
  }
  p03 <- make_plot(
    x_lab = "ICC",
    which_plot = which_plot,
    df = data.frame(
      x = rep(iccs, length(ow_vec)),
      y = v_ratios_3,
      n_wo = rep(c(1:length(ow_vec))-1, each=length(iccs))
    )
  )
  
  # Plot component 4: Clusters per sequence
  n_clust_per_seqs <- c(1:8)
  v_ratios_4 <- c()
  for (ow in ow_vec) {
    v_ratios_4 <- c(v_ratios_4, sapply(n_clust_per_seqs, function(x) {
      tryCatch(
        expr = {
          return(ss_ratio(
            power = 0.9,
            n_sequences = 6,
            n_clust_per_seq = x,
            effect_size = 0.1,
            icc = 0.1,
            n_omit = ow[1],
            n_wash = ow[2]
          ))
        },
        error = function(e) { return(NA) }
      )
    }))
  }
  p04 <- make_plot(
    x_lab = "# clusters per sequence",
    which_plot = which_plot,
    df = data.frame(
      x = rep(n_clust_per_seqs, length(ow_vec)),
      y = v_ratios_4,
      n_wo = rep(c(1:length(ow_vec))-1, each=length(n_clust_per_seqs))
    )
  )
  
  # Create combined plot
  plot <- ggpubr::ggarrange(p01, p02, p03, p04, ncol=2, nrow=2)
  if (!cfg2$suppress_title) {
    plot <- annotate_figure(
      plot,
      top = text_grob("Sample size ratio required for 90% power, TATE", size=14)
    )
  }
  ggsave(
    filename = paste0("../Figures + Tables/", cfg2$d, " fig_SSR_", which_plot,
                      ".pdf"),
    plot=plot, device="pdf", width=10, height=8
  )
  
}
