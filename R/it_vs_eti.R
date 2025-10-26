####################################.
##### Declare helper functions #####
####################################.

# Setup
cfg2 <- list(
  suppress_title = T,
  d = format(Sys.time(), "%Y-%m-%d")
)
# n_sequences <- c(3:8) # !!!!!
n_sequences <- seq(4,20,2)

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
ss_ratio <- function(power, n_sequences, n_clust_per_seq, effect_size, icc, cac,
                     n_omit, n_wash, staircase=F, n_before_and_after=NA) {
  
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
  
  return(ss_eti$n/ss_it$n)
  
}

# Plotting function
make_plot <- function(x_lab, which_plot, df) {
  
  if (which_plot=="n_omit") {
    # lab_col <- "# time points omitted from end"
    lab_col <- "Estimand"
    scale_y <- scale_y_continuous(breaks=seq(1,3,0.2), limits=c(1,3))
    theme_ <- theme(legend.position="bottom")
    df %<>% dplyr::mutate(
      n_wo = ifelse(n_wo==0, "TATE(0,S)", paste0("TATE(0,S-", n_wo, ")")),
      n_wo = factor(n_wo, levels=c("TATE(0,S)", "TATE(0,S-1)", "TATE(0,S-2)", "TATE(0,S-3)")),
      grp = x
    )
    # } else if (which_plot=="n_wash") {
  #   lab_col <- "# washout time points"
  #   scale_y <- scale_y_continuous(breaks=seq(1,7,0.5), limits=c(1,7.1))
  #   theme_ <- theme(legend.position="bottom")
  } else if (which_plot=="basic") {
    lab_col <- ""
    scale_y <- scale_y_continuous(breaks=seq(1,3,0.2), limits=c(1,3))
    theme_ <- theme(legend.position="none")
    df %<>% dplyr::mutate(grp=n_wo)
  # } else if (which_plot=="pte") {
  #   lab_col <- "s, PTE(s)"
  #   scale_y <- scale_y_continuous(breaks=c(1:12), limits=c(1,12.5))
  #   theme_ <- theme(legend.position="bottom")
  #   df %<>% dplyr::mutate(
  #     n_wo = pmax(2*n_wo,1) # Converts from c(0,1,2,3) to PTE scale c(1,2,4,6)
  #   )
  }
  
  if (x_lab=="# sequences") {
    scale_x_cts <- scale_x_continuous(breaks=n_sequences, minor_breaks=NULL)
  } else if (x_lab=="ICC") {
    scale_x_cts <- scale_x_continuous(minor_breaks=seq(0,0.2,0.01))
  } else if (x_lab=="CAC") {
    scale_x_cts <- scale_x_continuous(minor_breaks=seq(0.3,1,0.1))
  }
  plot <- ggplot(df, aes(x=x, y=y, group=grp, color=factor(n_wo),
                         shape=factor(n_wo))) +
    geom_line(color="darkgrey", linewidth=0.4) +
    geom_point() +
    scale_color_manual(values=c("#009E73", "#56B4E9", "#CC79A7", "#E69F00")) +
    scale_shape_manual(values=c(16,15,18,17)) +
    scale_x_cts +
    scale_y +
    labs(x=x_lab, y="ETI/IT sample size ratio", color=lab_col, shape=lab_col) +
    theme_
  
  return(plot)
  
}



#####################.
##### Plot: SSR #####
#####################.

for (which_plot in c("basic", "n_omit")) {
# for (which_plot in c("basic")) {
  
  if (which_plot=="basic") {
    ow_vec <- list(c(0,0))
  } else if (which_plot=="n_omit") {
    ow_vec <- list(c(0,0), c(1,0), c(2,0), c(3,0))
  # } else if (which_plot=="n_wash") {
  #   ow_vec <- list(c(0,0), c(0,1), c(0,2), c(0,3))
  # } else if (which_plot=="pte") {
  #   ow_vec <- list(c(5,0), c(4,1), c(2,3), c(0,5))
  }
  
  # Plot component 1: Sequences
  if (cfg$regen_objs) {
    
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
              icc = 0.05,
              cac = 1,
              n_omit = ow[1],
              n_wash = ow[2]
            ))
          },
          error = function(e) { return(NA) }
        )
      }))
    }
    saveRDS(v_ratios_1, paste0("objs/v_ratios_1_", which_plot, ".rds"))
    
  } else {
    
    v_ratios_1 <- readRDS(paste0("objs/v_ratios_1_", which_plot, ".rds"))
    
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
  
  # Plot component 3: ICCs
  iccs <- c(c(0,0.005),seq(0.01,0.2,0.01))
  if (cfg$regen_objs) {
    
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
              cac = 1,
              n_omit = ow[1],
              n_wash = ow[2]
            ))
          },
          error = function(e) { return(NA) }
        )
      }))
    }
    saveRDS(v_ratios_3, paste0("objs/v_ratios_3_", which_plot, ".rds"))
    
  } else {
    
    v_ratios_3 <- readRDS(paste0("objs/v_ratios_3_", which_plot, ".rds"))
    
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
  
  # Create combined plot
  # plot <- ggpubr::ggarrange(p01, p02, p03, p04, ncol=2, nrow=2)
  if (which_plot=="basic") {
    plot <- ggpubr::ggarrange(p01, p03, ncol=2)
  } else {
    plot <- ggpubr::ggarrange(p01, p03, ncol=2, legend="bottom", common.legend=T)
  }
  if (!cfg2$suppress_title) {
    plot <- annotate_figure(
      plot,
      top = text_grob("Sample size ratio required for 90% power, TATE", size=14)
    )
  }
  ggsave(
    filename = paste0("../Figures + Tables/", cfg2$d, " fig_SSR_", which_plot,
                      ".pdf"),
    plot=plot, device="pdf", width=8, height=4
  )
  
}



##################################.
##### Plot: Staircase design #####
##################################.

if (cfg$regen_objs) {
  
  df_res <- data.frame(
    nba = integer(),
    icc = double(),
    n_seq = integer(),
    n_clust_per_seq = integer(),
    ssr = double()
  )
  
  for (icc in c(0.001,0.01,0.05,0.1)) {
    for (n_clust_per_seq in c(2)) {
      for (n_seq in c(4,6,8)) {
        for (nba in c(1:5)) {
          
          if (nba==1) {
            ssr <- 1
          } else {
            ssr <- suppressMessages({
              ss_ratio(
                power = 0.9,
                n_sequences = n_seq,
                n_clust_per_seq = n_clust_per_seq,
                effect_size = 0.1,
                icc = icc,
                cac = 1,
                n_omit = 0,
                n_wash = 0,
                staircase = T,
                n_before_and_after = nba
              )
            })
          }
          
          df_res[nrow(df_res)+1,] <- list(nba, icc, n_seq, n_clust_per_seq, ssr)
          
        }
      }
    }
  }
  
  df_res %<>% dplyr::mutate(
    n_seq = paste0(n_seq, " sequences"),
    nba = nba*2
  )
  
  saveRDS(df_res, "objs/df_res.rds")
  
} else {
  
  df_res <- readRDS("objs/df_res.rds")
  
}

plot <- ggplot(df_res, aes(x=nba, y=ssr, color=factor(icc), shape=factor(icc))) +
  facet_grid(cols=dplyr::vars(n_seq)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values=c("#009E73", "#56B4E9", "#CC79A7", "#E69F00")) +
  scale_y_continuous(breaks=seq(1,3,0.5)) +
  scale_shape_manual(values=c(16,15,18,17)) +
  labs(x="Number of time points observed per sequence",
       y="ETI/IT sample size ratio", color="ICC", shape="ICC")
ggsave(
    filename = paste0("../Figures + Tables/", cfg2$d, " fig_staircase.pdf"),
  plot=plot, device="pdf", width=8, height=3
)
