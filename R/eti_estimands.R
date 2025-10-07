make_plot2 <- function(x_lab, which_plot, df) {
  
  if (which_plot=="n_omit") {
    df %<>% dplyr::mutate(
      n_wo = ifelse(n_wo==0, "TATE(0,S)", paste0("TATE(0,S-", n_wo, ")")),
      n_wo = factor(n_wo, levels=c("TATE(0,S)", "TATE(0,S-1)", "TATE(0,S-2)", "TATE(0,S-3)"))
    )
  } else if (which_plot=="n_wash") {
    df %<>% dplyr::mutate(
      n_wo = factor(paste0("TATE(", n_wo, ",S)"))
    )
  } else if (which_plot=="pte") {
    df %<>% dplyr::mutate(
      n_wo = pmax(2*n_wo,1) # Converts from c(0,1,2,3) to PTE scale c(1,2,4,6)
    )
  }
  
  plot <- ggplot(df, aes(x=n_wo, y=y, color=factor(x), group=factor(x))) +
    geom_line() +
    geom_point() +
    scale_color_manual(values=c("#009E73", "#56B4E9", "#CC79A7", "#E69F00")) +
    scale_y_log10(breaks=c(50,100,200,400,800,1600)) +
    labs(x="Estimand", y="Sample size required for 90% power (log scale)", color=x_lab) +
    theme(legend.position="right")
  
  if (which_plot=="pte") {
    plot <- plot + scale_x_continuous(
      breaks = c(1, 2, 4, 6),
      labels = c("PTE(1)", "PTE(2)", "PTE(4)", "PTE(6)"),
      minor_breaks = c(1:6)
    )
  }
  
  return(plot)
  
}



for (which_plot in c("n_omit", "n_wash", "pte")) {

  if (which_plot=="n_omit") {
    ow_vec <- list(c(0,0), c(1,0), c(2,0), c(3,0))
  } else if (which_plot=="n_wash") {
    ow_vec <- list(c(0,0), c(0,1), c(0,2), c(0,3))
  } else if (which_plot=="pte") {
    ow_vec <- list(c(5,0), c(4,1), c(2,3), c(0,5))
  }
  
  # Plot component 1: Sequences
  # n_sequences <- c(2:10)
  n_sequences <- c(4,6,8,10)
  ss_vec_1 <- c()
  for (ow in ow_vec) {
    ss_vec_1 <- c(ss_vec_1, sapply(n_sequences, function(x) {
      if (which_plot=="pte") { ow[1] <- x-(ow[2]+1) } # Hack to get PTE plot to work correctly
      tryCatch(
        expr = {
          return(calc_ss(
            power = 0.9,
            model = "ETI",
            n_sequences = x,
            n_clust_per_seq = 4,
            effect_size = 0.1,
            icc = 0.05,
            n_omit = ow[1],
            n_wash = ow[2]
          )$n)
        },
        error = function(e) { return(NA) }
      )
    }))
  }
  plot <- make_plot2(
    x_lab = "# sequences",
    which_plot = which_plot,
    df = data.frame(
      x = rep(n_sequences, length(ow_vec)),
      y = ss_vec_1,
      n_wo = rep(c(1:length(ow_vec))-1, each=length(n_sequences))
    )
  )
  
  ggsave(
    filename = paste0("../Figures + Tables/", cfg2$d, " fig_SS_ETI_", which_plot,
                      ".pdf"),
    plot=plot, device="pdf", width=6, height=4
  )
  
}



