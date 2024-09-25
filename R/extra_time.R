####################################.
##### Declare helper functions #####
####################################.

# Setup
cfg2 <- list(
  suppress_title = T,
  d = format(Sys.time(), "%Y-%m-%d")
)
n_extras <- list(c(3,0), c(2,0), c(1,0), c(0,0), c(0,0), c(0,1), c(0,2), c(0,3))
df_plot <- data.frame(
  "time" = integer(),
  "power" = double(),
  "which" = character(),
  "icc" = double()
)

# Populate dataframe
iccs <- c(0, 0.005, 0.01, 0.1, 0.2)

for (icc in iccs) {
  
  for (i in c(1:8)) {
    
    p <- calc_power(
      model = "ETI",
      n_sequences = 6,
      n_clust_per_seq = 4,
      n_ind_per_cell = 30,
      effect_size = 0.15,
      icc = icc,
      n_omit = 0,
      n_wash = 0,
      n_extra_c = n_extras[[i]][1],
      n_extra_t = n_extras[[i]][2]
    )
    
    df_plot[nrow(df_plot)+1,] <- list(
      sum(n_extras[[i]]),
      p,
      ifelse(i<=4, "Extra control time (start of design)",
             "Extra treatment time (end of design)"),
      paste0("ICC: ", icc)
    )
    
  }
  
}

# Generate plot
plot <- ggplot(df_plot, aes(x=time, y=power, color=which)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=c("#009E73", "#56B4E9")) +
  scale_y_continuous(labels=scales::percent) +
  facet_grid(cols=dplyr::vars(icc)) +
  labs(
    x = "# Extra time points",
    y = "Power",
    color = NULL
  ) +
  theme(legend.position="bottom")

ggsave(
  filename = paste0("../Figures + Tables/", cfg2$d,
                    " fig_power_extra_time_tate.pdf"),
  plot=plot, device="pdf", width=9, height=4
)
