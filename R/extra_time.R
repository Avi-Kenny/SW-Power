####################################.
##### Declare helper functions #####
####################################.

# Setup
cfg2 <- list(
  suppress_title = T,
  d = format(Sys.time(), "%Y-%m-%d")
)

create_df <- function(iccs, cacs, effect_sizes) {
  
  n_extras <- list(
    c(3,0), c(2,0), c(1,0), c(0,0), c(0,0), c(0,1), c(0,2), c(0,3)
  )
  df_plot <- data.frame(
    "time" = integer(),
    "power" = double(),
    "which" = character(),
    "icc" = character(),
    "cac" = character(),
    "effect_size" = double()
  )
  
  for (effect_size in effect_sizes) {
    for (icc in iccs) {
      for (cac in cacs) {
        for (i in c(1:8)) {
          
          p <- calc_power(
            model = "ETI",
            n_sequences = 6,
            n_clust_per_seq = 4,
            n_ind_per_cell = 30,
            effect_size = effect_size,
            icc = icc,
            cac = cac,
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
            paste0("ICC: ", icc),
            paste0("CAC: ", cac),
            effect_size
          )
          
          # This code was used to confirm that the power gain from increasing
          #     the sample size at baseline is equivalent to the gain from
          #     adding time points to the start of the study
          # if (i<=4) {
          #   
          #   p2 <- calc_power(
          #     model = "ETI",
          #     n_sequences = 6,
          #     n_clust_per_seq = 4,
          #     n_ind_per_cell = 30,
          #     effect_size = effect_size,
          #     icc = icc,
          #     cac = cac,
          #     n_omit = 0,
          #     n_wash = 0,
          #     n_extra_c = 0,
          #     n_extra_t = 0,
          #     n_baseline_scale = 1+n_extras[[i]][1]
          #   )
          #   
          #   df_plot[nrow(df_plot)+1,] <- list(
          #     sum(n_extras[[i]]),
          #     p2,
          #     "Extra observations in period 1",
          #     paste0("ICC: ", icc),
          #     paste0("CAC: ", cac),
          #     effect_size
          #   )
          #   
          # }
          
        }
      }
    }
  }
  
  return(df_plot)
  
}

# Generate plot 1 (CAC not varied)
df_plot_1 <- create_df(
  iccs = c(0, 0.005, 0.01, 0.05, 0.1),
  cacs = 1,
  effect_sizes = 0.15
)
plot_1 <- ggplot(df_plot_1, aes(x=time, y=power, color=which)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=c("#009E73", "#56B4E9")) +
  # scale_color_manual(values=c("#009E73", "#56B4E9", "#CC79A7")) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0.4,1),
    breaks = seq(0.4,1,0.1)
  ) +
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
  plot=plot_1, device="pdf", width=9, height=4
)

# Generate plot 2 (CAC, eff_size varied)
df_plot_2 <- create_df(
  iccs = c(0.05),
  cacs = c(0.5,0.75,1),
  effect_sizes = c(0.1,0.15,0.2,0.25)
)
df_plot_2 %<>% dplyr::mutate(
  effect_size = paste0("Eff. size: ", format(effect_size, nsmall=2))
)
plot_2 <- ggplot(df_plot_2, aes(x=time, y=power, color=which)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=c("#009E73", "#56B4E9")) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0.2,1),
    breaks = seq(0.2,1,0.2)
  ) +
  facet_grid(rows=dplyr::vars(cac), cols=dplyr::vars(effect_size)) +
  labs(
    x = "# Extra time points",
    y = "Power",
    color = NULL
  ) +
  theme(legend.position="bottom")

ggsave(
  filename = paste0("../Figures + Tables/", cfg2$d,
                    " fig_power_extra_time_tate_cac.pdf"),
  plot=plot_2, device="pdf", width=9, height=6
)

# Baseline power
# calc_power(model="ETI", n_sequences=6, n_clust_per_seq=4, n_ind_per_cell=30, effect_size=0.210, icc=0.05, cac=0.25, n_omit=0, n_wash=0, n_extra_c=0, n_extra_t=0)
df_1 <- create_df(iccs=0.05, cacs=1, effect_sizes=0.187)
df_2 <- create_df(iccs=0.05, cacs=0.75, effect_sizes=0.208)
df_3 <- create_df(iccs=0.05, cacs=0.5, effect_sizes=0.217)
df_4 <- create_df(iccs=0.05, cacs=0.25, effect_sizes=0.210)
df_plot_3 <- rbind(df_1,df_2,df_3,df_4)
df_plot_3 %<>% dplyr::mutate(
  effect_size = paste0("Eff. size: ", format(effect_size, nsmall=3)),
  scenario = factor(
    paste0(cac, "; ", effect_size),
    levels = c("CAC: 1; Eff. size: 0.187",
             "CAC: 0.75; Eff. size: 0.208",
             "CAC: 0.5; Eff. size: 0.217",
             "CAC: 0.25; Eff. size: 0.210")
  )
)
plot_3 <- ggplot(df_plot_3, aes(x=time, y=power, color=which)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=c("#009E73", "#56B4E9")) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0.6,1),
    breaks = seq(0.2,1,0.1)
  ) +
  facet_grid(cols=dplyr::vars(scenario)) +
  labs(
    x = "# Extra time points",
    y = "Power",
    color = NULL
  ) +
  theme(legend.position="bottom")

ggsave(
  filename = paste0("../Figures + Tables/", cfg2$d,
                    " fig_power_extra_time_tate_cac2.pdf"),
  plot=plot_3, device="pdf", width=9, height=4
)
