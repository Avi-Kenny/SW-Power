########################################.
##### Summarize simulation results #####
########################################.

sim <- readRDS("SimEngine.out/power_1_20240724.rds")

summ <- SimEngine::summarize(sim,
  list(stat="mean", x="reject", name="power")
  # list(stat="bias_pct", estimate="est", truth="true_tate"),
  # list(stat="sd", x="est"),
  # list(stat="coverage", estimate="est", se="se", truth="true_tate", name="cov")
)
if (sim$config$seed==764437398) { summ %<>% dplyr::filter(icc!=0.05) } # !!!!! TEMP

summ %<>% dplyr::mutate(
  # n_clust_per_seq = paste0(n_clust_per_seq, " clusters / seq"),
  icc = paste0("ICC: ", icc),
  estimand = factor(estimand, levels=c("TATE", "PTE-1", "PTE-S"))
)

plot <- ggplot(summ, aes(x=n_sequences, y=power, color=model)) +
  geom_line() +
  geom_point() +
  facet_grid(rows=dplyr::vars(icc), cols=dplyr::vars(estimand)) +
  scale_y_continuous(
    breaks = seq(0,1,0.2),
    limits = c(0,1),
    labels=scales::percent
  ) +
  labs(y="Power", x="Number of sequences", color="Model") +
  scale_color_manual(values=c("#009E73", "#56B4E9", "#CC79A7", "#E69F00")) +
  theme(legend.position="bottom")

ggsave(
  filename = paste0("../Figures + Tables/", cfg2$d,
                    " fig_models.pdf"),
  plot=plot, device="pdf", width=7, height=5
)



############################################.
##### Test dataset generating function #####
############################################.

if (F) {
  
  # Generate data
  dat <- generate_dataset(
    # data_type = "normal",
    # mu = 5,
    # sigma = 0.5,
    # tau = 0.2,
    # beta_j = c(0,2,4,6),
    # delta_s = c(5,5,5),
    # gamma_j = c(0,0,0,0),
    # n_sequences = 6,
    # n_clust_per_seq = 4,
    # n_ind_per_cluster = 10,
    # n_extra_time_points = 0
  )
  
  # Visualize cluster means
  dat_plot <- dat %>%
    dplyr::group_by(i,j,seq) %>%
    dplyr::summarize(y_bar=mean(y))
  ggplot(
    dat_plot,
    aes(x=j, y=y_bar, color=factor(seq), group=factor(i))
  ) +
    geom_line(alpha=0.5) +
    labs(color="Sequence")
  
}



############################.
##### VIZ: Scatterplot #####
############################.

# Figures produced: fig_1, fig_2

if (F) {
  
  # sim <- readRDS("SimEngine.out/estimation_1_20231112.rds")
  
  summ <- sim %>% SimEngine::summarize(
    list(stat="mean", x="est"),
    list(stat="bias", estimate="est", truth="true_tate")
  )
  
  print(summ)
  
}
