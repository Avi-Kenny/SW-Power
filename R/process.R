########################################.
##### Sim results: analysis models #####
########################################.

# Setup
cfg2 <- list(
  d = format(Sys.time(), "%Y-%m-%d")
)

sim <- readRDS("SimEngine.out/power_1_20250225.rds")

summ <- SimEngine::summarize(sim,
  list(stat="mean", x="reject", name="power")
  # list(stat="bias_pct", estimate="est", truth="true_tate"),
  # list(stat="sd", x="est"),
  # list(stat="coverage", estimate="est", se="se", truth="true_tate", name="cov")
)

summ %<>% dplyr::mutate(
  # n_clust_per_seq = paste0(n_clust_per_seq, " clusters / seq"),
  icc = paste0("ICC: ", icc),
  estimand = factor(estimand, levels=c("TATE", "PTE-1", "PTE-S"))
)

# Jitter "ETI cat cal" points
summ %<>% dplyr::mutate(
  power = ifelse(model=="ETI, cat cal time", power-0.01, power)
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
  scale_x_continuous(breaks=c(6,12,18,24)) +
  labs(y="Power", x="Number of sequences", color="Model") +
  scale_color_manual(values=c("#009E73", "#56B4E9", "#CC79A7", "#E69F00")) +
  theme(legend.position="bottom")

ggsave(
  filename = paste0("../Figures + Tables/", cfg2$d,
                    " fig_models.pdf"),
  plot=plot, device="pdf", width=7, height=5
)



#######################################.
##### Sim results: washout models #####
#######################################.

sim <- readRDS("SimEngine.out/washout_1_20250213.rds")

summ <- SimEngine::summarize(sim,
                             list(stat="mean", x="reject", name="power")
                             # list(stat="bias_pct", estimate="est", truth="true_tate"),
                             # list(stat="sd", x="est"),
                             # list(stat="coverage", estimate="est", se="se", truth="true_tate", name="cov")
)

summ %<>% dplyr::mutate(
  icc = paste0("ICC: ", icc),
  model = dplyr::case_when(
    model=="IT" ~ "IT*",
    model=="PIT" ~ "DCT",
    TRUE ~ model
  )
  # model = ifelse(model=="IT", "IT*", model)
)

# Jitter "ETI cat cal" points
summ %<>% dplyr::mutate(
  power = ifelse(model=="ETI, cat cal time", power-0.01, power)
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
                    " fig_washout.pdf"),
  plot=plot, device="pdf", width=7, height=5
)



