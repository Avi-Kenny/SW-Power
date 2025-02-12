
df_res <- data.frame(
  which = character(),
  icc = double(),
  n_seq = integer(),
  effect_size = double(),
  power = double()
)

# for (icc in c(0, 0.01, 0.05, 0.1)) {
for (icc in c(0.01)) {
  # for (n_seq in c(4,6,8,10)) {
  for (n_seq in c(10)) {
    # for (eff in c(0.15,0.185)) {
    for (eff in c(0.15)) {
      
      # Baseline power
      pwr_base <- calc_power(
        model = "ETI",
        n_sequences = n_seq,
        n_clust_per_seq = 2,
        n_ind_per_cell = 20,
        effect_size = eff,
        icc = icc,
        cac = 1,
        n_omit = 0,
        n_wash = 0
      )
      
      # Double n_ind_per_cell
      pwr_doubleInd <- calc_power(
        model = "ETI",
        n_sequences = n_seq,
        n_clust_per_seq = 2,
        n_ind_per_cell = 20*2,
        effect_size = eff,
        icc = icc,
        cac = 1,
        n_omit = 0,
        n_wash = 0
      )
      
      # Double n_clust_per_seq
      pwr_doubleClust <- calc_power(
        model = "ETI",
        n_sequences = n_seq,
        n_clust_per_seq = 2*2,
        n_ind_per_cell = 20,
        effect_size = eff,
        icc = icc,
        cac = 1,
        n_omit = 0,
        n_wash = 0
      )
      
      df_res[nrow(df_res)+1,] <- list("Baseline", icc, n_seq, eff, pwr_base)
      df_res[nrow(df_res)+1,] <- list("Double Ind", icc, n_seq, eff, pwr_doubleInd)
      df_res[nrow(df_res)+1,] <- list("Double Clust", icc, n_seq, eff, pwr_doubleClust)
      
    }
  }
}

df_res
# dplyr::filter(
#   df_res,
#   (icc==0.01 & effect_size==0.15) | (icc==0.05 & effect_size==0.185)
# )
