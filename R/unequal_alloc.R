##################################################.
##### Unequal allocation design optimization #####
##################################################.

create_heatmap <- function(n_mtx, title) {
  
  n_seq <- nrow(n_mtx)
  step_vec <- c(1:n_seq)
  
  data_table <- as.table(t(n_mtx))
  attr(data_table, "dimnames") %<>% lapply(., function(x) {
    as.numeric(as.factor(x))
  })
  data_df <- as.data.frame(data_table)
  data_df$modFreq <- data_df$Freq
  # data_df$modFreq <- log(data_df$Freq+1)
  heatmap_plot <- ggplot(data_df, aes(x = Var1, y = Var2, fill = modFreq)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low="beige", high="forestgreen") +
    labs(
      x = "Time",
      y = "Sequence",
      title = title
    ) +
    theme(legend.position="none") +
    scale_y_discrete(limits=rev) +
    geom_text(aes(label = Freq), color="#333", size=4) +
    geom_step(
      data = data.frame(
        x = c(step_vec+0.5,step_vec[n_seq]+0.51),
        y = c(rev(step_vec)+0.5, 0.5)
      ),
      mapping = aes(x=x, y=y),
      direction = "vh",
      inherit.aes = F,
      linewidth = 1,
      color = "brown"
    )
  
  return(heatmap_plot)  
  
}



estimands <- c("IT", "TATE(0,S)", "TATE(0,3)", "PTE(1)", "PTE(S)", "PTE(S-1)")

for (estimand in estimands) {
  
  set.seed(1)
  print(paste0("Start: ", estimand))
  print(Sys.time())
  n_clust_per_seq <- 4
  n_sequences <- 6
  n_per_cell <- 5
  n_reps <- 100000
  design <- swDsn(clusters=rep(n_clust_per_seq, n_sequences))
  power_max <- 0
  err_counter <- 0
  incr_counter <- 0
  
  if (estimand=="IT") {
    H <- NULL
  } else if (estimand=="TATE(0,S)") {
    H <- rep(1, n_col-1)
  } else if (estimand=="TATE(0,3)") {
    H <- c(rep(1,3), rep(0, n_col-4))
  } else if (estimand=="PTE(1)") {
    H <- c(1, rep(0, n_col-2))
  } else if (estimand=="PTE(S)") {
    H <- c(rep(0, n_col-2), 1)
  } else if (estimand=="PTE(S-1)") {
    H <- c(rep(0, n_col-3), 1, 0)
  } else {
    stop("`estimand` incorrectly specified")
  }
  
  # Create starting matrix
  n_matrix_old <- design$swDsn
  n_row <- nrow(n_matrix_old)
  n_col <- ncol(n_matrix_old)
  for (row in c(1:n_row)) {
    for (col in c(1:n_col)) {
      n_matrix_old[row,col] <- n_per_cell
    }
  }
  
  # Run algorithm
  for (i in c(1:n_reps)) {
    
    n_matrix_new <- n_matrix_old
    cell_1 <- c(sample(c(1:n_row), size=1), sample(c(1:n_col), size=1))
    cell_2 <- c(sample(c(1:n_row), size=1), sample(c(1:n_col), size=1))
    
    if (F) {
      # Option 1 (large jumps)
      n <- n_matrix_old[cell_1[1],cell_1[2]] + n_matrix_old[cell_2[1],cell_2[2]]
      n_1 <- sample(c(0:n), size=1)
      n_2 <- n - n_1
    } else {
      # Option 2 (small jumps)
      n_1 <- n_matrix_old[cell_1[1],cell_1[2]]
      n_2 <- n_matrix_old[cell_2[1],cell_2[2]]
      if (n_1!=0) {
        n_1 <- n_1 - 1
        n_2 <- n_2 + 1
      }
    }
    n_matrix_new[cell_1[1],cell_1[2]] <- n_1
    n_matrix_new[cell_2[1],cell_2[2]] <- n_2
    
    if (!identical(cell_1,cell_2)) {
      tryCatch(
        expr = {
          
          suppressMessages({
            power <- as.numeric(swPwr(
              design = design,
              distn = "gaussian",
              n = n_matrix_new,
              mu0 = 0,
              mu1 = 0.1,
              H = H,
              sigma = 1,
              icc = 0.1,
              cac = 1,
              alpha = 0.05
            ))
          })
          
          if (power>power_max) {
            incr_counter <<- incr_counter + 1
            power_max <- power
            n_matrix_old <- n_matrix_new
          }
          
        },
        error = function(e) {
          err_counter <<- err_counter + 1
        }
      )
    }
    
  }
  
  # print(n_matrix_new)
  # print(err_counter)
  # print(incr_counter)
  
  # Generate and save heatmap
  plot_title <- paste0("Quasi-optimal design for estimating ", estimand, " (",
                       format(n_reps, scientific=F, big.mark=","), " reps)")
  plot <- create_heatmap(n_mtx=n_matrix_new, title=plot_title)
  
  ggsave(
    filename = paste0("../Figures + Tables/fig_opt_design_",
                      fs::path_sanitize(estimand), ".pdf"),
    plot = plot,
    device = "pdf",
    width = 6,
    height = 4
  )
  
}



# as.numeric(swPwr(
#   design = design,
#   distn = "gaussian",
#   # n = n_matrix_old*75, # 58% power
#   n = n_matrix_new*75, # 90% power
#   mu0 = 0,
#   mu1 = 0.1,
#   H = H,
#   sigma = 1,
#   icc = 0.1,
#   cac = 1,
#   alpha = 0.05
# ))
