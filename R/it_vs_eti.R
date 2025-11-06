####################################.
##### Declare helper functions #####
####################################.

# Setup
cfg2 <- list(
  suppress_title = T,
  d = format(Sys.time(), "%Y-%m-%d")
)

# Plotting function
make_plot <- function(x_lab, which_plot, df, ymax, ylab, facets=T,
                      legend_pos="bottom") {
  
  df %<>% dplyr::mutate(
    n_seq = factor(
      paste(n_seq, "sequences"),
      levels = paste(sort(unique(df$n_seq)), "sequences")
    )
  )
  
  if (ymax==3) {
    scale_y <- scale_y_log10(breaks=seq(1,ymax,0.2), limits=c(1,ymax))
  } else if (ymax==25) {
    scale_y <- scale_y_log10(breaks=c(1,1.5,2,3,5,10,25), limits=c(1,ymax))
  }
  
  if (which_plot=="basic") {
    
    lab_col <- "CAC"
    theme_ <- theme(legend.position=legend_pos)
    scale_x_cts <- scale_x_continuous(minor_breaks=seq(0,0.2,0.01))
    df %<>% dplyr::mutate(grp=n_wo)
    
    plot <- ggplot(df, aes(x=icc, y=ssr, group=factor(cac), color=factor(cac),
                           shape=factor(cac))) +
      geom_line() +
      geom_point() +
      scale_color_manual(values=c("#009E73", "#56B4E9", "#CC79A7", "#E69F00")) +
      scale_shape_manual(values=c(16,15,18,17)) +
      scale_x_cts +
      scale_y +
      labs(x=x_lab, y=ylab, color=lab_col, shape=lab_col) +
      theme_
    
    if (facets) { plot <- plot + facet_wrap(~n_seq, ncol=2) }
    
  } else if (which_plot=="n_omit") {
    
    lab_col <- "Estimand"
    theme_ <- theme(legend.position=legend_pos)
    scale_x_cts <- scale_x_continuous(minor_breaks=seq(0,0.2,0.01))
    df %<>% dplyr::mutate(
      n_wo = ifelse(n_wo==0, "TATE(0,S)", paste0("TATE(0,S-", n_wo, ")")),
      n_wo = factor(n_wo, levels=c("TATE(0,S)", "TATE(0,S-1)", "TATE(0,S-2)", "TATE(0,S-3)")),
      grp = icc
    )
    
    plot <- ggplot(df, aes(x=icc, y=ssr, group=grp, color=factor(n_wo),
                           shape=factor(n_wo))) +
      geom_line(color="darkgrey", linewidth=0.4) +
      geom_point() +
      scale_color_manual(values=c("#009E73", "#56B4E9", "#CC79A7", "#E69F00")) +
      scale_shape_manual(values=c(16,15,18,17)) +
      scale_x_cts +
      scale_y +
      labs(x=x_lab, y=ylab, color=lab_col, shape=lab_col) +
      theme_
    
    if (facets) { plot <- plot + facet_wrap(~n_seq, ncol=2) }
    
  }
  
  return(plot)
  
}



##########################################.
##### Plot: SSR (basic; individuals) #####
##########################################.

if (cfg$regen_objs) {
  
  sim <- new_sim()
  sim %<>% set_levels(
    n_seq = c(6,9,12,15),
    # n_seq = c(6), # !!!!!
    icc = c(c(0,0.01,0.02),seq(0.03,0.15,0.03)),
    cac = c(0.6, 0.75, 0.9, 1),
    n_wo = list("0"=list(0,0))
  )
  sim %<>% set_config(
    parallel = T,
    n_cores = 8,
    packages = c("swCRTdesign")
  )
  sim %<>% set_script(function() {
    
    ss_ratio(
      power = 0.9,
      n_sequences = L$n_seq,
      n_clust_per_seq = 4,
      effect_size = 0.2,
      icc = L$icc,
      cac = L$cac,
      n_omit = L$n_wo[[1]],
      n_wash = L$n_wo[[2]]
    )
    
  })
  print("BASIC (ind): start")
  sim %<>% run()
  print("BASIC (ind): end")
  saveRDS(sim, "objs/sim_basic_ind.rds")
  
} else {
  
  sim <- readRDS("objs/sim_basic_ind.rds")
  
}

plot <- make_plot(
  x_lab = "ICC",
  which_plot = "basic",
  df = sim$results,
  ymax = 25,
  ylab = "Sample size ratio (individuals)"
)
ggsave(
  filename = paste0("../Figures + Tables/", cfg2$d, " fig_SSR_basic_ind.pdf"),
  plot=plot, device="pdf", width=8, height=5
)



#######################################.
##### Plot: SSR (basic; clusters) #####
#######################################.

if (cfg$regen_objs) {
  
  sim <- new_sim()
  sim %<>% set_levels(
    n_seq = c(6,9,12,15),
    icc = c(c(0,0.01,0.02),seq(0.03,0.15,0.03)),
    cac = c(0.6, 0.75, 0.9, 1),
    n_wo = list("0"=list(0,0))
  )
  sim %<>% set_config(
    parallel = T,
    n_cores = 8,
    packages = c("swCRTdesign")
  )
  sim %<>% set_script(function() {
    
    ss_ratio2(
      power = 0.9,
      n_sequences = L$n_seq,
      n_ind_per_cell = 5,
      effect_size = 0.2,
      icc = L$icc,
      cac = L$cac,
      n_omit = L$n_wo[[1]],
      n_wash = L$n_wo[[2]]
    )
    
  })
  print("BASIC (clust): start")
  sim %<>% run()
  print("BASIC (clust): end")
  saveRDS(sim, "objs/sim_basic_clust.rds")
  
} else {
  
  sim <- readRDS("objs/sim_basic_clust.rds")
  
}

plot <- make_plot(
  x_lab = "ICC",
  which_plot = "basic",
  df = sim$results,
  ymax = 3,
  ylab = "Sample size ratio (clusters)"
)
ggsave(
  filename = paste0("../Figures + Tables/", cfg2$d, " fig_SSR_basic_clust.pdf"),
  plot=plot, device="pdf", width=8, height=5
)



##############################.
##### Plot: SSR (n_omit) #####
##############################.

if (cfg$regen_objs) {
  
  sim <- new_sim()
  sim %<>% set_levels(
    n_seq = 6,
    icc = seq(0,0.15,0.01), # !!!!! New
    # icc = c(c(0,0.01,0.02),seq(0.03,0.15,0.03)),
    # icc = c(0,0.01,0.05,0.1,0.15), # !!!!!
    cac = 0.75,
    n_wo = list("0"=list(0,0), "1"=list(1,0), "3"=list(3,0))
  )
  sim %<>% set_config(
    parallel = T,
    n_cores = 8,
    packages = c("swCRTdesign")
  )
  sim %<>% set_script(function() {
    
    ss_ratio2(
      power = 0.9,
      n_sequences = L$n_seq,
      n_ind_per_cell = 5,
      effect_size = 0.2,
      icc = L$icc,
      cac = L$cac,
      n_omit = L$n_wo[[1]],
      n_wash = L$n_wo[[2]]
    )
    
  })
  print("N_OMIT: start")
  sim %<>% run()
  print("N_OMIT: end")
  saveRDS(sim, "objs/sim_n_omit.rds")
  
} else {
  
  sim <- readRDS("objs/sim_n_omit.rds")
  
}

plot <- make_plot(
  x_lab = "ICC",
  which_plot = "n_omit",
  df = sim$results,
  ymax = 3,
  ylab = "Sample size ratio (clusters)",
  facets = F,
  legend_pos = "right"
)
ggsave(
  filename = paste0("../Figures + Tables/", cfg2$d, " fig_SSR_n_omit.pdf"),
  plot=plot, device="pdf", width=6, height=4
)



# # Create combined plot
# # plot <- ggpubr::ggarrange(p01, p02, p03, p04, ncol=2, nrow=2)
# if (which_plot=="basic") {
#   plot <- ggpubr::ggarrange(p01, p03, ncol=2)
# } else {
#   plot <- ggpubr::ggarrange(p01, p03, ncol=2, legend="bottom", common.legend=T)
# }
# if (!cfg2$suppress_title) {
#   plot <- annotate_figure(
#     plot,
#     top = text_grob("Sample size ratio required for 90% power, TATE", size=14)
#   )
# }
# ggsave(
#   filename = paste0("../Figures + Tables/", cfg2$d, " fig_SSR_", which_plot,
#                     ".pdf"),
#   plot=plot, device="pdf", width=8, height=4
# )



##################################.
##### Plot: Staircase design #####
##################################.

if (cfg$regen_objs) {
  
  sim <- new_sim()
  sim %<>% set_levels(
    n_seq = c(4,6,8),
    icc = c(0.001,0.01,0.05,0.1),
    nba = c(1:5)
  )
  sim %<>% set_config(
    parallel = T,
    n_cores = 8, # Local
    # n_cores = 64, # Cluster
    packages = c("swCRTdesign")
  )
  sim %<>% set_script(function() {
    
    if (L$nba==1) {
      return(list(ss_eti=0, ss_it=0, ssr=1))
    } else {
      ss_ratio2(
        power = 0.9,
        n_sequences = L$n_seq,
        n_ind_per_cell = 5,
        effect_size = 0.2,
        icc = L$icc,
        cac = 0.75,
        n_omit = 0,
        n_wash = 0,
        staircase = T,
        n_before_and_after = L$nba,
        quiet = F
      )
    }
    
  })
  print("STAIRCASE: start")
  sim %<>% run()
  print("STAIRCASE: end")
  saveRDS(sim, "objs/sim_staircase.rds")
  
} else {
  
  sim <- readRDS("objs/sim_staircase.rds")
  
}

# plot
sim$results %<>% dplyr::mutate(
  n_seq = paste0(n_seq, " sequences"),
  nba = nba*2
)
plot <- ggplot(
  sim$results,
  aes(x=nba, y=ssr, color=factor(icc), shape=factor(icc))
) +
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
