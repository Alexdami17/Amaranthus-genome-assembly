#=========================================== Analyze Ks with GAM ========================================================================
#setwd("/home/tranel-lab-user/Alex_tmp/ka_ks_analysis/palmeri_kaks/palmer_mmseqs_out")
df <- read.delim("/home/tranel-lab-user/Alex_tmp/ka_ks_analysis/palmeri_kaks/palmer_mmseqs_out/palmeri.Chr03.kaks.out", header = TRUE)

# Convert bp → Mb if needed
df$start_mb <- df$start / 1e6

# Use to call a function that is placed within its own file 
source("/home/tranel-lab-user/Alex_tmp/ka_ks_analysis/scripts/fit_mcp_strata.R")

# Use to find potential outliers or extreme values dominating the dataset
df_tail <- df[order(df$Ks), c("gene_id", "start_mb", "Ks")]
tail(df_tail, 30)

# Run a quick QC using LOESS
qc <- analyze_ks(df, pos_col = "start_mb", ks_col = "Ks",
                 loess_span = 0.5, loess_degree = 1,
                 loess_se = TRUE, smooth_method = "loess",
                 k_outliers = 10, remove_outliers = TRUE,
                 x_label = "Genomic position (Mbp)",
                 y_max = 0.4,
                 tick_by = 10,
                 ylim = c(0, 0.4))
qc$plot 
#+ scale_x_continuous(limits = c(0, 20),breaks = seq(0, 100, by = 5))

# fit GAM. Although the code can toggle between GAM and LOESS
res_gam <- analyze_ks(
  df,
  pos_col       = "start_mb",
  ks_col        = "Ks",
  smooth_method = "gam",
  gam_family = "scat",
  gam_k         = 30,
  gam_method    = "REML",
  gam_se        = TRUE,
  k_outliers = 10,
  remove_outliers = TRUE,
  y_max = 0.4,
  ylim          = c(0, 0.4),
  tick_by       = 5,
  compute_derivatives = FALSE
)
res_gam$plot
summary(res_gam$gam_fit)

# Check fit
# <- gam(Ks ~ s(start_mb, k = 30), data = df, method = "REML")
#fit <- gam(Ks ~ s(start_mb, k = 100), data = df, method="GCV.Cp", select=FALSE)
#summary(fit)

# Save plots
ggsave("/home/tranel-lab-user/Alex_tmp/ka_ks_analysis/palmeri_kaks/palmer_mmseqs_out/qc_GAM_Chr03_plot.png", res_gam$plot, width = 8, height = 4, dpi = 300)
ggsave("/home/tranel-lab-user/Alex_tmp/ka_ks_analysis/palmeri_kaks/palmer_mmseqs_out/qc_GAM_Chr03_plot.pdf", res_gam$plot, width = 9, height = 4)
ggsave("/home/tranel-lab-user/Alex_tmp/ka_ks_analysis/palmeri_kaks/palmer_mmseqs_out/qc_GAM_Chr03_plot.svg", res_gam$plot, width = 9, height = 4)
#=========================================================================================================================================


#============================================== Option 1: fit a Bayesian MCP model =================================================================
# Unconstrained MCP (fastest/simplest)
set.seed(1234)
# Run strata analysis
res_mcp <- fit_mcp_strata(
  df,
  x_col  = "start_mb",
  y_col  = "Ks",
  max_cp = 1,
  iter   = 100000,
  adapt  = 20000,
  chains = 4,
  cores  = 4,
  y_max  = 0.4,
  ylim   = c(0, 0.4),
  rosner = TRUE,
  k_outliers = 15,
  remove_outliers = TRUE,
  tick_by = 10,
  x_label_custom = "Genomic position (Mbp)",
  enforce_min_segment = FALSE
)

# Best model info
res_mcp$best_n_changepoints   # best_n_changepoints is data-driven (0–max_cp)
res_mcp$cp_summary
res_mcp$change_points         # now should show two numeric cp positions
table(res_mcp$data$strata_best)

# Plot with plateaus + change points
res_mcp$plot                  # 3 colours + dashed cp lines

# Data with strata assignments
head(res_mcp$data[, c("gene_id", "start_mb", "Ks", "strata_best")])

# loo comparison table and weights
res_mcp$loo_compare
# Gives the pseudo-BMA weights reported in the paper
res_mcp$loo_weights

# Rebuild ggplot manually
ggplot(res_mcp$data_plot, res_mcp$layers$mapping) +
  res_mcp$layers$points +
  res_mcp$layers$segments +
  { if (!is.null(res_mcp$layers$vlines)) res_mcp$layers$vlines else NULL } +
  res_mcp$layers$color_scale +
  res_mcp$layers$base_theme +
  coord_cartesian(ylim = c(0, 0.3))

# Optionally get posterior cp panels later
cp_names <- paste0("cp_", seq_len(res_mcp$best_n_changepoints))
posterior_plots <- plot_mcp_cp_posterior(res_mcp, cp_pars = cp_names)
posterior_plots$fit_plot       # posterior fitted curves along x
# Make a better plot!
p <- posterior_plots$fit_plot

p2 <- p +
  scale_x_continuous(breaks = seq(0, 60, 5)) + # x ticks every 5 Mb
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 18)
  ) +
  labs(
    x = "Genomic position (Mbp)",
    y = expression(K[s]),
    title = "Posterior fit"
  ) + coord_cartesian(ylim = c(0, 0.4))
ggsave("/home/tranel-lab-user/Alex_tmp/ka_ks_analysis/palmeri_kaks/palmer_mmseqs_out/posterior_fit_Chr03_plot.svg", p2, width = 9, height = 4)
ggsave("/home/tranel-lab-user/Alex_tmp/ka_ks_analysis/palmeri_kaks/palmer_mmseqs_out/posterior_fit_Chr03_plot_2.svg", p2, width = 12, height = 5.3, units = "in")

# Check whether cp_1 is really stable (not drifting)
best_fit <- res_mcp$fits[[res_mcp$best_index]]
res_mcp$cp_summary
# The code below is plotting density. If the density is broad or multi-modal, 
# the "stratum boundary" is not stable (kernel density estimate?)
mcp::plot_pars(best_fit, pars = "cp_1")
# Look at segment means (plateaus) for interpretability
# The "y_mean" is the mean Ks value to report for the strata in the manuscript
res_mcp$seg_df

# Alternatively, you can compute the mean explicitly from the exact data used for plotting
res_mcp$data_plot %>%
  dplyr::group_by(strata_best) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_Ks = mean(Ks, na.rm=TRUE),
    sd_Ks = sd(Ks, na.rm=TRUE),
    .groups="drop"
  )

# Check Pareto-k diagnostics
# From Vehtari et al. (2017), the standard interpretation is k < 0.5 (excellent, LOO is reliable),
# 0.5 <- k < 0.7 (ok, usable), 0.7 <- k < 1.0 (warning, LOO becoming unstable), k >- 1.0 (problematic, LOO unreliable)
lapply(res_mcp$loo_list, function(x) max(x$diagnostics$pareto_k))

#==============================================================================================================


#=================================== Option 3: fit a Bayesian MCP model =======================================
# Constrained + CI (Option 1 + 3 combined; most reproducible)
res_mcp2 <- fit_mcp_strata(
  df, x_col="start_mb", y_col="Ks",
  max_cp=3,
  iter=100000, adapt=20000, chains=4, cores=4,
  y_max=0.3,
  rosner=TRUE, k_outliers=10, remove_outliers=TRUE,
  
  enforce_min_segment = TRUE,
  min_seg_n = 200,
  cp_buffer = 0.25,
  
  cp_ci_level = 0.95,
  cp_line_stat = "median",
  add_cp_ci_bars = TRUE
)

res_mcp2$loo_compare
res_mcp2$loo_weights

res_mcp2$best_n_changepoints
res_mcp2$change_points # # numeric cp position(s) for the best model; length = best_n_changepoints
table(res_mcp2$data$strata_best)

res_mcp2$plot         # dashed = chosen mean/median; dotted = CI bounds (optional)
res_mcp2$cp_summary   # table of cp mean/median/CI

