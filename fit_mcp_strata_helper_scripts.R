############################################################
## Ks (or dS) QC + Bayesian changepoint analysis helpers
## - analyze_ks(): Rosner outliers + LOESS
## - fit_mcp_strata(): Rosner outliers + mcp strata + loo
## - plot_mcp_cp_posterior(): helper for posterior / cp density
############################################################

# Required packages:
# install.packages(c("ggplot2", "dplyr", "EnvStats", "mcp", "loo", "rlang"))

library(ggplot2)
library(dplyr)
library(EnvStats)
library(mcp)
library(loo)
library(rlang)
library(mgcv)


############################################################
## 1. Quick QC plot with Rosner outliers + LOESS / GAM
############################################################

#' Analyze Ks (or dS) along a chromosome and detect outliers with Rosner's test
#'
#' Assumes the data frame has at least these columns:
#'   gene_id, start, Ka, Ks, Ka_Ks, Pvalue
#'
#' @param df              Data frame with Ks values and positions.
#' @param pos_col         Name of position column (e.g. "start" or "start_mb").
#' @param ks_col          Name of Ks column (e.g. "Ks" or "dS").
#' @param k_outliers      Fixed max number of outliers for rosnerTest.
#' @param k_prop          Optional proportion of points to allow as potential
#'                        outliers when k_outliers is NULL.
#' @param remove_outliers Logical; if TRUE, remove Rosner outliers from plot/data.
#' @param x_label         Label for x-axis.
#' @param ylim            Optional y-axis limits, e.g. c(0, 0.2).
#' @param tick_by         Optional spacing for X-axis ticks (e.g. 10 for 0,10,20,...).
#' @param smooth_method   "loess" (default) or "gam".
#'
#' LOESS-specific:
#' @param loess_span      LOESS span (fraction of data in each local fit).
#' @param loess_degree    LOESS polynomial degree (1 or 2).
#' @param loess_se        Logical; plot LOESS standard error ribbon?
#'
#' GAM-specific:
#' @param gam_k           Basis dimension k for s(x, k = gam_k).
#' @param gam_method      Smoothing parameter method, e.g. "REML" (default) or "GCV.Cp".
#' @param gam_se          Logical; plot GAM standard error ribbon?
#'
#' @return list with:
#'   - data:     data frame used for plotting (after optional outlier removal)
#'   - outliers: data frame of detected outliers
#'   - rosner:   raw rosnerTest result (or NULL if not run)
#'   - plot:     ggplot object

analyze_ks <- function(
    df,
    pos_col = "start",
    ks_col = "Ks",
    y_max = NULL,
    k_outliers = 15,
    k_prop = NULL,
    remove_outliers = TRUE,
    x_label = "Genomic position (Mbp)",
    ylim = NULL,
    tick_by = NULL,
    
    # ---- smoother choice ----
    smooth_method = c("loess", "gam"),
    
    # LOESS controls
    loess_span   = 0.4,
    loess_degree = 1,
    loess_se     = TRUE,
    
    # GAM controls
    gam_k       = 50,
    gam_method  = "REML",
    gam_se      = TRUE,
    
    # NEW: robust/heavy-tailed option (recommended when you see outlier sensitivity)
    gam_family  = c("scat", "gaussian"),
    
    # Prediction grid controls (GAM only)
    pred_n = 400,
    
    # Derivative controls (GAM only)
    compute_derivatives = TRUE,
    deriv_interval = "simultaneous",
    deriv_n = 200,
    
    # Rosner warning control
    cap_rosner_k_at_10 = FALSE
) {
  smooth_method <- match.arg(smooth_method)
  gam_family <- match.arg(gam_family)
  
  # ---- checks ----
  if (!pos_col %in% names(df)) stop("pos_col '", pos_col, "' not found in df.")
  if (!ks_col  %in% names(df)) stop("ks_col '", ks_col, "' not found in df.")
  
  # ---- helpers ----
  .pick_first_existing <- function(df_like, candidates) {
    hit <- intersect(candidates, names(df_like))
    if (length(hit) == 0) return(NULL)
    hit[1]
  }
  
  .coerce_numeric_vec <- function(x) {
    # Make a best effort to return an atomic numeric vector (or atomic character)
    if (is.data.frame(x)) {
      if (ncol(x) != 1) stop("Column resolves to multiple columns; provide a single column name.")
      x <- x[[1]]
    }
    if (is.list(x)) x <- unlist(x, use.names = FALSE)
    if (is.factor(x)) x <- as.character(x)
    
    # If character mostly numeric-like, coerce to numeric
    if (is.character(x)) {
      suppressWarnings(xn <- as.numeric(x))
      if (!all(is.na(xn)) && sum(!is.na(xn)) >= max(5, floor(0.5 * length(xn)))) {
        x <- xn
      }
    }
    
    if (!is.atomic(x)) stop("Column is not an atomic vector. Check column types.")
    x
  }
  
  .infer_xcol <- function(d, pos_col) {
    # 1) prefer exact pos_col
    if (!is.null(pos_col) && pos_col %in% names(d)) return(pos_col)
    
    # 2) common fallbacks
    candidates <- c("start_mb", "start", ".x", "x", "position")
    hit <- .pick_first_existing(d, candidates)
    if (!is.null(hit)) return(hit)
    
    # 3) fallback: a single non-dot numeric column (ignoring gratia metadata)
    nm <- names(d)
    keep <- nm[!grepl("^\\.", nm)]
    keep <- setdiff(keep, c("smooth", "by", "fs", "term", "component"))
    if (length(keep) > 0) {
      numish <- keep[vapply(d[keep], is.numeric, logical(1))]
      if (length(numish) >= 1) return(numish[1])
      return(keep[1])
    }
    NULL
  }
  
  .infer_dcol <- function(d) {
    dcol <- .pick_first_existing(d, c(".derivative", "derivative", "estimate", "est", "value", ".value"))
    if (!is.null(dcol)) return(dcol)
    nm <- names(d)
    hit <- nm[grepl("deriv", nm, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[1])
    NULL
  }
  
  .safe_gratia_derivatives <- function(fit, term = NULL, order = 1,
                                       interval = "simultaneous", n = 200, pos_col = NULL) {
    if (!requireNamespace("gratia", quietly = TRUE)) {
      return(list(ok = FALSE, msg = "Package 'gratia' not installed.", deriv = NULL, raw = NULL))
    }
    
    d <- tryCatch(
      gratia::derivatives(fit, term = term, order = order, interval = interval, n = n),
      error = function(e) e
    )
    if (inherits(d, "error")) {
      return(list(ok = FALSE, msg = conditionMessage(d), deriv = NULL, raw = NULL))
    }
    
    dcol <- .infer_dcol(d)
    xcol <- .infer_xcol(d, pos_col)
    
    lcol <- .pick_first_existing(d, c(".lower_ci", "lower_ci", "lower", "lo", "lwr"))
    ucol <- .pick_first_existing(d, c(".upper_ci", "upper_ci", "upper", "hi", "upr"))
    
    if (is.null(dcol) || is.null(xcol)) {
      return(list(
        ok = FALSE,
        msg = paste0(
          "Could not find derivative/x columns in gratia::derivatives() output. Found: ",
          paste(names(d), collapse = ", ")
        ),
        deriv = NULL,
        raw = d
      ))
    }
    
    out <- data.frame(
      x = d[[xcol]],
      derivative = d[[dcol]],
      lower = if (!is.null(lcol)) d[[lcol]] else NA_real_,
      upper = if (!is.null(ucol)) d[[ucol]] else NA_real_
    )
    
    list(ok = TRUE, msg = "ok", deriv = out, raw = d)
  }
  
  # ---- 1) Drop NA in BOTH position and Ks ----
  df2 <- df[!is.na(df[[ks_col]]) & !is.na(df[[pos_col]]), , drop = FALSE]
  if (nrow(df2) == 0) {
    stop("No non-NA values in columns '", pos_col, "' and '", ks_col, "' after filtering.")
  }
  
  # ---- 1b) Coerce pos and ks to numeric if possible ----
  df2[[pos_col]] <- .coerce_numeric_vec(df2[[pos_col]])
  df2[[ks_col]]  <- .coerce_numeric_vec(df2[[ks_col]])
  
  # Keep only finite points (prevents rank/qr/pathology later)
  df2 <- df2[is.finite(df2[[pos_col]]) & is.finite(df2[[ks_col]]), , drop = FALSE]
  if (nrow(df2) == 0) stop("No finite values left after coercion in '", pos_col, "' and '", ks_col, "'.")
  
  # ---- 1c) Optional y_max cutoff (TRUE exclusion from analysis; matches MCP behavior) ----
  if (!is.null(y_max)) {
    df2 <- df2[df2[[ks_col]] <= y_max, , drop = FALSE]
    if (nrow(df2) == 0) {
      stop("No points left after applying y_max cutoff (", y_max, ").")
    }
  }
  
  n <- nrow(df2)
  
  # ---- 2) Rosner outliers ----
  df2$rosner_outlier <- FALSE
  df_outliers <- df2[0, , drop = FALSE]
  ros <- NULL
  k_used <- 0
  
  max_theoretical <- n - 3
  if (max_theoretical >= 1) {
    if (!is.null(k_outliers)) {
      k <- min(k_outliers, max_theoretical)
    } else if (!is.null(k_prop)) {
      k <- floor(n * k_prop)
      if (k < 1) k <- 1
      if (k > max_theoretical) k <- max_theoretical
    } else {
      k <- min(10L, max_theoretical)
    }
    
    if (cap_rosner_k_at_10) k <- min(k, 10L)
    
    if (k > 0) {
      ros <- EnvStats::rosnerTest(df2[[ks_col]], k = k)
      outlier_idx <- which(ros$all.stats$Outlier)
      if (length(outlier_idx) > 0) df2$rosner_outlier[outlier_idx] <- TRUE
      df_outliers <- df2[df2$rosner_outlier, , drop = FALSE]
      k_used <- k
    }
  } else {
    warning("Not enough points (n = ", n, ") for rosnerTest (need at least k + 3). Skipping outlier test.")
  }
  
  message(
    "analyze_ks(): n = ", n,
    ", k_used = ", k_used,
    ", outliers_detected = ", nrow(df_outliers)
  )
  
  # ---- 3) remove outliers (optional) ----
  df_plot <- if (remove_outliers) df2[!df2$rosner_outlier, , drop = FALSE] else df2
  
  # Sanity: enough x-variation for any smoother
  if (nrow(df_plot) < 10) stop("Too few points after filtering/outlier removal to fit smoother.")
  if (length(unique(df_plot[[pos_col]])) < 2) stop("'", pos_col, "' has <2 unique values; cannot fit smoother.")
  
  # Sort by position
  df_plot <- df_plot[order(df_plot[[pos_col]], na.last = TRUE), , drop = FALSE]
  
  # ---- 4) base plot ----
  pos_sym <- rlang::sym(pos_col)
  ks_sym  <- rlang::sym(ks_col)
  
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = !!pos_sym, y = !!ks_sym)) +
    ggplot2::geom_point(size = 1)
  
  # ---- 5) smoother ----
  gam_fit <- NULL
  deriv1 <- NULL
  deriv2 <- NULL
  zc <- NULL
  high_curv <- NULL
  
  if (smooth_method == "loess") {
    
    p <- p +
      ggplot2::geom_smooth(
        method = "loess",
        span   = loess_span,
        se     = loess_se,
        method.args = list(degree = loess_degree)
      )
    
  } else {
    
    if (!requireNamespace("mgcv", quietly = TRUE)) {
      stop("smooth_method = 'gam' requires the 'mgcv' package. Please install it.")
    }
    
    # robust family choice
    fam <- switch(
      gam_family,
      scat = mgcv::scat(),
      gaussian = stats::gaussian()
    )
    
    # Fit GAM explicitly (no ggplot refit; plot always matches this model)
    # IMPORTANT: use s() (not mgcv::s()) in formula strings.
    # Also attach mgcv so s() resolves during formula evaluation.
    if (!"package:mgcv" %in% search()) {
      suppressPackageStartupMessages(library(mgcv))
    }
    
    gam_formula <- stats::as.formula(
      sprintf("%s ~ s(%s, k = %d)", ks_col, pos_col, as.integer(gam_k))
    )
    
    gam_fit <- mgcv::gam(
      formula = gam_formula,
      data    = df_plot,
      method  = gam_method,
      family  = fam
    )
    
    
    # Predict on a grid and plot from predictions (plot always matches the model)
    xseq <- seq(min(df_plot[[pos_col]], na.rm = TRUE),
                max(df_plot[[pos_col]], na.rm = TRUE),
                length.out = pred_n)
    newdata <- data.frame(xseq)
    names(newdata) <- pos_col
    
    pr <- mgcv::predict.gam(gam_fit, newdata = newdata, se.fit = gam_se, type = "response")
    
    smooth_df <- data.frame(
      x = xseq,
      fit = as.numeric(pr$fit),
      se  = if (gam_se) as.numeric(pr$se.fit) else NA_real_
    )
    names(smooth_df)[1] <- pos_col
    smooth_df <- smooth_df[order(smooth_df[[pos_col]]), , drop = FALSE]
    
    p <- p +
      ggplot2::geom_line(
        data = smooth_df,
        ggplot2::aes(x = !!rlang::sym(pos_col), y = .data$fit),
        inherit.aes = FALSE
      )
    
    if (gam_se) {
      p <- p +
        ggplot2::geom_ribbon(
          data = smooth_df,
          ggplot2::aes(
            x = !!rlang::sym(pos_col),
            ymin = .data$fit - 1.96 * .data$se,
            ymax = .data$fit + 1.96 * .data$se
          ),
          inherit.aes = FALSE,
          alpha = 0.2
        )
    }
    
    # ---- 5b) optional derivatives (never hard-stop) ----
    if (isTRUE(compute_derivatives)) {
      term_name <- paste0("s(", pos_col, ")")
      
      d1 <- .safe_gratia_derivatives(gam_fit, term = term_name, order = 1,
                                     interval = deriv_interval, n = deriv_n, pos_col = pos_col)
      if (d1$ok) {
        deriv1 <- d1$deriv
        
        # zero-crossings in first derivative
        zc <- deriv1 |>
          dplyr::mutate(sign = sign(derivative)) |>
          dplyr::mutate(change = sign != dplyr::lag(sign)) |>
          dplyr::filter(!is.na(change) & change)
        
        if (nrow(zc) > 0) {
          p <- p + ggplot2::geom_vline(
            data = zc,
            ggplot2::aes(xintercept = x),
            linetype = "dashed",
            alpha = 0.7,
            inherit.aes = FALSE
          )
        }
      } else {
        warning("GAM derivatives (order=1) skipped: ", d1$msg)
      }
      
      d2 <- .safe_gratia_derivatives(gam_fit, term = term_name, order = 2,
                                     interval = deriv_interval, n = deriv_n, pos_col = pos_col)
      if (d2$ok) {
        deriv2 <- d2$deriv
        thr <- stats::quantile(abs(deriv2$derivative), 0.95, na.rm = TRUE)
        high_curv <- deriv2 |>
          dplyr::filter(abs(derivative) >= thr)
      } else {
        warning("GAM derivatives (order=2) skipped: ", d2$msg)
      }
    }
  }
  
  # ---- 6) labels/theme ----
  p <- p +
    ggplot2::labs(x = x_label, y = expression(italic("K")[s])) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 16),
      axis.text.y = ggplot2::element_text(size = 16),
      axis.title.x = ggplot2::element_text(size = 18, vjust = -0.2),
      axis.title.y = ggplot2::element_text(size = 18, vjust = 2),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  if (!is.null(tick_by)) {
    rng <- range(df_plot[[pos_col]], na.rm = TRUE)
    if (is.finite(rng[1]) && is.finite(rng[2]) && rng[1] < rng[2]) {
      min_break <- floor(rng[1] / tick_by) * tick_by
      max_break <- ceiling(rng[2] / tick_by) * tick_by
      x_breaks  <- seq(min_break, max_break, by = tick_by)
      p <- p + ggplot2::scale_x_continuous(breaks = x_breaks)
    }
  }
  
  if (!is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(ylim = ylim)
  }
  
  list(
    data       = df_plot,
    outliers   = df_outliers,
    rosner     = ros,
    plot       = p,
    gam_fit    = gam_fit,
    deriv1     = deriv1,
    deriv2     = deriv2,
    zero_cross = zc,
    high_curv  = high_curv
  )
}


###########################################################
## fit_mcp_strata.R  (single clean file)
## - fit_mcp_strata() implementation
## - Robust Rosner outliers (optional)
## - Option 1: enforce minimum segment size + edge buffer
## - Option 3: extract cp posterior draws + credible intervals
## - Plot: points + per-stratum plateaus + cp lines (+ optional CI lines)
###########################################################

# Required packages:
# install.packages(c("ggplot2", "dplyr", "EnvStats", "mcp", "loo", "rlang"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(EnvStats)
  library(mcp)
  library(loo)
  library(rlang)
})

# ------------------------------------------------------------
# Helper: robustly extract cp posterior draws from an mcp fit
# ------------------------------------------------------------
.extract_cp_draws <- function(fit, n_cp) {
  if (n_cp <= 0) return(list(draws = NULL, names = character(0)))
  
  # ---- Route 0 (best): directly from coda mcmc.list stored by mcp ----
  if (!is.null(fit$mcmc_post) && inherits(fit$mcmc_post, "mcmc.list")) {
    # Safer than as.matrix(mcmc.list) which can yield NA colnames
    mats <- lapply(fit$mcmc_post, function(ch) {
      m <- as.matrix(ch)
      m
    })
    M <- do.call(rbind, mats)
    
    # keep only non-NA column names
    cn <- colnames(M)
    if (!is.null(cn)) {
      keep <- !is.na(cn) & nzchar(cn)
      M <- M[, keep, drop = FALSE]
      cn <- colnames(M)
    }
    
    cp_names <- grep("^cp_\\d+$", cn, value = TRUE)
    if (length(cp_names) >= 1) {
      cp_names <- cp_names[seq_len(min(length(cp_names), n_cp))]
      draws <- as.data.frame(M[, cp_names, drop = FALSE])
      return(list(draws = draws, names = cp_names))
    }
  }
  
  # ---- Route 1: as.data.frame(fit) works for some mcp versions ----
  draws_df <- NULL
  try(draws_df <- as.data.frame(fit), silent = TRUE)
  
  if (!is.null(draws_df) && nrow(draws_df) > 0) {
    cp_names <- grep("^cp_\\d+$", names(draws_df), value = TRUE)
    if (length(cp_names) >= 1) {
      cp_names <- cp_names[seq_len(min(length(cp_names), n_cp))]
      return(list(draws = draws_df[, cp_names, drop = FALSE], names = cp_names))
    }
  }
  
  # ---- Route 2: mcp::get_changepoints() (exists in some versions) ----
  if ("get_changepoints" %in% getNamespaceExports("mcp")) {
    out <- NULL
    try(out <- mcp::get_changepoints(fit), silent = TRUE)
    if (!is.null(out) && is.data.frame(out)) {
      cp_names <- grep("^cp_\\d+$", names(out), value = TRUE)
      if (length(cp_names) >= 1) {
        cp_names <- cp_names[seq_len(min(length(cp_names), n_cp))]
        return(list(draws = out[, cp_names, drop = FALSE], names = cp_names))
      }
    }
  }
  
  list(draws = NULL, names = character(0))
}

# ------------------------------------------------------------
# Helper: summarise cp draws into mean/median + CI
# ------------------------------------------------------------
.summarise_cp <- function(cp_draws, level = 0.95) {
  if (is.null(cp_draws) || ncol(cp_draws) == 0) return(NULL)
  
  alpha <- (1 - level) / 2
  out <- lapply(seq_len(ncol(cp_draws)), function(j) {
    v <- cp_draws[[j]]
    v <- v[is.finite(v)]
    if (!length(v)) return(NULL)
    data.frame(
      cp     = colnames(cp_draws)[j],
      mean   = mean(v),
      median = median(v),
      lower  = unname(quantile(v, probs = alpha)),
      upper  = unname(quantile(v, probs = 1 - alpha)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

# ------------------------------------------------------------
# Helper: build cp priors enforcing feasible cp regions
# - Ensures each segment has at least min_seg_n points
# - Optionally buffers away from chromosome edges by cp_buffer (x-units)
#
# NOTE:
#  - Truncation/ordering terms ("T(...)") were not added here,
#    because support varies across mcp backends/versions.
#  - In practice, the segment-size index bounds already reduce pathologies.
# ------------------------------------------------------------
.build_cp_priors <- function(x, n_cp, min_seg_n = 50, cp_buffer = NULL) {
  if (n_cp <= 0) return(list())
  
  x <- sort(as.numeric(x))
  x <- x[is.finite(x)]
  n <- length(x)
  
  if (n < (n_cp + 1) * min_seg_n) {
    stop(
      "min_seg_n too large: need at least (n_cp+1)*min_seg_n points.\n",
      "n = ", n, ", n_cp = ", n_cp, ", min_seg_n = ", min_seg_n
    )
  }
  
  buf <- if (is.null(cp_buffer)) 0 else as.numeric(cp_buffer)
  
  pri <- list()
  for (j in seq_len(n_cp)) {
    # Index-based feasible region for cp_j:
    # need >= min_seg_n points in each segment across ALL cps.
    lower_j <- x[min_seg_n * j + 1] + buf
    upper_j <- x[n - min_seg_n * (n_cp - j + 1)] - buf
    
    if (!is.finite(lower_j) || !is.finite(upper_j) || lower_j >= upper_j) {
      stop(
        "No feasible prior region for cp_", j,
        " after applying min_seg_n and cp_buffer.\n",
        "Try reducing min_seg_n and/or cp_buffer."
      )
    }
    
    pri[[paste0("cp_", j)]] <- sprintf("dunif(%f, %f)", lower_j, upper_j)
  }
  
  pri
}

# ------------------------------------------------------------
# MAIN: fit Bayesian changepoint models (mcp) and assign strata
# ------------------------------------------------------------
fit_mcp_strata <- function(df,
                           x_col = "start_mb",
                           y_col = "Ks",
                           max_cp = 3,
                           iter = 100000,
                           adapt = 20000,
                           chains = 4,
                           cores = 4,
                           y_max = NULL,            # TRUE exclusion from analysis
                           ylim = NULL,             # plotting only
                           rosner = FALSE,
                           k_outliers = 15,
                           k_prop = NULL,
                           remove_outliers = TRUE,
                           tick_by = NULL,
                           x_label_custom = NULL,
                           
                           # Option 1: constraints
                           enforce_min_segment = TRUE,
                           min_seg_n = 50,          # minimum points per stratum
                           cp_buffer = NULL,        # x-units (e.g., 0.25 Mb)
                           
                           # Option 3: cp credible intervals / plot controls
                           cp_ci_level = 0.95,
                           cp_line_stat = c("mean", "median"),
                           add_cp_ci_bars = FALSE
) {
  cp_line_stat <- match.arg(cp_line_stat)
  
  # ----------------------------
  # 1) Validate + clean
  # ----------------------------
  if (!x_col %in% names(df)) stop("Column '", x_col, "' not found in df.")
  if (!y_col %in% names(df)) stop("Column '", y_col, "' not found in df.")
  
  df2 <- df %>%
    dplyr::select(all_of(c(x_col, y_col)), dplyr::everything()) %>%
    dplyr::filter(!is.na(.data[[x_col]]), !is.na(.data[[y_col]]))
  
  # coerce numeric x
  if (!is.numeric(df2[[x_col]])) {
    suppressWarnings({ df2[[x_col]] <- as.numeric(df2[[x_col]]) })
    if (anyNA(df2[[x_col]])) stop("x_col could not be coerced to numeric cleanly.")
  }
  
  # optional saturation filter (TRUE exclusion from analysis)
  if (!is.null(y_max)) {
    df2 <- df2 %>% dplyr::filter(.data[[y_col]] <= y_max)
  }
  
  if (nrow(df2) < 10) stop("Not enough points after filtering to run changepoint analysis.")
  df2 <- df2 %>% arrange(.data[[x_col]])
  
  # ----------------------------
  # 2) Optional Rosner outliers
  # ----------------------------
  df2$rosner_outlier <- FALSE
  ros <- NULL
  k_used <- 0
  
  if (rosner) {
    n <- nrow(df2)
    max_theoretical <- n - 3
    
    if (max_theoretical >= 1) {
      if (!is.null(k_outliers)) {
        k <- min(k_outliers, max_theoretical)
      } else if (!is.null(k_prop)) {
        k <- floor(n * k_prop)
        if (k < 1) k <- 1
        if (k > max_theoretical) k <- max_theoretical
      } else {
        k <- min(10L, max_theoretical)
      }
      
      if (k > 0) {
        ros <- EnvStats::rosnerTest(df2[[y_col]], k = k)
        outlier_idx <- which(ros$all.stats$Outlier)
        if (length(outlier_idx) > 0) df2$rosner_outlier[outlier_idx] <- TRUE
        k_used <- k
      }
    }
    
    message(
      "Rosner outlier detection: n = ", nrow(df2),
      ", k_used = ", k_used,
      ", outliers_detected = ", sum(df2$rosner_outlier)
    )
  }
  
  df_fit <- if (rosner && remove_outliers) df2[!df2$rosner_outlier, , drop = FALSE] else df2
  if (nrow(df_fit) < 10) stop("Not enough points to fit mcp model after outlier handling.")
  
  # ----------------------------
  # 3) Build models: 0..max_cp
  # ----------------------------
  if (max_cp < 0 || max_cp > 3) stop("max_cp must be between 0 and 3.")
  
  build_model <- function(n_cp) {
    model_list <- vector("list", n_cp + 1)
    model_list[[1]] <- as.formula(paste(y_col, "~ 1"))
    if (n_cp > 0) for (i in 2:(n_cp + 1)) model_list[[i]] <- ~ 1
    model_list
  }
  
  cp_values <- 0:max_cp
  fits <- vector("list", length(cp_values))
  loo_list <- vector("list", length(cp_values))
  priors_used <- vector("list", length(cp_values))
  
  xvec_fit <- df_fit[[x_col]]
  
  for (i in seq_along(cp_values)) {
    n_cp <- cp_values[i]
    message("Fitting model with ", n_cp, " changepoint(s)...")
    
    model_i <- build_model(n_cp)
    
    # IMPORTANT: prior must ALWAYS be a list (never NULL)
    pri_i <- list()
    
    # Only add cp priors when we want constraints and n_cp > 0
    if (enforce_min_segment && n_cp > 0) {
      pri_i <- .build_cp_priors(
        x = xvec_fit,
        n_cp = n_cp,
        min_seg_n = min_seg_n,
        cp_buffer = cp_buffer
      )
    }
    
    priors_used[[i]] <- pri_i
    
    fit_i <- mcp::mcp(
      model  = model_i,
      data   = df_fit,
      par_x  = x_col,
      iter   = iter,
      adapt  = adapt,
      chains = chains,
      cores  = cores,
      prior  = pri_i
    )
    
    fits[[i]] <- fit_i
    loo_list[[i]] <- loo::loo(fit_i)
  }
  
  # ----------------------------
  # 4) Model comparison via LOO
  # ----------------------------
  loo_compare_tbl <- loo::loo_compare(loo_list)
  weights <- loo::loo_model_weights(loo_list, method = "pseudobma")
  
  wvec <- if (is.matrix(weights) || is.data.frame(weights)) as.numeric(weights[, 1]) else as.numeric(weights)
  
  best_index <- which.max(wvec)         # 1 = 0 cp, 2 = 1 cp, ...
  best_n_cp  <- cp_values[best_index]
  best_fit   <- fits[[best_index]]
  
  message("Best model has ", best_n_cp, " changepoint(s) (", best_n_cp + 1, " strata).")
  
  # ----------------------------
  # 5) Extract CPs (draws + summary)
  # ----------------------------
  cp_draws_obj <- .extract_cp_draws(best_fit, n_cp = best_n_cp)
  cp_draws <- cp_draws_obj$draws
  cp_summary <- .summarise_cp(cp_draws, level = cp_ci_level)
  
  change_points <- numeric(0)
  if (!is.null(cp_summary) && nrow(cp_summary) > 0) {
    change_points <- if (cp_line_stat == "mean") cp_summary$mean else cp_summary$median
    change_points <- sort(as.numeric(change_points))
  }
  
  message(
    "Extracted change_points (", cp_line_stat, "): ",
    if (length(change_points) == 0) "none" else paste(sprintf("%.3f", change_points), collapse = ", ")
  )
  
  # ----------------------------
  # 6) Assign strata to ALL points (df2)
  # ----------------------------
  if (length(change_points) > 0) {
    breaks <- c(-Inf, change_points, Inf)
    labels <- paste0("strata", seq_len(length(breaks) - 1))
    df2$strata_best <- cut(
      df2[[x_col]],
      breaks = breaks,
      labels = labels,
      include.lowest = TRUE,
      right = TRUE
    )
  } else {
    df2$strata_best <- factor("strata1")
  }
  
  # Plot data: may omit outliers if requested
  df_plot <- if (rosner && remove_outliers) df2[!df2$rosner_outlier, , drop = FALSE] else df2
  
  # ----------------------------
  # 7) Plateau segments (empirical mean per stratum)
  # ----------------------------
  seg_df <- df_plot %>%
    group_by(strata_best) %>%
    summarise(
      x_min  = min(.data[[x_col]]),
      x_max  = max(.data[[x_col]]),
      y_mean = mean(.data[[y_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  # ----------------------------
  # 8) Plot
  # ----------------------------
  mapping <- aes(x = !!sym(x_col), y = !!sym(y_col), color = strata_best)
  
  points_layer <- geom_point(size = 1)
  
  segments_layer <- geom_segment(
    data = seg_df,
    aes(x = x_min, xend = x_max, y = y_mean, yend = y_mean, color = strata_best),
    linewidth = 1,
    inherit.aes = FALSE
  )
  
  vlines_layer <- if (length(change_points) > 0) {
    geom_vline(xintercept = change_points, linetype = "dashed")
  } else {
    NULL
  }
  
  ci_vlines_layer <- NULL
  if (add_cp_ci_bars && !is.null(cp_summary) && nrow(cp_summary) > 0) {
    ci_vlines_layer <- list(
      geom_vline(xintercept = cp_summary$lower, linetype = "dotted", alpha = 0.35),
      geom_vline(xintercept = cp_summary$upper, linetype = "dotted", alpha = 0.35)
    )
  }
  
  color_scale <- scale_color_brewer(palette = "Set1")
  
  base_theme <- theme_bw() +
    theme(
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      panel.grid.minor = element_blank()
    )
  
  x_scale_layer <- NULL
  if (!is.null(tick_by)) {
    rng <- range(df_plot[[x_col]], na.rm = TRUE)
    if (is.finite(rng[1]) && is.finite(rng[2]) && rng[1] < rng[2]) {
      min_break <- floor(rng[1] / tick_by) * tick_by
      max_break <- ceiling(rng[2] / tick_by) * tick_by
      x_scale_layer <- scale_x_continuous(breaks = seq(min_break, max_break, by = tick_by))
    }
  }
  
  x_lab <- if (!is.null(x_label_custom)) x_label_custom else x_col
  
  p <- ggplot(df_plot, mapping) +
    points_layer +
    segments_layer +
    { if (!is.null(vlines_layer)) vlines_layer else NULL } +
    { if (!is.null(ci_vlines_layer)) ci_vlines_layer else NULL } +
    color_scale +
    { if (!is.null(x_scale_layer)) x_scale_layer else NULL } +
    labs(x = x_lab, y = y_col, color = "Stratum") +
    base_theme
  
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  
  # ----------------------------
  # 9) Return
  # ----------------------------
  list(
    data_full            = df2,          # all points with rosner_outlier + strata_best
    data                 = df2,          # alias
    data_plot            = df_plot,      # plotted points (may exclude outliers)
    seg_df               = seg_df,
    fits                 = fits,
    loo_list             = loo_list,
    loo_compare          = loo_compare_tbl,
    loo_weights          = weights,
    best_index           = best_index,
    best_n_changepoints  = best_n_cp,
    change_points        = change_points,
    cp_draws             = cp_draws,
    cp_summary           = cp_summary,
    priors_used          = priors_used,
    rosner_result        = ros,
    rosner_k_used        = k_used,
    plot                 = p,
    layers = list(
      mapping      = mapping,
      points       = points_layer,
      segments     = segments_layer,
      vlines       = vlines_layer,
      ci_vlines    = ci_vlines_layer,
      color_scale  = color_scale,
      base_theme   = base_theme,
      x_scale      = x_scale_layer
    )
  )
}

# ------------------------------------------------------------
# Posterior fit + changepoint density panels (best model)
# ------------------------------------------------------------
plot_mcp_cp_posterior <- function(res_mcp, cp_pars = "cp_1", ...) {
  if (is.null(res_mcp$fits) || is.null(res_mcp$best_index)) {
    stop("res_mcp must come from fit_mcp_strata() and contain $fits and $best_index.")
  }
  
  best_fit <- res_mcp$fits[[res_mcp$best_index]]
  
  # posterior fit along x (mcp plot method)
  fit_plot <- plot(best_fit, q_fit = TRUE)
  
  cp_pars <- as.character(cp_pars)
  cp_plots <- vector("list", length(cp_pars))
  names(cp_plots) <- cp_pars
  
  for (i in seq_along(cp_pars)) {
    cp_plots[[i]] <- mcp::plot_pars(best_fit, pars = cp_pars[i], ...)
  }
  
  list(fit_plot = fit_plot, cp_plots = cp_plots)
}

