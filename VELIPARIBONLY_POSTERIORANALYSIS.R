############################################################
# STAGE 3: CONTROL + DRUG-ONLY POSTERIOR ANALYSIS
############################################################

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(loo)
  library(bayesplot)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

set.seed(123)

stage3_root <- path.expand("~/Desktop/STADIO3/results")

model_dirs <- c(
  exponential = "stage3_control_drug_closed_form_exponential",
  logistic    = "stage3_control_drug_closed_form_logistic",
  gompertz    = "stage3_control_drug_closed_form_gompertz"
)

fit_names <- c(
  exponential = "fit_stage3_control_drug_exponential_full.rds",
  logistic    = "fit_stage3_control_drug_logistic_full.rds",
  gompertz    = "fit_stage3_control_drug_gompertz_full.rds"
)

fit_paths <- file.path(stage3_root, unname(model_dirs), unname(fit_names))
names(fit_paths) <- names(model_dirs)
obs_file <- file.path(
  stage3_root,
  model_dirs["exponential"],
  "control_drug_observations.csv"
)

stopifnot(all(file.exists(fit_paths)), file.exists(obs_file))

analysis_dir <- stage3_root
table_dir <- file.path(stage3_root, "tables")
figure_dir <- file.path(stage3_root, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

fits <- list(
  exponential = readRDS(fit_paths["exponential"]),
  logistic    = readRDS(fit_paths["logistic"]),
  gompertz    = readRDS(fit_paths["gompertz"])
)

stopifnot(all(vapply(fits, inherits, logical(1), what = "CmdStanMCMC")))

obs <- read.csv(obs_file, stringsAsFactors = FALSE) %>%
  arrange(mouse_id, time_days)

obs$group <- factor(obs$group, levels = c("Control", "Drug-only"))

stopifnot(
  nrow(obs) == 178L,
  length(unique(obs$mouse_id)) == 36L,
  sum(obs$is_drug == 0L) == 114L,
  sum(obs$is_drug == 1L) == 64L
)

model_label <- function(m) {
  switch(m,
         exponential = "Exponential",
         logistic = "Logistic",
         gompertz = "Gompertz")
}

key_vars <- function(m) {
  if (m == "exponential") {
    c("logV0_pop", "logr_pop", "sigma_logV0", "sigma_logr", "sigma")
  } else {
    c("logV0_pop", "logr_pop", "logK_pop",
      "sigma_logV0", "sigma_logr", "sigma_logK", "sigma")
  }
}

save_plot <- function(p, name, width = 8, height = 6) {
  ggsave(file.path(figure_dir, paste0(name, ".pdf")), p,
         width = width, height = height)
  ggsave(file.path(figure_dir, paste0(name, ".png")), p,
         width = width, height = height, dpi = 300)
}

error_summary <- function(d) {
  d %>% summarise(
    n = n(),
    mean_log_residual = mean(residual_log),
    sd_log_residual = sd(residual_log),
    mae_log_residual = mean(abs(residual_log)),
    rmse_log_residual = sqrt(mean(residual_log^2)),
    max_abs_log_residual = max(abs(residual_log)),
    mean_standardised_residual = mean(bayes_std_residual),
    sd_standardised_residual = sd(bayes_std_residual),
    predictive_coverage_95 = mean(pp_covered_95)
  )
}

hmc <- bind_rows(lapply(names(fits), function(m) {
  d <- fits[[m]]$diagnostic_summary()
  data.frame(
    model = model_label(m),
    chain = seq_along(d$num_divergent),
    divergences = d$num_divergent,
    max_treedepth_hits = d$num_max_treedepth,
    ebfmi = d$ebfmi
  )
}))

write.csv(hmc, file.path(table_dir, "hmc_diagnostics_by_chain.csv"), row.names = FALSE)

hmc_summary <- hmc %>%
  group_by(model) %>%
  summarise(
    total_divergences = sum(divergences),
    total_max_treedepth_hits = sum(max_treedepth_hits),
    min_ebfmi = min(ebfmi),
    max_ebfmi = max(ebfmi),
    .groups = "drop"
  )
write.csv(hmc_summary, file.path(table_dir, "hmc_diagnostics_model_summary.csv"), row.names = FALSE)

conv <- bind_rows(lapply(names(fits), function(m) {
  as.data.frame(fits[[m]]$summary(variables = key_vars(m))) %>%
    mutate(model = model_label(m), .before = 1)
}))

write.csv(conv, file.path(table_dir, "posterior_key_parameter_summary_all_models.csv"), row.names = FALSE)

conv_summary <- conv %>%
  group_by(model) %>%
  summarise(
    max_rhat = max(rhat, na.rm = TRUE),
    min_bulk_ess = min(ess_bulk, na.rm = TRUE),
    min_tail_ess = min(ess_tail, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(conv_summary, file.path(table_dir, "convergence_summary_all_models.csv"), row.names = FALSE)


for (m in names(fits)) {
  vars <- key_vars(m)
  darr <- posterior::as_draws_array(fits[[m]]$draws(variables = vars))
  
  p_trace <- bayesplot::mcmc_trace(darr, pars = vars,
                                   facet_args = list(ncol = 1, strip.position = "left")) +
    ggtitle(paste0("Stage 3 ", model_label(m), ": trace plots"))
  save_plot(p_trace, paste0("trace_", m), 10, ifelse(m == "exponential", 9, 12))
  
  p_post <- bayesplot::mcmc_areas(darr, pars = vars, prob = 0.8, prob_outer = 0.95) +
    ggtitle(paste0("Stage 3 ", model_label(m), ": posterior distributions"))
  save_plot(p_post, paste0("posterior_parameters_", m), 9, ifelse(m == "exponential", 6, 8))
}

n_prior <- 50000

prior_draws <- function(m) {
  r0 <- switch(m, exponential = 0.3, logistic = 0.4, gompertz = 0.2)
  x <- data.frame(
    logV0_pop = rnorm(n_prior, log(30), 0.5),
    logr_pop = rnorm(n_prior, log(r0), 0.5),
    sigma_logV0 = abs(rnorm(n_prior, 0, 0.5)),
    sigma_logr = abs(rnorm(n_prior, 0, 0.5)),
    sigma = abs(rnorm(n_prior, 0, 0.5))
  )
  if (m != "exponential") {
    x$logK_pop <- rnorm(n_prior, log(500), 0.7)
    x$sigma_logK <- abs(rnorm(n_prior, 0, 0.5))
  }
  x
}

for (m in names(fits)) {
  vars <- key_vars(m)
  
  pr <- prior_draws(m) %>%
    select(all_of(vars)) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    mutate(distribution = "Prior")
  
  po <- as.data.frame(posterior::as_draws_matrix(fits[[m]]$draws(variables = vars))) %>%
    select(all_of(vars)) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    mutate(distribution = "Posterior")
  
  pp <- bind_rows(pr, po)
  
  p <- ggplot(pp, aes(x = value, linetype = distribution, group = distribution)) +
    geom_density(linewidth = 0.85, adjust = 1.1) +
    facet_wrap(~ parameter, scales = "free", ncol = 2) +
    labs(
      title = paste0("Stage 3 ", model_label(m), ": prior vs posterior"),
      subtitle = "Population and between-mouse parameters",
      x = "Parameter value", y = "Density", linetype = NULL
    ) +
    theme_bw() +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
  
  save_plot(p, paste0("prior_vs_posterior_", m), 10, ifelse(m == "exponential", 8, 10))
}


analysis_list <- list()

for (m in names(fits)) {
  fit <- fits[[m]]
  
  mu <- posterior::as_draws_matrix(fit$draws("mu_out"))
  mu_log <- posterior::as_draws_matrix(fit$draws("mu_log_out"))
  yrep <- posterior::as_draws_matrix(fit$draws("y_rep"))
  sig <- as.numeric(posterior::as_draws_matrix(fit$draws("sigma"))[, 1])
  
  stopifnot(ncol(mu) == nrow(obs), ncol(mu_log) == nrow(obs), ncol(yrep) == nrow(obs))
  
  fitted_mean <- colMeans(mu)
  fitted_q025 <- apply(mu, 2, quantile, 0.025)
  fitted_q975 <- apply(mu, 2, quantile, 0.975)
  
  pp_q025 <- apply(yrep, 2, quantile, 0.025)
  pp_q50  <- apply(yrep, 2, quantile, 0.50)
  pp_q975 <- apply(yrep, 2, quantile, 0.975)
  
  numerator <- matrix(log(obs$volume), nrow = nrow(mu_log), ncol = ncol(mu_log), byrow = TRUE) - mu_log
  z_draws <- sweep(numerator, 1, sig, "/")
  
  d <- obs %>% mutate(
    model = model_label(m),
    fitted_mean = fitted_mean,
    fitted_q025 = fitted_q025,
    fitted_q975 = fitted_q975,
    residual_log = log(volume) - log(fitted_mean),
    bayes_std_residual = colMeans(z_draws),
    pp_q025 = pp_q025,
    pp_median = pp_q50,
    pp_q975 = pp_q975,
    pp_covered_95 = volume >= pp_q025 & volume <= pp_q975
  )
  
  analysis_list[[m]] <- d
  write.csv(d, file.path(table_dir, paste0("observation_diagnostics_", m, ".csv")), row.names = FALSE)
}

analysis_all <- bind_rows(analysis_list)
write.csv(analysis_all, file.path(table_dir, "observation_diagnostics_all_models.csv"), row.names = FALSE)


residual_overall <- analysis_all %>% group_by(model) %>% error_summary()
residual_by_group <- analysis_all %>% group_by(model, group) %>% error_summary()

write.csv(residual_overall, file.path(table_dir, "residual_summary_all_models.csv"), row.names = FALSE)
write.csv(residual_by_group, file.path(table_dir, "residual_summary_by_group.csv"), row.names = FALSE)


for (m in names(analysis_list)) {
  d <- analysis_list[[m]]
  ttl <- paste0("Stage 3 ", model_label(m), ": ")
  
  p1 <- ggplot(d, aes(fitted_mean, volume)) +
    geom_point(alpha = 0.7) + geom_abline(slope = 1, intercept = 0, linetype = 2) +
    scale_x_log10() + scale_y_log10() + theme_bw() +
    labs(title = paste0(ttl, "observed vs fitted"),
         x = "Posterior mean fitted volume (mm³, log scale)",
         y = "Observed volume (mm³, log scale)")
  save_plot(p1, paste0("observed_vs_fitted_", m))
  
  p2 <- ggplot(d, aes(fitted_mean, bayes_std_residual)) +
    geom_point(alpha = 0.7) + geom_hline(yintercept = 0, linetype = 2) +
    scale_x_log10() + theme_bw() +
    labs(title = paste0(ttl, "Bayesian standardised residuals vs fitted"),
         x = "Posterior mean fitted volume (mm³, log scale)", y = "Standardised residual")
  save_plot(p2, paste0("standardised_residuals_vs_fitted_", m))
  
  p3 <- ggplot(d, aes(time_days, bayes_std_residual)) +
    geom_point(alpha = 0.7) + geom_hline(yintercept = 0, linetype = 2) + theme_bw() +
    labs(title = paste0(ttl, "Bayesian standardised residuals vs time"),
         x = "Time (days)", y = "Standardised residual")
  save_plot(p3, paste0("standardised_residuals_vs_time_", m))
  
  p4 <- ggplot(d, aes(bayes_std_residual)) +
    geom_histogram(bins = 25, boundary = 0) + theme_bw() +
    labs(title = paste0(ttl, "standardised residual histogram"),
         x = "Standardised residual", y = "Count")
  save_plot(p4, paste0("residual_histogram_", m))
  
  p5 <- ggplot(d, aes(sample = bayes_std_residual)) +
    stat_qq() + stat_qq_line() + theme_bw() +
    labs(title = paste0(ttl, "Normal Q-Q plot"),
         x = "Theoretical Normal quantiles", y = "Observed residual quantiles")
  save_plot(p5, paste0("residual_qq_", m))
  
  p6 <- ggplot(d, aes(group, bayes_std_residual)) +
    geom_violin(trim = FALSE, alpha = 0.3) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = 2) + theme_bw() +
    labs(title = paste0(ttl, "residuals by cohort"),
         subtitle = "Consistency diagnostic; not an estimate of a direct drug effect",
         x = NULL, y = "Standardised residual")
  save_plot(p6, paste0("residual_distribution_by_group_", m))
}


for (m in names(analysis_list)) {
  d <- analysis_list[[m]]
  p <- ggplot(d, aes(time_days, volume)) +
    geom_ribbon(aes(ymin = fitted_q025, ymax = fitted_q975, group = mouse_id), alpha = 0.2) +
    geom_line(aes(y = fitted_mean, group = mouse_id), linewidth = 0.6) +
    geom_point(size = 1.1) +
    facet_wrap(~ mouse, scales = "free_y", ncol = 6) + theme_bw() +
    labs(title = paste0("Stage 3 ", model_label(m), ": mouse-specific fitted trajectories"),
         subtitle = "Ribbon = 95% interval for latent fitted mean",
         x = "Time (days)", y = "Tumour volume (mm³)") +
    theme(strip.text = element_text(size = 7))
  save_plot(p, paste0("mouse_specific_fits_", m), 14, 11)
}

n_ppc <- 50
for (m in names(fits)) {
  yrep <- posterior::as_draws_matrix(fits[[m]]$draws("y_rep"))
  draw_id <- sample(seq_len(nrow(yrep)), min(n_ppc, nrow(yrep)))
  
  p <- bayesplot::ppc_dens_overlay(log(obs$volume), log(yrep[draw_id, , drop = FALSE])) +
    ggtitle(paste0("Stage 3 ", model_label(m), ": posterior predictive density check"),
            subtitle = "Log tumour volume; replicated observations from y_rep")
  save_plot(p, paste0("ppc_density_pooled_", m))
  
  for (g in levels(obs$group)) {
    j <- which(obs$group == g)
    pg <- bayesplot::ppc_dens_overlay(log(obs$volume[j]), log(yrep[draw_id, j, drop = FALSE])) +
      ggtitle(paste0("Stage 3 ", model_label(m), ": PPC — ", g),
              subtitle = "Group-wise consistency check")
    save_plot(pg, paste0("ppc_density_", gsub("-", "_", tolower(g)), "_", m))
  }
}

coverage_overall <- analysis_all %>%
  group_by(model) %>%
  summarise(n = n(), n_covered = sum(pp_covered_95), coverage = mean(pp_covered_95), .groups = "drop")

coverage_by_group <- analysis_all %>%
  group_by(model, group) %>%
  summarise(n = n(), n_covered = sum(pp_covered_95), coverage = mean(pp_covered_95), .groups = "drop")

outside_ppi <- analysis_all %>%
  filter(!pp_covered_95) %>%
  select(model, mouse, mouse_id, group, time_days, volume, pp_q025, pp_median, pp_q975)

write.csv(coverage_overall, file.path(table_dir, "posterior_predictive_coverage_overall.csv"), row.names = FALSE)
write.csv(coverage_by_group, file.path(table_dir, "posterior_predictive_coverage_by_group.csv"), row.names = FALSE)
write.csv(outside_ppi, file.path(table_dir, "observations_outside_95pct_predictive_interval.csv"), row.names = FALSE)


breaks <- unique(quantile(obs$time_days, seq(0, 1, length.out = 5), na.rm = TRUE))
if (length(breaks) < 5) breaks <- seq(min(obs$time_days), max(obs$time_days), length.out = 5)
obs$time_bin <- cut(obs$time_days, breaks = breaks, include.lowest = TRUE, ordered_result = TRUE)

timebin_list <- list()
for (m in names(fits)) {
  yrep <- posterior::as_draws_matrix(fits[[m]]$draws("y_rep"))
  
  tab <- bind_rows(lapply(levels(obs$time_bin), function(b) {
    j <- which(obs$time_bin == b)
    rep_stat <- rowMeans(log(yrep[, j, drop = FALSE]))
    data.frame(
      model = model_label(m), time_bin = b, n_obs = length(j),
      observed_mean_log_volume = mean(log(obs$volume[j])),
      predictive_q025 = quantile(rep_stat, 0.025),
      predictive_median = quantile(rep_stat, 0.50),
      predictive_q975 = quantile(rep_stat, 0.975)
    )
  }))
  timebin_list[[m]] <- tab
  
  p <- ggplot(tab, aes(time_bin, predictive_median)) +
    geom_errorbar(aes(ymin = predictive_q025, ymax = predictive_q975), width = 0.15) +
    geom_point(size = 2) +
    geom_point(aes(y = observed_mean_log_volume), shape = 4, size = 3, stroke = 1) +
    theme_bw() +
    labs(title = paste0("Stage 3 ", model_label(m), ": time-binned PPC"),
         subtitle = "Cross = observed mean; interval = posterior predictive distribution",
         x = "Time bin (days)", y = "Mean log tumour volume") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_plot(p, paste0("ppc_time_binned_", m))
}

write.csv(bind_rows(timebin_list), file.path(table_dir, "time_binned_ppc_summary.csv"), row.names = FALSE)


loos <- list()
pareto_list <- list()

for (m in names(fits)) {
  ll <- posterior::as_draws_matrix(fits[[m]]$draws("log_lik"))
  loos[[m]] <- loo::loo(ll)
  saveRDS(loos[[m]], file.path(table_dir, paste0("loo_", m, ".rds")))
  
  k <- loo::pareto_k_values(loos[[m]])
  pareto_list[[m]] <- obs %>% mutate(
    model = model_label(m), obs_index = seq_len(n()), pareto_k = k,
    pareto_flag = case_when(pareto_k > 1 ~ ">1",
                            pareto_k > 0.7 ~ "0.7-1",
                            TRUE ~ "<=0.7")
  )
}

loo_summary <- bind_rows(lapply(names(loos), function(m) {
  x <- loos[[m]]$estimates
  data.frame(
    model = model_label(m),
    elpd_loo = x["elpd_loo", "Estimate"],
    se_elpd_loo = x["elpd_loo", "SE"],
    p_loo = x["p_loo", "Estimate"],
    looic = x["looic", "Estimate"]
  )
}))
write.csv(loo_summary, file.path(table_dir, "loo_summary_all_models.csv"), row.names = FALSE)

loo_cmp <- loo::loo_compare(list(
  exponential = loos$exponential,
  logistic = loos$logistic,
  gompertz = loos$gompertz
))
loo_cmp_df <- data.frame(model = rownames(loo_cmp), as.data.frame(loo_cmp), row.names = NULL)
write.csv(loo_cmp_df, file.path(table_dir, "loo_paired_comparison.csv"), row.names = FALSE)


pareto_all <- bind_rows(pareto_list)
problematic <- pareto_all %>% filter(pareto_k > 0.7) %>% arrange(model, desc(pareto_k))
pareto_summary <- pareto_all %>%
  group_by(model) %>%
  summarise(n_observations = n(), n_k_gt_0_7 = sum(pareto_k > 0.7),
            n_k_gt_1 = sum(pareto_k > 1), max_pareto_k = max(pareto_k), .groups = "drop")

write.csv(pareto_all, file.path(table_dir, "pareto_k_all_observations.csv"), row.names = FALSE)
write.csv(problematic, file.path(table_dir, "pareto_k_problematic_observations.csv"), row.names = FALSE)
write.csv(pareto_summary, file.path(table_dir, "pareto_k_summary.csv"), row.names = FALSE)

for (m in names(pareto_list)) {
  d <- pareto_list[[m]]
  p <- ggplot(d, aes(obs_index, pareto_k)) +
    geom_point(alpha = 0.75) +
    geom_hline(yintercept = 0.7, linetype = 2) +
    geom_hline(yintercept = 1.0, linetype = 3) +
    theme_bw() +
    labs(title = paste0("Stage 3 ", model_label(m), ": observation-level Pareto-k"),
         subtitle = "Reference lines at k = 0.7 and k = 1",
         x = "Observation index", y = "Pareto-k")
  save_plot(p, paste0("pareto_k_", m), 9, 6)
}

############################################################
# FOR Mechanistic consistency record
############################################################

mechanistic_summary <- data.frame(
  item = c(
    "Control mice", "Drug-only mice", "RT administrations in Stage 3",
    "Initial DNA damage D(0)", "Model-implied DNA damage without RT",
    "Direct drug effect estimated in Stage 3"
  ),
  value = c("23", "13", "None", "0", "D(t) = 0 exactly", "No")
)
write.csv(mechanistic_summary, file.path(table_dir, "stage3_mechanistic_consistency_summary.csv"), row.names = FALSE)

############################################################
# Stage 1 vs Stage 3 comparison
############################################################

stage1_fit_paths <- c(
  exponential = path.expand("~/Desktop/dissertation/results/fit_exp_control.rds"),
  logistic    = path.expand("~/Desktop/dissertation/results/fit_logistic_control.rds"),
  gompertz    = path.expand("~/Desktop/dissertation/results/fit_gompertz_control.rds")
)

stopifnot(all(file.exists(stage1_fit_paths)))

stage1_fits <- lapply(stage1_fit_paths, readRDS)
stopifnot(all(vapply(stage1_fits, inherits, logical(1), what = "CmdStanMCMC")))

for (m in c("logistic", "gompertz")) {
  mp <- stage1_fits[[m]]$metadata()$model_params
  if (!("logK_pop" %in% mp) || ("logK_extra_pop" %in% mp)) {
    stop(paste0(
      "Stage 1 ", m,
      " fit does not match the final K hierarchy (expected logK_pop and no logK_extra_pop)."
    ))
  }
}

stage13_vars <- function(m) {
  if (m == "exponential") {
    c("logV0_pop", "logr_pop", "sigma_logV0", "sigma_logr", "sigma")
  } else {
    c("logV0_pop", "logr_pop", "logK_pop",
      "sigma_logV0", "sigma_logr", "sigma_logK", "sigma")
  }
}

stage13_natural_vars <- function(m) {
  if (m == "exponential") c("V0_pop", "r_pop") else c("V0_pop", "r_pop", "K_pop")
}

extract_stage13_draws <- function(fit, m, stage_label) {
  vars <- stage13_vars(m)
  d <- as.data.frame(posterior::as_draws_matrix(fit$draws(variables = vars))) %>%
    select(all_of(vars))
  
  d$V0_pop <- exp(d$logV0_pop)
  d$r_pop  <- exp(d$logr_pop)
  if (m != "exponential") d$K_pop <- exp(d$logK_pop)
  d$stage <- stage_label
  d$model <- model_label(m)
  d
}

stage13_draws <- lapply(names(fits), function(m) {
  bind_rows(
    extract_stage13_draws(stage1_fits[[m]], m, "Stage 1: Control only"),
    extract_stage13_draws(fits[[m]],        m, "Stage 3: Control + Drug-only")
  )
})
names(stage13_draws) <- names(fits)

stage13_model_scale_summary <- bind_rows(lapply(names(fits), function(m) {
  d <- stage13_draws[[m]]
  vars <- stage13_vars(m)
  bind_rows(lapply(vars, function(v) {
    d %>%
      group_by(stage) %>%
      summarise(
        model = model_label(m),
        parameter = v,
        mean = mean(.data[[v]]),
        sd = sd(.data[[v]]),
        median = median(.data[[v]]),
        q025 = quantile(.data[[v]], 0.025),
        q975 = quantile(.data[[v]], 0.975),
        .groups = "drop"
      ) %>%
      select(model, parameter, stage, everything())
  }))
}))
write.csv(
  stage13_model_scale_summary,
  file.path(table_dir, "stage1_vs_stage3_posterior_summary_model_scale.csv"),
  row.names = FALSE
)

stage13_natural_summary <- bind_rows(lapply(names(fits), function(m) {
  d <- stage13_draws[[m]]
  vars <- stage13_natural_vars(m)
  bind_rows(lapply(vars, function(v) {
    d %>%
      group_by(stage) %>%
      summarise(
        model = model_label(m),
        parameter = v,
        mean = mean(.data[[v]]),
        sd = sd(.data[[v]]),
        median = median(.data[[v]]),
        q025 = quantile(.data[[v]], 0.025),
        q975 = quantile(.data[[v]], 0.975),
        .groups = "drop"
      ) %>%
      select(model, parameter, stage, everything())
  }))
}))
write.csv(
  stage13_natural_summary,
  file.path(table_dir, "stage1_vs_stage3_posterior_summary_natural_scale.csv"),
  row.names = FALSE
)

stage13_shift_summary <- stage13_natural_summary %>%
  select(model, parameter, stage, mean, median, q025, q975) %>%
  pivot_wider(
    names_from = stage,
    values_from = c(mean, median, q025, q975),
    names_sep = "__"
  ) %>%
  mutate(
    mean_shift_stage3_minus_stage1 =
      .data[["mean__Stage 3: Control + Drug-only"]] - .data[["mean__Stage 1: Control only"]],
    median_shift_stage3_minus_stage1 =
      .data[["median__Stage 3: Control + Drug-only"]] - .data[["median__Stage 1: Control only"]],
    relative_mean_shift_percent = 100 * mean_shift_stage3_minus_stage1 /
      .data[["mean__Stage 1: Control only"]]
  )
write.csv(
  stage13_shift_summary,
  file.path(table_dir, "stage1_vs_stage3_natural_scale_shift_summary.csv"),
  row.names = FALSE
)

for (m in names(fits)) {
  nat_vars <- stage13_natural_vars(m)
  plot_dat <- stage13_draws[[m]] %>%
    select(stage, all_of(nat_vars)) %>%
    pivot_longer(all_of(nat_vars), names_to = "parameter", values_to = "value")
  
  plot_dat$stage <- factor(
    plot_dat$stage,
    levels = c("Stage 1: Control only", "Stage 3: Control + Drug-only")
  )
  plot_dat$parameter <- factor(plot_dat$parameter, levels = nat_vars)
  
  subtitle_text <- if (m == "exponential") {
    "Population baseline volume and growth-rate distributions"
  } else {
    "Population baseline volume, growth-rate and carrying-capacity distributions"
  }
  
  p <- ggplot(plot_dat, aes(x = value, linetype = stage)) +
    geom_density(linewidth = 0.9, adjust = 1) +
    facet_wrap(~ parameter, scales = "free", ncol = 1) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    theme_bw() +
    labs(
      title = paste0(model_label(m), ": Stage 1 vs Stage 3 on natural scale"),
      subtitle = subtitle_text,
      x = "Parameter value",
      y = "Posterior density",
      linetype = NULL
    ) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey85"),
      strip.text = element_text(face = "plain")
    )
  
  save_plot(
    p,
    paste0("stage1_vs_stage3_natural_scale_", m),
    width = 10,
    height = ifelse(m == "exponential", 7, 9)
  )
}

for (m in names(fits)) {
  vars <- stage13_vars(m)
  plot_dat <- stage13_draws[[m]] %>%
    select(stage, all_of(vars)) %>%
    pivot_longer(all_of(vars), names_to = "parameter", values_to = "value")
  
  plot_dat$stage <- factor(
    plot_dat$stage,
    levels = c("Stage 1: Control only", "Stage 3: Control + Drug-only")
  )
  plot_dat$parameter <- factor(plot_dat$parameter, levels = vars)
  
  p <- ggplot(plot_dat, aes(x = value, linetype = stage)) +
    geom_density(linewidth = 0.8, adjust = 1) +
    facet_wrap(~ parameter, scales = "free", ncol = 2) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    theme_bw() +
    labs(
      title = paste0(model_label(m), ": Stage 1 vs Stage 3 posterior comparison"),
      subtitle = "Common fitted parameters",
      x = "Parameter value",
      y = "Posterior density",
      linetype = NULL
    ) +
    theme(legend.position = "bottom")
  
  save_plot(
    p,
    paste0("stage1_vs_stage3_model_scale_", m),
    width = 11,
    height = ifelse(m == "exponential", 8, 10)
  )
}


capture.output({
  cat("===== HMC SUMMARY =====\n"); print(hmc_summary)
  cat("\n===== RHAT / ESS SUMMARY =====\n"); print(conv_summary)
  cat("\n===== RESIDUAL SUMMARY =====\n"); print(residual_overall)
  cat("\n===== RESIDUAL SUMMARY BY GROUP =====\n"); print(residual_by_group)
  cat("\n===== 95% PREDICTIVE COVERAGE =====\n"); print(coverage_overall)
  cat("\n===== 95% PREDICTIVE COVERAGE BY GROUP =====\n"); print(coverage_by_group)
  cat("\n===== PAIRED LOO COMPARISON =====\n"); print(loo_cmp)
  cat("\n===== PARETO-k SUMMARY =====\n"); print(pareto_summary)
  cat("\n===== STAGE 1 vs STAGE 3: NATURAL-SCALE POSTERIOR SUMMARY =====\n"); print(stage13_natural_summary)
  cat("\n===== STAGE 1 vs STAGE 3: NATURAL-SCALE SHIFT SUMMARY =====\n"); print(stage13_shift_summary)
}, file = file.path(analysis_dir, "STAGE3_MASTER_SUMMARY.txt"))

capture.output(sessionInfo(), file = file.path(analysis_dir, "R_sessionInfo.txt"))

cat("\n====================================================\n")
cat("STAGE 3 POSTERIOR ANALYSIS COMPLETE\n")
cat("====================================================\n")
cat("Tables:  ", table_dir, "\n")
cat("Figures: ", figure_dir, "\n")
cat("Master summary: ", file.path(analysis_dir, "STAGE3_MASTER_SUMMARY.txt"), "\n")
cat("NO STAN MODEL WAS REFITTED.\n")
cat("====================================================\n")

drug_obs <- obs %>%
  filter(group == "Drug-only")

stopifnot(
  nrow(drug_obs) == 64L,
  length(unique(drug_obs$mouse_id)) == 13L
)

p_raw_drug <- ggplot(
  drug_obs,
  aes(x = time_days, y = volume)
) +
  geom_line(
    aes(group = mouse),
    colour = "grey50",
    linewidth = 0.5
  ) +
  geom_point(
    size = 1.8
  ) +
  facet_wrap(
    ~ mouse,
    scales = "free_y"
  ) +
  labs(
    x = "Time (days)",
    y = expression("Tumour volume (mm"^3*")")
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(figure_dir, "raw_drug_only_trajectories.pdf"),
  p_raw_drug,
  width = 10,
  height = 7
)

ggsave(
  file.path(figure_dir, "raw_drug_only_trajectories.png"),
  p_raw_drug,
  width = 10,
  height = 7,
  dpi = 300
)


for (m in names(fits)) {
  
  # Extract posterior draws from the completed Stage 3 fit
  d <- as.data.frame(
    posterior::as_draws_matrix(
      fits[[m]]$draws(
        variables = if (m == "exponential") {
          c(
            "logV0_pop",
            "logr_pop",
            "sigma_logV0",
            "sigma_logr",
            "sigma"
          )
        } else {
          c(
            "logV0_pop",
            "logr_pop",
            "logK_pop",
            "sigma_logV0",
            "sigma_logr",
            "sigma_logK",
            "sigma"
          )
        }
      )
    )
  )
  
  d$V0_pop <- exp(d$logV0_pop)
  d$r_pop  <- exp(d$logr_pop)
  
  if (m != "exponential") {
    d$K_pop <- exp(d$logK_pop)
  }
  
  if (m == "exponential") {
    
    plot_dat <- d %>%
      select(
        V0_pop,
        r_pop,
        sigma_logV0,
        sigma_logr,
        sigma
      ) %>%
      pivot_longer(
        cols = everything(),
        names_to = "parameter",
        values_to = "value"
      )
    
    parameter_levels <- c(
      "V0_pop",
      "r_pop",
      "sigma_logV0",
      "sigma_logr",
      "sigma"
    )
    
    parameter_labels <- c(
      V0_pop      = expression(V[0*",pop"]),
      r_pop       = expression(r["pop"]),
      sigma_logV0 = expression(sigma["log "*V[0]]),
      sigma_logr  = expression(sigma["log "*r]),
      sigma       = expression(sigma)
    )
    
  } else {
    
    plot_dat <- d %>%
      select(
        V0_pop,
        r_pop,
        K_pop,
        sigma_logV0,
        sigma_logr,
        sigma_logK,
        sigma
      ) %>%
      pivot_longer(
        cols = everything(),
        names_to = "parameter",
        values_to = "value"
      )
    
    parameter_levels <- c(
      "V0_pop",
      "r_pop",
      "K_pop",
      "sigma_logV0",
      "sigma_logr",
      "sigma_logK",
      "sigma"
    )
    
    parameter_labels <- c(
      V0_pop      = expression(V[0*",pop"]),
      r_pop       = expression(r["pop"]),
      K_pop       = expression(K["pop"]),
      sigma_logV0 = expression(sigma["log "*V[0]]),
      sigma_logr  = expression(sigma["log "*r]),
      sigma_logK  = expression(sigma["log "*K]),
      sigma       = expression(sigma)
    )
  }
  
  plot_dat$parameter <- factor(
    plot_dat$parameter,
    levels = parameter_levels
  )
  
  p_natural <- ggplot(
    plot_dat,
    aes(x = value)
  ) +
    geom_density(
      linewidth = 0.85,
      adjust = 1
    ) +
    facet_wrap(
      ~ parameter,
      scales = "free",
      ncol = 2,
      labeller = as_labeller(
        parameter_labels,
        default = label_parsed
      )
    ) +
    labs(
      x = "Parameter value",
      y = "Posterior density"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "plain"),
      panel.grid.minor = element_blank()
    )
  
  save_plot(
    p_natural,
    paste0(
      "posterior_parameters_natural_scale_",
      m
    ),
    width = 9,
    height = ifelse(
      m == "exponential",
      6.5,
      8
    )
  )
}



observed_fitted_panels <- lapply(names(analysis_list), function(m) {
  
  d <- analysis_list[[m]]
  
  ggplot(d, aes(x = fitted_mean, y = volume)) +
    geom_point(alpha = 0.65, size = 1.5) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = 2
    ) +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw() +
    labs(
      title = model_label(m),
      x = "Posterior mean fitted volume (mm³, log scale)",
      y = "Observed volume (mm³, log scale)"
    )
})

combined_observed_fitted <-
  wrap_plots(
    observed_fitted_panels,
    nrow = 1
  ) +
  plot_annotation(
    title = "Observed versus fitted tumour volumes"
  )

save_plot(
  combined_observed_fitted,
  "combined_observed_vs_fitted_stage3",
  width = 13,
  height = 4.5
)


residual_group_panels <- lapply(names(analysis_list), function(m) {
  
  d <- analysis_list[[m]]
  
  ggplot(
    d,
    aes(
      x = group,
      y = bayes_std_residual
    )
  ) +
    geom_violin(
      trim = FALSE,
      alpha = 0.30
    ) +
    geom_boxplot(
      width = 0.18,
      outlier.shape = NA
    ) +
    geom_hline(
      yintercept = 0,
      linetype = 2
    ) +
    theme_bw() +
    labs(
      title = model_label(m),
      x = NULL,
      y = "Standardised residual"
    )
})

combined_residuals_by_group <-
  wrap_plots(
    residual_group_panels,
    nrow = 1
  ) +
  plot_annotation(
    title = "Residual distributions by treatment group",
    subtitle = "Control and veliparib-only observations"
  )

save_plot(
  combined_residuals_by_group,
  "combined_residual_distribution_by_group_stage3",
  width = 13,
  height = 4.5
)

cat("\nCombined Stage 3 diagnostic figures saved.\n")


ppc_group_panels <- list()

for (m in names(fits)) {
  
  yrep <- posterior::as_draws_matrix(
    fits[[m]]$draws("y_rep")
  )
  
  # Fixed subset of posterior predictive draws for display
  set.seed(123)
  draw_id <- sample(
    seq_len(nrow(yrep)),
    min(n_ppc, nrow(yrep))
  )
  
  for (g in levels(obs$group)) {
    
    j <- which(obs$group == g)
    
    p <- bayesplot::ppc_dens_overlay(
      log(obs$volume[j]),
      log(yrep[draw_id, j, drop = FALSE])
    ) +
      ggtitle(
        paste0(model_label(m), " — ", g)
      ) +
      labs(
        x = "Log tumour volume",
        y = "Density"
      ) +
      theme_bw(base_size = 10) +
      theme(
        plot.title = element_text(
          size = 10,
          face = "plain"
        )
      )
    
    ppc_group_panels[[paste(m, g, sep = "_")]] <- p
  }
}

combined_ppc_density_by_group <-
  patchwork::wrap_plots(
    ppc_group_panels[
      c(
        "exponential_Control",
        "exponential_Drug-only",
        "logistic_Control",
        "logistic_Drug-only",
        "gompertz_Control",
        "gompertz_Drug-only"
      )
    ],
    ncol = 2
  ) +
  patchwork::plot_annotation(
    title = "Posterior predictive density checks by treatment group"
  )

save_plot(
  combined_ppc_density_by_group,
  "combined_ppc_density_by_group_stage3",
  width = 11,
  height = 11
)



time_ppc_panels <- lapply(
  names(timebin_list),
  function(m) {
    
    tab <- timebin_list[[m]]
    
    ggplot(
      tab,
      aes(
        x = time_bin,
        y = predictive_median
      )
    ) +
      geom_errorbar(
        aes(
          ymin = predictive_q025,
          ymax = predictive_q975
        ),
        width = 0.15
      ) +
      geom_point(
        size = 2
      ) +
      geom_point(
        aes(
          y = observed_mean_log_volume
        ),
        shape = 4,
        size = 3,
        stroke = 1
      ) +
      theme_bw(base_size = 10) +
      labs(
        title = model_label(m),
        x = "Time bin (days)",
        y = "Mean log tumour volume"
      ) +
      theme(
        axis.text.x = element_text(
          angle = 30,
          hjust = 1
        )
      )
  }
)

combined_ppc_time_binned <-
  patchwork::wrap_plots(
    time_ppc_panels,
    nrow = 1
  ) +
  patchwork::plot_annotation(
    title = "Time-binned posterior predictive checks",
    subtitle = "Cross = observed mean; interval = 95% posterior predictive interval"
  )

save_plot(
  combined_ppc_time_binned,
  "combined_ppc_time_binned_stage3",
  width = 13,
  height = 4.8
)


pareto_panels <- lapply(
  names(pareto_list),
  function(m) {
    
    d <- pareto_list[[m]]
    
    ggplot(
      d,
      aes(
        x = obs_index,
        y = pareto_k
      )
    ) +
      geom_point(
        alpha = 0.70,
        size = 1.5
      ) +
      geom_hline(
        yintercept = 0.7,
        linetype = 2
      ) +
      geom_hline(
        yintercept = 1.0,
        linetype = 3
      ) +
      theme_bw(base_size = 10) +
      labs(
        title = model_label(m),
        x = "Observation index",
        y = expression("Pareto-" * k)
      )
  }
)

combined_pareto_k <-
  patchwork::wrap_plots(
    pareto_panels,
    nrow = 1
  ) +
  patchwork::plot_annotation(
    title = "Observation-level PSIS-LOO diagnostics"
  )

save_plot(
  combined_pareto_k,
  "combined_pareto_k_stage3",
  width = 13,
  height = 4.5
)


stage13_combined_panels <- list()

for (m in names(fits)) {
  
  nat_vars <- stage13_natural_vars(m)
  
  plot_dat <- stage13_draws[[m]] %>%
    select(stage, all_of(nat_vars)) %>%
    pivot_longer(
      all_of(nat_vars),
      names_to = "parameter",
      values_to = "value"
    )
  
  plot_dat$stage <- factor(
    plot_dat$stage,
    levels = c(
      "Stage 1: Control only",
      "Stage 3: Control + Drug-only"
    )
  )
  
  parameter_labels <- c(
    V0_pop = expression(V[0*",pop"]),
    r_pop  = expression(r["pop"]),
    K_pop  = expression(K["pop"])
  )
  
  plot_dat$parameter <- factor(
    plot_dat$parameter,
    levels = nat_vars
  )
  
  p <- ggplot(
    plot_dat,
    aes(
      x = value,
      linetype = stage
    )
  ) +
    geom_density(
      linewidth = 0.8,
      adjust = 1
    ) +
    facet_wrap(
      ~ parameter,
      scales = "free",
      nrow = 1,
      labeller = as_labeller(
        parameter_labels,
        default = label_parsed
      )
    ) +
    scale_linetype_manual(
      values = c("solid", "dashed")
    ) +
    theme_bw(base_size = 10) +
    labs(
      title = model_label(m),
      x = "Parameter value",
      y = "Posterior density",
      linetype = NULL
    ) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "plain")
    )
  
  stage13_combined_panels[[m]] <- p
}

combined_stage1_vs_stage3 <-
  patchwork::wrap_plots(
    stage13_combined_panels,
    ncol = 1,
    guides = "collect"
  ) +
  patchwork::plot_annotation(
    title = "Stage 1 versus Stage 3 population-level posterior distributions",
    subtitle = "Control-only inference compared with joint control and veliparib-only inference"
  ) &
  theme(
    legend.position = "bottom"
  )

save_plot(
  combined_stage1_vs_stage3,
  "combined_stage1_vs_stage3_natural_scale",
  width = 12,
  height = 11
)


raw_drug_only_plot <-
  drug_obs %>%
  ggplot(
    aes(
      x = time_days,
      y = volume,
      group = mouse
    )
  ) +
  geom_line(
    linewidth = 0.55
  ) +
  geom_point(
    size = 1.7
  ) +
  facet_wrap(
    ~ mouse,
    scales = "free_y",
    ncol = 5
  ) +
  labs(
    x = "Time (days)",
    y = expression("Tumour volume (mm"^3 * ")")
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 9)
  )

raw_drug_only_plot

ggsave(
  "raw_drug_only_trajectories.pdf",
  raw_drug_only_plot,
  width = 10,
  height = 6.5
)

############################################################
# STAGE 3 POST-HOC:
# K_typical for consistent reporting with Stage 1
############################################################

library(posterior)

K_log_stage3 <- posterior::as_draws_matrix(
  fits$logistic$draws("K_mouse")
)

K_typical_log_stage3 <- apply(
  K_log_stage3,
  1,
  median
)

K_typical_log_summary <- quantile(
  K_typical_log_stage3,
  probs = c(0.025, 0.5, 0.975)
)


print(K_typical_log_summary)



K_gomp_stage3 <- posterior::as_draws_matrix(
  fits$gompertz$draws("K_mouse")
)

K_typical_gomp_stage3 <- apply(
  K_gomp_stage3,
  1,
  median
)

K_typical_gomp_summary <- quantile(
  K_typical_gomp_stage3,
  probs = c(0.025, 0.5, 0.975)
)


print(K_typical_gomp_summary)

fits$gompertz$summary() |>
  dplyr::arrange(ess_tail) |>
  dplyr::select(variable, rhat, ess_bulk, ess_tail) |>
  head(10)
fits$gompertz$summary() |>
  dplyr::arrange(ess_bulk) |>
  dplyr::select(variable, rhat, ess_bulk, ess_tail) |>
  head(10)

############################################################
# CORRECTED STAGE 1 vs STAGE 3 COMPARISON
# V0_pop, r_pop and K_typical
############################################################

library(dplyr)
library(ggplot2)

extract_scalar_natural <- function(fit, variable) {
  
  x <- posterior::as_draws_df(
    fit$draws(variables = variable)
  )[[variable]]
  
  x <- exp(x)
  
  stopifnot(
    length(x) > 0,
    all(is.finite(x))
  )
  
  x
}


extract_K_typical <- function(fit) {
  
  K <- posterior::as_draws_matrix(
    fit$draws(variables = "K_mouse")
  )
  
  # Keep only actual K_mouse columns
  K <- K[
    ,
    grepl("^K_mouse\\[", colnames(K)),
    drop = FALSE
  ]
  
  stopifnot(
    ncol(K) > 0,
    all(is.finite(K))
  )
  
  apply(
    K,
    1,
    median
  )
}


make_stage_summary <- function(fit, model, stage) {
  
  V0 <- extract_scalar_natural(
    fit,
    "logV0_pop"
  )
  
  r <- extract_scalar_natural(
    fit,
    "logr_pop"
  )
  
  out <- bind_rows(
    data.frame(
      parameter = "V0_pop",
      value = V0
    ),
    data.frame(
      parameter = "r_pop",
      value = r
    )
  )
  
  if (model %in% c("logistic", "gompertz")) {
    
    K_typical <- extract_K_typical(fit)
    
    out <- bind_rows(
      out,
      data.frame(
        parameter = "K_typical",
        value = K_typical
      )
    )
  }
  
  out |>
    mutate(
      model = tools::toTitleCase(model),
      stage = stage
    )
}


stage13_corrected_draws <- bind_rows(
  
  make_stage_summary(
    stage1_fits$exponential,
    "exponential",
    "Stage 1: Control only"
  ),
  
  make_stage_summary(
    fits$exponential,
    "exponential",
    "Stage 3: Control + Veliparib-only"
  ),
  
  make_stage_summary(
    stage1_fits$logistic,
    "logistic",
    "Stage 1: Control only"
  ),
  
  make_stage_summary(
    fits$logistic,
    "logistic",
    "Stage 3: Control + Veliparib-only"
  ),
  
  make_stage_summary(
    stage1_fits$gompertz,
    "gompertz",
    "Stage 1: Control only"
  ),
  
  make_stage_summary(
    fits$gompertz,
    "gompertz",
    "Stage 3: Control + Veliparib-only"
  )
)


stage13_corrected_summary <-
  stage13_corrected_draws |>
  group_by(
    model,
    parameter,
    stage
  ) |>
  summarise(
    median = median(value),
    q025 = quantile(value, 0.025),
    q975 = quantile(value, 0.975),
    .groups = "drop"
  )


print(
  as.data.frame(stage13_corrected_summary),
  row.names = FALSE
)



stage13_corrected_plot <-
  ggplot(
    stage13_corrected_draws,
    aes(
      x = value,
      linetype = stage
    )
  ) +
  geom_density(
    linewidth = 0.8
  ) +
  facet_grid(
    model ~ parameter,
    scales = "free",
    space = "free"
  ) +
  labs(
    x = "Parameter value",
    y = "Posterior density",
    linetype = NULL
  ) +
  theme_bw(
    base_size = 11
  )


ggsave(
  "combined_stage1_vs_stage3_natural_scale_corrected.pdf",
  stage13_corrected_plot,
  width = 10,
  height = 7
)



stage13_corrected_draws <-
  stage13_corrected_draws |>
  mutate(
    parameter_label = case_when(
      parameter == "V0_pop" ~ "V0_pop",
      parameter == "r_pop" ~ "r_pop",
      parameter == "K_typical" ~ "K_typical"
    ),
    panel = paste(model, parameter_label, sep = ": ")
  )

stage13_corrected_plot <-
  ggplot(
    stage13_corrected_draws,
    aes(
      x = value,
      linetype = stage
    )
  ) +
  geom_density(
    linewidth = 0.8
  ) +
  facet_wrap(
    ~ panel,
    scales = "free",
    ncol = 3
  ) +
  labs(
    x = "Parameter value",
    y = "Posterior density",
    linetype = NULL
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 9)
  )

ggsave(
  "combined_stage1_vs_stage3_natural_scale_corrected.pdf",
  stage13_corrected_plot,
  width = 10,
  height = 7
)

stage13_corrected_plot