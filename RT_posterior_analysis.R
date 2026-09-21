############################################################
# STAGE 2: CONTROL + RT-ONLY POSTERIOR ANALYSIS 
############################################################

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(loo)
  library(bayesplot)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

set.seed(123)
options(mc.cores = parallel::detectCores())

stage2_root <- path.expand("~/Desktop/STADIO2")

model_dirs <- c(
  exponential = "stage2_mfull_exp_ctrl_rt",
  logistic    = "stage2_mfull_logistic_ctrl_rt",
  gompertz    = "stage2_mfull_gompertz_ctrl_rt"
)

model_paths <- file.path(stage2_root, "results", model_dirs)
names(model_paths) <- names(model_dirs)

analysis_dir <- file.path(stage2_root, "posterior_analysis")
table_dir    <- file.path(stage2_root, "tables")
figure_dir   <- file.path(stage2_root, "figures")
object_dir   <- file.path(analysis_dir, "objects")

for (d in c(analysis_dir, table_dir, figure_dir, object_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

stopifnot(all(dir.exists(model_paths)))


# Logistic/Gompertz were saved as CmdStanR fit objects.
log_fit_file <- file.path(
  model_paths[["logistic"]],
  "fit_rt_logistic_MFULL_ctrl_rt_full.rds"
)

gomp_fit_file <- file.path(
  model_paths[["gompertz"]],
  "fit_rt_gompertz_MFULL_ctrl_rt_full.rds"
)

stopifnot(file.exists(log_fit_file), file.exists(gomp_fit_file))

fit_logistic <- readRDS(log_fit_file)
fit_gompertz <- readRDS(gomp_fit_file)

# The Exponential full run saved four CmdStan CSV files rather than a
# CmdStanMCMC RDS object  so we reconstruct the completed fit locally from those CSVs.
find_exp_full_csv <- function() {
  csv_dir <- file.path(model_paths[["exponential"]], "cmdstan_csv")
  stopifnot(dir.exists(csv_dir))
  
  x <- list.files(
    csv_dir,
    pattern = "stage2_mfull_exp_ctrl_rt_RK45_FULL.*\\.csv$",
    full.names = TRUE
  )
  
  x <- x[!grepl("fixed|validate|short|intermediate", basename(x), ignore.case = TRUE)]
  sort(x)
}

exp_csv_files <- find_exp_full_csv()

if (length(exp_csv_files) != 4L) {
  stop(
    "Could not identify exactly four Exponential full-fit CmdStan CSV files. Found: ",
    paste(exp_csv_files, collapse = ", ")
  )
}

fit_exponential <- cmdstanr::as_cmdstan_fit(exp_csv_files)

fits <- list(
  exponential = fit_exponential,
  logistic    = fit_logistic,
  gompertz    = fit_gompertz
)

stopifnot(all(vapply(fits, inherits, logical(1), what = "CmdStanMCMC")))

try(
  fit_exponential$save_object(
    file.path(object_dir, "fit_rt_exponential_MFULL_ctrl_rt_reconstructed.rds")
  ),
  silent = TRUE
)

writeLines(
  c(
    paste("Exponential CSV:", exp_csv_files),
    paste("Logistic fit:", log_fit_file),
    paste("Gompertz fit:", gomp_fit_file)
  ),
  con = file.path(analysis_dir, "completed_fit_provenance.txt")
)

obs_file <- file.path(
  model_paths[["exponential"]],
  "ctrl_rt_observations.csv"
)
rt_event_file <- file.path(
  model_paths[["exponential"]],
  "rt_events.csv"
)

stopifnot(file.exists(obs_file), file.exists(rt_event_file))

obs_all <- read.csv(obs_file, stringsAsFactors = FALSE) %>%
  arrange(mouse_id, time_days)

rt_events <- read.csv(rt_event_file, stringsAsFactors = FALSE) %>%
  arrange(mouse_id, time_days)

core_obs_cols <- intersect(
  c("mouse", "mouse_id", "time_days", "volume", "is_rt", "rt_index", "study"),
  names(obs_all)
)

for (m in c("logistic", "gompertz")) {
  f <- file.path(model_paths[[m]], "ctrl_rt_observations.csv")
  stopifnot(file.exists(f))
  z <- read.csv(f, stringsAsFactors = FALSE) %>% arrange(mouse_id, time_days)
  stopifnot(
    nrow(z) == nrow(obs_all),
    isTRUE(all.equal(z[, core_obs_cols, drop = FALSE],
                     obs_all[, core_obs_cols, drop = FALSE],
                     check.attributes = FALSE))
  )
}

stopifnot(
  nrow(obs_all) == 436L,
  length(unique(obs_all$mouse_id)) == 56L,
  nrow(rt_events) == 228L
)

if (!"is_rt" %in% names(obs_all)) {
  stop("ctrl_rt_observations.csv is missing is_rt.")
}

obs_all <- obs_all %>%
  mutate(
    is_rt = as.logical(is_rt),
    group = factor(ifelse(is_rt, "RT-only", "Control"),
                   levels = c("Control", "RT-only"))
  )

stopifnot(
  length(unique(obs_all$mouse_id[!obs_all$is_rt])) == 23L,
  length(unique(obs_all$mouse_id[obs_all$is_rt])) == 33L
)

baseline <- obs_all %>%
  group_by(mouse_id) %>%
  slice_min(time_days, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(mouse_id) %>%
  mutate(lik_index = row_number(), observation_part = "Baseline")

baseline_times <- baseline %>% select(mouse_id, baseline_time = time_days)
postbaseline <- obs_all %>%
  left_join(baseline_times, by = "mouse_id") %>%
  filter(time_days > baseline_time) %>%
  select(-baseline_time) %>%
  arrange(mouse_id, time_days) %>%
  mutate(
    lik_index = nrow(baseline) + row_number(),
    observation_part = "Post-baseline"
  )

obs_lik_order <- bind_rows(baseline, postbaseline)

stopifnot(
  nrow(baseline) == 56L,
  nrow(postbaseline) == 380L,
  nrow(obs_lik_order) == 436L,
  identical(obs_lik_order$lik_index, seq_len(436L))
)

write.csv(obs_all, file.path(table_dir, "stage2_observations_chronological.csv"), row.names = FALSE)
write.csv(obs_lik_order, file.path(table_dir, "stage2_observations_likelihood_order.csv"), row.names = FALSE)
write.csv(rt_events, file.path(table_dir, "stage2_rt_events.csv"), row.names = FALSE)


model_label <- function(m) {
  switch(m,
         exponential = "Exponential",
         logistic    = "Logistic",
         gompertz    = "Gompertz")
}

growth_prior_centre <- c(
  exponential = 0.3,
  logistic    = 0.4,
  gompertz    = 0.2
)

population_log_vars <- function(m) {
  x <- c(
    "logV0_pop", "logr_pop",
    "logalpha_pop", "logdelta_pop", "logbeta_a_pop",
    "logbeta_q_pop", "logk_rep_pop"
  )
  if (m != "exponential") x <- append(x, "logK_pop", after = 2)
  x
}

hierarchy_vars <- function(m) {
  x <- c(
    "sigma_logV0", "sigma_logr",
    "sigma_logalpha", "sigma_logdelta", "sigma_logbeta_a",
    "sigma_logbeta_q", "sigma_logk_rep", "sigma"
  )
  if (m != "exponential") x <- append(x, "sigma_logK", after = 2)
  x
}

key_vars <- function(m) c(population_log_vars(m), hierarchy_vars(m))

growth_vars <- function(m) {
  x <- c("logV0_pop", "logr_pop", "sigma_logV0", "sigma_logr", "sigma")
  if (m != "exponential") {
    x <- c("logV0_pop", "logr_pop", "logK_pop",
           "sigma_logV0", "sigma_logr", "sigma_logK", "sigma")
  }
  x
}

rt_vars <- function() c(
  "logalpha_pop", "logdelta_pop", "logbeta_a_pop",
  "logbeta_q_pop", "logk_rep_pop",
  "sigma_logalpha", "sigma_logdelta", "sigma_logbeta_a",
  "sigma_logbeta_q", "sigma_logk_rep"
)

natural_population_vars <- function(m) {
  x <- c("V0_pop", "r_pop", "alpha_pop", "delta_pop",
         "beta_a_pop", "beta_q_pop", "k_rep_pop")
  if (m != "exponential") x <- append(x, "K_pop", after = 2)
  x
}

save_plot <- function(p, name, width = 9, height = 6) {
  ggsave(file.path(figure_dir, paste0(name, ".pdf")), p,
         width = width, height = height, bg = "white")
  ggsave(file.path(figure_dir, paste0(name, ".png")), p,
         width = width, height = height, dpi = 300, bg = "white")
}

qfun <- function(x) {
  c(
    mean = mean(x),
    sd = sd(x),
    q2.5 = unname(quantile(x, 0.025)),
    median = unname(quantile(x, 0.5)),
    q97.5 = unname(quantile(x, 0.975))
  )
}

safe_draws <- function(fit, vars) {
  posterior::as_draws_matrix(fit$draws(variables = vars))
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
    predictive_coverage_95 = mean(pp_covered_95),
    .groups = "drop"
  )
}

rhalf_normal <- function(n, sd = 0.5) abs(rnorm(n, 0, sd))

########################
# HMC DIAGNOSTICS
########################

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

############################################################
#  RHAT / ESS / POSTERIOR SUMMARIES
############################################################

posterior_key_summary <- bind_rows(lapply(names(fits), function(m) {
  fits[[m]]$summary(variables = key_vars(m)) %>%
    as.data.frame() %>%
    mutate(model = model_label(m), .before = 1)
}))

write.csv(
  posterior_key_summary,
  file.path(table_dir, "posterior_key_parameter_summary_all_models.csv"),
  row.names = FALSE
)

convergence_summary <- posterior_key_summary %>%
  group_by(model) %>%
  summarise(
    max_rhat = max(rhat, na.rm = TRUE),
    min_bulk_ess = min(ess_bulk, na.rm = TRUE),
    min_tail_ess = min(ess_tail, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  convergence_summary,
  file.path(table_dir, "convergence_summary_all_models.csv"),
  row.names = FALSE
)

for (m in names(fits)) {
  sm <- fits[[m]]$summary() %>% as.data.frame()
  sm_non_z <- sm %>% filter(!grepl("^z_", variable))
  write.csv(
    sm_non_z,
    file.path(table_dir, paste0("posterior_all_non_z_parameters_", m, ".csv")),
    row.names = FALSE
  )
}

############################################################
# 6. TRACE / RANK / POSTERIOR AREA PLOTS
############################################################

for (m in names(fits)) {
  vars <- key_vars(m)
  arr <- posterior::as_draws_array(fits[[m]]$draws(variables = vars))
  
  p_trace <- bayesplot::mcmc_trace(
    arr,
    pars = vars,
    facet_args = list(ncol = 1, strip.position = "left")
  ) +
    ggtitle(paste0("Stage 2 ", model_label(m), ": trace plots"))
  
  save_plot(
    p_trace,
    paste0("trace_key_parameters_", m),
    width = 10,
    height = ifelse(m == "exponential", 20, 22)
  )
  
  p_rank <- bayesplot::mcmc_rank_overlay(arr, pars = vars) +
    ggtitle(paste0("Stage 2 ", model_label(m), ": rank plots"))
  save_plot(
    p_rank,
    paste0("rank_key_parameters_", m),
    width = 10,
    height = ifelse(m == "exponential", 14, 16)
  )
  
  p_area <- bayesplot::mcmc_areas(
    arr,
    pars = vars,
    prob = 0.8,
    prob_outer = 0.95
  ) +
    ggtitle(paste0("Stage 2 ", model_label(m), ": posterior distributions"))
  save_plot(
    p_area,
    paste0("posterior_key_parameters_", m),
    width = 10,
    height = ifelse(m == "exponential", 11, 12)
  )
}

############################################################
#  PRIOR vs POSTERIOR
############################################################

n_prior <- 50000L

prior_draws_stage2 <- function(m, n = n_prior) {
  x <- data.frame(
    logV0_pop = rnorm(n, log(30), 0.5),
    logr_pop = rnorm(n, log(growth_prior_centre[[m]]), 0.5),
    logalpha_pop = rnorm(n, log(0.9), 0.5),
    logdelta_pop = rnorm(n, log(1.2), 0.5),
    logbeta_a_pop = rnorm(n, log(26.4), 0.5),
    logbeta_q_pop = rnorm(n, log(4.8), 0.5),
    logk_rep_pop = rnorm(n, log(86.4), 0.5),
    sigma_logV0 = rhalf_normal(n),
    sigma_logr = rhalf_normal(n),
    sigma_logalpha = rhalf_normal(n),
    sigma_logdelta = rhalf_normal(n),
    sigma_logbeta_a = rhalf_normal(n),
    sigma_logbeta_q = rhalf_normal(n),
    sigma_logk_rep = rhalf_normal(n),
    sigma = rhalf_normal(n)
  )
  
  if (m != "exponential") {
    x$logK_pop <- rnorm(n, log(500), 0.7)
    x$sigma_logK <- rhalf_normal(n)
  }
  
  x
}

for (m in names(fits)) {
  prior <- prior_draws_stage2(m)
  
  groups <- list(
    growth = growth_vars(m),
    rt_mechanism = rt_vars(),
    population_location = population_log_vars(m)
  )
  
  # Avoid sigma appearing in the population-location group.
  for (g in names(groups)) {
    vars <- groups[[g]]
    
    pr <- prior %>%
      select(all_of(vars)) %>%
      pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
      mutate(distribution = "Prior")
    
    po <- as.data.frame(safe_draws(fits[[m]], vars)) %>%
      select(all_of(vars)) %>%
      pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
      mutate(distribution = "Posterior")
    
    pp <- bind_rows(pr, po)
    
    p <- ggplot(pp, aes(value, linetype = distribution, group = distribution)) +
      geom_density(linewidth = 0.8, adjust = 1.05) +
      facet_wrap(~ parameter, scales = "free", ncol = 2) +
      labs(
        title = paste0("Stage 2 ", model_label(m), ": prior vs posterior"),
        subtitle = switch(
          g,
          growth = "Intrinsic growth and observation parameters",
          rt_mechanism = "RT mechanism: population and between-RT-mouse parameters",
          population_location = "Population location parameters"
        ),
        x = "Parameter value",
        y = "Density",
        linetype = NULL
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        strip.text = element_text(face = "bold")
      )
    
    save_plot(
      p,
      paste0("prior_vs_posterior_", g, "_", m),
      width = 10,
      height = ifelse(length(vars) <= 7, 9, 12)
    )
  }
}

prior_posterior_natural <- list()

for (m in names(fits)) {
  pr <- prior_draws_stage2(m)
  po <- as.data.frame(safe_draws(fits[[m]], population_log_vars(m)))
  
  transform_population <- function(d, m) {
    out <- data.frame(
      V0_pop = exp(d$logV0_pop),
      r_pop = exp(d$logr_pop),
      alpha_pop = exp(d$logalpha_pop),
      delta_pop = exp(d$logdelta_pop),
      beta_a_pop = exp(d$logbeta_a_pop),
      beta_q_pop = exp(d$logbeta_q_pop),
      k_rep_pop = exp(d$logk_rep_pop)
    )
    if (m != "exponential") out$K_pop <- exp(d$logK_pop)
    out
  }
  
  pr_nat <- transform_population(pr, m)
  po_nat <- transform_population(po, m)
  
  summarize_nat <- function(d, distribution) {
    d %>%
      pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
      group_by(parameter) %>%
      summarise(
        mean = mean(value),
        sd = sd(value),
        q025 = quantile(value, 0.025),
        median = quantile(value, 0.5),
        q975 = quantile(value, 0.975),
        .groups = "drop"
      ) %>%
      mutate(distribution = distribution)
  }
  
  s <- bind_rows(
    summarize_nat(pr_nat, "Prior"),
    summarize_nat(po_nat, "Posterior")
  ) %>% mutate(model = model_label(m), .before = 1)
  
  prior_posterior_natural[[m]] <- s
  
  p <- ggplot(s, aes(median, parameter, shape = distribution)) +
    geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0.16,
                   position = position_dodge(width = 0.45)) +
    geom_point(size = 2.3, position = position_dodge(width = 0.45)) +
    facet_wrap(~ parameter, scales = "free_x", ncol = 2) +
    labs(
      title = paste0("Stage 2 ", model_label(m), ": prior vs posterior on natural scales"),
      subtitle = "Median and 95% interval; each facet has its own x-scale",
      x = "Natural-scale parameter value", y = NULL, shape = NULL
    ) +
    theme_bw() +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
  
  save_plot(p, paste0("prior_vs_posterior_natural_scale_", m), 10, 9)
}

prior_posterior_natural_all <- bind_rows(prior_posterior_natural)
write.csv(
  prior_posterior_natural_all,
  file.path(table_dir, "prior_vs_posterior_natural_scale_summary.csv"),
  row.names = FALSE
)




natural_summary <- list()

for (m in names(fits)) {
  logvars <- population_log_vars(m)
  d <- as.data.frame(safe_draws(fits[[m]], logvars))
  
  nat <- data.frame(
    V0_pop = exp(d$logV0_pop),
    r_pop = exp(d$logr_pop),
    alpha_pop = exp(d$logalpha_pop),
    delta_pop = exp(d$logdelta_pop),
    beta_a_pop = exp(d$logbeta_a_pop),
    beta_q_pop = exp(d$logbeta_q_pop),
    k_rep_pop = exp(d$logk_rep_pop)
  )
  if (m != "exponential") nat$K_pop <- exp(d$logK_pop)
  
  ss <- nat %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    group_by(parameter) %>%
    summarise(
      mean = mean(value), sd = sd(value),
      q025 = quantile(value, 0.025),
      median = quantile(value, 0.5),
      q975 = quantile(value, 0.975),
      .groups = "drop"
    ) %>%
    mutate(model = model_label(m), .before = 1)
  
  natural_summary[[m]] <- ss
}

natural_summary_all <- bind_rows(natural_summary)
write.csv(
  natural_summary_all,
  file.path(table_dir, "population_parameters_natural_scale_summary.csv"),
  row.names = FALSE
)



mouse_lookup <- obs_all %>%
  arrange(mouse_id, time_days) %>%
  group_by(mouse_id) %>%
  summarise(
    mouse = first(mouse),
    is_rt = first(is_rt),
    group = first(group),
    study = if ("study" %in% names(obs_all)) first(study) else NA,
    .groups = "drop"
  ) %>%
  arrange(mouse_id)

stopifnot(
  nrow(mouse_lookup) == 56L,
  !anyDuplicated(mouse_lookup$mouse_id)
)

rt_lookup <- mouse_lookup %>%
  filter(is_rt) %>%
  arrange(mouse_id) %>%
  mutate(rt_index = row_number())

stopifnot(
  nrow(rt_lookup) == 33L
)

summarise_indexed_draws <- function(fit, prefix, lookup, index_name) {
  mat <- safe_draws(fit, prefix)
  stopifnot(ncol(mat) == nrow(lookup))
  
  bind_rows(lapply(seq_len(ncol(mat)), function(j) {
    z <- mat[, j]
    cbind(
      lookup[j, , drop = FALSE],
      parameter = prefix,
      index = j,
      mean = mean(z), sd = sd(z),
      q025 = quantile(z, 0.025), median = quantile(z, 0.5), q975 = quantile(z, 0.975)
    )
  }))
}

mouse_parameter_summaries <- list()

for (m in names(fits)) {
  fit <- fits[[m]]
  pieces <- list(
    summarise_indexed_draws(fit, "V0_mouse", mouse_lookup, "mouse_id"),
    summarise_indexed_draws(fit, "r_mouse", mouse_lookup, "mouse_id")
  )
  
  if (m != "exponential") {
    pieces <- c(pieces, list(
      summarise_indexed_draws(fit, "K_mouse", mouse_lookup, "mouse_id")
    ))
  }
  
  for (p in c("alpha_mouse", "delta_mouse", "beta_a_mouse", "beta_q_mouse", "k_rep_mouse")) {
    pieces <- c(pieces, list(
      summarise_indexed_draws(fit, p, rt_lookup, "rt_index")
    ))
  }
  
  out <- bind_rows(pieces) %>% mutate(model = model_label(m), .before = 1)
  mouse_parameter_summaries[[m]] <- out
  
  write.csv(
    out,
    file.path(table_dir, paste0("mouse_specific_parameter_summary_", m, ".csv")),
    row.names = FALSE
  )
  
  # Forest plots
  for (pname in unique(out$parameter)) {
    z <- out %>% filter(parameter == pname)
    p <- ggplot(z, aes(median, factor(mouse, levels = rev(unique(mouse))))) +
      geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0.15) +
      geom_point(size = 1.6) +
      labs(
        title = paste0("Stage 2 ", model_label(m), ": ", pname),
        subtitle = ifelse(grepl("alpha|delta|beta|k_rep", pname),
                          "RT-only mice", "All 56 mice"),
        x = "Posterior median and 95% credible interval",
        y = "Mouse"
      ) + theme_bw()
    
    save_plot(
      p,
      paste0("mouse_forest_", pname, "_", m),
      width = 8,
      height = ifelse(nrow(z) > 40, 12, 9)
    )
  }
}

mouse_parameter_summary_all <- bind_rows(mouse_parameter_summaries)
write.csv(
  mouse_parameter_summary_all,
  file.path(table_dir, "mouse_specific_parameter_summary_all_models.csv"),
  row.names = FALSE
)

############################################################
#  CONTROLvsRT 
# This is a diagnostic check. It asks whether the
# posterior mouse-specific untreated-growth quantities show systematic
# separation by treatment cohort or by study after the joint hierarchy is
# fitted. 
############################################################

intrinsic_mouse <- mouse_parameter_summary_all %>%
  filter(parameter %in% c("V0_mouse", "r_mouse", "K_mouse"))

write.csv(
  intrinsic_mouse,
  file.path(table_dir, "intrinsic_growth_mouse_parameters_by_group_study.csv"),
  row.names = FALSE
)

for (m in names(fits)) {
  for (pname in c("V0_mouse", "r_mouse", if (m != "exponential") "K_mouse" else NULL)) {
    z <- intrinsic_mouse %>%
      filter(model == model_label(m), parameter == pname)
    
    p_group <- ggplot(z, aes(group, median)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.12, height = 0, size = 1.7) +
      labs(
        title = paste0(model_label(m), ": mouse-specific ", pname, " by cohort"),
        subtitle = "Descriptive check of the joint Control + RT hierarchy",
        x = NULL, y = "Posterior median"
      ) + theme_bw()
    save_plot(p_group, paste0("intrinsic_", pname, "_by_group_", m), 7, 6)
    
    if ("study" %in% names(z) && any(!is.na(z$study))) {
      p_study <- ggplot(z, aes(factor(study), median)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(aes(shape = group), width = 0.12, height = 0, size = 1.7) +
        labs(
          title = paste0(model_label(m), ": mouse-specific ", pname, " by study"),
          subtitle = "Points identify Control vs RT-only cohort",
          x = "Study", y = "Posterior median", shape = "Cohort"
        ) + theme_bw() + theme(legend.position = "bottom")
      save_plot(p_study, paste0("intrinsic_", pname, "_by_study_", m), 8, 6)
    }
  }
}


analysis_list <- list()

for (m in names(fits)) {
  fit <- fits[[m]]
  
  V0 <- safe_draws(fit, "V0_mouse")
  mu_post <- safe_draws(fit, "mu_out")
  mu_log_post <- safe_draws(fit, "mu_log_out")
  yrep_base <- safe_draws(fit, "y_rep_baseline")
  yrep_post <- safe_draws(fit, "y_rep")
  sigma_draw <- as.numeric(safe_draws(fit, "sigma")[, 1])
  
  stopifnot(
    ncol(V0) == 56L,
    ncol(mu_post) == 380L,
    ncol(mu_log_post) == 380L,
    ncol(yrep_base) == 56L,
    ncol(yrep_post) == 380L
  )
  
  mu_full <- cbind(V0, mu_post)
  mu_log_full <- cbind(log(V0), mu_log_post)
  yrep_full <- cbind(yrep_base, yrep_post)
  
  stopifnot(ncol(mu_full) == nrow(obs_lik_order))
  
  fitted_mean <- colMeans(mu_full)
  fitted_q025 <- apply(mu_full, 2, quantile, 0.025)
  fitted_q975 <- apply(mu_full, 2, quantile, 0.975)
  
  pp_q025 <- apply(yrep_full, 2, quantile, 0.025)
  pp_q50 <- apply(yrep_full, 2, quantile, 0.50)
  pp_q975 <- apply(yrep_full, 2, quantile, 0.975)
  
  numerator <- matrix(
    log(obs_lik_order$volume),
    nrow = nrow(mu_log_full),
    ncol = ncol(mu_log_full),
    byrow = TRUE
  ) - mu_log_full
  
  z_draws <- sweep(numerator, 1, sigma_draw, "/")
  
  d <- obs_lik_order %>% mutate(
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
write.csv(
  analysis_all,
  file.path(table_dir, "observation_diagnostics_all_models.csv"),
  row.names = FALSE
)

############################################################
#  NUMERICAL RESIDUAL - COVERAGE SUMMARIES
############################################################

residual_overall <- analysis_all %>% group_by(model) %>% error_summary()
residual_by_group <- analysis_all %>% group_by(model, group) %>% error_summary()
residual_by_part <- analysis_all %>% group_by(model, observation_part) %>% error_summary()
residual_by_group_part <- analysis_all %>%
  group_by(model, group, observation_part) %>% error_summary()

write.csv(residual_overall, file.path(table_dir, "residual_summary_all_models.csv"), row.names = FALSE)
write.csv(residual_by_group, file.path(table_dir, "residual_summary_by_group.csv"), row.names = FALSE)
write.csv(residual_by_part, file.path(table_dir, "residual_summary_baseline_postbaseline.csv"), row.names = FALSE)
write.csv(residual_by_group_part, file.path(table_dir, "residual_summary_by_group_and_part.csv"), row.names = FALSE)

############################################################
#  OBSERVEDvsFITTED and RESIDUAL DIAGNOSTIC FIGURES
############################################################

for (m in names(analysis_list)) {
  d <- analysis_list[[m]]
  ttl <- paste0("Stage 2 ", model_label(m), ": ")
  
  p1 <- ggplot(d, aes(fitted_mean, volume, shape = group)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point(alpha = 0.8) +
    scale_x_log10() + scale_y_log10() +
    labs(
      title = paste0(ttl, "observed vs fitted"),
      x = "Posterior mean fitted volume (mm3, log scale)",
      y = "Observed volume (mm3, log scale)", shape = "Cohort"
    ) + theme_bw() + theme(legend.position = "bottom")
  save_plot(p1, paste0("observed_vs_fitted_", m), 8, 6)
  
  p2 <- ggplot(d, aes(fitted_mean, bayes_std_residual, shape = group)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(alpha = 0.75) + scale_x_log10() +
    labs(
      title = paste0(ttl, "Bayesian standardised residuals vs fitted"),
      x = "Posterior mean fitted volume (mm3, log scale)",
      y = "Standardised residual", shape = "Cohort"
    ) + theme_bw() + theme(legend.position = "bottom")
  save_plot(p2, paste0("standardised_residuals_vs_fitted_", m), 8, 6)
  
  p3 <- ggplot(d, aes(time_days, bayes_std_residual, shape = group)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(alpha = 0.75) +
    labs(
      title = paste0(ttl, "Bayesian standardised residuals vs time"),
      x = "Time (days)", y = "Standardised residual", shape = "Cohort"
    ) + theme_bw() + theme(legend.position = "bottom")
  save_plot(p3, paste0("standardised_residuals_vs_time_", m), 8, 6)
  
  p4 <- ggplot(d, aes(bayes_std_residual)) +
    geom_histogram(bins = 30) +
    labs(title = paste0(ttl, "standardised residual histogram"),
         x = "Standardised residual", y = "Count") + theme_bw()
  save_plot(p4, paste0("residual_histogram_", m), 7, 6)
  
  p5 <- ggplot(d, aes(sample = bayes_std_residual)) +
    stat_qq() + stat_qq_line() +
    labs(title = paste0(ttl, "standardised residual Q-Q plot"),
         x = "Theoretical Normal quantiles", y = "Observed residual quantiles") +
    theme_bw()
  save_plot(p5, paste0("residual_qq_", m), 7, 6)
  
  p6 <- ggplot(d, aes(group, bayes_std_residual)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.12, alpha = 0.45) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(title = paste0(ttl, "residuals by cohort"), x = NULL,
         y = "Standardised residual") + theme_bw()
  save_plot(p6, paste0("residual_distribution_by_group_", m), 7, 6)
}

############################################################
#  MOUSE-SPECIFIC FITTED TRAJECTORIES
############################################################

for (m in names(analysis_list)) {
  d <- analysis_list[[m]] %>%
    arrange(mouse_id, time_days)
  
  e <- rt_events %>% distinct(mouse, time_days)
  
  p <- ggplot(d, aes(time_days, volume)) +
    geom_ribbon(aes(ymin = fitted_q025, ymax = fitted_q975), alpha = 0.2) +
    geom_line(aes(y = fitted_mean), linewidth = 0.65) +
    geom_point(size = 1.2) +
    geom_vline(
      data = e,
      aes(xintercept = time_days),
      linetype = 3,
      linewidth = 0.25,
      inherit.aes = FALSE
    ) +
    facet_wrap(~ mouse, scales = "free_y", ncol = 7) +
    labs(
      title = paste0("Stage 2 ", model_label(m), ": fitted trajectories by mouse"),
      subtitle = "Points: observed; line/ribbon: latent fitted mean and 95% CrI; vertical lines: RT fractions",
      x = "Time (days)", y = "Tumour volume (mm3)"
    ) + theme_bw() + theme(strip.text = element_text(size = 7))
  
  save_plot(p, paste0("mouse_fitted_trajectories_", m), 16, 14)
}

############################################################
#  POSTERIOR PREDICTIVE DENSITY CHECKS USING y_rep
############################################################

for (m in names(fits)) {
  fit <- fits[[m]]
  yrep_full <- cbind(
    safe_draws(fit, "y_rep_baseline"),
    safe_draws(fit, "y_rep")
  )
  
  # Thin draws only for plotting clarity; inference above uses all draws.
  take <- unique(round(seq(1, nrow(yrep_full), length.out = min(100, nrow(yrep_full)))))
  
  subsets <- list(
    all = seq_len(nrow(obs_lik_order)),
    control = which(obs_lik_order$group == "Control"),
    rt_only = which(obs_lik_order$group == "RT-only")
  )
  
  for (s in names(subsets)) {
    idx <- subsets[[s]]
    p <- bayesplot::ppc_dens_overlay(
      y = log(obs_lik_order$volume[idx]),
      yrep = log(yrep_full[take, idx, drop = FALSE])
    ) +
      ggtitle(
        paste0("Stage 2 ", model_label(m), ": posterior predictive density — ", s),
        subtitle = "Log tumour volume; replicated observations from y_rep"
      )
    save_plot(p, paste0("ppc_density_", s, "_", m), 9, 6)
  }
}


#  95% POSTERIOR PREDICTIVE 


coverage_summary <- analysis_all %>%
  group_by(model) %>%
  summarise(n = n(), coverage_95 = mean(pp_covered_95), .groups = "drop")

coverage_by_group <- analysis_all %>%
  group_by(model, group) %>%
  summarise(n = n(), coverage_95 = mean(pp_covered_95), .groups = "drop")

coverage_by_group_part <- analysis_all %>%
  group_by(model, group, observation_part) %>%
  summarise(n = n(), coverage_95 = mean(pp_covered_95), .groups = "drop")

write.csv(coverage_summary, file.path(table_dir, "posterior_predictive_coverage.csv"), row.names = FALSE)
write.csv(coverage_by_group, file.path(table_dir, "posterior_predictive_coverage_by_group.csv"), row.names = FALSE)
write.csv(coverage_by_group_part, file.path(table_dir, "posterior_predictive_coverage_by_group_and_part.csv"), row.names = FALSE)

############################################################
#  TIME-BINNED POSTERIOR PREDICTIVE CHECKS
############################################################

# Fixed bins make cross-model comparison straightforward.
breaks <- c(-Inf, 2, 5, 10, 15, 20, Inf)
labels <- c("0-2", "2-5", "5-10", "10-15", "15-20", ">20")

obs_binned <- obs_lik_order %>%
  mutate(time_bin = cut(time_days, breaks = breaks, labels = labels, right = TRUE))

time_ppc_all <- list()

for (m in names(fits)) {
  yrep_full <- cbind(
    safe_draws(fits[[m]], "y_rep_baseline"),
    safe_draws(fits[[m]], "y_rep")
  )
  
  rows <- list()
  rr <- 0L
  
  for (g in levels(obs_binned$group)) {
    for (b in levels(obs_binned$time_bin)) {
      idx <- which(obs_binned$group == g & obs_binned$time_bin == b)
      if (length(idx) == 0L) next
      
      rep_stat <- rowMeans(log(yrep_full[, idx, drop = FALSE]))
      obs_stat <- mean(log(obs_binned$volume[idx]))
      
      rr <- rr + 1L
      rows[[rr]] <- data.frame(
        model = model_label(m), group = g, time_bin = b,
        n_observations = length(idx), observed_mean_log_volume = obs_stat,
        predictive_q025 = quantile(rep_stat, 0.025),
        predictive_median = quantile(rep_stat, 0.5),
        predictive_q975 = quantile(rep_stat, 0.975),
        covered_95 = obs_stat >= quantile(rep_stat, 0.025) &
          obs_stat <= quantile(rep_stat, 0.975)
      )
    }
  }
  
  z <- bind_rows(rows)
  time_ppc_all[[m]] <- z
  
  p <- ggplot(z, aes(time_bin, observed_mean_log_volume, shape = group)) +
    geom_errorbar(aes(ymin = predictive_q025, ymax = predictive_q975),
                  position = position_dodge(width = 0.4), width = 0.15) +
    geom_point(position = position_dodge(width = 0.4), size = 2.2) +
    labs(
      title = paste0("Stage 2 ", model_label(m), ": time-binned posterior predictive check"),
      subtitle = "Point: observed mean log-volume; interval: posterior predictive 95% interval",
      x = "Time bin (days)", y = "Mean log tumour volume", shape = "Cohort"
    ) + theme_bw() + theme(legend.position = "bottom")
  
  save_plot(p, paste0("ppc_time_binned_", m), 9, 6)
}

write.csv(
  bind_rows(time_ppc_all),
  file.path(table_dir, "time_binned_posterior_predictive_checks.csv"),
  row.names = FALSE
)

############################################################
#  OBSERVATION-LEVEL PSIS-LOO + PAIRED MODEL COMPARISON
#
# 
# The completed Exponential fit stored log_lik on the
# log-volume Normal scale:
#
#   normal_lpdf(log(volume) | mu_log, sigma)
#
# but the corrected Logistic and Gompertz fits stored
# log_lik on the equivalent raw-volume Lognormal scale:
#
#   lognormal_lpdf(volume | mu_log, sigma)
#
# Logistic and Gompertz log_lik
# values are therefore transformed to the same log-volume
# Normal scale as Exponential by adding log(volume).
#.
############################################################

log_y <- log(obs_lik_order$volume)

log_lik_corrected <- list()
loos <- list()
pareto_list <- list()

for (m in names(fits)) {
  
  ll <- safe_draws(fits[[m]], "log_lik")
  
  stopifnot(
    ncol(ll) == nrow(obs_lik_order),
    nrow(ll) == 4000L
  )
  
  if (m == "exponential") {
        ll_corrected <- ll
    
  } else {
      ll_corrected <- sweep(
      ll,
      MARGIN = 2,
      STATS = log_y,
      FUN = "+"
    )
  }
  
  log_lik_corrected[[m]] <- ll_corrected
  loos[[m]] <- loo::loo(ll_corrected)
  
  saveRDS(
    loos[[m]],
    file.path(
      object_dir,
      paste0("loo_corrected_", m, ".rds")
    )
  )
  
  k <- loo::pareto_k_values(loos[[m]])
  
  pareto_list[[m]] <- obs_lik_order %>%
    mutate(
      model = model_label(m),
      obs_index = seq_len(n()),
      pareto_k = k,
      pareto_flag = case_when(
        pareto_k > 1   ~ ">1",
        pareto_k > 0.7 ~ "0.7-1",
        pareto_k > 0.5 ~ "0.5-0.7",
        TRUE           ~ "<=0.5"
      )
    )
}


loo_summary <- bind_rows(
  lapply(names(loos), function(m) {
    
    x <- loos[[m]]$estimates
    
    data.frame(
      model = model_label(m),
      elpd_loo = x["elpd_loo", "Estimate"],
      se_elpd_loo = x["elpd_loo", "SE"],
      p_loo = x["p_loo", "Estimate"],
      looic = x["looic", "Estimate"]
    )
  })
)

write.csv(
  loo_summary,
  file.path(
    table_dir,
    "loo_summary_all_models.csv"
  ),
  row.names = FALSE
)


loo_cmp <- loo::loo_compare(
  list(
    exponential = loos$exponential,
    logistic = loos$logistic,
    gompertz = loos$gompertz
  )
)

loo_cmp_df <- data.frame(
  model = rownames(loo_cmp),
  as.data.frame(loo_cmp),
  row.names = NULL
)

write.csv(
  loo_cmp_df,
  file.path(
    table_dir,
    "loo_paired_comparison.csv"
  ),
  row.names = FALSE
)

likelihood_scale_record <- data.frame(
  model = c(
    "Exponential",
    "Logistic",
    "Gompertz"
  ),
  saved_log_lik_scale = c(
    "Normal likelihood for log(volume)",
    "Lognormal likelihood for volume",
    "Lognormal likelihood for volume"
  ),
  loo_comparison_scale = rep(
    "Normal likelihood for log(volume)",
    3
  ),
  correction_applied = c(
    "None",
    "+ log(volume)",
    "+ log(volume)"
  )
)

write.csv(
  likelihood_scale_record,
  file.path(
    table_dir,
    "loo_likelihood_scale_correction.csv"
  ),
  row.names = FALSE
)


print(loo_summary)

print(loo_cmp)

############################################################
#  PARETO-k TABLES / FIGURES etc.
############################################################

pareto_all <- bind_rows(pareto_list)
problematic <- pareto_all %>%
  filter(pareto_k > 0.7) %>%
  arrange(model, desc(pareto_k))

pareto_summary <- pareto_all %>%
  group_by(model) %>%
  summarise(
    n_observations = n(),
    n_k_gt_0_5 = sum(pareto_k > 0.5),
    n_k_gt_0_7 = sum(pareto_k > 0.7),
    n_k_gt_1 = sum(pareto_k > 1),
    max_pareto_k = max(pareto_k),
    .groups = "drop"
  )

pareto_by_mouse <- pareto_all %>%
  group_by(model, mouse, mouse_id, group) %>%
  summarise(
    n = n(), max_pareto_k = max(pareto_k),
    n_k_gt_0_7 = sum(pareto_k > 0.7),
    .groups = "drop"
  ) %>% arrange(model, desc(max_pareto_k))

write.csv(pareto_all, file.path(table_dir, "pareto_k_all_observations.csv"), row.names = FALSE)
write.csv(problematic, file.path(table_dir, "pareto_k_problematic_observations.csv"), row.names = FALSE)
write.csv(pareto_summary, file.path(table_dir, "pareto_k_summary.csv"), row.names = FALSE)
write.csv(pareto_by_mouse, file.path(table_dir, "pareto_k_by_mouse.csv"), row.names = FALSE)

for (m in names(pareto_list)) {
  d <- pareto_list[[m]]
  p <- ggplot(d, aes(obs_index, pareto_k, shape = group)) +
    geom_hline(yintercept = c(0.5, 0.7, 1), linetype = c(3, 2, 1)) +
    geom_point(alpha = 0.8) +
    labs(
      title = paste0("Stage 2 ", model_label(m), ": PSIS-LOO Pareto-k"),
      x = "Likelihood contribution index (baseline first)",
      y = "Pareto k", shape = "Cohort"
    ) + theme_bw() + theme(legend.position = "bottom")
  save_plot(p, paste0("pareto_k_", m), 9, 6)
}

############################################################
# RT-MECHANISM POSTERIOR DEPENDENCE - IDENTIFIABILITY
############################################################

rt_population_vars <- c(
  "logalpha_pop", "logdelta_pop", "logbeta_a_pop",
  "logbeta_q_pop", "logk_rep_pop"
)

rt_corr_all <- list()

for (m in names(fits)) {
  d <- as.data.frame(safe_draws(fits[[m]], rt_population_vars))
  cc <- cor(d)
  long <- as.data.frame(as.table(cc), stringsAsFactors = FALSE)
  names(long) <- c("parameter_1", "parameter_2", "correlation")
  long$model <- model_label(m)
  rt_corr_all[[m]] <- long
  
  write.csv(
    long,
    file.path(table_dir, paste0("rt_population_parameter_correlations_", m, ".csv")),
    row.names = FALSE
  )
  
  arr <- fits[[m]]$draws(variables = rt_population_vars)
  p <- bayesplot::mcmc_pairs(
    arr,
    pars = rt_population_vars,
    off_diag_args = list(size = 0.5, alpha = 0.25)
  )
  ggsave(
    file.path(figure_dir, paste0("rt_population_pairs_", m, ".pdf")),
    p, width = 11, height = 11, bg = "white"
  )
  ggsave(
    file.path(figure_dir, paste0("rt_population_pairs_", m, ".png")),
    p, width = 11, height = 11, dpi = 300, bg = "white"
  )
}

write.csv(
  bind_rows(rt_corr_all),
  file.path(table_dir, "rt_population_parameter_correlations_all_models.csv"),
  row.names = FALSE
)

rt_nat <- natural_summary_all %>%
  filter(parameter %in% c("alpha_pop", "delta_pop", "beta_a_pop", "beta_q_pop", "k_rep_pop"))

p_rt_cross <- ggplot(rt_nat, aes(median, model)) +
  geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0.15) +
  geom_point(size = 2) +
  facet_wrap(~ parameter, scales = "free_x", ncol = 2) +
  labs(
    title = "Stage 2 RT mechanism: cross-growth-model posterior comparison",
    subtitle = "Natural-scale posterior medians and 95% credible intervals",
    x = "Parameter value", y = NULL
  ) + theme_bw()

save_plot(p_rt_cross, "rt_mechanism_cross_model_intervals", 10, 8)


mechanistic_summary <- data.frame(
  item = c(
    "Control mice", "RT-only mice", "RT administrations",
    "Control RT-specific random effects",
    "RT damage update at fraction",
    "Damage between fractions",
    "Likelihood contributions",
    "Posterior predictive quantities"
  ),
  value = c(
    "23", "33", "228",
    "None: controls have rt_index = 0",
    "D receives an alpha_mouse x dose jump",
    "D decays according to k_rep",
    "56 baseline + 380 post-baseline = 436",
    "y_rep_baseline + y_rep"
  )
)

write.csv(
  mechanistic_summary,
  file.path(table_dir, "stage2_mechanistic_consistency_summary.csv"),
  row.names = FALSE
)

############################################################
# STAGE 1 vs STAGE 2 POSTERIOR COMPARISON
#to assess whether adding the RT cohort gives broadly compatible
# inference for the underlying untreated-growth distribution.
############################################################

stage1_fit_paths <- c(
  exponential = path.expand("~/Desktop/STADIO1/results/fit_exp_control.rds"),
  logistic    = path.expand("~/Desktop/STADIO1/results/fit_logistic_control.rds"),
  gompertz    = path.expand("~/Desktop/STADIO1/results/fit_gompertz_control.rds")
)

stage1_available <- all(file.exists(stage1_fit_paths))

if (stage1_available) {
  stage1_fits <- list(
    exponential = readRDS(stage1_fit_paths[["exponential"]]),
    logistic    = readRDS(stage1_fit_paths[["logistic"]]),
    gompertz    = readRDS(stage1_fit_paths[["gompertz"]])
  )
  
  stopifnot(all(vapply(stage1_fits, inherits, logical(1), what = "CmdStanMCMC")))
  
  stage_vars <- list(
    exponential = c("logV0_pop", "logr_pop", "sigma_logV0", "sigma_logr", "sigma"),
    logistic = c("logV0_pop", "logr_pop", "logK_pop",
                 "sigma_logV0", "sigma_logr", "sigma_logK", "sigma"),
    gompertz = c("logV0_pop", "logr_pop", "logK_pop",
                 "sigma_logV0", "sigma_logr", "sigma_logK", "sigma")
  )
  
  stage_comp <- list()
  stage_diff <- list()
  stage_nat <- list()
  
  for (m in names(fits)) {
    vars <- stage_vars[[m]]
    
    s1 <- stage1_fits[[m]]$summary(variables = vars) %>%
      as.data.frame() %>%
      mutate(stage = "Stage 1: Control only", model = model_label(m), .before = 1)
    
    s2 <- fits[[m]]$summary(variables = vars) %>%
      as.data.frame() %>%
      mutate(stage = "Stage 2: Control + RT-only", model = model_label(m), .before = 1)
    
    stage_comp[[m]] <- bind_rows(s1, s2)
    
    d1 <- as.data.frame(safe_draws(stage1_fits[[m]], vars)) %>%
      pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
      mutate(stage = "Stage 1: Control only")
    d2 <- as.data.frame(safe_draws(fits[[m]], vars)) %>%
      pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
      mutate(stage = "Stage 2: Control + RT-only")
    
    pd <- bind_rows(d1, d2)
    p <- ggplot(pd, aes(value, linetype = stage, group = stage)) +
      geom_density(linewidth = 0.8) +
      facet_wrap(~ parameter, scales = "free", ncol = 2) +
      labs(
        title = paste0(model_label(m), ": Stage 1 vs Stage 2 growth posterior"),
        subtitle = "Control-only fit versus joint Control + RT-only fit",
        x = "Parameter value", y = "Density", linetype = NULL
      ) + theme_bw() + theme(legend.position = "bottom")
    save_plot(p, paste0("stage1_vs_stage2_growth_posterior_", m), 10, 9)
    
    location_vars <- c("logV0_pop", "logr_pop", if (m != "exponential") "logK_pop" else NULL)
    a <- as.data.frame(safe_draws(stage1_fits[[m]], location_vars))
    b <- as.data.frame(safe_draws(fits[[m]], location_vars))
    
    make_nat <- function(x) {
      out <- data.frame(V0_pop = exp(x$logV0_pop), r_pop = exp(x$logr_pop))
      if (m != "exponential") out$K_pop <- exp(x$logK_pop)
      out
    }
    
    an <- make_nat(a) %>%
      pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
      mutate(stage = "Stage 1: Control only")
    bn <- make_nat(b) %>%
      pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
      mutate(stage = "Stage 2: Control + RT-only")
    
    sn <- bind_rows(an, bn) %>%
      group_by(stage, parameter) %>%
      summarise(
        mean = mean(value), sd = sd(value),
        q025 = quantile(value, .025), median = quantile(value, .5),
        q975 = quantile(value, .975), .groups = "drop"
      ) %>% mutate(model = model_label(m), .before = 1)
    stage_nat[[m]] <- sn
    
    pnat <- ggplot(sn, aes(median, stage)) +
      geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0.15) +
      geom_point(size = 2) +
      facet_wrap(~ parameter, scales = "free_x") +
      labs(
        title = paste0(model_label(m), ": Stage 1 vs Stage 2 on natural scales"),
        x = "Posterior median and 95% credible interval", y = NULL
      ) + theme_bw()
    save_plot(pnat, paste0("stage1_vs_stage2_growth_natural_scale_", m), 9, 6)
    
    nd <- min(nrow(a), nrow(b))
    ia <- sample(seq_len(nrow(a)), nd, replace = FALSE)
    ib <- sample(seq_len(nrow(b)), nd, replace = FALSE)
    
    diffs <- bind_rows(lapply(vars, function(v) {
      delta <- b[ib, v] - a[ia, v]
      data.frame(
        model = model_label(m), parameter = v,
        mean_difference_stage2_minus_stage1 = mean(delta),
        median_difference = median(delta),
        q025_difference = quantile(delta, .025),
        q975_difference = quantile(delta, .975),
        probability_stage2_greater_stage1 = mean(delta > 0)
      )
    }))
    stage_diff[[m]] <- diffs
  }
  
  write.csv(bind_rows(stage_comp),
            file.path(table_dir, "stage1_vs_stage2_posterior_summary.csv"), row.names = FALSE)
  write.csv(bind_rows(stage_nat),
            file.path(table_dir, "stage1_vs_stage2_natural_scale_summary.csv"), row.names = FALSE)
  write.csv(bind_rows(stage_diff),
            file.path(table_dir, "stage1_vs_stage2_posterior_difference_summary.csv"), row.names = FALSE)
  
} else {
  warning(
    "Stage-1 fit objects were not all found, so the optional Stage 1 vs Stage 2 comparison was skipped. Expected:\n",
    paste(stage1_fit_paths, collapse = "\n")
  )
  writeLines(
    c("Stage 1 vs Stage 2 comparison skipped: Stage-1 fit object(s) missing.",
      paste(stage1_fit_paths, file.exists(stage1_fit_paths))),
    file.path(analysis_dir, "stage1_comparison_SKIPPED.txt")
  )
}

############################################################
#  CROSS-MODEL FIT / PPC COMPARISON TABLE
############################################################

master_model_comparison <- loo_summary %>%
  left_join(convergence_summary, by = "model") %>%
  left_join(hmc_summary, by = "model") %>%
  left_join(residual_overall %>% select(model, rmse_log_residual, predictive_coverage_95),
            by = "model") %>%
  left_join(pareto_summary %>% select(model, n_k_gt_0_7, n_k_gt_1, max_pareto_k),
            by = "model") %>%
  arrange(desc(elpd_loo))

write.csv(
  master_model_comparison,
  file.path(table_dir, "master_model_comparison.csv"),
  row.names = FALSE
)

p_model <- ggplot(master_model_comparison,
                  aes(reorder(model, elpd_loo), elpd_loo)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = elpd_loo - se_elpd_loo,
                    ymax = elpd_loo + se_elpd_loo), width = 0.12) +
  coord_flip() +
  labs(
    title = "Stage 2 model comparison by PSIS-LOO",
    subtitle = "Point estimate ± one reported standard error of elpd_loo",
    x = NULL, y = "elpd_loo"
  ) + theme_bw()

save_plot(p_model, "model_comparison_elpd_loo", 8, 5.5)


integrity <- data.frame(
  check = c(
    "number_of_mice", "number_of_control_mice", "number_of_rt_mice",
    "number_of_observations", "number_of_baselines",
    "number_of_postbaseline_observations", "number_of_rt_events",
    "same_observation_order_all_models",
    "all_three_fit_objects_loaded",
    "analysis_refits_stan"
  ),
  value = c(
    56, 23, 33, 436, 56, 380, 228,
    TRUE, TRUE, FALSE
  )
)

write.csv(integrity, file.path(table_dir, "analysis_integrity_checks.csv"), row.names = FALSE)

# Raw RT-only tumour volume trajectories and RT events


rt_event_times <-
  rt_events %>%
  distinct(
    mouse,
    time_days
  )

raw_rt_only_plot <-
  obs_all %>%
  filter(is_rt) %>%
  ggplot(
    aes(
      x = time_days,
      y = volume,
      group = mouse
    )
  ) +
  geom_vline(
    data = rt_event_times,
    aes(xintercept = time_days),
    linetype = "dashed",
    linewidth = 0.3,
    inherit.aes = FALSE
  ) +
  geom_line(
    linewidth = 0.5
  ) +
  geom_point(
    size = 1.5
  ) +
  facet_wrap(
    ~ mouse,
    scales = "free_y",
    ncol = 6
  ) +
  labs(
    x = "Time (days)",
    y = expression("Tumour volume (mm"^3 * ")")
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    strip.text = element_text(size = 8)
  )

raw_rt_only_plot

ggsave(
  filename = "raw_rt_only_trajectories.pdf",
  plot = raw_rt_only_plot,
  width = 10,
  height = 8
)


############################################################
# Posterior parameter dependence / geometry
############################################################

library(posterior)
library(dplyr)
library(ggplot2)

dir.create(
  "results/stage2_parameter_geometry",
  recursive = TRUE,
  showWarnings = FALSE
)

exp_geom <- as.data.frame(
  posterior::as_draws_df(
    fits$exponential$draws(
      variables = c(
        "logr_pop",
        "logalpha_pop",
        "logbeta_a_pop",
        "logbeta_q_pop",
        "logk_rep_pop"
      )
    )
  )
)

log_geom <- as.data.frame(
  posterior::as_draws_df(
    fits$logistic$draws(
      variables = c(
        "logr_pop",
        "logK_pop",
        "logalpha_pop",
        "logbeta_a_pop",
        "logbeta_q_pop",
        "logk_rep_pop"
      )
    )
  )
)

gomp_geom <- as.data.frame(
  posterior::as_draws_df(
    fits$gompertz$draws(
      variables = c(
        "logr_pop",
        "logK_pop",
        "logalpha_pop",
        "logbeta_a_pop",
        "logbeta_q_pop",
        "logk_rep_pop"
      )
    )
  )
)


############################################################
# r-K posterior dependence Logistic and Gompertz 
############################################################

log_rK_cor <- cor(
  log_geom$logr_pop,
  log_geom$logK_pop
)

gomp_rK_cor <- cor(
  gomp_geom$logr_pop,
  gomp_geom$logK_pop
)

cat(
  "Logistic: cor(log r_pop, log K_pop) =",
  round(log_rK_cor, 3),
  "\n"
)

cat(
  "Gompertz: cor(log r_pop, log K_pop) =",
  round(gomp_rK_cor, 3),
  "\n"
)

rK_correlations <- data.frame(
  model = c("Logistic", "Gompertz"),
  correlation = c(
    log_rK_cor,
    gomp_rK_cor
  )
)

write.csv(
  rK_correlations,
  "results/stage2_parameter_geometry/rK_correlations.csv",
  row.names = FALSE
)



#  Plot r-K posterior geometry


rK_plot_data <- bind_rows(
  data.frame(
    model = "Logistic",
    log_K_pop = log_geom$logK_pop,
    log_r_pop = log_geom$logr_pop
  ),
  data.frame(
    model = "Gompertz",
    log_K_pop = gomp_geom$logK_pop,
    log_r_pop = gomp_geom$logr_pop
  )
)

rK_plot <- ggplot(
  rK_plot_data,
  aes(
    x = log_K_pop,
    y = log_r_pop
  )
) +
  geom_point(
    alpha = 0.15,
    size = 0.7
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.8
  ) +
  facet_wrap(
    ~model,
    scales = "free"
  ) +
  labs(
    x = expression(log(K[pop])),
    y = expression(log(r[pop]))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank()
  )

print(rK_plot)

ggsave(
  "results/stage2_parameter_geometry/stage2_rK_geometry.pdf",
  rK_plot,
  width = 8,
  height = 4.5
)



# beta_a - beta_q posterior dependence

beta_correlations <- data.frame(
  model = c(
    "Exponential",
    "Logistic",
    "Gompertz"
  ),
  correlation = c(
    cor(
      exp_geom$logbeta_a_pop,
      exp_geom$logbeta_q_pop
    ),
    cor(
      log_geom$logbeta_a_pop,
      log_geom$logbeta_q_pop
    ),
    cor(
      gomp_geom$logbeta_a_pop,
      gomp_geom$logbeta_q_pop
    )
  )
)


print(
  beta_correlations,
  row.names = FALSE
)

write.csv(
  beta_correlations,
  "results/stage2_parameter_geometry/beta_a_beta_q_correlations.csv",
  row.names = FALSE
)


beta_plot_data <- bind_rows(
  data.frame(
    model = "Exponential",
    log_beta_a = exp_geom$logbeta_a_pop,
    log_beta_q = exp_geom$logbeta_q_pop
  ),
  data.frame(
    model = "Logistic",
    log_beta_a = log_geom$logbeta_a_pop,
    log_beta_q = log_geom$logbeta_q_pop
  ),
  data.frame(
    model = "Gompertz",
    log_beta_a = gomp_geom$logbeta_a_pop,
    log_beta_q = gomp_geom$logbeta_q_pop
  )
)

beta_plot <- ggplot(
  beta_plot_data,
  aes(
    x = log_beta_a,
    y = log_beta_q
  )
) +
  geom_point(
    alpha = 0.15,
    size = 0.7
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.8
  ) +
  facet_wrap(
    ~model,
    scales = "free"
  ) +
  labs(
    x = expression(log(beta[a*",pop"])),
    y = expression(log(beta[q*",pop"]))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank()
  )

print(beta_plot)

ggsave(
  "results/stage2_parameter_geometry/stage2_beta_a_beta_q_geometry.pdf",
  beta_plot,
  width = 9,
  height = 4
)

#  alpha-k_rep posterior dependence


alpha_krep_correlations <- data.frame(
  model = c(
    "Exponential",
    "Logistic",
    "Gompertz"
  ),
  correlation = c(
    cor(
      exp_geom$logalpha_pop,
      exp_geom$logk_rep_pop
    ),
    cor(
      log_geom$logalpha_pop,
      log_geom$logk_rep_pop
    ),
    cor(
      gomp_geom$logalpha_pop,
      gomp_geom$logk_rep_pop
    )
  )
)

print(
  alpha_krep_correlations,
  row.names = FALSE
)

write.csv(
  alpha_krep_correlations,
  "results/stage2_parameter_geometry/alpha_krep_correlations.csv",
  row.names = FALSE
)

alpha_krep_plot_data <- bind_rows(
  data.frame(
    model = "Exponential",
    log_alpha = exp_geom$logalpha_pop,
    log_krep = exp_geom$logk_rep_pop
  ),
  data.frame(
    model = "Logistic",
    log_alpha = log_geom$logalpha_pop,
    log_krep = log_geom$logk_rep_pop
  ),
  data.frame(
    model = "Gompertz",
    log_alpha = gomp_geom$logalpha_pop,
    log_krep = gomp_geom$logk_rep_pop
  )
)

alpha_krep_plot <- ggplot(
  alpha_krep_plot_data,
  aes(
    x = log_krep,
    y = log_alpha
  )
) +
  geom_point(
    alpha = 0.15,
    size = 0.7
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.8
  ) +
  facet_wrap(
    ~model,
    scales = "free"
  ) +
  labs(
    x = expression(log(k[rep*",pop"])),
    y = expression(log(alpha[pop]))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank()
  )

print(alpha_krep_plot)

ggsave(
  "results/stage2_parameter_geometry/stage2_alpha_krep_geometry.pdf",
  alpha_krep_plot,
  width = 9,
  height = 4
)

#  Combined numerical summary

geometry_summary <- bind_rows(
  data.frame(
    relationship = "log(r_pop) vs log(K_pop)",
    model = rK_correlations$model,
    correlation = rK_correlations$correlation
  ),
  data.frame(
    relationship = "log(beta_a_pop) vs log(beta_q_pop)",
    model = beta_correlations$model,
    correlation = beta_correlations$correlation
  ),
  data.frame(
    relationship = "log(alpha_pop) vs log(k_rep_pop)",
    model = alpha_krep_correlations$model,
    correlation = alpha_krep_correlations$correlation
  )
)

print(
  geometry_summary,
  row.names = FALSE
)

write.csv(
  geometry_summary,
  "results/stage2_parameter_geometry/posterior_dependence_summary.csv",
  row.names = FALSE
)



#  MASTER CONSOLE SUMMARY


master_summary_file <- file.path(analysis_dir, "STAGE2_MASTER_SUMMARY.txt")

capture.output({
  cat("===== STAGE 2 CONTROL + RT POSTERIOR ANALYSIS =====\n\n")
  cat("All new outputs saved under:\n", analysis_dir, "\n\n")
  
  cat("===== DATA =====\n")
  print(integrity)
  
  cat("\n===== HMC SUMMARY =====\n")
  print(hmc_summary)
  
  cat("\n===== CONVERGENCE SUMMARY =====\n")
  print(convergence_summary)
  
  cat("\n===== NATURAL-SCALE POPULATION PARAMETERS =====\n")
  print(natural_summary_all)
  
  cat("\n===== RESIDUAL SUMMARY =====\n")
  print(residual_overall)
  
  cat("\n===== RESIDUAL SUMMARY BY COHORT =====\n")
  print(residual_by_group)
  
  cat("\n===== POSTERIOR PREDICTIVE COVERAGE BY COHORT =====\n")
  print(coverage_by_group)
  
  cat("\n===== LOO SUMMARY =====\n")
  print(loo_summary)
  
  cat("\n===== PAIRED LOO COMPARISON =====\n")
  print(loo_cmp)
  
  cat("\n===== PARETO-k SUMMARY =====\n")
  print(pareto_summary)
  
  cat("\n===== MASTER MODEL COMPARISON =====\n")
  print(master_model_comparison)
  
  cat("\n===== JAMIE PRIOR-vs-POSTERIOR =====\n")
  cat("Completed for every model: growth, RT mechanism, and natural-scale population parameters.\n")
  
  cat("\n===== STAGE 1 vs STAGE 2 =====\n")
  cat("Stage-1 comparison available:", stage1_available, "\n")
}, file = master_summary_file)

capture.output(
  sessionInfo(),
  file = file.path(analysis_dir, "R_sessionInfo.txt")
)
