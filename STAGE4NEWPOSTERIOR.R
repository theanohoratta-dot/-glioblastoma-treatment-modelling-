############################################################
# STAGE4NEWPOSTERIOR — FINAL STAGE 4 POSTERIOR ANALYSIS
############################################################
suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(tidyr)
  library(posterior)
  library(ggplot2)
  library(bayesplot)
  library(loo)
})

set.seed(123)
n_prior_stage4 <- 50000

fit_root <- path.expand(Sys.getenv(
  "STAGE4_FIT_ROOT",
  unset = "~/Desktop/STADIO4/FINAL_FITS"
))

analysis_dir <- path.expand(Sys.getenv(
  "STAGE4_NEWRESULTS",
  unset = "~/Desktop/STADIO4/NEWRESULTS"
))


data_env <- Sys.getenv("STAGE4_DATA", unset = "")
data_candidates <- unique(path.expand(c(
  data_env,
  "~/Desktop/STADIO4/lemassonParpi.csv",
  "~/Desktop/STADIO4/lemassonParpi (2)(20260913-202754).csv",
  "~/Desktop/dissertation/data/lemassonParpi.csv",
  "~/Desktop/lemassonParpi.csv"
)))
data_candidates <- data_candidates[nzchar(data_candidates)]
data_hits <- data_candidates[file.exists(data_candidates)]
if (length(data_hits) < 1L) {
  stop(
    "Could not find lemassonParpi.csv. Set STAGE4_DATA or copy the CSV to ",
    "~/Desktop/STADIO4/lemassonParpi.csv"
  )
}
data_path <- normalizePath(data_hits[1], mustWork = TRUE)

# Helper: locate exactly one completed file anywhere below fit_root.
find_unique_stage4_file <- function(filename) {
  hits <- list.files(
    fit_root,
    pattern = paste0("^", filename, "$"),
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(hits) != 1L) {
    stop(
      "Expected exactly one file named '", filename,
      "' below ", fit_root,
      ", but found ", length(hits), ".\n",
      if (length(hits) > 0) paste(hits, collapse = "\n") else ""
    )
  }
  normalizePath(hits, mustWork = TRUE)
}

fit_exp_path <- find_unique_stage4_file(
  "fit_stage4NEW_exp_reducesum.rds"
)
fit_log_path <- find_unique_stage4_file(
  "fit_stage4NEW_log_reducesum.rds"
)
fit_gomp_path <- find_unique_stage4_file(
  "fit_stage4NEW_gomp_reducesum.rds"
)

figure_dir   <- file.path(analysis_dir, "figures")
table_dir    <- file.path(analysis_dir, "tables")
object_dir   <- file.path(analysis_dir, "objects")

for (d in c(analysis_dir, figure_dir, table_dir, object_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("\n===== STAGE 4 FINAL POSTERIOR ANALYSIS =====\n")
cat("Fit root     :", fit_root, "\n")
cat("NEWRESULTS   :", analysis_dir, "\n")
cat("Data         :", data_path, "\n\n")



required_files <- c(
  data_path,
  fit_exp_path,
  fit_log_path,
  fit_gomp_path
)

file_check <- data.frame(
  file = required_files,
  exists = file.exists(required_files)
)

print(
  file_check,
  row.names = FALSE
)

if (!all(file_check$exists)) {
  stop(
    "At least one required Stage 4 analysis file is missing. ",
    "Check the paths printed above before continuing."
  )
}

fit_exp_stage4 <- readRDS(
  fit_exp_path
)

fit_log_stage4 <- readRDS(
  fit_log_path
)

fit_gomp_stage4 <- readRDS(
  fit_gomp_path
)


fit_classes <- data.frame(
  model = c(
    "Exponential",
    "Logistic",
    "Gompertz"
  ),
  class = c(
    paste(class(fit_exp_stage4), collapse = " / "),
    paste(class(fit_log_stage4), collapse = " / "),
    paste(class(fit_gomp_stage4), collapse = " / ")
  )
)

print(
  fit_classes,
  row.names = FALSE
)

cat("\nDraw dimensions:\n")

cat(
  "Exponential:",
  paste(dim(fit_exp_stage4$draws()), collapse = " x "),
  "\n"
)

cat(
  "Logistic:",
  paste(dim(fit_log_stage4$draws()), collapse = " x "),
  "\n"
)

cat(
  "Gompertz:",
  paste(dim(fit_gomp_stage4$draws()), collapse = " x "),
  "\n"
)


data_raw <- read.csv(
  data_path
)

required_columns <- c(
  "mouse",
  "time",
  "volume",
  "observation",
  "dose",
  "modality",
  "study"
)

if (!all(required_columns %in% names(data_raw))) {
  stop(
    "The dataset does not contain all columns expected by ",
    "the fitted Stage 4 model."
  )
}

cat(
  "\nDataset loaded successfully:",
  nrow(data_raw),
  "rows and",
  ncol(data_raw),
  "columns.\n"
)

cat(
  "Distinct mice:",
  dplyr::n_distinct(data_raw$mouse),
  "\n"
)


data_stage4 <- data_raw %>%
  mutate(
    time_days = time / 24
  )

stopifnot(
  all(is.finite(data_stage4$time_days)),
  all(data_stage4$time_days >= 0)
)


mouse_treatment <- data_stage4 %>%
  group_by(mouse) %>%
  summarise(
    has_drug = any(modality == 1, na.rm = TRUE),
    has_rt   = any(modality == 2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    group = case_when(
      !has_drug & !has_rt ~ 1L,
      has_drug & !has_rt ~ 2L,
      !has_drug &  has_rt ~ 3L,
      has_drug &  has_rt ~ 4L
    ),
    group_name = c(
      "control",
      "drug_only",
      "rt_only",
      "combination"
    )[group]
  ) %>%
  arrange(mouse) %>%
  mutate(
    mouse_id = row_number()
  )

stopifnot(
  nrow(mouse_treatment) == 78,
  identical(
    mouse_treatment$mouse_id,
    seq_len(nrow(mouse_treatment))
  )
)


group_counts <- mouse_treatment %>%
  count(
    group,
    group_name,
    name = "n_mice"
  )

print(
  group_counts,
  n = Inf
)

stopifnot(
  group_counts$n_mice[group_counts$group == 1] == 23,
  group_counts$n_mice[group_counts$group == 2] == 13,
  group_counts$n_mice[group_counts$group == 3] == 33,
  group_counts$n_mice[group_counts$group == 4] == 9
)


obs_stage4 <- data_stage4 %>%
  filter(
    observation == 1
  ) %>%
  left_join(
    mouse_treatment %>%
      select(
        mouse,
        mouse_id,
        group,
        group_name
      ),
    by = "mouse"
  ) %>%
  arrange(
    mouse_id,
    time_days
  )

stopifnot(
  all(!is.na(obs_stage4$mouse_id)),
  all(is.finite(obs_stage4$volume)),
  all(obs_stage4$volume > 0)
)

cat(
  "\nNumber of tumour-volume observations:",
  nrow(obs_stage4),
  "\n"
)

cat(
  "Observation time range (days):",
  min(obs_stage4$time_days),
  "to",
  max(obs_stage4$time_days),
  "\n"
)

cat(
  "Observed volume range (mm^3):",
  min(obs_stage4$volume),
  "to",
  max(obs_stage4$volume),
  "\n"
)

obs_group_counts <- obs_stage4 %>%
  count(
    group,
    group_name,
    name = "n_observations"
  )

print(
  as.data.frame(obs_group_counts),
  row.names = FALSE
)

write.csv(
  group_counts,
  file.path(
    table_dir,
    "stage4_mouse_counts_by_group.csv"
  ),
  row.names = FALSE
)

write.csv(
  obs_group_counts,
  file.path(
    table_dir,
    "stage4_observation_counts_by_group.csv"
  ),
  row.names = FALSE
)

# Combination mice only
combination_obs <- obs_stage4 %>%
  filter(group == 4)

# RT administration times for combination mice
combination_rt_events <- data_stage4 %>%
  filter(
    mouse %in% combination_obs$mouse,
    observation == 0,
    modality == 2
  ) %>%
  distinct(
    mouse,
    time_days
  )

# Raw observed trajectories
raw_combination_plot <- ggplot(
  combination_obs,
  aes(
    x = time_days,
    y = volume,
    group = mouse
  )
) +
  geom_vline(
    data = combination_rt_events,
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
    ncol = 3
  ) +
  labs(
    x = "Time (days)",
    y = expression("Tumour volume (mm"^3 * ")")
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    strip.text = element_text(size = 9)
  )

print(raw_combination_plot)

ggsave(
  filename = file.path(
    figure_dir,
    "raw_combination_trajectories.pdf"
  ),
  plot = raw_combination_plot,
  width = 9,
  height = 7
)

ggsave(
  filename = file.path(
    figure_dir,
    "raw_combination_trajectories.png"
  ),
  plot = raw_combination_plot,
  width = 9,
  height = 7,
  dpi = 300
)

fits_stage4 <- list(
  Exponential = fit_exp_stage4,
  Logistic    = fit_log_stage4,
  Gompertz    = fit_gomp_stage4
)


hmc_summary_list <- lapply(
  names(fits_stage4),
  function(model_name) {
    
    fit <- fits_stage4[[model_name]]
    
    diag <- fit$diagnostic_summary()
    
    data.frame(
      model = model_name,
      chain = seq_along(diag$num_divergent),
      divergences = diag$num_divergent,
      max_treedepth = diag$num_max_treedepth,
      ebfmi = diag$ebfmi
    )
  }
)

hmc_summary <- bind_rows(
  hmc_summary_list
)

print(
  as.data.frame(hmc_summary),
  row.names = FALSE
)

write.csv(
  hmc_summary,
  file.path(
    table_dir,
    "stage4_hmc_diagnostics.csv"
  ),
  row.names = FALSE
)

hmc_model_summary <- hmc_summary %>%
  group_by(model) %>%
  summarise(
    total_divergences = sum(divergences),
    total_max_treedepth = sum(max_treedepth),
    min_ebfmi = min(ebfmi),
    .groups = "drop"
  )

print(
  as.data.frame(hmc_model_summary),
  row.names = FALSE
)

write.csv(
  hmc_model_summary,
  file.path(
    table_dir,
    "stage4_hmc_model_summary.csv"
  ),
  row.names = FALSE
)
############################################################
# PARAMETER CONVERGENCE — EXPONENTIAL
############################################################

params_exp <- c(
  "logV0_pop",
  "logr_pop",
  
  "sigma_logV0",
  "sigma_logr",
  
  "logdelta_pop",
  "logalpha_pop",
  "logbeta_a_pop",
  "logbeta_q_pop",
  "logk_rep_pop",
  
  "sigma_logdelta",
  "sigma_logalpha",
  "sigma_logbeta_a",
  "sigma_logbeta_q",
  "sigma_logk_rep",
  
  "logk_abs_pop",
  "logk_cl_pop",
  "logk_repi_pop",
  "sigma_logk_abs",
  "sigma_logk_cl",
  "sigma_logk_repi",
  
  "sigma"
)

conv_exp <- fit_exp_stage4$summary(
  variables = params_exp
) %>%
  mutate(
    model = "Exponential",
    .before = 1
  )

print(
  as.data.frame(
    conv_exp[, c(
      "model",
      "variable",
      "rhat",
      "ess_bulk",
      "ess_tail"
    )]
  ),
  row.names = FALSE
)
############################################################
# KEY PARAMETER CONVERGENCE — LOGISTIC
############################################################

params_log <- c(
  "logV0_pop",
  "logr_pop",
  "logK_pop",
  
  "sigma_logV0",
  "sigma_logr",
  "sigma_logK",
  
  "logdelta_pop",
  "logalpha_pop",
  "logbeta_a_pop",
  "logbeta_q_pop",
  "logk_rep_pop",
  
  "sigma_logdelta",
  "sigma_logalpha",
  "sigma_logbeta_a",
  "sigma_logbeta_q",
  "sigma_logk_rep",
  
  "logk_abs_pop",
  "logk_cl_pop",
  "logk_repi_pop",
  "sigma_logk_abs",
  "sigma_logk_cl",
  "sigma_logk_repi",
  
  "sigma"
)

conv_log <- fit_log_stage4$summary(
  variables = params_log
) %>%
  mutate(
    model = "Logistic",
    .before = 1
  )

print(
  as.data.frame(
    conv_log[, c(
      "model",
      "variable",
      "rhat",
      "ess_bulk",
      "ess_tail"
    )]
  ),
  row.names = FALSE
)
############################################################
#KEY PARAMETER CONVERGENCE — GOMPERTZ
############################################################

params_gomp <- c(
  "logV0_pop",
  "logr_pop",
  "logK_pop",
  
  "sigma_logV0",
  "sigma_logr",
  "sigma_logK",
  
  "logdelta_pop",
  "logalpha_pop",
  "logbeta_a_pop",
  "logbeta_q_pop",
  "logk_rep_pop",
  
  "sigma_logdelta",
  "sigma_logalpha",
  "sigma_logbeta_a",
  "sigma_logbeta_q",
  "sigma_logk_rep",
  
  "logk_abs_pop",
  "logk_cl_pop",
  "logk_repi_pop",
  "sigma_logk_abs",
  "sigma_logk_cl",
  "sigma_logk_repi",
  
  "sigma"
)

conv_gomp <- fit_gomp_stage4$summary(
  variables = params_gomp
) %>%
  mutate(
    model = "Gompertz",
    .before = 1
  )

print(
  as.data.frame(
    conv_gomp[, c(
      "model",
      "variable",
      "rhat",
      "ess_bulk",
      "ess_tail"
    )]
  ),
  row.names = FALSE
)

convergence_stage4 <- bind_rows(
  conv_exp,
  conv_log,
  conv_gomp
)

write.csv(
  convergence_stage4,
  file.path(
    table_dir,
    "stage4_convergence_all_models.csv"
  ),
  row.names = FALSE
)

convergence_summary <- convergence_stage4 %>%
  group_by(model) %>%
  summarise(
    max_rhat = max(rhat, na.rm = TRUE),
    min_ess_bulk = min(ess_bulk, na.rm = TRUE),
    min_ess_tail = min(ess_tail, na.rm = TRUE),
    .groups = "drop"
  )

print(
  as.data.frame(convergence_summary),
  row.names = FALSE
)

write.csv(
  convergence_summary,
  file.path(
    table_dir,
    "stage4_convergence_summary.csv"
  ),
  row.names = FALSE
)




log_check_params <- c(
  "sigma_logalpha",
  "sigma_logk_rep",
  "sigma_logr",
  "sigma_logbeta_a",
  "sigma_logK"
)

log_check_draws <- fit_log_stage4$draws(
  variables = log_check_params
)

# Trace plots
p_log_trace <- bayesplot::mcmc_trace(
  log_check_draws
)

ggsave(
  filename = file.path(
    figure_dir,
    "stage4_logistic_targeted_traceplots.pdf"
  ),
  plot = p_log_trace,
  width = 10,
  height = 8
)

print(p_log_trace)




p_log_rank <- bayesplot::mcmc_rank_overlay(
  posterior::as_draws_array(
    log_check_draws
  )
)

ggsave(
  filename = file.path(
    figure_dir,
    "stage4_logistic_targeted_rankplots.pdf"
  ),
  plot = p_log_rank,
  width = 10,
  height = 8
)

print(p_log_rank)



sigma_alpha_draws <- posterior::as_draws_df(
  fit_log_stage4$draws(
    variables = "sigma_logalpha"
  )
)

sigma_alpha_chain_summary <- sigma_alpha_draws %>%
  group_by(.chain) %>%
  summarise(
    mean = mean(sigma_logalpha),
    median = median(sigma_logalpha),
    sd = sd(sigma_logalpha),
    q05 = quantile(sigma_logalpha, 0.05),
    q95 = quantile(sigma_logalpha, 0.95),
    .groups = "drop"
  )

print(
  as.data.frame(sigma_alpha_chain_summary),
  row.names = FALSE
)
write.csv(
  sigma_alpha_chain_summary,
  file.path(
    table_dir,
    "stage4_logistic_sigma_logalpha_chain_summary.csv"
  ),
  row.names = FALSE
)
############################################################
# POSTERIOR PARAMETER SUMMARIES
############################################################


make_posterior_summary <- function(
    fit,
    variables,
    model_name
) {
  
  draws <- fit$draws(
    variables = variables
  )
  
  out <- posterior::summarise_draws(
    draws,
    mean,
    median,
    sd,
    q2.5 = function(x) {
      unname(quantile(x, 0.025))
    },
    q97.5 = function(x) {
      unname(quantile(x, 0.975))
    },
    rhat,
    ess_bulk,
    ess_tail
  )
  
  out %>%
    as.data.frame() %>%
    mutate(
      model = model_name,
      .before = 1
    )
}


posterior_summary_exp <- make_posterior_summary(
  fit = fit_exp_stage4,
  variables = params_exp,
  model_name = "Exponential"
)

posterior_summary_log <- make_posterior_summary(
  fit = fit_log_stage4,
  variables = params_log,
  model_name = "Logistic"
)

posterior_summary_gomp <- make_posterior_summary(
  fit = fit_gomp_stage4,
  variables = params_gomp,
  model_name = "Gompertz"
)


posterior_summary_all <- bind_rows(
  posterior_summary_exp,
  posterior_summary_log,
  posterior_summary_gomp
)

write.csv(
  posterior_summary_all,
  file.path(
    table_dir,
    "stage4_posterior_summary_log_scale.csv"
  ),
  row.names = FALSE
)
############################################################
# PRINT  POPULATION PARAMETERS
############################################################

key_population_parameters <- c(
  "logV0_pop",
  "logr_pop",
  "logK_pop",
  "logdelta_pop",
  "logalpha_pop",
  "logbeta_a_pop",
  "logbeta_q_pop",
  "logk_rep_pop",
  "logk_abs_pop",
  "logk_cl_pop",
  "logk_repi_pop",
  "sigma_logk_abs",
  "sigma_logk_cl",
  "sigma_logk_repi",
  "sigma"
)

posterior_key <- posterior_summary_all %>%
  filter(
    variable %in% key_population_parameters
  ) %>%
  select(
    model,
    variable,
    mean,
    median,
    sd,
    q2.5,
    q97.5,
    rhat,
    ess_bulk,
    ess_tail
  )

print(
  as.data.frame(posterior_key),
  row.names = FALSE
)
############################################################
#  NATURAL-SCALE POSTERIOR SUMMARIES
############################################################

natural_scale_summary <- function(
    fit,
    model_name,
    include_K = FALSE
) {
  
  parameters <- c(
    "logV0_pop",
    "logr_pop",
    "logdelta_pop",
    "logalpha_pop",
    "logbeta_a_pop",
    "logbeta_q_pop",
    "logk_rep_pop",
    "logk_abs_pop",
    "logk_cl_pop",
    "logk_repi_pop"
  )
  
  if (include_K) {
    parameters <- c(
      parameters,
      "logK_pop"
    )
  }
  
  draws <- posterior::as_draws_df(
    fit$draws(
      variables = parameters
    )
  )
  
  result <- lapply(
    parameters,
    function(par) {
      
      x <- exp(draws[[par]])
      
      data.frame(
        model = model_name,
        parameter = sub("^log", "", par),
        mean = mean(x),
        median = median(x),
        sd = sd(x),
        q2.5 = quantile(x, 0.025),
        q97.5 = quantile(x, 0.975)
      )
    }
  )
  
  bind_rows(result)
}
natural_exp <- natural_scale_summary(
  fit_exp_stage4,
  "Exponential",
  include_K = FALSE
)

natural_log <- natural_scale_summary(
  fit_log_stage4,
  "Logistic",
  include_K = TRUE
)

natural_gomp <- natural_scale_summary(
  fit_gomp_stage4,
  "Gompertz",
  include_K = TRUE
)

natural_stage4 <- bind_rows(
  natural_exp,
  natural_log,
  natural_gomp
)

print(
  as.data.frame(natural_stage4),
  row.names = FALSE
)

write.csv(
  natural_stage4,
  file.path(
    table_dir,
    "stage4_posterior_summary_natural_scale.csv"
  ),
  row.names = FALSE
)
############################################################
# 4. PRIOR VS POSTERIOR
############################################################

make_prior_posterior_data <- function(
    fit,
    model_name,
    parameter,
    prior_mean,
    prior_sd,
    n_prior = 4000
) {
  
  posterior_draws <- posterior::as_draws_df(
    fit$draws(
      variables = parameter
    )
  )[[parameter]]
  
  prior_draws <- rnorm(
    n_prior,
    mean = prior_mean,
    sd = prior_sd
  )
  
  bind_rows(
    data.frame(
      model = model_name,
      parameter = parameter,
      distribution = "Prior",
      value = prior_draws
    ),
    data.frame(
      model = model_name,
      parameter = parameter,
      distribution = "Posterior",
      value = posterior_draws
    )
  )
}
############################################################
#  PRIOR SPECIFICATIONS
############################################################

prior_spec_exp <- data.frame(
  parameter = c(
    "logV0_pop",
    "logr_pop",
    "logdelta_pop",
    "logalpha_pop",
    "logbeta_a_pop",
    "logbeta_q_pop",
    "logk_rep_pop",
    "logk_abs_pop",
    "logk_cl_pop",
    "logk_repi_pop"
  ),
  prior_mean = c(
    log(30),
    log(0.3),
    log(1.2),
    log(0.9),
    log(26.4),
    log(4.8),
    log(86.4),
    log(96),
    log(19.2),
    log(0.04)
  ),
  prior_sd = c(
    0.5,
    0.5,
    0.5,
    0.5,
    0.5,
    0.5,
    0.5,
    0.35,
    0.35,
    0.35
  )
)

prior_spec_log <- data.frame(
  parameter = c(
    "logV0_pop",
    "logr_pop",
    "logK_pop",
    "logdelta_pop",
    "logalpha_pop",
    "logbeta_a_pop",
    "logbeta_q_pop",
    "logk_rep_pop",
    "logk_abs_pop",
    "logk_cl_pop",
    "logk_repi_pop"
  ),
  prior_mean = c(
    log(30),
    log(0.4),
    log(500),
    log(1.2),
    log(0.9),
    log(26.4),
    log(4.8),
    log(86.4),
    log(96),
    log(19.2),
    log(0.04)
  ),
  prior_sd = c(
    0.5,
    0.5,
    0.7,
    0.5,
    0.5,
    0.5,
    0.5,
    0.5,
    0.35,
    0.35,
    0.35
  )
)

prior_spec_gomp <- prior_spec_log

prior_spec_gomp$prior_mean[
  prior_spec_gomp$parameter == "logr_pop"
] <- log(0.2)
############################################################
#  BUILD PRIOR VS POSTERIOR DATA
############################################################

build_prior_post_model <- function(
    fit,
    model_name,
    prior_spec
) {
  
  bind_rows(
    lapply(
      seq_len(nrow(prior_spec)),
      function(i) {
        
        make_prior_posterior_data(
          fit = fit,
          model_name = model_name,
          parameter = prior_spec$parameter[i],
          prior_mean = prior_spec$prior_mean[i],
          prior_sd = prior_spec$prior_sd[i]
        )
      }
    )
  )
}

prior_post_exp <- build_prior_post_model(
  fit_exp_stage4,
  "Exponential",
  prior_spec_exp
)

prior_post_log <- build_prior_post_model(
  fit_log_stage4,
  "Logistic",
  prior_spec_log
)

prior_post_gomp <- build_prior_post_model(
  fit_gomp_stage4,
  "Gompertz",
  prior_spec_gomp
)

prior_post_stage4 <- bind_rows(
  prior_post_exp,
  prior_post_log,
  prior_post_gomp
)
############################################################
#  PRIOR VS POSTERIOR FIGURES
############################################################

make_prior_post_plot <- function(
    dat,
    model_name
) {
  
  ggplot(
    dat,
    aes(
      x = value,
      linetype = distribution
    )
  ) +
    geom_density(
      linewidth = 0.8,
      adjust = 1
    ) +
    facet_wrap(
      ~ parameter,
      scales = "free",
      ncol = 3
    ) +
    labs(
      title = paste0(
        model_name,
        " model: Full joint treatment model — Prior vs posterior"
      ),
      x = "Parameter value (log scale)",
      y = "Density",
      linetype = NULL
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom"
    )
}


p_prior_post_exp <- make_prior_post_plot(
  prior_post_exp,
  "Exponential"
)

p_prior_post_log <- make_prior_post_plot(
  prior_post_log,
  "Logistic"
)

p_prior_post_gomp <- make_prior_post_plot(
  prior_post_gomp,
  "Gompertz"
)


ggsave(
  file.path(
    figure_dir,
    "stage4_prior_vs_posterior_exponential.pdf"
  ),
  p_prior_post_exp,
  width = 11,
  height = 9
)

ggsave(
  file.path(
    figure_dir,
    "stage4_prior_vs_posterior_logistic.pdf"
  ),
  p_prior_post_log,
  width = 11,
  height = 10
)

ggsave(
  file.path(
    figure_dir,
    "stage4_prior_vs_posterior_gompertz.pdf"
  ),
  p_prior_post_gomp,
  width = 11,
  height = 10
)

print(p_prior_post_exp)
print(p_prior_post_log)
print(p_prior_post_gomp)

############################################################
#  CLEAN PRIOR VS POSTERIOR FIGURES
############################################################

parameter_labels <- c(
  logV0_pop       = "log(V0_pop)",
  logr_pop        = "log(r_pop)",
  logK_pop  = "log(K_pop)",
  logdelta_pop    = "log(delta_pop)",
  logalpha_pop    = "log(alpha_pop)",
  logbeta_a_pop   = "log(beta_a_pop)",
  logbeta_q_pop   = "log(beta_q_pop)",
  logk_rep_pop    = "log(k_rep_pop)",
  logk_abs_pop    = "log(k_abs,pop)",
  logk_cl_pop     = "log(k_cl,pop)",
  logk_repi_pop   = "log(k_repi,pop)"
)

growth_parameters <- c(
  "logV0_pop",
  "logr_pop",
  "logK_pop"
)

treatment_parameters <- c(
  "logdelta_pop",
  "logalpha_pop",
  "logbeta_a_pop",
  "logbeta_q_pop",
  "logk_rep_pop",
  "logk_abs_pop",
  "logk_cl_pop",
  "logk_repi_pop"
)


make_clean_prior_post_plot <- function(
    dat,
    model_name,
    parameters,
    plot_title,
    ncol = 2
) {
  
  plot_dat <- dat %>%
    filter(parameter %in% parameters)
  
  ggplot(
    plot_dat,
    aes(
      x = value,
      linetype = distribution
    )
  ) +
    geom_density(
      linewidth = 0.9,
      adjust = 1
    ) +
    facet_wrap(
      ~ parameter,
      scales = "free",
      ncol = ncol,
      labeller = as_labeller(parameter_labels)
    ) +
    scale_linetype_manual(
      values = c(
        "Prior" = "dashed",
        "Posterior" = "solid"
      )
    ) +
    labs(
      title = paste0(
        model_name,
        " model: Full joint treatment model — ",
        plot_title
      ),
      x = "Parameter value (log scale)",
      y = "Density",
      linetype = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 11),
      plot.title = element_text(
        size = 14,
        face = "bold"
      ),
      panel.spacing = grid::unit(1.2, "lines")
    )
}
############################################################
#  GROWTH-PARAMETER PRIOR VS POSTERIOR
############################################################

p_growth_exp <- make_clean_prior_post_plot(
  prior_post_exp,
  "Exponential",
  intersect(
    growth_parameters,
    unique(prior_post_exp$parameter)
  ),
  "Growth parameters",
  ncol = 2
)

p_growth_log <- make_clean_prior_post_plot(
  prior_post_log,
  "Logistic",
  intersect(
    growth_parameters,
    unique(prior_post_log$parameter)
  ),
  "Growth parameters",
  ncol = 3
)

p_growth_gomp <- make_clean_prior_post_plot(
  prior_post_gomp,
  "Gompertz",
  intersect(
    growth_parameters,
    unique(prior_post_gomp$parameter)
  ),
  "Growth parameters",
  ncol = 3
)

ggsave(
  file.path(
    figure_dir,
    "stage4_prior_vs_posterior_growth_exponential.pdf"
  ),
  p_growth_exp,
  width = 9,
  height = 5
)

ggsave(
  file.path(
    figure_dir,
    "stage4_prior_vs_posterior_growth_logistic.pdf"
  ),
  p_growth_log,
  width = 11,
  height = 4.5
)

ggsave(
  file.path(
    figure_dir,
    "stage4_prior_vs_posterior_growth_gompertz.pdf"
  ),
  p_growth_gomp,
  width = 11,
  height = 4.5
)
############################################################
# TREATMENT-PARAMETER PRIOR VS POSTERIOR
############################################################

p_treat_exp <- make_clean_prior_post_plot(
  prior_post_exp,
  "Exponential",
  treatment_parameters,
  "Treatment-mechanism parameters",
  ncol = 2
)

p_treat_log <- make_clean_prior_post_plot(
  prior_post_log,
  "Logistic",
  treatment_parameters,
  "Treatment-mechanism parameters",
  ncol = 2
)

p_treat_gomp <- make_clean_prior_post_plot(
  prior_post_gomp,
  "Gompertz",
  treatment_parameters,
  "Treatment-mechanism parameters",
  ncol = 2
)

ggsave(
  file.path(
    figure_dir,
    "stage4_prior_vs_posterior_treatment_exponential.pdf"
  ),
  p_treat_exp,
  width = 10,
  height = 11
)

ggsave(
  file.path(
    figure_dir,
    "stage4_prior_vs_posterior_treatment_logistic.pdf"
  ),
  p_treat_log,
  width = 10,
  height = 11
)

ggsave(
  file.path(
    figure_dir,
    "stage4_prior_vs_posterior_treatment_gompertz.pdf"
  ),
  p_treat_gomp,
  width = 10,
  height = 11
)

print(p_growth_exp)
print(p_growth_log)
print(p_growth_gomp)

print(p_treat_exp)
print(p_treat_log)
print(p_treat_gomp)

############################################################
# HIERARCHICAL SD PRIOR VS POSTERIOR
############################################################

make_sd_prior_post_data <- function(
    fit,
    model_name,
    parameters,
    prior_sd = 0.5,
    n_prior = 4000
) {
  
  bind_rows(
    lapply(
      parameters,
      function(par) {
        
        posterior_draws <- posterior::as_draws_df(
          fit$draws(
            variables = par
          )
        )[[par]]
        
        # Half-normal N(0, prior_sd), constrained > 0
        prior_draws <- abs(
          rnorm(
            n_prior,
            mean = 0,
            sd = prior_sd
          )
        )
        
        bind_rows(
          data.frame(
            model = model_name,
            parameter = par,
            distribution = "Prior",
            value = prior_draws
          ),
          data.frame(
            model = model_name,
            parameter = par,
            distribution = "Posterior",
            value = posterior_draws
          )
        )
      }
    )
  )
}
############################################################
# HIERARCHICAL SD PARAMETER SETS
############################################################

sd_params_exp <- c(
  "sigma_logV0",
  "sigma_logr",
  "sigma_logdelta",
  "sigma_logalpha",
  "sigma_logbeta_a",
  "sigma_logbeta_q",
  "sigma_logk_rep"
)

sd_params_log <- c(
  "sigma_logV0",
  "sigma_logr",
  "sigma_logK",
  "sigma_logdelta",
  "sigma_logalpha",
  "sigma_logbeta_a",
  "sigma_logbeta_q",
  "sigma_logk_rep"
)

sd_params_gomp <- sd_params_log



sd_prior_post_exp <- make_sd_prior_post_data(
  fit_exp_stage4,
  "Exponential",
  sd_params_exp
)

sd_prior_post_log <- make_sd_prior_post_data(
  fit_log_stage4,
  "Logistic",
  sd_params_log
)

sd_prior_post_gomp <- make_sd_prior_post_data(
  fit_gomp_stage4,
  "Gompertz",
  sd_params_gomp
)
############################################################
#  HIERARCHICAL SD PRIOR VS POSTERIOR FIGURES
############################################################

sd_parameter_labels <- c(
  sigma_logV0     = "sigma_logV0",
  sigma_logr      = "sigma_logr",
  sigma_logK      = "sigma_logK",
  sigma_logdelta  = "sigma_logdelta",
  sigma_logalpha  = "sigma_logalpha",
  sigma_logbeta_a = "sigma_logbeta_a",
  sigma_logbeta_q = "sigma_logbeta_q",
  sigma_logk_rep  = "sigma_logk_rep"
)

make_sd_plot <- function(
    dat,
    model_name
) {
  
  ggplot(
    dat,
    aes(
      x = value,
      linetype = distribution
    )
  ) +
    geom_density(
      linewidth = 0.9
    ) +
    facet_wrap(
      ~ parameter,
      scales = "free",
      ncol = 2,
      labeller = as_labeller(sd_parameter_labels)
    ) +
    scale_linetype_manual(
      values = c(
        "Prior" = "dashed",
        "Posterior" = "solid"
      )
    ) +
    labs(
      title = paste0(
        model_name,
        " model: Full joint treatment model — Hierarchical variability"
      ),
      x = "Hierarchical standard deviation",
      y = "Density",
      linetype = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 11),
      plot.title = element_text(
        size = 14,
        face = "bold"
      ),
      panel.spacing = grid::unit(1.2, "lines")
    )
}


p_sd_exp <- make_sd_plot(
  sd_prior_post_exp,
  "Exponential"
)

p_sd_log <- make_sd_plot(
  sd_prior_post_log,
  "Logistic"
)

p_sd_gomp <- make_sd_plot(
  sd_prior_post_gomp,
  "Gompertz"
)


ggsave(
  file.path(
    figure_dir,
    "stage4_prior_vs_posterior_sd_exponential.pdf"
  ),
  p_sd_exp,
  width = 10,
  height = 10
)

ggsave(
  file.path(
    figure_dir,
    "stage4_prior_vs_posterior_sd_logistic.pdf"
  ),
  p_sd_log,
  width = 10,
  height = 11
)

ggsave(
  file.path(
    figure_dir,
    "stage4_prior_vs_posterior_sd_gompertz.pdf"
  ),
  p_sd_gomp,
  width = 10,
  height = 11
)

print(p_sd_exp)
print(p_sd_log)
print(p_sd_gomp)
############################################################
#  POSTERIOR CHECK OF THE EFFECTIVE DNA-REPAIR MULTIPLIER
#      R(t) = k_repi * C2(t)
############################################################


pk_grid_step_days <- 0.02

combo_drug_events <- data_stage4 %>%
  filter(
    modality == 1,
    mouse %in% mouse_treatment$mouse[mouse_treatment$group == 4]
  ) %>%
  transmute(
    mouse,
    time_days,
    drug_dose = dose
  ) %>%
  group_by(
    mouse,
    time_days
  ) %>%
  summarise(
    drug_dose = sum(drug_dose),
    .groups = "drop"
  ) %>%
  arrange(
    mouse,
    time_days
  )

stopifnot(
  n_distinct(combo_drug_events$mouse) == 9,
  all(is.finite(combo_drug_events$time_days)),
  all(combo_drug_events$time_days >= 0),
  all(is.finite(combo_drug_events$drug_dose)),
  all(combo_drug_events$drug_dose > 0)
)

combo_followup <- obs_stage4 %>%
  filter(
    group == 4
  ) %>%
  group_by(
    mouse,
    mouse_id
  ) %>%
  summarise(
    final_time = max(time_days),
    .groups = "drop"
  ) %>%
  arrange(
    mouse_id
  )

stopifnot(
  nrow(combo_followup) == 9,
  all(combo_followup$final_time > 0)
)

propagate_pk_vectorised <- function(
    C1,
    C2,
    dt,
    k_abs,
    k_cl
) {
  if (dt < 0) {
    stop("Negative PK propagation interval encountered.")
  }

  if (dt == 0) {
    return(list(C1 = C1, C2 = C2))
  }

  e_abs <- exp(-k_abs * dt)
  e_cl  <- exp(-k_cl  * dt)

  C1_new <- C1 * e_abs
  C2_new <- numeric(length(C2))

  near_equal <- abs(k_cl - k_abs) < 1e-8

  if (any(!near_equal)) {
    ii <- !near_equal
    C2_new[ii] <-
      C2[ii] * e_cl[ii] +
      k_abs[ii] * C1[ii] *
      (e_abs[ii] - e_cl[ii]) /
      (k_cl[ii] - k_abs[ii])
  }

  if (any(near_equal)) {
    ii <- near_equal
    k_mid <- 0.5 * (k_abs[ii] + k_cl[ii])
    C2_new[ii] <-
      exp(-k_mid * dt) *
      (C2[ii] + k_mid * C1[ii] * dt)
  }

  C1_new <- pmax(C1_new, 0)
  C2_new <- pmax(C2_new, 0)

  list(
    C1 = C1_new,
    C2 = C2_new
  )
}

check_repair_multiplier_one_model <- function(
    fit,
    model_name,
    grid_step_days = pk_grid_step_days
) {

  n_comb <- nrow(combo_followup)
  stopifnot(n_comb == 9L)

  mouse_draw_results <- vector("list", n_comb)
  n_draws <- NULL
  overall_max_ratio <- NULL
  n_grid_above_one <- NULL
  n_grid_total <- NULL

  for (mm in seq_len(n_comb)) {
    mouse_now <- combo_followup$mouse[mm]
    mouse_id_now <- combo_followup$mouse_id[mm]
    final_time_now <- combo_followup$final_time[mm]

    c_index <- mm

    pk_draws <- posterior::as_draws_df(
      fit$draws(
        variables = c(
          paste0("k_abs_comb[", c_index, "]"),
          paste0("k_cl_comb[", c_index, "]"),
          paste0("k_repi_comb[", c_index, "]")
        )
      )
    )

    k_abs <- pk_draws[[paste0("k_abs_comb[", c_index, "]")]]
    k_cl <- pk_draws[[paste0("k_cl_comb[", c_index, "]")]]
    k_repi <- pk_draws[[paste0("k_repi_comb[", c_index, "]")]]

    if (is.null(n_draws)) {
      n_draws <- length(k_repi)
      overall_max_ratio <- rep(0, n_draws)
      n_grid_above_one <- integer(n_draws)
      n_grid_total <- integer(n_draws)
    }

    stopifnot(
      length(k_abs) == n_draws,
      length(k_cl) == n_draws,
      length(k_repi) == n_draws,
      all(is.finite(k_abs)),
      all(is.finite(k_cl)),
      all(is.finite(k_repi)),
      all(k_abs > 0),
      all(k_cl > 0),
      all(k_repi > 0)
    )

    events_now <- combo_drug_events %>%
      filter(
        mouse == mouse_now,
        time_days <= final_time_now
      ) %>%
      arrange(time_days)

    regular_grid <- seq(
      from = 0,
      to = final_time_now,
      by = grid_step_days
    )

    timeline <- sort(unique(c(
      0,
      regular_grid,
      final_time_now,
      events_now$time_days
    )))

    C1 <- rep(0, n_draws)
    C2 <- rep(0, n_draws)
    max_ratio_mouse <- rep(0, n_draws)
    n_above_mouse <- integer(n_draws)

    previous_time <- timeline[1]

    for (tt in seq_along(timeline)) {
      current_time <- timeline[tt]

      if (tt > 1) {
        propagated <- propagate_pk_vectorised(
          C1 = C1,
          C2 = C2,
          dt = current_time - previous_time,
          k_abs = k_abs,
          k_cl = k_cl
        )
        C1 <- propagated$C1
        C2 <- propagated$C2
      }

      # Drug dose is added to C1 at the event/segment start.
      dose_here <- events_now$drug_dose[
        abs(events_now$time_days - current_time) < 1e-10
      ]
      if (length(dose_here) > 0) {
        C1 <- C1 + sum(dose_here)
      }

      repair_ratio <- k_repi * C2

      max_ratio_mouse <- pmax(max_ratio_mouse, repair_ratio)
      n_above_mouse <- n_above_mouse + as.integer(repair_ratio > 1)

      previous_time <- current_time
    }

    overall_max_ratio <- pmax(overall_max_ratio, max_ratio_mouse)
    n_grid_above_one <- n_grid_above_one + n_above_mouse
    n_grid_total <- n_grid_total + length(timeline)

    mouse_draw_results[[mm]] <- data.frame(
      model = model_name,
      draw = seq_len(n_draws),
      comb_index = c_index,
      mouse = mouse_now,
      mouse_id = mouse_id_now,
      max_krepi_C2 = max_ratio_mouse,
      ever_above_1 = max_ratio_mouse > 1,
      fraction_grid_above_1 = n_above_mouse / length(timeline)
    )
  }

  mouse_draw_results <- bind_rows(mouse_draw_results)

  overall_draw_results <- data.frame(
    model = model_name,
    draw = seq_len(n_draws),
    max_krepi_C2 = overall_max_ratio,
    ever_above_1 = overall_max_ratio > 1,
    fraction_grid_above_1 = n_grid_above_one / n_grid_total,
    min_repair_multiplier = 1 - overall_max_ratio
  )

  list(
    overall = overall_draw_results,
    by_mouse = mouse_draw_results
  )
}

repair_check_exp <- check_repair_multiplier_one_model(
  fit = fit_exp_stage4,
  model_name = "Exponential"
)

repair_check_log <- check_repair_multiplier_one_model(
  fit = fit_log_stage4,
  model_name = "Logistic"
)

repair_check_gomp <- check_repair_multiplier_one_model(
  fit = fit_gomp_stage4,
  model_name = "Gompertz"
)

repair_check_overall <- bind_rows(
  repair_check_exp$overall,
  repair_check_log$overall,
  repair_check_gomp$overall
)

repair_check_by_mouse <- bind_rows(
  repair_check_exp$by_mouse,
  repair_check_log$by_mouse,
  repair_check_gomp$by_mouse
)

repair_check_summary <- repair_check_overall %>%
  group_by(
    model
  ) %>%
  summarise(
    median_max_krepi_C2 = median(max_krepi_C2),
    q025_max_krepi_C2 = quantile(max_krepi_C2, 0.025),
    q975_max_krepi_C2 = quantile(max_krepi_C2, 0.975),
    probability_ever_above_1 = mean(ever_above_1),
    median_fraction_grid_above_1 = median(fraction_grid_above_1),
    q975_fraction_grid_above_1 = quantile(fraction_grid_above_1, 0.975),
    .groups = "drop"
  )

repair_check_mouse_summary <- repair_check_by_mouse %>%
  group_by(
    model,
    mouse,
    mouse_id
  ) %>%
  summarise(
    median_max_krepi_C2 = median(max_krepi_C2),
    q025_max_krepi_C2 = quantile(max_krepi_C2, 0.025),
    q975_max_krepi_C2 = quantile(max_krepi_C2, 0.975),
    probability_ever_above_1 = mean(ever_above_1),
    median_fraction_grid_above_1 = median(fraction_grid_above_1),
    .groups = "drop"
  )

cat(
  "\nPosterior diagnostic summary for k_repi * C2(t):\n"
)

print(
  as.data.frame(repair_check_summary),
  row.names = FALSE
)

cat(
  "\nMouse-specific posterior diagnostic summary:\n"
)

print(
  as.data.frame(repair_check_mouse_summary),
  row.names = FALSE
)

write.csv(
  repair_check_summary,
  file.path(
    table_dir,
    "stage4_repair_multiplier_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  repair_check_mouse_summary,
  file.path(
    table_dir,
    "stage4_repair_multiplier_by_combination_mouse.csv"
  ),
  row.names = FALSE
)

write.csv(
  repair_check_overall,
  file.path(
    table_dir,
    "stage4_repair_multiplier_draw_level.csv"
  ),
  row.names = FALSE
)

# Diagnostic figure
p_repair_multiplier <- ggplot(
  repair_check_overall,
  aes(
    x = max_krepi_C2
  )
) +
  geom_histogram(
    bins = 60,
    boundary = 0
  ) +
  geom_vline(
    xintercept = 1,
    linetype = 2,
    linewidth = 0.8
  ) +
  facet_wrap(
    ~ model,
    scales = "free_y"
  ) +
  labs(
    title = "Stage 4 posterior check of the effective DNA-repair multiplier",
    subtitle = "Maximum of k_repi * C2(t) across all combination mice in each posterior draw",
    x = "Maximum k_repi * C2(t)",
    y = "Posterior draws"
  ) +
  theme_bw(base_size = 12)

print(
  p_repair_multiplier
)

ggsave(
  filename = file.path(
    figure_dir,
    "stage4_repair_multiplier_posterior_check.pdf"
  ),
  plot = p_repair_multiplier,
  width = 9,
  height = 5.5
)

ggsave(
  filename = file.path(
    figure_dir,
    "stage4_repair_multiplier_posterior_check.png"
  ),
  plot = p_repair_multiplier,
  width = 9,
  height = 5.5,
  dpi = 300
)

cat(
  paste0(
    "\nInterpretation:\n",
    "  probability_ever_above_1 is the posterior probability that at least one ",
    "combination-mouse trajectory enters k_repi*C2(t) > 1 during observed follow-up.\n",
    "  Values near zero mean the fitted model almost always remains in the positive-repair regime.\n",
    "  Large values mean the posterior frequently uses the regime where D can grow between RT events.\n"
  )
)

############################################################
# 5. POSTERIOR PREDICTIVE AND RESIDUAL DIAGNOSTICS
############################################################


fit_gq_exp  <- fit_exp_stage4
fit_gq_log  <- fit_log_stage4
fit_gq_gomp <- fit_gomp_stage4

required_gq_variables <- c(
  "mu_log_baseline",
  "mu_baseline",
  "y_rep_baseline",
  "log_lik_baseline",
  "mu_log_post",
  "mu_post",
  "y_rep_post",
  "log_lik_post",
  "min_repair_factor_checked",
  "negative_repair_factor_count"
)

check_gq_presence <- function(fit, model_name) {
  for (v in required_gq_variables) {
    ok <- tryCatch({
      x <- fit$draws(variables = v)
      length(x) > 0L
    }, error = function(e) FALSE)
    if (!ok) {
      stop(model_name, " final fit is missing built-in GQ variable: ", v)
    }
  }
  cat(model_name, ": built-in generated quantities verified.\n")
}

check_gq_presence(fit_gq_exp,  "Exponential")
check_gq_presence(fit_gq_log,  "Logistic")
check_gq_presence(fit_gq_gomp, "Gompertz")


obs_stage4 <- obs_stage4 %>%
  arrange(
    mouse_id,
    time_days
  ) %>%
  mutate(
    obs_order = row_number()
  )


baseline_stage4 <- obs_stage4 %>%
  group_by(
    mouse,
    mouse_id,
    group,
    group_name
  ) %>%
  slice_min(
    time_days,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(
    mouse_id
  )


obs_post_stage4 <- obs_stage4 %>%
  left_join(
    baseline_stage4 %>%
      select(
        mouse,
        t0 = time_days
      ),
    by = "mouse"
  ) %>%
  filter(
    time_days > t0
  ) %>%
  select(
    -t0
  ) %>%
  arrange(
    mouse_id,
    time_days
  ) %>%
  mutate(
    obs_index = row_number()
  )


stopifnot(
  nrow(baseline_stage4) == 78,
  nrow(obs_post_stage4) == 530,
  nrow(obs_stage4) == 608
)

cat(
  "\nBaseline observations:",
  nrow(baseline_stage4),
  "\n"
)

cat(
  "Post-baseline observations:",
  nrow(obs_post_stage4),
  "\n"
)

cat(
  "Total observations:",
  nrow(obs_stage4),
  "\n"
)


extract_stage4_gq <- function(
    fit_gq,
    model_name
) {
  
  mu_log_baseline <-
    posterior::as_draws_matrix(
      fit_gq$draws(
        variables = "mu_log_baseline"
      )
    )
  
  mu_baseline <-
    posterior::as_draws_matrix(
      fit_gq$draws(
        variables = "mu_baseline"
      )
    )
  
  y_rep_baseline <-
    posterior::as_draws_matrix(
      fit_gq$draws(
        variables = "y_rep_baseline"
      )
    )
  
  log_lik_baseline <-
    posterior::as_draws_matrix(
      fit_gq$draws(
        variables = "log_lik_baseline"
      )
    )
  
  
  mu_log_post <-
    posterior::as_draws_matrix(
      fit_gq$draws(
        variables = "mu_log_post"
      )
    )
  
  mu_post <-
    posterior::as_draws_matrix(
      fit_gq$draws(
        variables = "mu_post"
      )
    )
  
  y_rep_post <-
    posterior::as_draws_matrix(
      fit_gq$draws(
        variables = "y_rep_post"
      )
    )
  
  log_lik_post <-
    posterior::as_draws_matrix(
      fit_gq$draws(
        variables = "log_lik_post"
      )
    )
  
  
  stopifnot(
    nrow(mu_log_baseline) == 4000,
    ncol(mu_log_baseline) == 78,
    
    nrow(mu_log_post) == 4000,
    ncol(mu_log_post) == 530
  )
  
  
  n_draws <- nrow(mu_log_baseline)
  n_total <- nrow(obs_stage4)
  
  
  mu_log_all <-
    matrix(
      NA_real_,
      nrow = n_draws,
      ncol = n_total
    )
  
  mu_all <- mu_log_all
  y_rep_all <- mu_log_all
  log_lik_all <- mu_log_all
  
  
  baseline_positions <-
    baseline_stage4$obs_order
  
  post_positions <-
    obs_post_stage4$obs_order
  
  
  mu_log_all[, baseline_positions] <-
    mu_log_baseline
  
  mu_log_all[, post_positions] <-
    mu_log_post
  
  
  mu_all[, baseline_positions] <-
    mu_baseline
  
  mu_all[, post_positions] <-
    mu_post
  
  
  y_rep_all[, baseline_positions] <-
    y_rep_baseline
  
  y_rep_all[, post_positions] <-
    y_rep_post
  
  
  log_lik_all[, baseline_positions] <-
    log_lik_baseline
  
  log_lik_all[, post_positions] <-
    log_lik_post
  
  
  stopifnot(
    all(is.finite(mu_log_all)),
    all(is.finite(mu_all)),
    all(mu_all > 0),
    all(is.finite(y_rep_all)),
    all(y_rep_all > 0),
    all(is.finite(log_lik_all))
  )
  
  
  cat(
    model_name,
    ": reconstructed",
    nrow(mu_all),
    "draws x",
    ncol(mu_all),
    "observations.\n"
  )
  
  
  list(
    mu_log = mu_log_all,
    mu = mu_all,
    y_rep = y_rep_all,
    log_lik = log_lik_all
  )
}


gq_exp <- extract_stage4_gq(
  fit_gq_exp,
  "Exponential"
)

gq_log <- extract_stage4_gq(
  fit_gq_log,
  "Logistic"
)

gq_gomp <- extract_stage4_gq(
  fit_gq_gomp,
  "Gompertz"
)
############################################################
#  FITTED AND POSTERIOR PREDICTIVE SUMMARIES
############################################################

make_prediction_summary <- function(
    gq,
    model_name
) {
  
  fitted_mean <-
    colMeans(
      gq$mu
    )
  
  fitted_median <-
    apply(
      gq$mu,
      2,
      median
    )
  
  fitted_q025 <-
    apply(
      gq$mu,
      2,
      quantile,
      probs = 0.025
    )
  
  fitted_q975 <-
    apply(
      gq$mu,
      2,
      quantile,
      probs = 0.975
    )
  
  
  pred_median <-
    apply(
      gq$y_rep,
      2,
      median
    )
  
  pred_q025 <-
    apply(
      gq$y_rep,
      2,
      quantile,
      probs = 0.025
    )
  
  pred_q975 <-
    apply(
      gq$y_rep,
      2,
      quantile,
      probs = 0.975
    )
  
  
  mean_mu_log <-
    colMeans(
      gq$mu_log
    )
  
  
  obs_stage4 %>%
    mutate(
      model = model_name,
      
      fitted_mean = fitted_mean,
      fitted_median = fitted_median,
      fitted_q025 = fitted_q025,
      fitted_q975 = fitted_q975,
      
      pred_median = pred_median,
      pred_q025 = pred_q025,
      pred_q975 = pred_q975,
      
      log_residual =
        log(volume) - mean_mu_log,
      
      abs_log_residual =
        abs(log_residual),
      
      sq_log_residual =
        log_residual^2,
      
      covered_95 =
        volume >= pred_q025 &
        volume <= pred_q975
    )
}
pred_exp <- make_prediction_summary(
  gq_exp,
  "Exponential"
)

pred_log <- make_prediction_summary(
  gq_log,
  "Logistic"
)

pred_gomp <- make_prediction_summary(
  gq_gomp,
  "Gompertz"
)


prediction_stage4 <- bind_rows(
  pred_exp,
  pred_log,
  pred_gomp
)


write.csv(
  prediction_stage4,
  file.path(
    table_dir,
    "stage4_observation_level_predictions.csv"
  ),
  row.names = FALSE
)
############################################################
#  RESIDUAL AND PPC NUMERICAL SUMMARY
############################################################

residual_ppc_summary <- prediction_stage4 %>%
  group_by(
    model
  ) %>%
  summarise(
    n = n(),
    
    mean_log_residual =
      mean(log_residual),
    
    sd_log_residual =
      sd(log_residual),
    
    mae_log =
      mean(abs_log_residual),
    
    rmse_log =
      sqrt(mean(sq_log_residual)),
    
    max_abs_log_residual =
      max(abs_log_residual),
    
    coverage_95 =
      mean(covered_95),
    
    .groups = "drop"
  )


print(
  as.data.frame(
    residual_ppc_summary
  ),
  row.names = FALSE
)


write.csv(
  residual_ppc_summary,
  file.path(
    table_dir,
    "stage4_residual_ppc_summary.csv"
  ),
  row.names = FALSE
)
############################################################
#  RESIDUAL AND PPC SUMMARY BY TREATMENT GROUP
############################################################

residual_ppc_by_group <- prediction_stage4 %>%
  group_by(
    model,
    group,
    group_name
  ) %>%
  summarise(
    n = n(),
    
    mean_log_residual =
      mean(log_residual),
    
    sd_log_residual =
      sd(log_residual),
    
    mae_log =
      mean(abs_log_residual),
    
    rmse_log =
      sqrt(mean(sq_log_residual)),
    
    coverage_95 =
      mean(covered_95),
    
    .groups = "drop"
  )


print(
  as.data.frame(
    residual_ppc_by_group
  ),
  row.names = FALSE
)


write.csv(
  residual_ppc_by_group,
  file.path(
    table_dir,
    "stage4_residual_ppc_by_group.csv"
  ),
  row.names = FALSE
)
############################################################
#  PSIS-LOO MODEL COMPARISON
############################################################

loo_exp_stage4 <- loo::loo(
  gq_exp$log_lik
)

loo_log_stage4 <- loo::loo(
  gq_log$log_lik
)

loo_gomp_stage4 <- loo::loo(
  gq_gomp$log_lik
)


print(
  loo_exp_stage4
)

print(
  loo_log_stage4
)

print(
  loo_gomp_stage4
)


loo_compare_stage4 <- loo::loo_compare(
  list(
    Exponential = loo_exp_stage4,
    Logistic = loo_log_stage4,
    Gompertz = loo_gomp_stage4
  )
)

print(loo_compare_stage4)

loo_compare_table <- as.data.frame(
  loo_compare_stage4
)

loo_compare_table$model <-
  rownames(
    loo_compare_table
  )

rownames(
  loo_compare_table
) <- NULL


loo_compare_table <- loo_compare_table %>%
  select(
    model,
    everything()
  )


write.csv(
  loo_compare_table,
  file.path(
    table_dir,
    "stage4_loo_model_comparison.csv"
  ),
  row.names = FALSE
)
############################################################
#  PARETO-K DIAGNOSTICS
############################################################

extract_pareto <- function(
    loo_object,
    model_name
) {
  
  k <- loo::pareto_k_values(
    loo_object
  )
  
  obs_stage4 %>%
    mutate(
      model = model_name,
      pareto_k = k
    )
}


pareto_exp <- extract_pareto(
  loo_exp_stage4,
  "Exponential"
)

pareto_log <- extract_pareto(
  loo_log_stage4,
  "Logistic"
)

pareto_gomp <- extract_pareto(
  loo_gomp_stage4,
  "Gompertz"
)


pareto_stage4 <- bind_rows(
  pareto_exp,
  pareto_log,
  pareto_gomp
)


pareto_summary_stage4 <- pareto_stage4 %>%
  group_by(
    model
  ) %>%
  summarise(
    n_obs = n(),
    n_k_gt_0_5 = sum(pareto_k > 0.5),
    n_k_gt_0_7 = sum(pareto_k > 0.7),
    n_k_gt_1 = sum(pareto_k > 1),
    max_k = max(pareto_k),
    .groups = "drop"
  )


print(
  as.data.frame(
    pareto_summary_stage4
  ),
  row.names = FALSE
)


write.csv(
  pareto_stage4,
  file.path(
    table_dir,
    "stage4_pareto_k_observation_level.csv"
  ),
  row.names = FALSE
)

write.csv(
  pareto_summary_stage4,
  file.path(
    table_dir,
    "stage4_pareto_k_summary.csv"
  ),
  row.names = FALSE
)
############################################################
#  IDENTIFY PROBLEMATIC PARETO-K OBSERVATIONS
############################################################

problematic_pareto_stage4 <- pareto_stage4 %>%
  filter(
    pareto_k > 0.7
  ) %>%
  arrange(
    model,
    desc(pareto_k)
  ) %>%
  select(
    model,
    obs_order,
    mouse,
    mouse_id,
    group,
    group_name,
    time_days,
    volume,
    pareto_k
  )


print(
  as.data.frame(
    problematic_pareto_stage4
  ),
  row.names = FALSE
)


write.csv(
  problematic_pareto_stage4,
  file.path(
    table_dir,
    "stage4_problematic_pareto_observations.csv"
  ),
  row.names = FALSE
)
############################################################
#  INSPECT HIGH-PARETO-K OBSERVATIONS WITH RESIDUALS
############################################################

pareto_with_predictions <- pareto_stage4 %>%
  select(
    model,
    obs_order,
    pareto_k
  ) %>%
  left_join(
    prediction_stage4 %>%
      select(
        model,
        obs_order,
        mouse,
        mouse_id,
        group,
        group_name,
        time_days,
        volume,
        fitted_median,
        pred_q025,
        pred_q975,
        log_residual
      ),
    by = c(
      "model",
      "obs_order"
    )
  )

high_pareto_with_predictions <- pareto_with_predictions %>%
  filter(
    pareto_k > 0.7
  ) %>%
  arrange(
    model,
    desc(pareto_k)
  )

print(
  as.data.frame(
    high_pareto_with_predictions
  ),
  row.names = FALSE
)
############################################################
#  INSPECT MOUSE 77 TRAJECTORY
############################################################

mouse77_check <- prediction_stage4 %>%
  filter(
    mouse == 77
  ) %>%
  select(
    model,
    mouse,
    time_days,
    volume,
    fitted_median,
    pred_q025,
    pred_q975,
    log_residual
  ) %>%
  arrange(
    model,
    time_days
  )

print(
  as.data.frame(
    mouse77_check
  ),
  row.names = FALSE
)
############################################################
#  VERIFY MOUSE 77 DIRECTLY AGAINST RAW DATA
############################################################

mouse77_raw <- data_raw %>%
  filter(
    mouse == 77
  ) %>%
  mutate(
    time_days = time / 24
  ) %>%
  arrange(
    time,
    modality
  ) %>%
  select(
    mouse,
    time,
    time_days,
    observation,
    volume,
    modality,
    dose,
    study
  )

print(
  as.data.frame(mouse77_raw),
  row.names = FALSE
)
mouse77_raw_obs <- data_raw %>%
  filter(
    mouse == 77,
    observation == 1
  ) %>%
  mutate(
    time_days = time / 24
  ) %>%
  arrange(
    time_days
  ) %>%
  select(
    time_days,
    volume
  )

mouse77_analysis_obs <- obs_stage4 %>%
  filter(
    mouse == 77
  ) %>%
  arrange(
    time_days
  ) %>%
  select(
    time_days,
    volume
  )

print(mouse77_raw_obs)
print(mouse77_analysis_obs)

stopifnot(
  identical(
    mouse77_raw_obs$time_days,
    mouse77_analysis_obs$time_days
  ),
  identical(
    mouse77_raw_obs$volume,
    mouse77_analysis_obs$volume
  )
)

cat(
  "\nMouse 77 observations exactly match the original CSV.\n"
)
############################################################
#  FINAL OBSERVATION-ORDER CONSISTENCY CHECK
############################################################

baseline_positions <- sort(
  baseline_stage4$obs_order
)

post_positions <- sort(
  obs_post_stage4$obs_order
)

all_positions <- sort(
  c(
    baseline_positions,
    post_positions
  )
)

stopifnot(
  length(baseline_positions) == 78,
  length(post_positions) == 530,
  length(all_positions) == 608,
  length(unique(all_positions)) == 608,
  identical(
    all_positions,
    1:608
  ),
  length(
    intersect(
      baseline_positions,
      post_positions
    )
  ) == 0
)

cat(
  "\nPASS: baseline and post-baseline positions form an exact, non-overlapping partition of observations 1:608.\n"
)
stopifnot(
  identical(dim(gq_exp$mu), c(4000L, 608L)),
  identical(dim(gq_exp$y_rep), c(4000L, 608L)),
  identical(dim(gq_exp$log_lik), c(4000L, 608L)),
  
  identical(dim(gq_log$mu), c(4000L, 608L)),
  identical(dim(gq_log$y_rep), c(4000L, 608L)),
  identical(dim(gq_log$log_lik), c(4000L, 608L)),
  
  identical(dim(gq_gomp$mu), c(4000L, 608L)),
  identical(dim(gq_gomp$y_rep), c(4000L, 608L)),
  identical(dim(gq_gomp$log_lik), c(4000L, 608L))
)

cat(
  "PASS: all Stage 4 generated-quantity matrices are 4000 x 608.\n"
)
############################################################
#  STAGE 4 POSTERIOR PREDICTIVE TRAJECTORIES
############################################################

ppc_stage4_plot_data <- prediction_stage4 %>%
  select(
    model,
    mouse,
    mouse_id,
    group,
    group_name,
    time_days,
    volume,
    pred_median,
    pred_q025,
    pred_q975
  )

p_stage4_ppc <- ggplot(
  ppc_stage4_plot_data,
  aes(
    x = time_days
  )
) +
  geom_ribbon(
    aes(
      ymin = pred_q025,
      ymax = pred_q975,
      group = interaction(model, mouse)
    ),
    alpha = 0.20
  ) +
  geom_line(
    aes(
      y = pred_median,
      group = interaction(model, mouse)
    ),
    linewidth = 0.45
  ) +
  geom_point(
    aes(
      y = volume
    ),
    size = 0.9
  ) +
  facet_grid(
    model ~ group_name,
    scales = "free_y"
  ) +
  labs(
    title = "Full joint treatment model: posterior predictive tumour trajectories",
    subtitle = "Points: observed volumes; lines: posterior predictive medians; ribbons: 95% posterior predictive intervals",
    x = "Time (days)",
    y = expression("Tumour volume (mm"^3*")")
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 9),
    plot.title = element_text(face = "bold")
  )

print(
  p_stage4_ppc
)

ggsave(
  filename = file.path(
    figure_dir,
    "stage4_posterior_predictive_trajectories_by_group.pdf"
  ),
  plot = p_stage4_ppc,
  width = 13,
  height = 8
)

ggsave(
  filename = file.path(
    figure_dir,
    "stage4_posterior_predictive_trajectories_by_group.png"
  ),
  plot = p_stage4_ppc,
  width = 13,
  height = 8,
  dpi = 300
)
############################################################
#  DISSERTATION-QUALITY PPC FIGURES BY MODEL
############################################################
make_stage4_ppc_plot <- function(
    data,
    model_name
) {
  
  plot_data <- data %>%
    filter(
      model == model_name
    ) %>%
    mutate(
      group_name = factor(
        group_name,
        levels = c(
          "control",
          "drug_only",
          "rt_only",
          "combination"
        ),
        labels = c(
          "Control",
          "Drug only",
          "RT only",
          "Combination"
        )
      )
    )
  
  p <- ggplot(
    plot_data,
    aes(
      x = time_days
    )
  ) +
    geom_ribbon(
      aes(
        ymin = pred_q025,
        ymax = pred_q975,
        group = mouse
      ),
      alpha = 0.15
    ) +
    geom_line(
      aes(
        y = pred_median,
        group = mouse
      ),
      linewidth = 0.45
    ) +
    geom_point(
      aes(
        y = volume
      ),
      size = 1.1
    ) +
    facet_wrap(
      ~ group_name,
      nrow = 2,
      scales = "free_y"
    ) +
    labs(
      title = paste0(
        model_name,
        " model: Full joint treatment model"
      ),
      subtitle = "Posterior predictive trajectories",
      x = "Time (days)",
      y = expression("Tumour volume (mm"^3*")")
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      legend.position = "none",
      strip.text = element_text(
        face = "bold"
      ),
      plot.title = element_text(
        face = "bold"
      )
    )
  
  return(p)
}

p_ppc_exp <- make_stage4_ppc_plot(
  prediction_stage4,
  "Exponential"
)

p_ppc_log <- make_stage4_ppc_plot(
  prediction_stage4,
  "Logistic"
)

p_ppc_gomp <- make_stage4_ppc_plot(
  prediction_stage4,
  "Gompertz"
)


print(p_ppc_exp)
print(p_ppc_log)
print(p_ppc_gomp)


ggsave(
  file.path(
    figure_dir,
    "stage4_ppc_exponential.pdf"
  ),
  p_ppc_exp,
  width = 10,
  height = 7
)

ggsave(
  file.path(
    figure_dir,
    "stage4_ppc_logistic.pdf"
  ),
  p_ppc_log,
  width = 10,
  height = 7
)

ggsave(
  file.path(
    figure_dir,
    "stage4_ppc_gompertz.pdf"
  ),
  p_ppc_gomp,
  width = 10,
  height = 7
)
############################################################
#  VERIFY MODEL PREDICTIONS ARE DISTINCT
############################################################

prediction_difference_check <- data.frame(
  
  comparison = c(
    "Exponential vs Logistic",
    "Exponential vs Gompertz",
    "Logistic vs Gompertz"
  ),
  
  mean_abs_difference = c(
    
    mean(
      abs(
        pred_exp$fitted_median -
          pred_log$fitted_median
      )
    ),
    
    mean(
      abs(
        pred_exp$fitted_median -
          pred_gomp$fitted_median
      )
    ),
    
    mean(
      abs(
        pred_log$fitted_median -
          pred_gomp$fitted_median
      )
    )
  ),
  
  max_abs_difference = c(
    
    max(
      abs(
        pred_exp$fitted_median -
          pred_log$fitted_median
      )
    ),
    
    max(
      abs(
        pred_exp$fitted_median -
          pred_gomp$fitted_median
      )
    ),
    
    max(
      abs(
        pred_log$fitted_median -
          pred_gomp$fitted_median
      )
    )
  ),
  
  exactly_identical = c(
    
    identical(
      pred_exp$fitted_median,
      pred_log$fitted_median
    ),
    
    identical(
      pred_exp$fitted_median,
      pred_gomp$fitted_median
    ),
    
    identical(
      pred_log$fitted_median,
      pred_gomp$fitted_median
    )
  )
)

print(
  prediction_difference_check,
  row.names = FALSE
)
############################################################
#  RESIDUAL DIAGNOSTIC FIGURES
############################################################

residual_plot_data <- prediction_stage4 %>%
  mutate(
    group_name = factor(
      group_name,
      levels = c(
        "control",
        "drug_only",
        "rt_only",
        "combination"
      ),
      labels = c(
        "Control",
        "Drug only",
        "RT only",
        "Combination"
      )
    )
  )


############################################################
#  RESIDUALS VS FITTED VALUES
############################################################

p_resid_fitted <- ggplot(
  residual_plot_data,
  aes(
    x = fitted_median,
    y = log_residual
  )
) +
  geom_point(
    alpha = 0.45,
    size = 1.1
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.7
  ) +
  facet_grid(
    model ~ group_name,
    scales = "free_x"
  ) +
  labs(
    title = "Full joint treatment model: residuals versus fitted tumour volume",
    x = expression("Posterior median fitted volume (mm"^3*")"),
    y = "Log-scale residual"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

print(p_resid_fitted)


ggsave(
  file.path(
    figure_dir,
    "stage4_residuals_vs_fitted.pdf"
  ),
  p_resid_fitted,
  width = 12,
  height = 8
)


############################################################
#  RESIDUALS VS TIME
############################################################

p_resid_time <- ggplot(
  residual_plot_data,
  aes(
    x = time_days,
    y = log_residual
  )
) +
  geom_point(
    alpha = 0.45,
    size = 1.1
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.7
  ) +
  facet_grid(
    model ~ group_name
  ) +
  labs(
    title = "Full joint treatment model: residuals versus time",
    x = "Time (days)",
    y = "Log-scale residual"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

print(p_resid_time)


ggsave(
  file.path(
    figure_dir,
    "stage4_residuals_vs_time.pdf"
  ),
  p_resid_time,
  width = 12,
  height = 8
)
############################################################
#  PARETO-K DIAGNOSTIC FIGURE
############################################################

pareto_plot_data <- pareto_stage4 %>%
  mutate(
    group_name = factor(
      group_name,
      levels = c(
        "control",
        "drug_only",
        "rt_only",
        "combination"
      ),
      labels = c(
        "Control",
        "Drug only",
        "RT only",
        "Combination"
      )
    )
  )


p_pareto <- ggplot(
  pareto_plot_data,
  aes(
    x = obs_order,
    y = pareto_k
  )
) +
  geom_point(
    alpha = 0.65,
    size = 1.2
  ) +
  geom_hline(
    yintercept = 0.7,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_hline(
    yintercept = 1.0,
    linetype = "dotted",
    linewidth = 0.6
  ) +
  facet_grid(
    model ~ group_name,
    scales = "free_x"
  ) +
  labs(
    title = "Full joint treatment model: PSIS-LOO Pareto-k diagnostics",
    subtitle = "Dashed line: k = 0.7; dotted line: k = 1.0",
    x = "Observation index",
    y = "Pareto-k"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

print(p_pareto)


ggsave(
  file.path(
    figure_dir,
    "stage4_pareto_k_diagnostics.pdf"
  ),
  p_pareto,
  width = 12,
  height = 8
)
############################################################
#  PRIOR VS POSTERIOR
#     Stage 4 actual fitted priors
############################################################

n_prior_stage4 <- 50000

stage4_prior_draws <- function(
    model_name
) {
  
  r0 <- switch(
    model_name,
    Exponential = 0.3,
    Logistic = 0.4,
    Gompertz = 0.2
  )
  
  
  prior <- data.frame(
    
    logV0_pop =
      rnorm(
        n_prior_stage4,
        log(30),
        0.5
      ),
    
    logr_pop =
      rnorm(
        n_prior_stage4,
        log(r0),
        0.5
      ),
    
    logdelta_pop =
      rnorm(
        n_prior_stage4,
        log(1.2),
        0.5
      ),
    
    logalpha_pop =
      rnorm(
        n_prior_stage4,
        log(0.9),
        0.5
      ),
    
    logbeta_a_pop =
      rnorm(
        n_prior_stage4,
        log(26.4),
        0.5
      ),
    
    logbeta_q_pop =
      rnorm(
        n_prior_stage4,
        log(4.8),
        0.5
      ),
    
    logk_rep_pop =
      rnorm(
        n_prior_stage4,
        log(86.4),
        0.5
      ),
    
    logk_abs_pop =
      rnorm(
        n_prior_stage4,
        log(96),
        0.35
      ),
    
    logk_cl_pop =
      rnorm(
        n_prior_stage4,
        log(19.2),
        0.35
      ),
    
    logk_repi_pop =
      rnorm(
        n_prior_stage4,
        log(0.04),
        0.35
      ),
    
    sigma_logV0 =
      abs(
        rnorm(
          n_prior_stage4,
          0,
          0.5
        )
      ),
    
    sigma_logr =
      abs(
        rnorm(
          n_prior_stage4,
          0,
          0.5
        )
      ),
    
    sigma_logdelta =
      abs(
        rnorm(
          n_prior_stage4,
          0,
          0.5
        )
      ),
    
    sigma_logalpha =
      abs(
        rnorm(
          n_prior_stage4,
          0,
          0.5
        )
      ),
    
    sigma_logbeta_a =
      abs(
        rnorm(
          n_prior_stage4,
          0,
          0.5
        )
      ),
    
    sigma_logbeta_q =
      abs(
        rnorm(
          n_prior_stage4,
          0,
          0.5
        )
      ),
    
    sigma_logk_rep =
      abs(
        rnorm(
          n_prior_stage4,
          0,
          0.5
        )
      ),
    
    sigma =
      abs(
        rnorm(
          n_prior_stage4,
          0,
          0.5
        )
      )
  )
  
  
  if (
    model_name %in%
    c(
      "Logistic",
      "Gompertz"
    )
  ) {
    
    prior$logK_pop <-
      rnorm(
        n_prior_stage4,
        log(500),
        0.7
      )
    
    prior$sigma_logK <-
      abs(
        rnorm(
          n_prior_stage4,
          0,
          0.5
        )
      )
  }
  
  
  prior
}
############################################################
#  PRIOR VS POSTERIOR:
#  POPULATION + MECHANISTIC PARAMETERS
############################################################

stage4_pop_vars <- function(
    model_name
) {
  
  vars <- c(
    "logV0_pop",
    "logr_pop",
    "logdelta_pop",
    "logalpha_pop",
    "logbeta_a_pop",
    "logbeta_q_pop",
    "logk_rep_pop",
    "logk_abs_pop",
    "logk_cl_pop",
    "logk_repi_pop"
  )
  
  if (
    model_name %in%
    c(
      "Logistic",
      "Gompertz"
    )
  ) {
    
    vars <- c(
      vars,
      "logK_pop"
    )
  }
  
  vars
}


make_stage4_prior_posterior_population <- function(
    fit,
    model_name
) {
  
  vars <-
    stage4_pop_vars(
      model_name
    )
  
  
  prior_long <-
    stage4_prior_draws(
      model_name
    ) %>%
    select(
      all_of(vars)
    ) %>%
    pivot_longer(
      everything(),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    mutate(
      distribution = "Prior"
    )
  
  
  posterior_long <-
    as.data.frame(
      posterior::as_draws_matrix(
        fit$draws(
          variables = vars
        )
      )
    ) %>%
    select(
      all_of(vars)
    ) %>%
    pivot_longer(
      everything(),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    mutate(
      distribution = "Posterior"
    )
  
  
  plot_data <-
    bind_rows(
      prior_long,
      posterior_long
    ) %>%
    mutate(
      distribution = factor(
        distribution,
        levels = c(
          "Prior",
          "Posterior"
        )
      )
    )
  
  
  ggplot(
    plot_data,
    aes(
      x = value,
      linetype = distribution,
      group = distribution
    )
  ) +
    geom_density(
      linewidth = 0.85,
      adjust = 1.1
    ) +
    facet_wrap(
      ~ parameter,
      scales = "free",
      ncol = 3
    ) +
    labs(
      title = paste0(
        model_name,
        " model: Full joint treatment model — Prior vs posterior"
      ),
      subtitle =
        "Population growth, treatment-response and pharmacokinetic parameters",
      x =
        "Parameter value on fitted log scale",
      y =
        "Density",
      linetype =
        NULL
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      legend.position =
        "bottom",
      strip.text =
        element_text(
          face = "bold"
        ),
      plot.title =
        element_text(
          face = "bold"
        )
    )
}


p_priorpost_pop_exp <-
  make_stage4_prior_posterior_population(
    fit_exp_stage4,
    "Exponential"
  )

p_priorpost_pop_log <-
  make_stage4_prior_posterior_population(
    fit_log_stage4,
    "Logistic"
  )

p_priorpost_pop_gomp <-
  make_stage4_prior_posterior_population(
    fit_gomp_stage4,
    "Gompertz"
  )


print(
  p_priorpost_pop_exp
)

print(
  p_priorpost_pop_log
)

print(
  p_priorpost_pop_gomp
)


ggsave(
  file.path(
    figure_dir,
    "stage4_prior_vs_posterior_population_exponential.pdf"
  ),
  p_priorpost_pop_exp,
  width = 11,
  height = 8
)

ggsave(
  file.path(
    figure_dir,
    "stage4_prior_vs_posterior_population_logistic.pdf"
  ),
  p_priorpost_pop_log,
  width = 11,
  height = 8
)

ggsave(
  file.path(
    figure_dir,
    "stage4_prior_vs_posterior_population_gompertz.pdf"
  ),
  p_priorpost_pop_gomp,
  width = 11,
  height = 8
)


############################################################
#  PRIOR VS POSTERIOR:
#  BETWEEN-MOUSE VARIABILITY + RESIDUAL SD
############################################################

stage4_sd_vars <- function(
    model_name
) {
  
  vars <- c(
    "sigma_logV0",
    "sigma_logr",
    "sigma_logdelta",
    "sigma_logalpha",
    "sigma_logbeta_a",
    "sigma_logbeta_q",
    "sigma_logk_rep",
    "sigma"
  )
  
  if (
    model_name %in%
    c(
      "Logistic",
      "Gompertz"
    )
  ) {
    
    vars <- c(
      vars,
      "sigma_logK"
    )
  }
  
  vars
}


make_stage4_prior_posterior_sd <- function(
    fit,
    model_name
) {
  
  vars <-
    stage4_sd_vars(
      model_name
    )
  
  
  prior_long <-
    stage4_prior_draws(
      model_name
    ) %>%
    select(
      all_of(vars)
    ) %>%
    pivot_longer(
      everything(),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    mutate(
      distribution = "Prior"
    )
  
  
  posterior_long <-
    as.data.frame(
      posterior::as_draws_matrix(
        fit$draws(
          variables = vars
        )
      )
    ) %>%
    select(
      all_of(vars)
    ) %>%
    pivot_longer(
      everything(),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    mutate(
      distribution = "Posterior"
    )
  
  
  plot_data <-
    bind_rows(
      prior_long,
      posterior_long
    ) %>%
    mutate(
      distribution = factor(
        distribution,
        levels = c(
          "Prior",
          "Posterior"
        )
      )
    )
  
  
  ggplot(
    plot_data,
    aes(
      x = value,
      linetype = distribution,
      group = distribution
    )
  ) +
    geom_density(
      linewidth = 0.85,
      adjust = 1.1
    ) +
    facet_wrap(
      ~ parameter,
      scales = "free",
      ncol = 3
    ) +
    labs(
      title = paste0(
        model_name,
        " model: Full joint treatment model — Prior vs posterior"
      ),
      subtitle =
        "Between-mouse variability and residual observation variability",
      x =
        "Parameter value",
      y =
        "Density",
      linetype =
        NULL
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      legend.position =
        "bottom",
      strip.text =
        element_text(
          face = "bold"
        ),
      plot.title =
        element_text(
          face = "bold"
        )
    )
}


p_priorpost_sd_exp <-
  make_stage4_prior_posterior_sd(
    fit_exp_stage4,
    "Exponential"
  )

p_priorpost_sd_log <-
  make_stage4_prior_posterior_sd(
    fit_log_stage4,
    "Logistic"
  )

p_priorpost_sd_gomp <-
  make_stage4_prior_posterior_sd(
    fit_gomp_stage4,
    "Gompertz"
  )


print(
  p_priorpost_sd_exp
)

print(
  p_priorpost_sd_log
)

print(
  p_priorpost_sd_gomp
)


ggsave(
  file.path(
    figure_dir,
    "stage4_prior_vs_posterior_sd_exponential.pdf"
  ),
  p_priorpost_sd_exp,
  width = 11,
  height = 7
)

ggsave(
  file.path(
    figure_dir,
    "stage4_prior_vs_posterior_sd_logistic.pdf"
  ),
  p_priorpost_sd_log,
  width = 11,
  height = 7
)

ggsave(
  file.path(
    figure_dir,
    "stage4_prior_vs_posterior_sd_gompertz.pdf"
  ),
  p_priorpost_sd_gomp,
  width = 11,
  height = 7
)
############################################################
#  DISSERTATION PRIOR VS POSTERIOR:
#  MEDIAN + 95% INTERVALS ON NATURAL SCALE
############################################################

make_stage4_priorpost_interval_data <- function(
    fit,
    model_name
) {
  
  vars <-
    stage4_pop_vars(
      model_name
    )

  
  prior_draws <-
    stage4_prior_draws(
      model_name
    ) %>%
    select(
      all_of(vars)
    )
  
  posterior_draws <-
    as.data.frame(
      posterior::as_draws_matrix(
        fit$draws(
          variables = vars
        )
      )
    ) %>%
    select(
      all_of(vars)
    )

  
  prior_long <-
    prior_draws %>%
    pivot_longer(
      everything(),
      names_to = "stan_parameter",
      values_to = "log_value"
    ) %>%
    mutate(
      value = exp(log_value),
      distribution = "Prior"
    )
  
  
  posterior_long <-
    posterior_draws %>%
    pivot_longer(
      everything(),
      names_to = "stan_parameter",
      values_to = "log_value"
    ) %>%
    mutate(
      value = exp(log_value),
      distribution = "Posterior"
    )
  
  
  combined <-
    bind_rows(
      prior_long,
      posterior_long
    ) %>%
    mutate(
      
      parameter = case_when(
        
        stan_parameter == "logV0_pop" ~
          "V0 (mm³)",
        
        stan_parameter == "logr_pop" ~
          "r (day^-1)",
        
        stan_parameter == "logK_pop" ~
         "K (mm³)",
        
        stan_parameter == "logalpha_pop" ~
          "alpha",
        
        stan_parameter == "logbeta_a_pop" ~
          "beta_a (day^-1)",
        
        stan_parameter == "logbeta_q_pop" ~
          "beta_q (day^-1)",
        
        stan_parameter == "logdelta_pop" ~
          "delta (day^-1)",
        
        stan_parameter == "logk_rep_pop" ~
          "k_rep (day^-1)",
        
        stan_parameter == "logk_abs_pop" ~
          "k_abs (day^-1)",
        
        stan_parameter == "logk_cl_pop" ~
          "k_cl (day^-1)",
        
        stan_parameter == "logk_repi_pop" ~
          "k_repi",
        
        TRUE ~
          stan_parameter
      ),
      
      
      group = case_when(
        
        stan_parameter %in%
          c(
            "logV0_pop",
            "logr_pop",
            "logK_pop"
          ) ~
          "Growth",
        
        stan_parameter %in%
          c(
            "logalpha_pop",
            "logbeta_a_pop",
            "logbeta_q_pop",
            "logdelta_pop",
            "logk_rep_pop"
          ) ~
          "Treatment response",
        
        stan_parameter %in%
          c(
            "logk_abs_pop",
            "logk_cl_pop",
            "logk_repi_pop"
          ) ~
          "PK / PARP inhibition"
      )
    )
  
  summary_data <-
    combined %>%
    group_by(
      group,
      parameter,
      distribution
    ) %>%
    summarise(
      
      median =
        median(value),
      
      lower =
        quantile(
          value,
          0.025
        ),
      
      upper =
        quantile(
          value,
          0.975
        ),
      
      .groups = "drop"
    ) %>%
    mutate(
      
      distribution = factor(
        distribution,
        levels = c(
          "Prior",
          "Posterior"
        )
      )
    )
  
  
  summary_data
}

priorpost_interval_exp <-
  make_stage4_priorpost_interval_data(
    fit_exp_stage4,
    "Exponential"
  )

priorpost_interval_log <-
  make_stage4_priorpost_interval_data(
    fit_log_stage4,
    "Logistic"
  )

priorpost_interval_gomp <-
  make_stage4_priorpost_interval_data(
    fit_gomp_stage4,
    "Gompertz"
  )
############################################################
# DISSERTATION-QUALITY PLOT FUNCTION
############################################################

plot_stage4_priorpost_intervals <- function(
    summary_data,
    model_name,
    parameter_group
) {
  
  plot_data <-
    summary_data %>%
    filter(
      group == parameter_group
    ) %>%
    mutate(
      distribution = factor(
        distribution,
        levels = c(
          "Prior",
          "Posterior"
        )
      )
    )
  
  
  ggplot(
    plot_data,
    aes(
      x = median,
      y = distribution,
      xmin = lower,
      xmax = upper,
      shape = distribution
    )
  ) +
    
    geom_errorbarh(
      height = 0,
      linewidth = 0.8
    ) +
    
    geom_point(
      size = 3
    ) +
    
    facet_wrap(
      ~ parameter,
      scales = "free_x",
      ncol = 2
    ) +
    
    scale_shape_manual(
      values = c(
        "Prior" = 16,
        "Posterior" = 17
      )
    ) +
    
    labs(
      title = paste0(
        model_name,
        " model: Full joint treatment model"
      ),
      subtitle = paste0(
        parameter_group,
        ": prior and posterior summaries"
      ),
      x = "Median and 95% interval",
      y = NULL,
      shape = NULL
    ) +
    
    theme_bw(
      base_size = 12
    ) +
    
    theme(
      legend.position = "none",
      
      strip.text = element_text(
        face = "bold",
        size = 11
      ),
      
      plot.title = element_text(
        face = "bold",
        size = 15,
        margin = margin(
          b = 4
        )
      ),
      
      plot.subtitle = element_text(
        size = 12,
        margin = margin(
          b = 10
        )
      ),
      
      plot.title.position = "plot",
      
      panel.spacing = grid::unit(
        1.1,
        "lines"
      ),
      
      plot.margin = margin(
        t = 12,
        r = 15,
        b = 12,
        l = 12
      )
    )
}
p_priorpost_exp_growth <-
  plot_stage4_priorpost_intervals(
    priorpost_interval_exp,
    "Exponential",
    "Growth"
  )

p_priorpost_exp_treatment <-
  plot_stage4_priorpost_intervals(
    priorpost_interval_exp,
    "Exponential",
    "Treatment response"
  )

p_priorpost_exp_pk <-
  plot_stage4_priorpost_intervals(
    priorpost_interval_exp,
    "Exponential",
    "PK / PARP inhibition"
  )


print(p_priorpost_exp_growth)
print(p_priorpost_exp_treatment)
print(p_priorpost_exp_pk)
p_priorpost_log_growth <-
  plot_stage4_priorpost_intervals(
    priorpost_interval_log,
    "Logistic",
    "Growth"
  )

p_priorpost_log_treatment <-
  plot_stage4_priorpost_intervals(
    priorpost_interval_log,
    "Logistic",
    "Treatment response"
  )

p_priorpost_log_pk <-
  plot_stage4_priorpost_intervals(
    priorpost_interval_log,
    "Logistic",
    "PK / PARP inhibition"
  )


print(p_priorpost_log_growth)
print(p_priorpost_log_treatment)
print(p_priorpost_log_pk)
p_priorpost_gomp_growth <-
  plot_stage4_priorpost_intervals(
    priorpost_interval_gomp,
    "Gompertz",
    "Growth"
  )

p_priorpost_gomp_treatment <-
  plot_stage4_priorpost_intervals(
    priorpost_interval_gomp,
    "Gompertz",
    "Treatment response"
  )

p_priorpost_gomp_pk <-
  plot_stage4_priorpost_intervals(
    priorpost_interval_gomp,
    "Gompertz",
    "PK / PARP inhibition"
  )


print(p_priorpost_gomp_growth)
print(p_priorpost_gomp_treatment)
print(p_priorpost_gomp_pk)


ggsave(
  file.path(
    figure_dir,
    "stage4_priorpost_growth_exponential.pdf"
  ),
  p_priorpost_exp_growth,
  width = 8,
  height = 5.5
)

ggsave(
  file.path(
    figure_dir,
    "stage4_priorpost_treatment_exponential.pdf"
  ),
  p_priorpost_exp_treatment,
  width = 8,
  height = 8
)

ggsave(
  file.path(
    figure_dir,
    "stage4_priorpost_pk_exponential.pdf"
  ),
  p_priorpost_exp_pk,
  width = 8,
  height = 6
)


ggsave(
  file.path(
    figure_dir,
    "stage4_priorpost_growth_logistic.pdf"
  ),
  p_priorpost_log_growth,
  width = 8,
  height = 6
)

ggsave(
  file.path(
    figure_dir,
    "stage4_priorpost_treatment_logistic.pdf"
  ),
  p_priorpost_log_treatment,
  width = 8,
  height = 8
)

ggsave(
  file.path(
    figure_dir,
    "stage4_priorpost_pk_logistic.pdf"
  ),
  p_priorpost_log_pk,
  width = 8,
  height = 6
)


ggsave(
  file.path(
    figure_dir,
    "stage4_priorpost_growth_gompertz.pdf"
  ),
  p_priorpost_gomp_growth,
  width = 8,
  height = 6
)

ggsave(
  file.path(
    figure_dir,
    "stage4_priorpost_treatment_gompertz.pdf"
  ),
  p_priorpost_gomp_treatment,
  width = 8,
  height = 8
)

ggsave(
  file.path(
    figure_dir,
    "stage4_priorpost_pk_gompertz.pdf"
  ),
  p_priorpost_gomp_pk,
  width = 8,
  height = 6
)

############################################################
# RESIDUAL HISTOGRAMS
############################################################

p_resid_hist <- ggplot(
  residual_plot_data,
  aes(
    x = log_residual
  )
) +
  geom_histogram(
    bins = 30,
    fill = "grey80",
    colour = "black"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  facet_grid(
    model ~ group_name,
    scales = "free_y"
  ) +
  labs(
    title = "Full joint treatment model: residual distributions",
    x = "Log-scale residual",
    y = "Count"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    strip.text = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

print(
  p_resid_hist
)

ggsave(
  file.path(
    figure_dir,
    "stage4_residual_histograms.pdf"
  ),
  p_resid_hist,
  width = 12,
  height = 8
)

############################################################
# RESIDUAL NORMAL Q-Q PLOTS
############################################################

p_resid_qq <- ggplot(
  residual_plot_data,
  aes(
    sample = log_residual
  )
) +
  stat_qq(
    alpha = 0.55,
    size = 1
  ) +
  stat_qq_line(
    linewidth = 0.6
  ) +
  facet_grid(
    model ~ group_name
  ) +
  labs(
    title = "Full joint treatment model: residual Q-Q plots",
    x = "Theoretical quantiles",
    y = "Observed residual quantiles"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    strip.text = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

print(
  p_resid_qq
)

ggsave(
  file.path(
    figure_dir,
    "stage4_residual_qq_plots.pdf"
  ),
  p_resid_qq,
  width = 12,
  height = 8
)
############################################################
# POSTERIOR PREDICTIVE DENSITY OVERLAYS
############################################################

set.seed(123)

# Use a manageable subset of genuine posterior predictive draws
n_ppc_draws <- 50

ppc_draw_ids <- sample(
  seq_len(nrow(gq_exp$y_rep)),
  n_ppc_draws,
  replace = FALSE
)


make_ppc_density_plot <- function(
    gq,
    model_name
) {
  
  observed_log_volume <-
    log(obs_stage4$volume)
  
  yrep_log <-
    log(
      gq$y_rep[
        ppc_draw_ids,
        ,
        drop = FALSE
      ]
    )
  
  
  p <- bayesplot::ppc_dens_overlay(
    y = observed_log_volume,
    yrep = yrep_log
  ) +
    labs(
      title = paste0(
        model_name,
        " model: Full joint treatment model"
      ),
      subtitle =
        "Posterior predictive density check",
      x =
        "Log tumour volume",
      y =
        "Density"
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      plot.title = element_text(
        face = "bold"
      )
    )
  
  return(p)
}


p_ppc_density_exp <-
  make_ppc_density_plot(
    gq_exp,
    "Exponential"
  )

p_ppc_density_log <-
  make_ppc_density_plot(
    gq_log,
    "Logistic"
  )

p_ppc_density_gomp <-
  make_ppc_density_plot(
    gq_gomp,
    "Gompertz"
  )


print(p_ppc_density_exp)
print(p_ppc_density_log)
print(p_ppc_density_gomp)


ggsave(
  file.path(
    figure_dir,
    "stage4_ppc_density_exponential.pdf"
  ),
  p_ppc_density_exp,
  width = 8,
  height = 5.5
)

ggsave(
  file.path(
    figure_dir,
    "stage4_ppc_density_logistic.pdf"
  ),
  p_ppc_density_log,
  width = 8,
  height = 5.5
)

ggsave(
  file.path(
    figure_dir,
    "stage4_ppc_density_gompertz.pdf"
  ),
  p_ppc_density_gomp,
  width = 8,
  height = 5.5
)
############################################################
#  TIME-BINNED POSTERIOR PREDICTIVE CHECKS
############################################################

time_breaks <- seq(
  floor(min(obs_stage4$time_days)),
  ceiling(max(obs_stage4$time_days)) + 5,
  by = 5
)


time_bin_data <- obs_stage4 %>%
  mutate(
    time_bin = cut(
      time_days,
      breaks = time_breaks,
      include.lowest = TRUE,
      right = FALSE
    )
  )

time_bin_info <- time_bin_data %>%
  mutate(
    obs_index = row_number()
  ) %>%
  group_by(
    time_bin
  ) %>%
  summarise(
    time_mid =
      mean(time_days),
    n_obs =
      n(),
    observed_median =
      median(log(volume)),
    .groups = "drop"
  )


make_time_binned_ppc <- function(
    gq,
    model_name
) {
  
  bin_indices <-
    split(
      seq_len(nrow(time_bin_data)),
      time_bin_data$time_bin
    )
  
  
  bin_indices <-
    bin_indices[
      lengths(bin_indices) > 0
    ]
  
  
  predictive_bin_medians <-
    sapply(
      bin_indices,
      function(idx) {
        
        apply(
          log(
            gq$y_rep[
              ,
              idx,
              drop = FALSE
            ]
          ),
          1,
          median
        )
      }
    )
  
  
  if (
    is.vector(
      predictive_bin_medians
    )
  ) {
    
    predictive_bin_medians <-
      matrix(
        predictive_bin_medians,
        ncol = 1
      )
  }
  
  
  predictive_summary <-
    data.frame(
      
      time_bin =
        names(bin_indices),
      
      predictive_median =
        apply(
          predictive_bin_medians,
          2,
          median
        ),
      
      predictive_q025 =
        apply(
          predictive_bin_medians,
          2,
          quantile,
          probs = 0.025
        ),
      
      predictive_q975 =
        apply(
          predictive_bin_medians,
          2,
          quantile,
          probs = 0.975
        )
    )
  
  
  plot_data <-
    time_bin_info %>%
    mutate(
      time_bin =
        as.character(time_bin)
    ) %>%
    inner_join(
      predictive_summary,
      by = "time_bin"
    ) %>%
    
    filter(
      
      n_obs >= 5
      
    )
  
  
  p <- ggplot(
    plot_data,
    aes(
      x = time_mid
    )
  ) +
    geom_ribbon(
      aes(
        ymin = predictive_q025,
        ymax = predictive_q975
      ),
      alpha = 0.20
    ) +
    geom_line(
      aes(
        y = predictive_median
      ),
      linewidth = 0.7
    ) +
    geom_point(
      aes(
        y = observed_median
      ),
      size = 2.3
    ) +
    labs(
      title = paste0(
        model_name,
        " model: Full joint treatment model"
      ),
      subtitle =
        "Time-binned posterior predictive check",
      x =
        "Time (days)",
      y =
        "Median log tumour volume"
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      plot.title = element_text(
        face = "bold"
      )
    )
  
  
  list(
    plot = p,
    summary = plot_data
  )
}


time_ppc_exp <-
  make_time_binned_ppc(
    gq_exp,
    "Exponential"
  )

time_ppc_log <-
  make_time_binned_ppc(
    gq_log,
    "Logistic"
  )

time_ppc_gomp <-
  make_time_binned_ppc(
    gq_gomp,
    "Gompertz"
  )


print(time_ppc_exp$plot)
print(time_ppc_log$plot)
print(time_ppc_gomp$plot)

print(
  time_bin_info %>%
    select(
      time_bin,
      time_mid,
      n_obs,
      observed_median
    ),
  n = Inf
)

ggsave(
  file.path(
    figure_dir,
    "stage4_time_binned_ppc_exponential.pdf"
  ),
  time_ppc_exp$plot,
  width = 8,
  height = 5.5
)

ggsave(
  file.path(
    figure_dir,
    "stage4_time_binned_ppc_logistic.pdf"
  ),
  time_ppc_log$plot,
  width = 8,
  height = 5.5
)

ggsave(
  file.path(
    figure_dir,
    "stage4_time_binned_ppc_gompertz.pdf"
  ),
  time_ppc_gomp$plot,
  width = 8,
  height = 5.5
)


write.csv(
  time_ppc_exp$summary,
  file.path(
    table_dir,
    "stage4_time_binned_ppc_exponential.csv"
  ),
  row.names = FALSE
)

write.csv(
  time_ppc_log$summary,
  file.path(
    table_dir,
    "stage4_time_binned_ppc_logistic.csv"
  ),
  row.names = FALSE
)

write.csv(
  time_ppc_gomp$summary,
  file.path(
    table_dir,
    "stage4_time_binned_ppc_gompertz.csv"
  ),
  row.names = FALSE
)

############################################################
# 11.1 FULL SCIENTIFIC PARAMETER CONVERGENCE SUMMARY
############################################################

scientific_parameter_regex <- paste0(
  "^(logV0_pop|logr_pop|logK_pop|",
  "sigma_logV0|sigma_logr|sigma_logK|",
  "logdelta_pop|logalpha_pop|logbeta_a_pop|logbeta_q_pop|logk_rep_pop|",
  "sigma_logdelta|sigma_logalpha|sigma_logbeta_a|sigma_logbeta_q|sigma_logk_rep|",
  "logk_abs_pop|logk_cl_pop|logk_repi_pop|sigma|",
  "V0_mouse\\[|r_mouse\\[|K_mouse\\[|",
  "delta_rt\\[|alpha_rt\\[|beta_a_rt\\[|beta_q_rt\\[|k_rep_rt\\[)"
)

make_full_scientific_summary <- function(fit, model_name) {
  sm <- fit$summary()
  keep <- grepl(scientific_parameter_regex, sm$variable)
  sm <- sm[keep, , drop = FALSE]
  sm$model <- model_name
  sm %>%
    select(model, variable, everything())
}

full_conv_exp <- make_full_scientific_summary(fit_exp_stage4, "Exponential")
full_conv_log <- make_full_scientific_summary(fit_log_stage4, "Logistic")
full_conv_gomp <- make_full_scientific_summary(fit_gomp_stage4, "Gompertz")

full_convergence_stage4 <- bind_rows(
  full_conv_exp,
  full_conv_log,
  full_conv_gomp
)

write.csv(
  full_convergence_stage4,
  file.path(table_dir, "stage4_full_scientific_parameter_convergence.csv"),
  row.names = FALSE
)

full_convergence_compact <- full_convergence_stage4 %>%
  group_by(model) %>%
  summarise(
    n_parameters = n(),
    max_rhat = max(rhat, na.rm = TRUE),
    variable_max_rhat = variable[which.max(rhat)],
    min_bulk_ess = min(ess_bulk, na.rm = TRUE),
    variable_min_bulk_ess = variable[which.min(ess_bulk)],
    min_tail_ess = min(ess_tail, na.rm = TRUE),
    variable_min_tail_ess = variable[which.min(ess_tail)],
    n_rhat_gt_1_01 = sum(rhat > 1.01, na.rm = TRUE),
    n_rhat_gt_1_05 = sum(rhat > 1.05, na.rm = TRUE),
    .groups = "drop"
  )

print(full_convergence_compact)
write.csv(
  full_convergence_compact,
  file.path(table_dir, "stage4_full_convergence_compact.csv"),
  row.names = FALSE
)

############################################################
# STANDARDISED KEY TRACE AND RANK PLOTS
############################################################

key_trace_vars_common <- c(
  "logV0_pop", "logr_pop",
  "logdelta_pop", "logalpha_pop",
  "logbeta_a_pop", "logbeta_q_pop", "logk_rep_pop",
  "logk_abs_pop", "logk_cl_pop", "logk_repi_pop",
  "sigma_logV0", "sigma_logr",
  "sigma_logdelta", "sigma_logalpha",
  "sigma_logbeta_a", "sigma_logbeta_q", "sigma_logk_rep",
  "sigma"
)

save_key_trace_rank <- function(fit, model_name, include_K = FALSE) {
  vars <- key_trace_vars_common
  if (include_K) {
    vars <- c("logV0_pop", "logr_pop", "logK_pop", "sigma_logK",
              setdiff(vars, c("logV0_pop", "logr_pop")))
  }

  draws_arr <- fit$draws(variables = vars)

  p_trace <- bayesplot::mcmc_trace(draws_arr) +
    ggtitle(paste0(model_name, ": key-parameter trace plots"))
  p_rank <- bayesplot::mcmc_rank_overlay(draws_arr) +
    ggtitle(paste0(model_name, ": key-parameter rank plots"))

  ggsave(
    file.path(figure_dir, paste0("stage4_key_trace_", tolower(model_name), ".pdf")),
    p_trace, width = 12, height = 10
  )
  ggsave(
    file.path(figure_dir, paste0("stage4_key_rank_", tolower(model_name), ".pdf")),
    p_rank, width = 12, height = 10
  )
}

save_key_trace_rank(fit_exp_stage4, "Exponential", FALSE)
save_key_trace_rank(fit_log_stage4, "Logistic", TRUE)
save_key_trace_rank(fit_gomp_stage4, "Gompertz", TRUE)

############################################################
#  POSTERIOR PARAMETER DEPENDENCE / IDENTIFIABILITY
############################################################

correlation_parameters <- c(
  "logV0_pop", "logr_pop", "logK_pop",
  "logdelta_pop", "logalpha_pop",
  "logbeta_a_pop", "logbeta_q_pop", "logk_rep_pop",
  "logk_abs_pop", "logk_cl_pop", "logk_repi_pop"
)

make_posterior_correlation <- function(fit, model_name) {
  available <- fit$metadata()$model_params
  vars <- correlation_parameters[correlation_parameters %in% available]

  # Fallback for CmdStanR versions where metadata$model_params is absent
  if (length(vars) < 2L) {
    all_names <- colnames(posterior::as_draws_matrix(fit$draws()))
    vars <- correlation_parameters[correlation_parameters %in% all_names]
  }

  if (length(vars) < 2L) return(NULL)

  d <- posterior::as_draws_df(fit$draws(variables = vars))
  mat <- as.data.frame(d[, vars, drop = FALSE])
  cmat <- cor(mat)

  long <- as.data.frame(as.table(cmat), stringsAsFactors = FALSE)
  names(long) <- c("parameter_1", "parameter_2", "correlation")
  long$model <- model_name

  write.csv(
    long,
    file.path(table_dir, paste0("stage4_posterior_correlations_", tolower(model_name), ".csv")),
    row.names = FALSE
  )

  p <- ggplot(long, aes(parameter_1, parameter_2, fill = correlation)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", correlation)), size = 2.7) +
    scale_fill_gradient2(limits = c(-1, 1), midpoint = 0) +
    coord_equal() +
    labs(
      title = paste0(model_name, ": posterior parameter correlations"),
      x = NULL, y = NULL, fill = "Correlation"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )

  ggsave(
    file.path(figure_dir, paste0("stage4_posterior_correlations_", tolower(model_name), ".pdf")),
    p, width = 10, height = 8
  )

  long
}

corr_exp <- make_posterior_correlation(fit_exp_stage4, "Exponential")
corr_log <- make_posterior_correlation(fit_log_stage4, "Logistic")
corr_gomp <- make_posterior_correlation(fit_gomp_stage4, "Gompertz")

corr_all <- bind_rows(corr_exp, corr_log, corr_gomp)
if (nrow(corr_all) > 0L) {
  write.csv(
    corr_all,
    file.path(table_dir, "stage4_posterior_correlations_all_models.csv"),
    row.names = FALSE
  )
}

############################################################
#  BUILT-IN DNA-REPAIR MULTIPLIER GQ CHECK
############################################################

repair_gq_summary <- function(fit, model_name) {
  min_factor <- as.numeric(
    posterior::as_draws_matrix(
      fit$draws(variables = "min_repair_factor_checked")
    )[, 1]
  )
  negative_count <- as.numeric(
    posterior::as_draws_matrix(
      fit$draws(variables = "negative_repair_factor_count")
    )[, 1]
  )

  data.frame(
    model = model_name,
    median_min_repair_factor = median(min_factor),
    q025_min_repair_factor = quantile(min_factor, 0.025),
    q975_min_repair_factor = quantile(min_factor, 0.975),
    minimum_over_draws = min(min_factor),
    probability_any_negative_checked = mean(negative_count > 0),
    maximum_negative_count_checked = max(negative_count)
  )
}

repair_gq_all <- bind_rows(
  repair_gq_summary(fit_exp_stage4, "Exponential"),
  repair_gq_summary(fit_log_stage4, "Logistic"),
  repair_gq_summary(fit_gomp_stage4, "Gompertz")
)

print(repair_gq_all)
write.csv(
  repair_gq_all,
  file.path(table_dir, "stage4_builtin_repair_factor_gq_summary.csv"),
  row.names = FALSE
)

############################################################
#  POSTERIOR PREDICTIVE COVERAGE BY MOUSE
############################################################

coverage_by_mouse <- prediction_stage4 %>%
  mutate(
    covered_95 = volume >= pred_q025 & volume <= pred_q975
  ) %>%
  group_by(model, mouse, mouse_id, group, group_name) %>%
  summarise(
    n_obs = n(),
    coverage_95 = mean(covered_95),
    mean_abs_log_residual = mean(abs(log_residual)),
    rmse_log = sqrt(mean(log_residual^2)),
    .groups = "drop"
  )

write.csv(
  coverage_by_mouse,
  file.path(table_dir, "stage4_predictive_coverage_by_mouse.csv"),
  row.names = FALSE
)

p_coverage_mouse <- ggplot(
  coverage_by_mouse,
  aes(x = reorder(as.factor(mouse), coverage_95), y = coverage_95)
) +
  geom_point() +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  facet_grid(model ~ group_name, scales = "free_x", space = "free_x") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Stage 4: 95% posterior predictive coverage by mouse",
    x = "Mouse",
    y = "Empirical 95% predictive coverage"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(figure_dir, "stage4_predictive_coverage_by_mouse.pdf"),
  p_coverage_mouse,
  width = 14,
  height = 8
)

############################################################
#  TARGETED POSTERIOR PREDICTIVE TEST STATISTICS
############################################################

ppc_statistic_summary <- function(gq, model_name) {
  y_obs <- log(obs_stage4$volume)
  y_rep <- log(gq$y_rep)

  obs_stats <- c(
    mean_log_volume = mean(y_obs),
    sd_log_volume = sd(y_obs),
    median_log_volume = median(y_obs),
    max_log_volume = max(y_obs)
  )

  rep_stats <- cbind(
    mean_log_volume = rowMeans(y_rep),
    sd_log_volume = apply(y_rep, 1, sd),
    median_log_volume = apply(y_rep, 1, median),
    max_log_volume = apply(y_rep, 1, max)
  )

  bind_rows(lapply(seq_along(obs_stats), function(i) {
    nm <- names(obs_stats)[i]
    x <- rep_stats[, nm]
    data.frame(
      model = model_name,
      statistic = nm,
      observed = unname(obs_stats[i]),
      predictive_median = median(x),
      predictive_q025 = quantile(x, 0.025),
      predictive_q975 = quantile(x, 0.975),
      posterior_predictive_tail_probability = mean(x >= obs_stats[i])
    )
  }))
}

ppc_stats_all <- bind_rows(
  ppc_statistic_summary(gq_exp, "Exponential"),
  ppc_statistic_summary(gq_log, "Logistic"),
  ppc_statistic_summary(gq_gomp, "Gompertz")
)

print(ppc_stats_all)
write.csv(
  ppc_stats_all,
  file.path(table_dir, "stage4_targeted_ppc_statistics.csv"),
  row.names = FALSE
)

############################################################
#  PARETO-k INFLUENCE SUMMARISED BY MOUSE
############################################################

pareto_mouse_summary <- pareto_stage4 %>%
  group_by(model, mouse, mouse_id, group, group_name) %>%
  summarise(
    n_obs = n(),
    max_pareto_k = max(pareto_k),
    mean_pareto_k = mean(pareto_k),
    n_k_gt_0_7 = sum(pareto_k > 0.7),
    n_k_gt_1 = sum(pareto_k > 1),
    .groups = "drop"
  ) %>%
  arrange(model, desc(max_pareto_k))

write.csv(
  pareto_mouse_summary,
  file.path(table_dir, "stage4_pareto_k_by_mouse.csv"),
  row.names = FALSE
)

p_pareto_mouse <- ggplot(
  pareto_mouse_summary,
  aes(x = reorder(as.factor(mouse), max_pareto_k), y = max_pareto_k)
) +
  geom_point() +
  geom_hline(yintercept = 0.7, linetype = "dashed") +
  geom_hline(yintercept = 1.0, linetype = "dotted") +
  facet_grid(model ~ group_name, scales = "free_x", space = "free_x") +
  labs(
    title = "Stage 4: maximum observation-level Pareto-k by mouse",
    x = "Mouse",
    y = "Maximum Pareto-k"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(figure_dir, "stage4_pareto_k_by_mouse.pdf"),
  p_pareto_mouse,
  width = 14,
  height = 8
)

############################################################
#  OPTIONAL LOO-PIT CALIBRATION CHECK
############################################################


compute_loo_pit_manual <- function(gq, model_name) {
  out <- tryCatch({
    ps <- loo::psis(-gq$log_lik)
    w <- loo::weights(ps, normalize = TRUE, log = FALSE)
    y <- obs_stage4$volume
    yrep <- gq$y_rep

    pit <- vapply(seq_along(y), function(j) {
      sum(w[, j] * (yrep[, j] <= y[j]))
    }, numeric(1))

    d <- obs_stage4 %>%
      mutate(model = model_name, loo_pit = pit)

    list(ok = TRUE, data = d, message = NA_character_)
  }, error = function(e) {
    list(ok = FALSE, data = NULL, message = conditionMessage(e))
  })
  out
}

pit_results <- list(
  Exponential = compute_loo_pit_manual(gq_exp, "Exponential"),
  Logistic = compute_loo_pit_manual(gq_log, "Logistic"),
  Gompertz = compute_loo_pit_manual(gq_gomp, "Gompertz")
)

pit_ok <- vapply(pit_results, `[[`, logical(1), "ok")
if (all(pit_ok)) {
  loo_pit_stage4 <- bind_rows(lapply(pit_results, `[[`, "data"))
  write.csv(
    loo_pit_stage4,
    file.path(table_dir, "stage4_loo_pit_values.csv"),
    row.names = FALSE
  )

  p_loo_pit <- ggplot(loo_pit_stage4, aes(x = loo_pit)) +
    geom_histogram(bins = 10, boundary = 0, closed = "left") +
    facet_wrap(~ model, ncol = 1) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(
      title = "Stage 4: LOO-PIT calibration diagnostic",
      x = "LOO-PIT",
      y = "Count"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(
    file.path(figure_dir, "stage4_loo_pit_histograms.pdf"),
    p_loo_pit,
    width = 8,
    height = 9
  )
} else {
  pit_status <- data.frame(
    model = names(pit_results),
    success = pit_ok,
    message = vapply(pit_results, function(x) {
      if (is.null(x$message) || is.na(x$message)) "" else x$message
    }, character(1))
  )
  write.csv(
    pit_status,
    file.path(table_dir, "stage4_loo_pit_status.csv"),
    row.names = FALSE
  )
  warning("LOO-PIT could not be computed for at least one model; see status CSV.")
}

############################################################
#  MASTER SUMMARY / SESSION INFORMATION
############################################################

capture.output(
  {
    cat("===== STAGE 4 FINAL POSTERIOR ANALYSIS =====\n\n")
    cat("===== HMC MODEL SUMMARY =====\n")
    print(hmc_model_summary)
    cat("\n===== FULL CONVERGENCE SUMMARY =====\n")
    print(full_convergence_compact)
    cat("\n===== RESIDUAL / PPC SUMMARY =====\n")
    print(residual_ppc_summary)
    cat("\n===== RESIDUAL / PPC BY GROUP =====\n")
    print(residual_ppc_by_group)
    cat("\n===== LOO COMPARISON =====\n")
    print(loo_compare_stage4)
    cat("\n===== PARETO-k SUMMARY =====\n")
    print(pareto_summary_stage4)
    cat("\n===== BUILT-IN REPAIR FACTOR CHECK =====\n")
    print(repair_gq_all)
    cat("\n===== TARGETED PPC STATISTICS =====\n")
    print(ppc_stats_all)
  },
  file = file.path(analysis_dir, "stage4_master_summary.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(analysis_dir, "stage4_session_info.txt")
)

saveRDS(
  list(
    hmc = hmc_summary,
    convergence = full_convergence_stage4,
    natural_scale = natural_stage4,
    predictions = prediction_stage4,
    residual_ppc = residual_ppc_summary,
    residual_ppc_by_group = residual_ppc_by_group,
    loo = list(
      exponential = loo_exp_stage4,
      logistic = loo_log_stage4,
      gompertz = loo_gomp_stage4,
      comparison = loo_compare_stage4
    ),
    pareto = pareto_stage4,
    repair_factor = repair_gq_all,
    ppc_statistics = ppc_stats_all
  ),
  file = file.path(object_dir, "stage4_final_analysis_objects.rds")
)

cat("\n============================================================\n")
cat("STAGE 4 FINAL POSTERIOR ANALYSIS COMPLETE\n")
cat("Tables :", table_dir, "\n")
cat("Figures:", figure_dir, "\n")
cat("Objects:", object_dir, "\n")
cat("============================================================\n")

############################################################
# STAGE 4 K_typical
############################################################

library(posterior)

# Logistic
K_log_stage4 <- posterior::as_draws_matrix(
  fit_log_stage4$draws(variables = "K_mouse")
)

K_log_stage4 <- K_log_stage4[
  ,
  grepl("^K_mouse\\[", colnames(K_log_stage4)),
  drop = FALSE
]

K_typical_log_stage4 <- apply(
  K_log_stage4,
  1,
  median
)

K_typical_log_stage4_summary <- quantile(
  K_typical_log_stage4,
  probs = c(0.025, 0.5, 0.975)
)

cat("\n===== STAGE 4 LOGISTIC K_typical =====\n")
print(K_typical_log_stage4_summary)


# Gompertz
K_gomp_stage4 <- posterior::as_draws_matrix(
  fit_gomp_stage4$draws(variables = "K_mouse")
)

K_gomp_stage4 <- K_gomp_stage4[
  ,
  grepl("^K_mouse\\[", colnames(K_gomp_stage4)),
  drop = FALSE
]

K_typical_gomp_stage4 <- apply(
  K_gomp_stage4,
  1,
  median
)

K_typical_gomp_stage4_summary <- quantile(
  K_typical_gomp_stage4,
  probs = c(0.025, 0.5, 0.975)
)

cat("\n===== STAGE 4 GOMPERTZ K_typical =====\n")
print(K_typical_gomp_stage4_summary)

############################################################
# STAGE 4: INDIVIDUAL COMBINATION-MOUSE TRAJECTORIES
# POSTERIOR ANALYSIS ONLY 
############################################################

library(dplyr)
library(ggplot2)

combination_trajectory_data <- prediction_stage4 %>%
  filter(
    group_name == "combination"
  )

p_combination_trajectories <- ggplot(
  combination_trajectory_data,
  aes(x = time_days)
) +
  geom_ribbon(
    aes(
      ymin = pred_q025,
      ymax = pred_q975
    ),
    alpha = 0.20
  ) +
  geom_line(
    aes(y = pred_median),
    linewidth = 0.55
  ) +
  geom_point(
    aes(y = volume),
    size = 1.3
  ) +
  facet_grid(
    model ~ mouse,
    scales = "free_y"
  ) +
  labs(
    x = "Time (days)",
    y = expression("Tumour volume (mm"^3*")")
  ) +
  theme_bw(
    base_size = 10
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8)
  )

print(p_combination_trajectories)

ggsave(
  file.path(
    figure_dir,
    "stage4_combination_mouse_trajectories.pdf"
  ),
  p_combination_trajectories,
  width = 14,
  height = 8
)