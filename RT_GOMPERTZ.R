##################################################
#RT-stage: CONTROL + RT-ONLY MICE — GOMPERTZ
##################################################

library(cmdstanr)
cmdstanr::set_cmdstan_path("/home/ucahhor/Scratch/cmdstan-2.35.0")
library(dplyr)
library(tidyr)
library(posterior)
library(loo)
library(bayesplot)
library(ggplot2)
library(patchwork)
library(deSolve)

project_dir <- if (dir.exists(path.expand("~/Scratch/dissertation"))) {
  path.expand("~/Scratch/dissertation")
} else {
  path.expand("~/Desktop/dissertation")
}
setwd(project_dir)

data <- read.csv(file.path(project_dir, "data", "lemassonParpi.csv"))

results_dir <- "results/stage2_mfull_gompertz_ctrl_rt"
figures_dir <- "figures/stage2_mfull_gompertz_ctrl_rt"
stan_dir <- "stan"
csv_dir <- file.path(results_dir, "cmdstan_csv")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stan_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

RUN_INTERMEDIATE_DIAGNOSTICS <- FALSE
RUN_FULL_FITS <- TRUE
RUN_DOWNSTREAM_ANALYSIS <- FALSE

mouse_groups <- data %>%
  group_by(mouse) %>%
  summarise(
    has_drug = any(modality == 1, na.rm = TRUE),
    has_rt = any(modality == 2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    group = case_when(
      !has_drug & !has_rt ~ "control",
      !has_drug &  has_rt ~ "rt_only",
      TRUE                ~ "exclude"
    )
  )

control_mice <- mouse_groups %>%
  filter(group == "control") %>%
  pull(mouse) %>%
  sort()

rt_only_mice <- mouse_groups %>%
  filter(group == "rt_only") %>%
  pull(mouse) %>%
  sort()

included_mice <- sort(c(control_mice, rt_only_mice))

cat("Control mice:", length(control_mice), "\n")
cat("RT-only mice:", length(rt_only_mice), "\n")
cat("Total Stage 2 mice:", length(included_mice), "\n")

stopifnot(
  length(control_mice) == 23L,
  length(rt_only_mice) == 33L,
  length(included_mice) == 56L,
  length(intersect(control_mice, rt_only_mice)) == 0L
)

rt_obs <- data %>%
  filter(
    mouse %in% included_mice,
    observation == 1
  ) %>%
  mutate(time_days = time / 24) %>%
  arrange(mouse, time_days)

mouse_lookup <- tibble(
  mouse = sort(unique(rt_obs$mouse))
) %>%
  mutate(
    mouse_id = row_number(),
    is_rt = mouse %in% rt_only_mice
  )

rt_mouse_lookup <- mouse_lookup %>%
  filter(is_rt) %>%
  mutate(rt_index = row_number()) %>%
  select(mouse, mouse_id, rt_index)

mouse_lookup <- mouse_lookup %>%
  left_join(
    rt_mouse_lookup %>% select(mouse, rt_index),
    by = "mouse"
  ) %>%
  mutate(rt_index = ifelse(is.na(rt_index), 0L, as.integer(rt_index)))

rt_obs <- rt_obs %>%
  left_join(mouse_lookup, by = "mouse") %>%
  arrange(mouse_id, time_days)

rt_events <- data %>%
  filter(
    mouse %in% rt_only_mice,
    modality == 2
  ) %>%
  mutate(time_days = time / 24) %>%
  left_join(
    mouse_lookup %>% select(mouse, mouse_id, rt_index),
    by = "mouse"
  ) %>%
  arrange(mouse_id, time_days)

stopifnot(
  nrow(rt_obs) > 0,
  nrow(rt_events) > 0,
  all(!is.na(rt_obs$mouse)),
  all(!is.na(rt_obs$mouse_id)),
  all(is.finite(rt_obs$time_days)),
  all(is.finite(rt_obs$volume)),
  all(rt_obs$time_days >= 0),
  all(rt_obs$volume > 0),
  all(!is.na(rt_events$mouse)),
  all(!is.na(rt_events$mouse_id)),
  all(!is.na(rt_events$rt_index)),
  all(rt_events$rt_index > 0),
  all(is.finite(rt_events$time_days)),
  all(rt_events$time_days >= 0),
  all(!is.na(rt_events$dose)),
  all(is.finite(rt_events$dose)),
  all(rt_events$dose > 0),
  setequal(unique(rt_obs$mouse), included_mice),
  setequal(unique(rt_events$mouse), rt_only_mice),
  nrow(mouse_lookup) == 56L,
  nrow(rt_mouse_lookup) == 33L,
  sum(mouse_lookup$is_rt) == 33L,
  sum(mouse_lookup$rt_index == 0L) == 23L,
  setequal(rt_mouse_lookup$rt_index, seq_len(33L))
)

duplicate_rt_observations <- rt_obs %>%
  count(mouse, time_days) %>%
  filter(n > 1)
print(duplicate_rt_observations)
stopifnot(nrow(duplicate_rt_observations) == 0)

duplicate_rt_events <- rt_events %>%
  count(mouse, time_days, dose) %>%
  filter(n > 1)
print(duplicate_rt_events)

rt_dataset_summary <- data.frame(
  number_of_mice = n_distinct(rt_obs$mouse),
  number_of_control_mice = length(control_mice),
  number_of_rt_mice = length(rt_only_mice),
  number_of_observations = nrow(rt_obs),
  number_of_rt_events = nrow(rt_events),
  minimum_observation_time_days = min(rt_obs$time_days),
  maximum_observation_time_days = max(rt_obs$time_days),
  minimum_volume = min(rt_obs$volume),
  maximum_volume = max(rt_obs$volume),
  minimum_rt_dose = min(rt_events$dose),
  maximum_rt_dose = max(rt_events$dose)
)
print(rt_dataset_summary)
write.csv(rt_dataset_summary,
          file.path(results_dir, "ctrl_rt_dataset_summary.csv"),
          row.names = FALSE)

rt_schedule_by_mouse <- rt_events %>%
  group_by(mouse, mouse_id) %>%
  summarise(
    number_of_rt_fractions = n(),
    first_rt_time_days = min(time_days),
    last_rt_time_days = max(time_days),
    total_rt_dose = sum(dose),
    minimum_fraction_dose = min(dose),
    maximum_fraction_dose = max(dose),
    .groups = "drop"
  )
print(rt_schedule_by_mouse, n = Inf)
write.csv(rt_schedule_by_mouse,
          file.path(results_dir, "rt_schedule_by_mouse.csv"),
          row.names = FALSE)

first_rt_by_mouse <- rt_events %>%
  group_by(mouse) %>%
  summarise(first_rt_time_days = min(time_days), .groups = "drop")

pretreatment_summary <- rt_obs %>%
  left_join(first_rt_by_mouse, by = "mouse") %>%
  mutate(is_rt = mouse %in% rt_only_mice) %>%
  group_by(mouse, is_rt) %>%
  summarise(
    n_observations = n(),
    n_before_first_rt = if (first(is_rt)) {
      sum(time_days < first(first_rt_time_days))
    } else {
      NA_integer_
    },
    n_at_or_after_first_rt = if (first(is_rt)) {
      sum(time_days >= first(first_rt_time_days))
    } else {
      NA_integer_
    },
    .groups = "drop"
  )
print(pretreatment_summary, n = Inf)
write.csv(pretreatment_summary,
          file.path(results_dir, "pretreatment_observation_summary.csv"),
          row.names = FALSE)

write.csv(rt_obs, file.path(results_dir, "ctrl_rt_observations.csv"), row.names = FALSE)
write.csv(rt_events, file.path(results_dir, "rt_events.csv"), row.names = FALSE)

rt_baseline <- rt_obs %>%
  group_by(mouse, mouse_id) %>%
  slice_min(order_by = time_days, n = 1, with_ties = FALSE) %>%
  transmute(mouse, mouse_id, t0 = time_days, V_t0 = volume) %>%
  ungroup() %>%
  arrange(mouse_id)

stopifnot(
  nrow(rt_baseline) == nrow(mouse_lookup),
  all(rt_baseline$mouse_id == mouse_lookup$mouse_id),
  all(is.finite(rt_baseline$t0)),
  all(is.finite(rt_baseline$V_t0)),
  all(rt_baseline$V_t0 > 0)
)

rt_baseline_check <- rt_baseline %>%
  left_join(first_rt_by_mouse, by = "mouse") %>%
  mutate(is_rt = mouse %in% rt_only_mice)

stopifnot(
  all(
    rt_baseline_check$t0[rt_baseline_check$is_rt] <=
      rt_baseline_check$first_rt_time_days[rt_baseline_check$is_rt]
  ),
  all(is.na(rt_baseline_check$first_rt_time_days[!rt_baseline_check$is_rt]))
)

nonzero_baseline_times <- rt_baseline %>% filter(t0 != 0)
print(nonzero_baseline_times)

rt_obs_postbaseline <- rt_obs %>%
  left_join(rt_baseline %>% select(mouse, t0), by = "mouse") %>%
  filter(time_days > t0) %>%
  select(-t0) %>%
  arrange(mouse_id, time_days) %>%
  mutate(obs_index = row_number())

stopifnot(
  nrow(rt_obs_postbaseline) > 0,
  all(rt_obs_postbaseline$volume > 0),
  all(is.finite(rt_obs_postbaseline$time_days))
)
cat("Post-baseline observations:", nrow(rt_obs_postbaseline), "\n")

rt_events_collapsed <- rt_events %>%
  group_by(mouse, mouse_id, time_days) %>%
  summarise(dose = sum(dose), .groups = "drop") %>%
  arrange(mouse_id, time_days)

same_time_rt_and_observation <- inner_join(
  rt_obs_postbaseline %>% select(mouse, observation_time = time_days),
  rt_events_collapsed %>% select(mouse, rt_time = time_days, dose),
  by = "mouse",
  relationship = "many-to-many"
) %>%
  filter(observation_time == rt_time)
print(same_time_rt_and_observation)
stopifnot(nrow(same_time_rt_and_observation) == 0)


segment_list <- list()
solve_time_list <- list()
segment_counter <- 0L
solve_counter <- 0L

for (m in seq_len(nrow(rt_baseline))) {
  
  mouse_current <- rt_baseline$mouse[m]
  t0_current <- rt_baseline$t0[m]
  
  obs_current <- rt_obs_postbaseline %>%
    filter(mouse_id == m) %>%
    arrange(time_days)
  stopifnot(nrow(obs_current) > 0)
  
  final_time_current <- max(obs_current$time_days)
  
  rt_current <- rt_events_collapsed %>%
    filter(mouse_id == m, time_days <= final_time_current) %>%
    arrange(time_days)
  
  segment_start_times <- sort(unique(c(
    t0_current,
    rt_current$time_days[rt_current$time_days < final_time_current]
  )))
  
  for (s_local in seq_along(segment_start_times)) {
    
    segment_start <- segment_start_times[s_local]
    
    if (s_local < length(segment_start_times)) {
      segment_end <- segment_start_times[s_local + 1L]
    } else {
      segment_end <- final_time_current
    }
    
    if (segment_end <= segment_start) next
    
    dose_at_start <- sum(
      rt_current$dose[rt_current$time_days == segment_start]
    )
    
    obs_segment <- obs_current %>%
      filter(time_days > segment_start, time_days <= segment_end)
    
    output_times <- sort(unique(c(obs_segment$time_days, segment_end)))
    
    output_obs_index <- match(output_times, obs_segment$time_days)
    output_obs_index[is.na(output_obs_index)] <- 0L
    positive_index <- output_obs_index > 0
    output_obs_index[positive_index] <-
      obs_segment$obs_index[output_obs_index[positive_index]]
    
    segment_counter <- segment_counter + 1L
    solve_start_current <- solve_counter + 1L
    
    solve_time_list[[segment_counter]] <- data.frame(
      segment_id = segment_counter,
      time_days = output_times,
      obs_index = as.integer(output_obs_index)
    )
    
    solve_counter <- solve_counter + length(output_times)
    solve_end_current <- solve_counter
    
    segment_list[[segment_counter]] <- data.frame(
      segment_id = segment_counter,
      mouse = mouse_current,
      mouse_id = m,
      segment_start_time = segment_start,
      segment_end_time = segment_end,
      rt_dose_at_start = dose_at_start,
      solve_start = solve_start_current,
      solve_end = solve_end_current
    )
  }
}

stan_segment_table <- bind_rows(segment_list)
stan_solve_table <- bind_rows(solve_time_list)

segment_index_by_mouse <- stan_segment_table %>%
  group_by(mouse, mouse_id) %>%
  summarise(
    segment_start = min(segment_id),
    segment_end = max(segment_id),
    n_segments = n(),
    .groups = "drop"
  ) %>%
  arrange(mouse_id)

mapped_observation_indices <-
  stan_solve_table$obs_index[stan_solve_table$obs_index > 0]

stopifnot(
  length(mapped_observation_indices) == nrow(rt_obs_postbaseline),
  length(unique(mapped_observation_indices)) == nrow(rt_obs_postbaseline),
  setequal(mapped_observation_indices, rt_obs_postbaseline$obs_index)
)

segment_timing_check <- stan_solve_table %>%
  left_join(
    stan_segment_table %>%
      select(segment_id, segment_start_time, segment_end_time),
    by = "segment_id"
  )

stopifnot(
  all(segment_timing_check$time_days > segment_timing_check$segment_start_time),
  all(segment_timing_check$time_days <= segment_timing_check$segment_end_time),
  all(is.finite(segment_timing_check$time_days)),
  nrow(segment_index_by_mouse) == nrow(rt_baseline),
  all(segment_index_by_mouse$mouse_id == rt_baseline$mouse_id),
  all(stan_segment_table$segment_end_time > stan_segment_table$segment_start_time)
)

write.csv(stan_segment_table,
          file.path(results_dir, "stan_solver_segments.csv"),
          row.names = FALSE)
write.csv(stan_solve_table,
          file.path(results_dir, "stan_solver_output_times.csv"),
          row.names = FALSE)

M <- nrow(rt_baseline)
M_rt <- nrow(rt_mouse_lookup)
N_obs <- nrow(rt_obs_postbaseline)
N_segments <- nrow(stan_segment_table)
N_solve <- nrow(stan_solve_table)

cat(
  "\nNumber of mice:", M,
  "\nNumber of RT mice:", M_rt,
  "\nNumber of control mice:", M - M_rt,
  "\nNumber of post-baseline observations:", N_obs,
  "\nNumber of ODE segments:", N_segments,
  "\nNumber of solver output times:", N_solve,
  "\n"
)
ode_settings_full <- list(
  rel_tol = 1e-6,
  abs_tol = 1e-6,
  max_num_steps = 1000000L
)

ode_settings_validation <- list(
  rel_tol = 1e-4,
  abs_tol = 1e-4,
  max_num_steps = 50000L
)

rt_common_stan_data <- list(
  M = M,
  M_rt = M_rt,
  rt_index = mouse_lookup$rt_index,
  N_obs = N_obs,
  N_segments = N_segments,
  N_solve = N_solve,
  segment_start_time = stan_segment_table$segment_start_time,
  segment_rt_dose = stan_segment_table$rt_dose_at_start,
  solve_start = stan_segment_table$solve_start,
  solve_end = stan_segment_table$solve_end,
  solve_time = stan_solve_table$time_days,
  solve_obs_index = stan_solve_table$obs_index,
  mouse_segment_start = segment_index_by_mouse$segment_start,
  mouse_segment_end = segment_index_by_mouse$segment_end,
  baseline_volume = rt_baseline$V_t0,
  volume = rt_obs_postbaseline$volume
)

rt_exp_stan_data <- c(rt_common_stan_data, ode_settings_full)
rt_exp_stan_data_validation <- c(rt_common_stan_data, ode_settings_validation)

rt_logistic_stan_data <- c(
  rt_common_stan_data,
  ode_settings_full
)
rt_logistic_stan_data_validation <- c(
  rt_common_stan_data,
  ode_settings_validation
)

rt_gompertz_stan_data <- c(
  rt_common_stan_data,
  ode_settings_full
)

rt_gompertz_stan_data_validation <- c(
  rt_common_stan_data,
  ode_settings_validation
)


##################################################
#  Stage-2 Stan models
##################################################

make_rt_exponential_stan_model <- function() {
  paste0(
    "functions {

  vector rt_exponential_reduced_ode(
      real t,
      vector y,
      vector theta,
      real P_start,
      real D_start,
      real t_start
  ) {
    real r = theta[1];
    real delta = theta[2];
    real beta_a = theta[3];
    real beta_q = theta[4];
    real k_rep = theta[5];

    real A = y[1];
    real dt = t - t_start;
    real repair_factor = exp(-k_rep * dt);
    real D = D_start * repair_factor;
    real log_P_ratio =
      r * dt
      - ((beta_a + beta_q) * D_start / k_rep)
        * (1 - repair_factor);
    real P = P_start * exp(log_P_ratio);

    vector[2] dydt;
    dydt[1] = beta_a * D * P - delta * A;
    dydt[2] = r * P - delta * A;
    return dydt;
  }
}

data {
  int<lower=1> M;
  int<lower=1> M_rt;
  array[M] int<lower=0, upper=M_rt> rt_index;
  int<lower=1> N_obs;
  int<lower=1> N_segments;
  int<lower=1> N_solve;

  array[N_segments] real<lower=0> segment_start_time;
  array[N_segments] real<lower=0> segment_rt_dose;
  array[N_segments] int<lower=1, upper=N_solve> solve_start;
  array[N_segments] int<lower=1, upper=N_solve> solve_end;
  array[N_solve] real<lower=0> solve_time;
  array[N_solve] int<lower=0, upper=N_obs> solve_obs_index;
  array[M] int<lower=1, upper=N_segments> mouse_segment_start;
  array[M] int<lower=1, upper=N_segments> mouse_segment_end;

  vector<lower=0>[M] baseline_volume;
  vector<lower=0>[N_obs] volume;

  real<lower=0> rel_tol;
  real<lower=0> abs_tol;
  int<lower=1> max_num_steps;
}

parameters {
  real logV0_pop;
  real logr_pop;
  real logalpha_pop;

  real<lower=0> sigma_logV0;
  real<lower=0> sigma_logr;
  real<lower=0> sigma_logalpha;

  vector[M] z_logV0;
  vector[M] z_logr;
  vector[M_rt] z_logalpha;

  real logdelta_pop;
  real logbeta_a_pop;
  real logbeta_q_pop;
  real logk_rep_pop;

  real<lower=0> sigma;
}

transformed parameters {
  real V0_pop = exp(logV0_pop);
  real r_pop = exp(logr_pop);
  real alpha_pop = exp(logalpha_pop);
  real delta_pop = exp(logdelta_pop);
  real beta_a_pop = exp(logbeta_a_pop);
  real beta_q_pop = exp(logbeta_q_pop);
  real k_rep_pop = exp(logk_rep_pop);

  vector[M] V0_mouse;
  vector[M] r_mouse;
  vector[M_rt] alpha_mouse;
  vector[N_obs] mu_out;
  vector[N_obs] mu_log_out;

  for (m in 1:M) {
    V0_mouse[m] = exp(logV0_pop + sigma_logV0 * z_logV0[m]);
    r_mouse[m] = exp(logr_pop + sigma_logr * z_logr[m]);
  }
  for (j in 1:M_rt) {
    alpha_mouse[j] = exp(logalpha_pop + sigma_logalpha * z_logalpha[j]);
  }

  for (m in 1:M) {
    vector[2] y_current;
    real P_current = V0_mouse[m];
    real D_current = 0;
    vector[5] theta;

    y_current[1] = 0;              // A
    y_current[2] = V0_mouse[m];    // total V

    theta[1] = r_mouse[m];
    theta[2] = delta_pop;
    theta[3] = beta_a_pop;
    theta[4] = beta_q_pop;
    theta[5] = k_rep_pop;

    for (s in mouse_segment_start[m]:mouse_segment_end[m]) {
      real interval_start = segment_start_time[s];
      if (rt_index[m] > 0)
        D_current += alpha_mouse[rt_index[m]] * segment_rt_dose[s];

      int n_times = solve_end[s] - solve_start[s] + 1;
      array[n_times] real solve_times;
      array[n_times] vector[2] y_solution;

      for (j in 1:n_times)
        solve_times[j] = solve_time[solve_start[s] + j - 1];

      y_solution = ode_rk45_tol(
        rt_exponential_reduced_ode,
        y_current,
        interval_start,
        solve_times,
        rel_tol,
        abs_tol,
        max_num_steps,
        theta,
        P_current,
        D_current,
        interval_start
      );

      for (j in 1:n_times) {
        int g = solve_start[s] + j - 1;
        int obs = solve_obs_index[g];
        if (obs > 0) {
          real V_pred = fmax(y_solution[j][2], 1e-10);
          mu_out[obs] = V_pred;
          mu_log_out[obs] = log(V_pred);
        }
      }

      {
        real dt_end = solve_times[n_times] - interval_start;
        real repair_factor_end = exp(-k_rep_pop * dt_end);
        real log_P_ratio_end =
          r_mouse[m] * dt_end
          - ((beta_a_pop + beta_q_pop) * D_current / k_rep_pop)
            * (1 - repair_factor_end);

        P_current *= exp(log_P_ratio_end);
        y_current = y_solution[n_times];
        D_current *= repair_factor_end;
      }
    }
  }
}

model {
  logV0_pop ~ normal(log(30), 0.5);
  logr_pop ~ normal(log(0.3), 0.5);
  sigma_logV0 ~ normal(0, 0.5);
  sigma_logr ~ normal(0, 0.5);
  z_logV0 ~ std_normal();
  z_logr ~ std_normal();

  logalpha_pop ~ normal(log(0.9), 0.5);
  sigma_logalpha ~ normal(0, 0.5);
  z_logalpha ~ std_normal();

  logdelta_pop ~ normal(log(1.2), 0.5);
  logbeta_a_pop ~ normal(log(26.4), 0.5);
  logbeta_q_pop ~ normal(log(4.8), 0.5);
  logk_rep_pop ~ normal(log(86.4), 0.5);

  sigma ~ normal(0, 0.5);

  log(baseline_volume) ~ normal(log(V0_mouse), sigma);
  log(volume) ~ normal(mu_log_out, sigma);
}

generated quantities {
  vector[M + N_obs] log_lik;
  vector[M] y_rep_baseline;
  vector[N_obs] y_rep;

  for (m in 1:M) {
    log_lik[m] = lognormal_lpdf(
  baseline_volume[m] |
  log(V0_mouse[m]),
  sigma
);
    y_rep_baseline[m] = lognormal_rng(log(V0_mouse[m]), sigma);
  }
  for (n in 1:N_obs) {
    log_lik[M + n] = lognormal_lpdf(
  volume[n] |
  mu_log_out[n],
  sigma
);
    y_rep[n] = lognormal_rng(mu_log_out[n], sigma);
  }
}")
}

make_rt_bounded_stan_model <- function(growth_model = c("logistic", "gompertz")) {
  growth_model <- match.arg(growth_model)
  function_name <- paste0("rt_", growth_model, "_ode")
  growth_rhs <- switch(
    growth_model,
    logistic = "r * P * (1 - V / K)",
    gompertz = "r * P * log(K / V)"
  )
  growth_prior <- if (growth_model == "logistic") "0.4" else "0.2"
  
  paste0(
    "functions {
  vector ", function_name, "(
      real t,
      vector y,
      vector theta,
      real D_start,
      real t_start
  ) {
    real r = theta[1];
    real delta = theta[2];
    real beta_a = theta[3];
    real beta_q = theta[4];
    real k_rep = theta[5];
    real K = theta[6];

    real P = y[1];
    real Q = y[2];
    real A = y[3];
    real V = fmax(P + Q + A, 1e-10);
    real D = D_start * exp(-k_rep * (t - t_start));

    vector[3] dydt;
    dydt[1] = ", growth_rhs, " - beta_a * D * P - beta_q * D * P;
    dydt[2] = beta_q * D * P;
    dydt[3] = beta_a * D * P - delta * A;
    return dydt;
  }
}

data {
  int<lower=1> M;
  int<lower=1> M_rt;
  array[M] int<lower=0, upper=M_rt> rt_index;
  int<lower=1> N_obs;
  int<lower=1> N_segments;
  int<lower=1> N_solve;

  array[N_segments] real<lower=0> segment_start_time;
  array[N_segments] real<lower=0> segment_rt_dose;
  array[N_segments] int<lower=1, upper=N_solve> solve_start;
  array[N_segments] int<lower=1, upper=N_solve> solve_end;
  array[N_solve] real<lower=0> solve_time;
  array[N_solve] int<lower=0, upper=N_obs> solve_obs_index;
  array[M] int<lower=1, upper=N_segments> mouse_segment_start;
  array[M] int<lower=1, upper=N_segments> mouse_segment_end;

  vector<lower=0>[M] baseline_volume;
  vector<lower=0>[N_obs] volume;

  real<lower=0> rel_tol;
  real<lower=0> abs_tol;
  int<lower=1> max_num_steps;
}

parameters {
  real logV0_pop;
  real logr_pop;
  real logalpha_pop;
  real logK_pop;

  real<lower=0> sigma_logV0;
  real<lower=0> sigma_logr;
  real<lower=0> sigma_logalpha;
  real<lower=0> sigma_logK;

  vector[M] z_logV0;
  vector[M] z_logr;
  vector[M_rt] z_logalpha;
  vector[M] z_logK;

  real logdelta_pop;
  real logbeta_a_pop;
  real logbeta_q_pop;
  real logk_rep_pop;

  real<lower=0> sigma_logdelta;
  real<lower=0> sigma_logbeta_a;
  real<lower=0> sigma_logbeta_q;
  real<lower=0> sigma_logk_rep;

  vector[M_rt] z_logdelta;
  vector[M_rt] z_logbeta_a;
  vector[M_rt] z_logbeta_q;
  vector[M_rt] z_logk_rep;

  real<lower=0> sigma;
}

transformed parameters {
  real V0_pop = exp(logV0_pop);
  real r_pop = exp(logr_pop);
  real alpha_pop = exp(logalpha_pop);
  real delta_pop = exp(logdelta_pop);
  real beta_a_pop = exp(logbeta_a_pop);
  real beta_q_pop = exp(logbeta_q_pop);
  real k_rep_pop = exp(logk_rep_pop);

  vector[M] V0_mouse;
  vector[M] r_mouse;
  vector[M_rt] alpha_mouse;
  vector[M] K_mouse;
  vector[M_rt] delta_mouse;
  vector[M_rt] beta_a_mouse;
  vector[M_rt] beta_q_mouse;
  vector[M_rt] k_rep_mouse;
  vector[N_obs] mu_out;
  vector[N_obs] mu_log_out;

  for (m in 1:M) {
    V0_mouse[m] = exp(logV0_pop + sigma_logV0 * z_logV0[m]);
    r_mouse[m] = exp(logr_pop + sigma_logr * z_logr[m]);
    K_mouse[m] =
  exp(
    logK_pop +
    sigma_logK * z_logK[m]
  );
  }

  for (j in 1:M_rt) {
    alpha_mouse[j] = exp(logalpha_pop + sigma_logalpha * z_logalpha[j]);
    delta_mouse[j] = exp(logdelta_pop + sigma_logdelta * z_logdelta[j]);
    beta_a_mouse[j] = exp(logbeta_a_pop + sigma_logbeta_a * z_logbeta_a[j]);
    beta_q_mouse[j] = exp(logbeta_q_pop + sigma_logbeta_q * z_logbeta_q[j]);
    k_rep_mouse[j] = exp(logk_rep_pop + sigma_logk_rep * z_logk_rep[j]);
  }

  for (m in 1:M) {
    vector[3] y_current;
    real D_current = 0;
    vector[6] theta;

    y_current[1] = V0_mouse[m];
    y_current[2] = 0;
    y_current[3] = 0;

    theta[1] = r_mouse[m];
    if (rt_index[m] > 0) {
      theta[2] = delta_mouse[rt_index[m]];
      theta[3] = beta_a_mouse[rt_index[m]];
      theta[4] = beta_q_mouse[rt_index[m]];
      theta[5] = k_rep_mouse[rt_index[m]];
    } else {
      theta[2] = delta_pop;
      theta[3] = beta_a_pop;
      theta[4] = beta_q_pop;
      theta[5] = k_rep_pop;
    }
    theta[6] = K_mouse[m];

    for (s in mouse_segment_start[m]:mouse_segment_end[m]) {
      real interval_start = segment_start_time[s];
      if (rt_index[m] > 0)
        D_current += alpha_mouse[rt_index[m]] * segment_rt_dose[s];

      int n_times = solve_end[s] - solve_start[s] + 1;
      array[n_times] real solve_times;
      array[n_times] vector[3] y_solution;

      for (j in 1:n_times)
        solve_times[j] = solve_time[solve_start[s] + j - 1];

      y_solution = ode_rk45_tol(
        ", function_name, ",
        y_current,
        interval_start,
        solve_times,
        rel_tol,
        abs_tol,
        max_num_steps,
        theta,
        D_current,
        interval_start
      );

      for (j in 1:n_times) {
        int g = solve_start[s] + j - 1;
        int obs = solve_obs_index[g];
        if (obs > 0) {
          real V_pred = fmax(
            y_solution[j][1] + y_solution[j][2] + y_solution[j][3],
            1e-10
          );
          mu_out[obs] = V_pred;
          mu_log_out[obs] = log(V_pred);
        }
      }

      y_current = y_solution[n_times];
      D_current *= exp(-theta[5] * (solve_times[n_times] - interval_start));
    }
  }
}

model {
  logV0_pop ~ normal(log(30), 0.5);
  logr_pop ~ normal(log(", growth_prior, "), 0.5);
  logK_pop ~ normal(log(500), 0.7);

  sigma_logV0 ~ normal(0, 0.5);
  sigma_logr ~ normal(0, 0.5);
  sigma_logK ~ normal(0, 0.5);
  z_logV0 ~ std_normal();
  z_logr ~ std_normal();
  z_logK ~ std_normal();

  logalpha_pop ~ normal(log(0.9), 0.5);
  sigma_logalpha ~ normal(0, 0.5);
  z_logalpha ~ std_normal();

  logdelta_pop ~ normal(log(1.2), 0.5);
  logbeta_a_pop ~ normal(log(26.4), 0.5);
  logbeta_q_pop ~ normal(log(4.8), 0.5);
  logk_rep_pop ~ normal(log(86.4), 0.5);

  sigma_logdelta ~ normal(0, 0.5);
  sigma_logbeta_a ~ normal(0, 0.5);
  sigma_logbeta_q ~ normal(0, 0.5);
  sigma_logk_rep ~ normal(0, 0.5);

  z_logdelta ~ std_normal();
  z_logbeta_a ~ std_normal();
  z_logbeta_q ~ std_normal();
  z_logk_rep ~ std_normal();

  sigma ~ normal(0, 0.5);

  log(baseline_volume) ~ normal(log(V0_mouse), sigma);
  log(volume) ~ normal(mu_log_out, sigma);
}

generated quantities {
  vector[M + N_obs] log_lik;
  vector[M] y_rep_baseline;
  vector[N_obs] y_rep;

  for (m in 1:M) {
    log_lik[m] = lognormal_lpdf(
  baseline_volume[m] |
  log(V0_mouse[m]),
  sigma
);
    y_rep_baseline[m] = lognormal_rng(log(V0_mouse[m]), sigma);
  }
  for (n in 1:N_obs) {
    log_lik[M + n] = lognormal_lpdf(
  volume[n] |
  mu_log_out[n],
  sigma
);
    y_rep[n] = lognormal_rng(mu_log_out[n], sigma);
  }
}")
}

stan_files <- c(
  exponential = file.path(stan_dir, "exponential_rt_MFULL_gompertz_build.stan"),
  logistic = file.path(stan_dir, "logistic_rt_MFULL_gompertz_build.stan"),
  gompertz = file.path(stan_dir, "gompertz_rt_MFULL.stan")
)

writeLines(make_rt_exponential_stan_model(), stan_files[["exponential"]])
writeLines(make_rt_bounded_stan_model("logistic"), stan_files[["logistic"]])
writeLines(make_rt_bounded_stan_model("gompertz"), stan_files[["gompertz"]])


mod_rt_gompertz <- cmdstan_model(stan_files[["gompertz"]])


base_init <- function(include_K = FALSE, growth_prior = 0.3) {
  
  x <- list(
    logV0_pop = log(30),
    sigma_logV0 = 0.15,
    z_logV0 = rep(0, M),
    
    logr_pop = log(growth_prior),
    sigma_logr = 0.15,
    z_logr = rep(0, M),
    
    logalpha_pop = log(0.9),
    sigma_logalpha = 0.15,
    z_logalpha = rep(0, M_rt),
    
    logdelta_pop = log(1.2),
    logbeta_a_pop = log(26.4),
    logbeta_q_pop = log(4.8),
    logk_rep_pop = log(86.4),
    
    sigma_logdelta = 0.15,
    sigma_logbeta_a = 0.15,
    sigma_logbeta_q = 0.15,
    sigma_logk_rep = 0.15,
    
    z_logdelta = rep(0, M_rt),
    z_logbeta_a = rep(0, M_rt),
    z_logbeta_q = rep(0, M_rt),
    z_logk_rep = rep(0, M_rt),
    
    sigma = 0.3
  )
  
  if (include_K) {
    x$logK_pop <- log(500)
    x$sigma_logK <- 0.15
    x$z_logK <- rep(0, M)
  }
  
  x
}

rt_exp_init <- function() base_init(FALSE, 0.3)
rt_logistic_init <- function() base_init(TRUE, 0.4)
rt_gompertz_init <- function() base_init(TRUE, 0.2)


screen_fit <- function(fit, treedepth_limit = 12) {
  
  sd <- fit$sampler_diagnostics(format = "draws_df")
  ds <- fit$diagnostic_summary()
  
  out <- data.frame(
    sampling_draws = nrow(sd),
    divergences = sum(sd$divergent__),
    max_treedepth_hits = sum(sd$treedepth__ >= treedepth_limit),
    max_treedepth_percent = 100 * mean(sd$treedepth__ >= treedepth_limit),
    ebfmi = ds$ebfmi
  )
  
  out$passes <-
    out$divergences == 0 &&
    out$max_treedepth_hits == 0 &&
    is.finite(out$ebfmi) &&
    out$ebfmi > 0.3
  
  out
}

run_short_validation <- function(model, stan_data, init_fun, basename) {
  
  system.time({
    fit <- model$sample(
      data = stan_data,
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 50,
      iter_sampling = 50,
      adapt_delta = 0.95,
      max_treedepth = 12,
      seed = 123,
      init = init_fun,
      output_dir = csv_dir,
      output_basename = paste0(basename, "_short_validation"),
      refresh = 10
    )
  }) -> timing
  
  list(
    fit = fit,
    timing = timing,
    screen = screen_fit(fit, 12)
  )
}

run_intermediate_validation <- function(model, stan_data, init_fun, basename) {
  
  system.time({
    fit <- model$sample(
      data = stan_data,
      chains = 1,
      parallel_chains = 1,
      iter_warmup = 200,
      iter_sampling = 100,
      adapt_delta = 0.95,
      max_treedepth = 12,
      seed = 123,
      init = init_fun,
      output_dir = csv_dir,
      output_basename = paste0(basename, "_intermediate_validation"),
      refresh = 25
    )
  }) -> timing
  
  list(
    fit = fit,
    timing = timing,
    screen = screen_fit(fit, 12)
  )
}

##################################################
# full fit: GOMPERTZ
##################################################

system.time({
  fit_rt_gompertz <- mod_rt_gompertz$sample(
    data = rt_gompertz_stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 1000,
    adapt_delta = 0.99,
    max_treedepth = 12,
    seed = 123,
    init = rt_gompertz_init,
    output_dir = csv_dir,
    output_basename = "gompertz_MFULL_ctrl_rt_full",
    refresh = 50
  )
}) -> full_timing

print(full_timing)
print(fit_rt_gompertz$diagnostic_summary())

fit_rt_gompertz$save_object(
  file.path(results_dir, "fit_rt_gompertz_MFULL_ctrl_rt_full.rds")
)

cat("\nMFULL GOMPERTZ CONTROL + RT FULL FIT COMPLETE\n")

quit(save = "no", status = 0)


load_fit_if_available <- function(existing_name, pattern) {
  if (exists(existing_name, inherits = TRUE)) {
    return(get(existing_name, inherits = TRUE))
  }
  
  csv <- list.files(
    csv_dir,
    pattern = pattern,
    full.names = TRUE
  )
  
  if (length(csv) == 4) {
    return(cmdstanr::read_cmdstan_csv(csv))
  }
  
  NULL
}

fit_rt_exp <- load_fit_if_available(
  "fit_rt_exp",
  "^exponential_full.*\\.csv$"
)
fit_rt_logistic <- load_fit_if_available(
  "fit_rt_logistic",
  "^logistic_full.*\\.csv$"
)
fit_rt_gompertz <- load_fit_if_available(
  "fit_rt_gompertz",
  "^gompertz_full.*\\.csv$"
)


get_draws <- function(fit, variables = NULL) {
  
  if (inherits(fit, "CmdStanMCMC")) {
    return(fit$draws(variables = variables))
  }
  
  d <- fit$post_warmup_draws
  
  if (!is.null(variables)) {
    all_vars <- dimnames(d)$variable
    
    keep <- unique(unlist(lapply(
      variables,
      function(v) {
        which(
          all_vars == v |
            startsWith(all_vars, paste0(v, "["))
        )
      }
    )))
    
    d <- d[, , keep, drop = FALSE]
  }
  
  d
}


get_summary <- function(fit, variables = NULL) {
  
  if (inherits(fit, "CmdStanMCMC")) {
    return(fit$summary(variables = variables))
  }
  
  posterior::summarise_draws(
    get_draws(fit, variables),
    mean,
    median,
    sd,
    mad,
    ~quantile(.x, 0.05),
    ~quantile(.x, 0.95),
    posterior::rhat,
    posterior::ess_bulk,
    posterior::ess_tail
  )
}


have_all_full_fits <-
  exists("fit_rt_exp") &&
  !is.null(fit_rt_exp) &&
  exists("fit_rt_logistic") &&
  !is.null(fit_rt_logistic) &&
  exists("fit_rt_gompertz") &&
  !is.null(fit_rt_gompertz)

#  Population and mouse-specific summaries


population_vars_common <- c(
  "V0_pop", "sigma_logV0",
  "r_pop", "sigma_logr",
  "alpha_pop", "sigma_logalpha",
  "delta_pop",
  "beta_a_pop",
  "beta_q_pop",
  "k_rep_pop",
  "sigma"
)

write_model_summaries <- function(fit, model_name, include_K = FALSE) {
  
  vars <- population_vars_common
  if (include_K) {
    vars <- c(
      vars,
      "K_pop",
      "logK_pop",
      "sigma_logK"
    )
  }
  
  pop_summary <- get_summary(fit, vars)
  write.csv(
    pop_summary,
    file.path(results_dir, paste0(model_name, "_population_parameter_summary.csv")),
    row.names = FALSE
  )
  
  mouse_vars <- c(
    "V0_mouse", "r_mouse", "alpha_mouse"
  )
  if (include_K) mouse_vars <- c(mouse_vars, "K_mouse")
  
  mouse_summary <- get_summary(fit, mouse_vars)
  write.csv(
    mouse_summary,
    file.path(results_dir, paste0(model_name, "_mouse_parameter_summary.csv")),
    row.names = FALSE
  )
  
  invisible(list(population = pop_summary, mouse = mouse_summary))
}

if (have_all_full_fits && RUN_DOWNSTREAM_ANALYSIS) {
  
  exp_summaries <- write_model_summaries(fit_rt_exp, "exponential", FALSE)
  logistic_summaries <- write_model_summaries(fit_rt_logistic, "logistic", TRUE)
  gompertz_summaries <- write_model_summaries(fit_rt_gompertz, "gompertz", TRUE)
}


#  Posterior predictive checks


make_ppc <- function(fit, model_name, n_draws = 50) {
  
  yrep <- as_draws_matrix(get_draws(fit,"y_rep"))
  set.seed(123)
  ids <- sample(seq_len(nrow(yrep)), min(n_draws, nrow(yrep)))
  
  p_raw <- ppc_dens_overlay(
    y = rt_obs_postbaseline$volume,
    yrep = yrep[ids, , drop = FALSE]
  ) +
    ggtitle(paste("Posterior predictive check:", model_name))
  
  p_log <- ppc_dens_overlay(
    y = log(rt_obs_postbaseline$volume),
    yrep = log(yrep[ids, , drop = FALSE])
  ) +
    ggtitle(paste("Posterior predictive check (log scale):", model_name))
  
  ggsave(
    file.path(figures_dir, paste0(model_name, "_ppc_raw.pdf")),
    p_raw,
    width = 7,
    height = 5
  )
  
  ggsave(
    file.path(figures_dir, paste0(model_name, "_ppc_log.pdf")),
    p_log,
    width = 7,
    height = 5
  )
  
  invisible(list(raw = p_raw, log = p_log))
}

if (have_all_full_fits && RUN_DOWNSTREAM_ANALYSIS) {
  ppc_exp <- make_ppc(fit_rt_exp, "exponential")
  ppc_logistic <- make_ppc(fit_rt_logistic, "logistic")
  ppc_gompertz <- make_ppc(fit_rt_gompertz, "gompertz")
}

#  PSIS-LOO model comparison

if (have_all_full_fits && RUN_DOWNSTREAM_ANALYSIS) {
  
  loo_exp <- loo(get_draws(fit_rt_exp, "log_lik"))
  loo_logistic <- loo(get_draws(fit_rt_logistic, "log_lik"))
  loo_gompertz <- loo(get_draws(fit_rt_gompertz, "log_lik"))
  
  loo_comparison <- loo_compare(
    list(
      exponential = loo_exp,
      logistic = loo_logistic,
      gompertz = loo_gompertz
    )
  )
  
  print(loo_comparison)
  
  write.csv(
    as.data.frame(loo_comparison),
    file.path(results_dir, "psis_loo_model_comparison.csv")
  )
  
  pareto_summary <- data.frame(
    model = c("Exponential", "Logistic", "Gompertz"),
    n_k_gt_0_7 = c(
      sum(pareto_k_values(loo_exp) > 0.7),
      sum(pareto_k_values(loo_logistic) > 0.7),
      sum(pareto_k_values(loo_gompertz) > 0.7)
    ),
    n_k_ge_1 = c(
      sum(pareto_k_values(loo_exp) >= 1),
      sum(pareto_k_values(loo_logistic) >= 1),
      sum(pareto_k_values(loo_gompertz) >= 1)
    ),
    max_pareto_k = c(
      max(pareto_k_values(loo_exp)),
      max(pareto_k_values(loo_logistic)),
      max(pareto_k_values(loo_gompertz))
    )
  )
  
  print(pareto_summary)
  write.csv(
    pareto_summary,
    file.path(results_dir, "psis_loo_pareto_k_summary.csv"),
    row.names = FALSE
  )
}

#  Posterior correlation / identifiability check


make_population_correlation <- function(fit, model_name) {
  
  mat <- as_draws_matrix(
    get_draws(
      fit,
      c(
        "alpha_pop",
        "beta_a_pop",
        "beta_q_pop",
        "k_rep_pop",
        "delta_pop",
        "r_pop"
      )
    )
  )
  
  cor_mat <- cor(mat)
  
  write.csv(
    cor_mat,
    file.path(results_dir, paste0(model_name, "_population_posterior_correlations.csv"))
  )
  
  cor_mat
}

if (have_all_full_fits && RUN_DOWNSTREAM_ANALYSIS) {
  cor_exp <- make_population_correlation(fit_rt_exp, "exponential")
  cor_logistic <- make_population_correlation(fit_rt_logistic, "logistic")
  cor_gompertz <- make_population_correlation(fit_rt_gompertz, "gompertz")
}


#  Observed versus posterior-median predictions

make_prediction_table <- function(fit, model_name) {
  
  mu_draws <- as_draws_matrix(get_draws(fit, "mu_out"))
  
  q <- apply(
    mu_draws,
    2,
    quantile,
    probs = c(0.05, 0.5, 0.95)
  )
  
  out <- rt_obs_postbaseline %>%
    select(mouse, mouse_id, time_days, volume, obs_index) %>%
    mutate(
      pred_q05 = q[1, ],
      pred_median = q[2, ],
      pred_q95 = q[3, ]
    )
  
  write.csv(
    out,
    file.path(results_dir, paste0(model_name, "_observation_predictions.csv")),
    row.names = FALSE
  )
  
  p <- ggplot(out, aes(x = time_days)) +
    geom_ribbon(
      aes(ymin = pred_q05, ymax = pred_q95),
      alpha = 0.2
    ) +
    geom_line(aes(y = pred_median)) +
    geom_point(aes(y = volume), size = 0.8) +
    facet_wrap(~ mouse, scales = "free_y") +
    labs(
      x = "Time (days)",
      y = "Tumour volume",
      title = paste("Observed and fitted RT-only trajectories:", model_name)
    ) +
    theme_bw()
  
  ggsave(
    file.path(figures_dir, paste0(model_name, "_trajectories_by_mouse.pdf")),
    p,
    width = 12,
    height = 10
  )
  
  invisible(out)
}

if (have_all_full_fits && RUN_DOWNSTREAM_ANALYSIS) {
  pred_exp <- make_prediction_table(fit_rt_exp, "exponential")
  pred_logistic <- make_prediction_table(fit_rt_logistic, "logistic")
  pred_gompertz <- make_prediction_table(fit_rt_gompertz, "gompertz")
}

#  Controlled RT-schedule simulations

simulate_rt_schedule <- function(
    growth_model = c("exponential", "logistic", "gompertz"),
    times,
    rt_times,
    rt_doses,
    V0,
    r,
    alpha,
    delta,
    beta_a,
    beta_q,
    k_rep,
    K = NULL
) {
  
  growth_model <- match.arg(growth_model)
  
  state <- c(P = V0, Q = 0, A = 0)
  D_current <- 0
  current_time <- min(times)
  
  out <- data.frame(
    time = current_time,
    P = state["P"],
    Q = state["Q"],
    A = state["A"],
    D = D_current
  )
  
  event_times <- sort(unique(c(times[times > current_time], rt_times)))
  
  rhs <- function(t, y, parms) {
    
    P <- y[1]
    Q <- y[2]
    A <- y[3]
    D <- parms$D_start * exp(-k_rep * (t - parms$t_start))
    V <- max(P + Q + A, 1e-10)
    
    growth <- switch(
      growth_model,
      exponential = r * P,
      logistic = r * P * (1 - V / K),
      gompertz = r * P * log(K / V)
    )
    
    list(c(
      growth - beta_a * D * P - beta_q * D * P,
      beta_q * D * P,
      beta_a * D * P - delta * A
    ))
  }
  
  for (next_time in event_times) {
    
    # Apply RT exactly at the start of an interval.
    dose_now <- sum(rt_doses[rt_times == current_time])
    if (length(dose_now) == 0 || is.na(dose_now)) dose_now <- 0
    D_current <- D_current + alpha * dose_now
    
    if (next_time > current_time) {
      
      sol <- deSolve::ode(
        y = state,
        times = c(current_time, next_time),
        func = rhs,
        parms = list(D_start = D_current, t_start = current_time),
        method = "lsoda"
      )
      
      state <- as.numeric(sol[nrow(sol), c("P", "Q", "A")])
      names(state) <- c("P", "Q", "A")
      
      D_current <- D_current * exp(-k_rep * (next_time - current_time))
      current_time <- next_time
    }
    
    if (current_time %in% times) {
      out <- rbind(
        out,
        data.frame(
          time = current_time,
          P = state["P"],
          Q = state["Q"],
          A = state["A"],
          D = D_current
        )
      )
    }
  }
  
  out <- out %>%
    distinct(time, .keep_all = TRUE) %>%
    arrange(time) %>%
    mutate(total_volume = P + Q + A)
  
  out
}

extract_population_medians <- function(fit, include_K = FALSE) {
  
  vars <- c(
    "r_pop", "alpha_pop", "delta_pop",
    "beta_a_pop", "beta_q_pop", "k_rep_pop"
  )
  
  d <- as_draws_matrix(get_draws(fit, vars))
  
  out <- apply(d, 2, median)
  
  if (include_K) {
    
    kdraw <- as_draws_matrix(
      get_draws(
        fit,
        "logK_pop"
      )
    )[, 1]
    
    out <- c(
      out,
      K = median(
        exp(kdraw)
      )
    )
  }
  
  out
}

if (have_all_full_fits && RUN_DOWNSTREAM_ANALYSIS) {
  
  schedule_list <- list(
    `1Gy_x5` = data.frame(time = 1:5 + 1/12, dose = 1),
    `2Gy_x5` = data.frame(time = 1:5 + 1/12, dose = 2),
    `4Gy_x5` = data.frame(time = 1:5 + 1/12, dose = 4),
    `1Gy_x10` = data.frame(
      time = c(1:5, 8:12) + 1/12,
      dose = 1
    )
  )
  
  simulation_times <- seq(0, 25, by = 0.25)
  V0_typical <- median(rt_baseline$V_t0)
  
  simulation_results <- list()
  
  fits_for_sim <- list(
    exponential = fit_rt_exp,
    logistic = fit_rt_logistic,
    gompertz = fit_rt_gompertz
  )
  
  for (model_name in names(fits_for_sim)) {
    
    include_K <- model_name != "exponential"
    pars <- extract_population_medians(fits_for_sim[[model_name]], include_K)
    
    for (schedule_name in names(schedule_list)) {
      
      sch <- schedule_list[[schedule_name]]
      
      sim <- simulate_rt_schedule(
        growth_model = model_name,
        times = simulation_times,
        rt_times = sch$time,
        rt_doses = sch$dose,
        V0 = V0_typical,
        r = pars["r_pop"],
        alpha = pars["alpha_pop"],
        delta = pars["delta_pop"],
        beta_a = pars["beta_a_pop"],
        beta_q = pars["beta_q_pop"],
        k_rep = pars["k_rep_pop"],
        K = if (include_K) pars["K"] else NULL
      ) %>%
        mutate(
          model = model_name,
          schedule = schedule_name
        )
      
      simulation_results[[paste(model_name, schedule_name, sep = "_")]] <- sim
    }
  }
  
  controlled_simulations <- bind_rows(simulation_results)
  
  write.csv(
    controlled_simulations,
    file.path(results_dir, "controlled_rt_schedule_simulations.csv"),
    row.names = FALSE
  )
  
  p_sim <- ggplot(
    controlled_simulations,
    aes(x = time, y = total_volume, linetype = schedule)
  ) +
    geom_line() +
    facet_wrap(~ model, scales = "free_y") +
    labs(
      x = "Time (days)",
      y = "Simulated tumour volume",
      title = "Controlled RT schedule simulations",
      linetype = "RT schedule"
    ) +
    theme_bw()
  
  ggsave(
    file.path(figures_dir, "controlled_rt_schedule_simulations.pdf"),
    p_sim,
    width = 9,
    height = 6
  )
}

#  Residual diagnostics

if (have_all_full_fits && RUN_DOWNSTREAM_ANALYSIS) {
  
  make_residual_diagnostics <- function(pred_table, model_name) {
    
    residual_data <- pred_table %>%
      mutate(
        log_observed = log(volume),
        log_predicted = log(pred_median),
        residual = log_observed - log_predicted,
        model = model_name
      )
    
    write.csv(
      residual_data,
      file.path(
        results_dir,
        paste0(model_name, "_residual_diagnostics.csv")
      ),
      row.names = FALSE
    )
    
    residual_data
  }
  
  
  resid_exp <- make_residual_diagnostics(
    pred_exp,
    "Exponential"
  )
  
  resid_logistic <- make_residual_diagnostics(
    pred_logistic,
    "Logistic"
  )
  
  resid_gompertz <- make_residual_diagnostics(
    pred_gompertz,
    "Gompertz"
  )
  
  
  all_residuals <- bind_rows(
    resid_exp,
    resid_logistic,
    resid_gompertz
  )
  
  
  # Residuals against time
  p_resid_time <- ggplot(
    all_residuals,
    aes(x = time_days, y = residual)
  ) +
    geom_hline(
      yintercept = 0,
      linetype = 2
    ) +
    geom_point(
      alpha = 0.55,
      size = 1
    ) +
    geom_smooth(
      method = "loess",
      se = FALSE
    ) +
    facet_wrap(~ model) +
    labs(
      x = "Time (days)",
      y = "Log-scale residual",
      title = "Residuals versus time"
    ) +
    theme_bw()
  
  ggsave(
    file.path(
      figures_dir,
      "residuals_vs_time_all_models.pdf"
    ),
    p_resid_time,
    width = 9,
    height = 5
  )
  
  
  # Residual distributions
  p_resid_density <- ggplot(
    all_residuals,
    aes(x = residual)
  ) +
    geom_density() +
    geom_vline(
      xintercept = 0,
      linetype = 2
    ) +
    facet_wrap(~ model) +
    labs(
      x = "Log-scale residual",
      y = "Density",
      title = "Residual distributions"
    ) +
    theme_bw()
  
  ggsave(
    file.path(
      figures_dir,
      "residual_distributions_all_models.pdf"
    ),
    p_resid_density,
    width = 9,
    height = 5
  )
  
  
  # Numerical residual summaries
  residual_summary <- all_residuals %>%
    group_by(model) %>%
    summarise(
      mean_residual = mean(residual),
      sd_residual = sd(residual),
      median_absolute_residual = median(abs(residual)),
      rmse_log = sqrt(mean(residual^2)),
      .groups = "drop"
    )
  
  print(residual_summary)
  
  write.csv(
    residual_summary,
    file.path(
      results_dir,
      "residual_summary_all_models.csv"
    ),
    row.names = FALSE
  )
  

  #  Master model comparison and ELPD plot
    loo_table <- data.frame(
    model = c("Exponential", "Logistic", "Gompertz"),
    elpd_loo = c(
      loo_exp$estimates["elpd_loo", "Estimate"],
      loo_logistic$estimates["elpd_loo", "Estimate"],
      loo_gompertz$estimates["elpd_loo", "Estimate"]
    ),
    se_elpd_loo = c(
      loo_exp$estimates["elpd_loo", "SE"],
      loo_logistic$estimates["elpd_loo", "SE"],
      loo_gompertz$estimates["elpd_loo", "SE"]
    ),
    looic = c(
      loo_exp$estimates["looic", "Estimate"],
      loo_logistic$estimates["looic", "Estimate"],
      loo_gompertz$estimates["looic", "Estimate"]
    ),
    max_pareto_k = c(
      max(pareto_k_values(loo_exp)),
      max(pareto_k_values(loo_logistic)),
      max(pareto_k_values(loo_gompertz))
    ),
    n_pareto_k_gt_0_7 = c(
      sum(pareto_k_values(loo_exp) > 0.7),
      sum(pareto_k_values(loo_logistic) > 0.7),
      sum(pareto_k_values(loo_gompertz) > 0.7)
    )
  )
  
  master_model_comparison <- loo_table %>%
    left_join(
      residual_summary %>%
        select(model, rmse_log, mean_residual),
      by = "model"
    )
  
  print(master_model_comparison)
  
  write.csv(
    master_model_comparison,
    file.path(results_dir, "master_model_comparison.csv"),
    row.names = FALSE
  )
  

  # ELPD differences relative to best model
  elpd_df <- as.data.frame(loo_comparison)
  elpd_df$model <- rownames(elpd_df)
  rownames(elpd_df) <- NULL
  
  elpd_df$model <- factor(
    elpd_df$model,
    levels = rev(elpd_df$model)
  )
  
  p_elpd <- ggplot(
    elpd_df,
    aes(x = elpd_diff, y = model)
  ) +
    geom_vline(
      xintercept = 0,
      linetype = 2
    ) +
    geom_point(size = 2.5) +
    geom_errorbarh(
      aes(
        xmin = elpd_diff - se_diff,
        xmax = elpd_diff + se_diff
      ),
      height = 0.15
    ) +
    labs(
      x = "ELPD difference relative to best model",
      y = NULL,
      title = "PSIS-LOO model comparison"
    ) +
    theme_bw()
  
  ggsave(
    file.path(figures_dir, "psis_loo_elpd_comparison.pdf"),
    p_elpd,
    width = 7,
    height = 4
  )
  
  #  Mechanistic parameter comparison

  make_parameter_interval_table <- function(fit, model_name) {
    
    vars <- c(
      "alpha_pop",
      "beta_a_pop",
      "beta_q_pop",
      "k_rep_pop",
      "delta_pop"
    )
    
    d <- as_draws_matrix(get_draws(fit, vars))
    
    out <- data.frame(
      parameter = colnames(d),
      median = apply(d, 2, median),
      q05 = apply(d, 2, quantile, probs = 0.05),
      q95 = apply(d, 2, quantile, probs = 0.95),
      model = model_name
    )
    
    out
  }
  
  
  param_exp <- make_parameter_interval_table(
    fit_rt_exp,
    "Exponential"
  )
  
  param_logistic <- make_parameter_interval_table(
    fit_rt_logistic,
    "Logistic"
  )
  
  param_gompertz <- make_parameter_interval_table(
    fit_rt_gompertz,
    "Gompertz"
  )
  
  
  parameter_intervals <- bind_rows(
    param_exp,
    param_logistic,
    param_gompertz
  )
  
  write.csv(
    parameter_intervals,
    file.path(
      results_dir,
      "mechanistic_parameter_intervals_all_models.csv"
    ),
    row.names = FALSE
  )
  
  
  p_parameters <- ggplot(
    parameter_intervals,
    aes(
      x = median,
      y = model
    )
  ) +
    geom_point(size = 2) +
    geom_errorbar(
      aes(
        xmin = q05,
        xmax = q95
      ),
      orientation = "y",
      width = 0.15
    ) +
    facet_wrap(
      ~ parameter,
      scales = "free_x"
    ) +
    labs(
      x = "Posterior median and 90% credible interval",
      y = NULL,
      title = "Mechanistic parameter estimates across growth models"
    ) +
    theme_bw()
  
  ggsave(
    file.path(
      figures_dir,
      "mechanistic_parameter_intervals_all_models.pdf"
    ),
    p_parameters,
    width = 10,
    height = 6
  )
  
  print(parameter_intervals)
  

  #  Mechanistic compartment and DNA-damage plots

  compartment_long <- controlled_simulations %>%
    select(
      time,
      P,
      Q,
      A,
      D,
      model,
      schedule
    ) %>%
    tidyr::pivot_longer(
      cols = c(P, Q, A),
      names_to = "compartment",
      values_to = "volume"
    )
  
  # Tumour compartments
  p_compartments <- ggplot(
    compartment_long,
    aes(
      x = time,
      y = volume,
      linetype = compartment
    )
  ) +
    geom_line() +
    facet_grid(
      model ~ schedule,
      scales = "free_y"
    ) +
    labs(
      x = "Time (days)",
      y = "Simulated compartment volume",
      title = "Mechanistic tumour compartments under RT schedules",
      linetype = "Compartment"
    ) +
    theme_bw()
  
  ggsave(
    file.path(
      figures_dir,
      "controlled_rt_compartment_simulations.pdf"
    ),
    p_compartments,
    width = 12,
    height = 8
  )
  
  

  #  DNA-damage trajectories

  
  damage_results <- list()
  
  for (model_name in names(fits_for_sim)) {
    
    pars <- extract_population_medians(
      fits_for_sim[[model_name]],
      include_K = model_name != "exponential"
    )
    
    alpha_now <- pars["alpha_pop"]
    krep_now <- pars["k_rep_pop"]
    
    for (schedule_name in names(schedule_list)) {
      
      sch <- schedule_list[[schedule_name]]
      
      t_dense <- sort(unique(c(
        seq(0, 25, by = 0.005),
        sch$time
      )))
      
      D_dense <- sapply(
        t_dense,
        function(t) {
          
          past_events <- sch$time <= t
          
          if (!any(past_events)) {
            return(0)
          }
          
          sum(
            alpha_now *
              sch$dose[past_events] *
              exp(
                -krep_now *
                  (t - sch$time[past_events])
              )
          )
        }
      )
      
      damage_results[[
        paste(model_name, schedule_name, sep = "_")
      ]] <- data.frame(
        time = t_dense,
        D = D_dense,
        model = model_name,
        schedule = schedule_name
      )
    }
  }
  
  damage_simulations <- bind_rows(damage_results)
  
  write.csv(
    damage_simulations,
    file.path(
      results_dir,
      "controlled_rt_DNA_damage_high_resolution.csv"
    ),
    row.names = FALSE
  )
  
  p_damage <- ggplot(
    damage_simulations,
    aes(
      x = time,
      y = D,
      linetype = schedule
    )
  ) +
    geom_line() +
    facet_wrap(~ model, scales = "free_y") +
    labs(
      x = "Time (days)",
      y = "DNA damage state D",
      title = "DNA-damage dynamics under RT schedules",
      linetype = "RT schedule"
    ) +
    theme_bw()
  
  ggsave(
    file.path(
      figures_dir,
      "controlled_rt_DNA_damage_simulations.pdf"
    ),
    p_damage,
    width = 9,
    height = 6
  )
  
  #  Quantitative controlled-schedule comparison

  
  schedule_summary <- controlled_simulations %>%
    group_by(model, schedule) %>%
    summarise(
      final_volume = total_volume[which.max(time)],
      minimum_volume = min(total_volume),
      time_of_minimum = time[which.min(total_volume)],
      auc_volume = sum(
        diff(time) *
          (head(total_volume, -1) + tail(total_volume, -1)) / 2
      ),
      .groups = "drop"
    )
  
  schedule_summary <- schedule_summary %>%
    group_by(model) %>%
    mutate(
      reference_final_volume = max(final_volume),
      final_volume_reduction_percent =
        100 * (
          1 - final_volume / reference_final_volume
        )
    ) %>%
    ungroup()
  
  print(schedule_summary)
  
  write.csv(
    schedule_summary,
    file.path(
      results_dir,
      "controlled_rt_schedule_quantitative_summary.csv"
    ),
    row.names = FALSE
  )
  
    p_final_schedule <- ggplot(
    schedule_summary,
    aes(
      x = schedule,
      y = final_volume
    )
  ) +
    geom_col() +
    facet_wrap(
      ~ model,
      scales = "free_y"
    ) +
    labs(
      x = "RT schedule",
      y = "Tumour volume at day 25",
      title = "Final tumour volume under controlled RT schedules"
    ) +
    theme_bw()
  
  ggsave(
    file.path(
      figures_dir,
      "controlled_rt_schedule_final_volume.pdf"
    ),
    p_final_schedule,
    width = 9,
    height = 5
  )
  
 