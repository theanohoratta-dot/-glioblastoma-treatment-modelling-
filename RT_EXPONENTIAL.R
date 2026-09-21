##################################################
#Rt-stage: CONTROL + RT-ONLY MICE — EXPONENTIAL
##################################################

.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))

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

data <- read.csv(
  file.path(project_dir, "data", "lemassonParpi.csv")
)

results_dir <- file.path("results", "stage2_mfull_exp_ctrl_rt")
figures_dir <- "figuresstage3"
stan_dir <- "stan"
csv_dir <- file.path(results_dir, "cmdstan_csv")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stan_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

RUN_INTERMEDIATE_DIAGNOSTICS <- FALSE
RUN_FULL_FITS <- TRUE
RUN_DOWNSTREAM_ANALYSIS <- TRUE


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
cat("All post-baseline observations are mapped exactly once.\n")

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
cat("Solver segmentation checks passed.\n")

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

rt_max_volume_mouse <- rt_obs %>%
  group_by(mouse_id) %>%
  summarise(max_volume = max(volume), .groups = "drop") %>%
  arrange(mouse_id)

stopifnot(
  nrow(rt_max_volume_mouse) == M,
  all(rt_max_volume_mouse$mouse_id == seq_len(M)),
  all(is.finite(rt_max_volume_mouse$max_volume)),
  all(rt_max_volume_mouse$max_volume > 0)
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
  list(max_volume_mouse = rt_max_volume_mouse$max_volume),
  ode_settings_full
)
rt_logistic_stan_data_validation <- c(
  rt_common_stan_data,
  list(max_volume_mouse = rt_max_volume_mouse$max_volume),
  ode_settings_validation
)

rt_gompertz_stan_data <- c(
  rt_common_stan_data,
  list(max_volume_mouse = rt_max_volume_mouse$max_volume),
  ode_settings_full
)
rt_gompertz_stan_data_validation <- c(
  rt_common_stan_data,
  list(max_volume_mouse = rt_max_volume_mouse$max_volume),
  ode_settings_validation
)

make_rt_stan_model <- function(growth_model = c("exponential")) {
  
  growth_model <- match.arg(growth_model)
  
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
  real logdelta_pop;
  real logalpha_pop;
  real logbeta_a_pop;
  real logbeta_q_pop;
  real logk_rep_pop;

  real<lower=0> sigma_logV0;
  real<lower=0> sigma_logr;
  real<lower=0> sigma_logdelta;
  real<lower=0> sigma_logalpha;
  real<lower=0> sigma_logbeta_a;
  real<lower=0> sigma_logbeta_q;
  real<lower=0> sigma_logk_rep;

  vector[M] z_logV0;
  vector[M] z_logr;
  vector[M_rt] z_logdelta;
  vector[M_rt] z_logalpha;
  vector[M_rt] z_logbeta_a;
  vector[M_rt] z_logbeta_q;
  vector[M_rt] z_logk_rep;

  real<lower=0> sigma;
}

transformed parameters {

  real V0_pop = exp(logV0_pop);
  real r_pop = exp(logr_pop);
  real delta_pop = exp(logdelta_pop);
  real alpha_pop = exp(logalpha_pop);
  real beta_a_pop = exp(logbeta_a_pop);
  real beta_q_pop = exp(logbeta_q_pop);
  real k_rep_pop = exp(logk_rep_pop);

  vector[M] V0_mouse;
  vector[M] r_mouse;
  vector[M_rt] delta_mouse;
  vector[M_rt] alpha_mouse;
  vector[M_rt] beta_a_mouse;
  vector[M_rt] beta_q_mouse;
  vector[M_rt] k_rep_mouse;

  vector[N_obs] mu_out;
  vector[N_obs] mu_log_out;

  for (m in 1:M) {
    V0_mouse[m] =
      exp(logV0_pop + sigma_logV0 * z_logV0[m]);

    r_mouse[m] =
      exp(logr_pop + sigma_logr * z_logr[m]);

  }

  for (j in 1:M_rt) {
    delta_mouse[j] =
      exp(logdelta_pop + sigma_logdelta * z_logdelta[j]);

    alpha_mouse[j] =
      exp(logalpha_pop + sigma_logalpha * z_logalpha[j]);

    beta_a_mouse[j] =
      exp(logbeta_a_pop + sigma_logbeta_a * z_logbeta_a[j]);

    beta_q_mouse[j] =
      exp(logbeta_q_pop + sigma_logbeta_q * z_logbeta_q[j]);

    k_rep_mouse[j] =
      exp(logk_rep_pop + sigma_logk_rep * z_logk_rep[j]);
  }

  for (m in 1:M) {
    vector[2] y_current;

    real P_current;
    real D_current;
    vector[5] theta;

    P_current = V0_mouse[m];

    y_current[1] = 0;
    y_current[2] = V0_mouse[m];

    D_current = 0;

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

    for (s in mouse_segment_start[m]:mouse_segment_end[m]) {

      real interval_start = segment_start_time[s];
      
      if (rt_index[m] > 0) {
        D_current =
          D_current + alpha_mouse[rt_index[m]] * segment_rt_dose[s];
      }

      int n_times = solve_end[s] - solve_start[s] + 1;
      array[n_times] real solve_times;
      array[n_times] vector[2] y_solution;

      for (j in 1:n_times) {
        solve_times[j] = solve_time[solve_start[s] + j - 1];
      }

      y_solution =
        ode_rk45_tol(
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

        int global_solve_index = solve_start[s] + j - 1;
        int observation_index = solve_obs_index[global_solve_index];

        if (observation_index > 0) {

          real V_pred = fmax(y_solution[j][2], 1e-10);

          mu_out[observation_index] = V_pred;
          mu_log_out[observation_index] = log(V_pred);
        }
      }

      {
        real dt_end = solve_times[n_times] - interval_start;
        real repair_factor_end =
          exp(-theta[5] * dt_end);

        real log_P_ratio_end =
          r_mouse[m] * dt_end
          - (
              (theta[3] + theta[4])
              * D_current
              / theta[5]
            )
            * (1 - repair_factor_end);
        P_current =
          P_current * exp(log_P_ratio_end);
        y_current = y_solution[n_times];
        D_current =
          D_current * repair_factor_end;
      }
    }
  }
}

model {

  // Same priors as the previous all-mouse-specific diagnostic.

  logV0_pop ~ normal(log(30), 0.5);
  logr_pop ~ normal(log(0.3), 0.5);

  // hour-to-day conversions:
  // delta: 5/100 per hour -> 1.2/day
  // alpha: 0.9 
  // beta_a: 1.1/hour -> 26.4/day
  // beta_q: 0.2/hour -> 4.8/day
  // k_rep: 3.6/hour -> 86.4/day
  logdelta_pop ~ normal(log(1.2), 0.5);
  logalpha_pop ~ normal(log(0.9), 0.5);
  logbeta_a_pop ~ normal(log(26.4), 0.5);
  logbeta_q_pop ~ normal(log(4.8), 0.5);
  logk_rep_pop ~ normal(log(86.4), 0.5);

  sigma_logV0 ~ normal(0, 0.5);
  sigma_logr ~ normal(0, 0.5);
  sigma_logdelta ~ normal(0, 0.5);
  sigma_logalpha ~ normal(0, 0.5);
  sigma_logbeta_a ~ normal(0, 0.5);
  sigma_logbeta_q ~ normal(0, 0.5);
  sigma_logk_rep ~ normal(0, 0.5);

  z_logV0 ~ std_normal();
  z_logr ~ std_normal();
  z_logdelta ~ std_normal();
  z_logalpha ~ std_normal();
  z_logbeta_a ~ std_normal();
  z_logbeta_q ~ std_normal();
  z_logk_rep ~ std_normal();

  sigma ~ normal(0, 0.5);

  log(baseline_volume) ~ normal(log(V0_mouse), sigma);
  log(volume) ~ normal(mu_log_out, sigma);
}

generated quantities {

  vector[M + N_obs] log_lik;
  vector[N_obs] y_rep;
  vector[M] y_rep_baseline;

  for (m in 1:M) {
    log_lik[m] =
      normal_lpdf(log(baseline_volume[m]) | log(V0_mouse[m]), sigma);

    y_rep_baseline[m] =
      lognormal_rng(log(V0_mouse[m]), sigma);
  }

  for (n in 1:N_obs) {
    log_lik[M + n] =
      normal_lpdf(log(volume[n]) | mu_log_out[n], sigma);

    y_rep[n] =
      lognormal_rng(mu_log_out[n], sigma);
  }
}
")
}

stan_files <- c(
  exponential = file.path(stan_dir, "stage2_mfull_exp_ctrl_rt.stan")
)

for (nm in names(stan_files)) {
  writeLines(make_rt_stan_model(nm), stan_files[[nm]])
}

##################################################
# 6. Syntax checks and compilation
##################################################

for (nm in names(stan_files)) {
  cat("\nChecking Stan syntax:", nm, "\n")
  cmdstan_model(stan_files[[nm]], compile = FALSE)$check_syntax()
}

mod_rt_exp <- cmdstan_model(stan_files[["exponential"]])

base_init <- function(include_K = FALSE, growth_prior = 0.3) {
  
  x <- list(
    logV0_pop = log(30),
    logr_pop = log(growth_prior),
    logdelta_pop = log(1.2),
    logalpha_pop = log(0.9),
    logbeta_a_pop = log(26.4),
    logbeta_q_pop = log(4.8),
    logk_rep_pop = log(86.4),
    
    sigma_logV0 = 0.15,
    sigma_logr = 0.15,
    sigma_logdelta = 0.15,
    sigma_logalpha = 0.15,
    sigma_logbeta_a = 0.15,
    sigma_logbeta_q = 0.15,
    sigma_logk_rep = 0.15,
    
    z_logV0 = rep(0, M),
    z_logr = rep(0, M),
    z_logdelta = rep(0, M_rt),
    z_logalpha = rep(0, M_rt),
    z_logbeta_a = rep(0, M_rt),
    z_logbeta_q = rep(0, M_rt),
    z_logk_rep = rep(0, M_rt),
    
    sigma = 0.2
  )
  
  if (include_K) {
    x$logK_extra_pop <- log(100)
    x$sigma_logK <- 0.15
    x$z_logK <- rep(0, M)
  }
  
  x
}

rt_exp_init <- function() base_init(FALSE, 0.3)
rt_logistic_init <- function() base_init(TRUE, 0.4)
rt_gompertz_init <- function() base_init(TRUE, 0.2)

fixed_timing <- system.time({
  exp_fixed <- mod_rt_exp$sample(
    data = rt_exp_stan_data_validation,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 0,
    iter_sampling = 1,
    fixed_param = TRUE,
    seed = 123,
    init = rt_exp_init,
    output_dir = csv_dir,
    output_basename = "stage2_mfull_exp_ctrl_rt_fixed_parameter_test",
    refresh = 1
  )
})

print(fixed_timing)

mu_fixed <- exp_fixed$draws(
  variables = "mu_out",
  format = "matrix"
)

mu_fixed <- as.numeric(mu_fixed[1, ])

fixed_check <- data.frame(
  n_predictions = length(mu_fixed),
  all_finite = all(is.finite(mu_fixed)),
  all_positive = all(mu_fixed > 0),
  min_prediction = min(mu_fixed),
  max_prediction = max(mu_fixed)
)

print(fixed_check)

stopifnot(
  fixed_check$all_finite,
  fixed_check$all_positive
)

write.csv(
  fixed_check,
  file.path(results_dir, "exponential_fixed_parameter_test.csv"),
  row.names = FALSE
)

cat("\n===== 56-MOUSE CONTROL + RT EXPONENTIAL MFULL HMC FIT =====\n")
cat("Mice in test:", M, "\n")
cat("Post-baseline observations:", N_obs, "\n")
cat("ODE segments:", N_segments, "\n")
cat("Solver output times:", N_solve, "\n")

hmc_timing <- system.time({
  
  exp_hmc_ctrl_rt <- mod_rt_exp$sample(
    data = rt_exp_stan_data,
    
    chains = 4,
    parallel_chains = 4,
    
    iter_warmup = 500,
    iter_sampling = 1000,
    
    adapt_delta = 0.99,
    max_treedepth = 12,
    
    seed = 123,
    init = rt_exp_init,
    
    output_dir = csv_dir,
    output_basename = "stage2_mfull_exp_ctrl_rt_RK45_FULL",
    
    refresh = 10
  )
})

print(hmc_timing)

saveRDS(
  exp_hmc_ctrl_rt$output_files(),
  file.path(results_dir, "cmdstan_output_files.rds")
)

hmc_diag <- exp_hmc_ctrl_rt$diagnostic_summary()
print(hmc_diag)
write.csv(
  as.data.frame(hmc_diag),
  file.path(results_dir, "hmc_diagnostic_summary.csv"),
  row.names = FALSE
)

summary_vars <- c(
  "logV0_pop",
  "logr_pop",
  "logdelta_pop",
  "logalpha_pop",
  "logbeta_a_pop",
  "logbeta_q_pop",
  "logk_rep_pop",
  "sigma_logV0",
  "sigma_logr",
  "sigma_logdelta",
  "sigma_logalpha",
  "sigma_logbeta_a",
  "sigma_logbeta_q",
  "sigma_logk_rep",
  "sigma"
)

param_summary <- exp_hmc_ctrl_rt$summary(variables = summary_vars)
print(param_summary)
write.csv(
  as.data.frame(param_summary),
  file.path(results_dir, "population_parameter_summary.csv"),
  row.names = FALSE
)

finite_rhat <- param_summary$rhat[is.finite(param_summary$rhat)]
finite_ess_bulk <- param_summary$ess_bulk[is.finite(param_summary$ess_bulk)]
finite_ess_tail <- param_summary$ess_tail[is.finite(param_summary$ess_tail)]
cat("Max R-hat among reported population/scale parameters:",
    if (length(finite_rhat)) max(finite_rhat) else NA_real_, "\n")
cat("Min bulk ESS among reported population/scale parameters:",
    if (length(finite_ess_bulk)) min(finite_ess_bulk) else NA_real_, "\n")
cat("Min tail ESS among reported population/scale parameters:",
    if (length(finite_ess_tail)) min(finite_ess_tail) else NA_real_, "\n")

mechanistic_draws <- exp_hmc_ctrl_rt$draws(
  variables = c("alpha_mouse", "beta_a_mouse", "beta_q_mouse", "k_rep_mouse"),
  format = "matrix"
)

ridge_rows <- vector("list", M_rt)

for (m in seq_len(M_rt)) {
  alpha_m <- mechanistic_draws[, sprintf("alpha_mouse[%d]", m)]
  beta_a_m <- mechanistic_draws[, sprintf("beta_a_mouse[%d]", m)]
  beta_q_m <- mechanistic_draws[, sprintf("beta_q_mouse[%d]", m)]
  krep_m <- mechanistic_draws[, sprintf("k_rep_mouse[%d]", m)]
  
  log_alpha <- log(alpha_m)
  log_beta_a <- log(beta_a_m)
  log_beta_q <- log(beta_q_m)
  log_krep <- log(krep_m)
  log_beta_sum <- log(beta_a_m + beta_q_m)
  
  effective_strength <- alpha_m * (beta_a_m + beta_q_m) / krep_m
  
  ridge_rows[[m]] <- data.frame(
    mouse_id = m,
    mouse = rt_mouse_lookup$mouse[m],
    cor_logalpha_logbeta_a = as.numeric(cor(log_alpha, log_beta_a)),
    cor_logalpha_logbeta_q = as.numeric(cor(log_alpha, log_beta_q)),
    cor_logalpha_logbeta_sum = as.numeric(cor(log_alpha, log_beta_sum)),
    cor_logalpha_logkrep = as.numeric(cor(log_alpha, log_krep)),
    effective_strength_mean = mean(effective_strength),
    effective_strength_sd = sd(effective_strength),
    effective_strength_q05 = unname(quantile(effective_strength, 0.05)),
    effective_strength_q50 = unname(quantile(effective_strength, 0.50)),
    effective_strength_q95 = unname(quantile(effective_strength, 0.95))
  )
}

ridge_diagnostics <- bind_rows(ridge_rows)
print(ridge_diagnostics)
write.csv(
  ridge_diagnostics,
  file.path(results_dir, "mouse_specific_ridge_diagnostics.csv"),
  row.names = FALSE
)

trace_vars <- c(
  "logV0_pop", "logr_pop",
  "logdelta_pop", "logalpha_pop",
  "logbeta_a_pop", "logbeta_q_pop", "logk_rep_pop",
  "sigma_logV0", "sigma_logr",
  "sigma_logdelta", "sigma_logalpha",
  "sigma_logbeta_a", "sigma_logbeta_q", "sigma_logk_rep",
  "sigma"
)

trace_plot <- bayesplot::mcmc_trace(
  exp_hmc_ctrl_rt$draws(variables = trace_vars)
)
ggsave(
  file.path(figures_dir, "stage2_exp_ctrl_rt_trace_population_hyperparameters.png"),
  trace_plot, width = 14, height = 12, dpi = 300
)

set.seed(123)
posterior_hyper <- as_draws_matrix(
  exp_hmc_ctrl_rt$draws(variables = trace_vars)
)
n_pp <- nrow(posterior_hyper)

prior_spec <- list(
  logV0_pop      = function(n) rnorm(n, log(30), 0.5),
  logr_pop       = function(n) rnorm(n, log(0.3), 0.5),
  logdelta_pop   = function(n) rnorm(n, log(1.2), 0.5),
  logalpha_pop   = function(n) rnorm(n, log(0.9), 0.5),
  logbeta_a_pop  = function(n) rnorm(n, log(26.4), 0.5),
  logbeta_q_pop  = function(n) rnorm(n, log(4.8), 0.5),
  logk_rep_pop   = function(n) rnorm(n, log(86.4), 0.5),
  sigma_logV0    = function(n) abs(rnorm(n, 0, 0.5)),
  sigma_logr     = function(n) abs(rnorm(n, 0, 0.5)),
  sigma_logdelta = function(n) abs(rnorm(n, 0, 0.5)),
  sigma_logalpha = function(n) abs(rnorm(n, 0, 0.5)),
  sigma_logbeta_a= function(n) abs(rnorm(n, 0, 0.5)),
  sigma_logbeta_q= function(n) abs(rnorm(n, 0, 0.5)),
  sigma_logk_rep = function(n) abs(rnorm(n, 0, 0.5)),
  sigma           = function(n) abs(rnorm(n, 0, 0.5))
)

parameter_labels <- c(
  logV0_pop = "log population V0",
  logr_pop = "log population r",
  logdelta_pop = "log population delta",
  logalpha_pop = "log population alpha",
  logbeta_a_pop = "log population beta_a",
  logbeta_q_pop = "log population beta_q",
  logk_rep_pop = "log population k_rep",
  sigma_logV0 = "between-mouse SD: log V0",
  sigma_logr = "between-mouse SD: log r",
  sigma_logdelta = "between-RT-mouse SD: log delta",
  sigma_logalpha = "between-RT-mouse SD: log alpha",
  sigma_logbeta_a = "between-RT-mouse SD: log beta_a",
  sigma_logbeta_q = "between-RT-mouse SD: log beta_q",
  sigma_logk_rep = "between-RT-mouse SD: log k_rep",
  sigma = "residual log-volume SD"
)

prior_posterior_data <- bind_rows(lapply(names(prior_spec), function(v) {
  bind_rows(
    data.frame(
      parameter = parameter_labels[[v]],
      value = prior_spec[[v]](n_pp),
      distribution = "Prior"
    ),
    data.frame(
      parameter = parameter_labels[[v]],
      value = posterior_hyper[, v],
      distribution = "Posterior"
    )
  )
}))

prior_posterior_data$distribution <- factor(
  prior_posterior_data$distribution,
  levels = c("Prior", "Posterior")
)

prior_posterior_plot <- ggplot(
  prior_posterior_data,
  aes(x = value, colour = distribution, fill = distribution)
) +
  geom_density(alpha = 0.25, adjust = 1) +
  facet_wrap(~ parameter, scales = "free", ncol = 3) +
  labs(
    title = "Stage 2 exponential MFULL: prior and posterior distributions",
    x = "Parameter value",
    y = "Density",
    colour = "Distribution",
    fill = "Distribution"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(
  file.path(figures_dir, "stage2_exp_ctrl_rt_prior_vs_posterior_population_hyperparameters.png"),
  prior_posterior_plot, width = 15, height = 15, dpi = 300
)


sample_half_normal <- function(n, sd = 0.5) abs(rnorm(n, 0, sd))

make_induced_prior_posterior <- function(
    fit, stan_variable, mouse_ids,
    log_pop_mean, log_pop_sd = 0.5, sigma_sd = 0.5,
    seed = 123) {
  
  post <- as_draws_matrix(fit$draws(variables = stan_variable))
  stopifnot(ncol(post) == length(mouse_ids))
  n_draws <- nrow(post)
  
  set.seed(seed)
  log_pop_prior <- rnorm(n_draws, log_pop_mean, log_pop_sd)
  sigma_prior <- sample_half_normal(n_draws, sigma_sd)
  
  out <- vector("list", length(mouse_ids))
  for (j in seq_along(mouse_ids)) {
    z_prior <- rnorm(n_draws, 0, 1)
    induced_prior <- exp(log_pop_prior + sigma_prior * z_prior)
    
    out[[j]] <- bind_rows(
      data.frame(
        mouse = factor(mouse_ids[j], levels = mouse_ids),
        value = induced_prior,
        distribution = "Induced prior"
      ),
      data.frame(
        mouse = factor(mouse_ids[j], levels = mouse_ids),
        value = post[, j],
        distribution = "Posterior"
      )
    )
  }
  
  bind_rows(out)
}

plot_induced_prior_posterior <- function(dat, title_text, file_name, ncol = 6) {
  dat$distribution <- factor(
    dat$distribution,
    levels = c("Induced prior", "Posterior")
  )
  
  p <- ggplot(
    dat,
    aes(x = value, colour = distribution, fill = distribution)
  ) +
    geom_density(alpha = 0.20, adjust = 1) +
    facet_wrap(~ mouse, scales = "free", ncol = ncol) +
    labs(
      title = title_text,
      subtitle = "Prior is induced by the fitted non-centred hierarchical model",
      x = "Natural-scale parameter value",
      y = "Density",
      colour = "Distribution",
      fill = "Distribution"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(file.path(figures_dir, file_name), p, width = 16, height = 14, dpi = 300)
  invisible(p)
}

induced_V0 <- make_induced_prior_posterior(
  exp_hmc_ctrl_rt, "V0_mouse", mouse_lookup$mouse,
  log_pop_mean = log(30), seed = 123
)
plot_induced_prior_posterior(
  induced_V0,
  "Stage 2 exponential MFULL: mouse-specific V0 prior vs posterior",
  "stage2_exp_ctrl_rt_induced_prior_vs_posterior_V0_by_mouse.png"
)

induced_r <- make_induced_prior_posterior(
  exp_hmc_ctrl_rt, "r_mouse", mouse_lookup$mouse,
  log_pop_mean = log(0.3), seed = 124
)
plot_induced_prior_posterior(
  induced_r,
  "Stage 2 exponential MFULL: mouse-specific r prior vs posterior",
  "stage2_exp_ctrl_rt_induced_prior_vs_posterior_r_by_mouse.png"
)

rt_induced_specs <- list(
  delta = list(stan = "delta_mouse", centre = log(1.2), seed = 125),
  alpha = list(stan = "alpha_mouse", centre = log(0.9), seed = 126),
  beta_a = list(stan = "beta_a_mouse", centre = log(26.4), seed = 127),
  beta_q = list(stan = "beta_q_mouse", centre = log(4.8), seed = 128),
  k_rep = list(stan = "k_rep_mouse", centre = log(86.4), seed = 129)
)

for (nm in names(rt_induced_specs)) {
  spec <- rt_induced_specs[[nm]]
  dat_induced <- make_induced_prior_posterior(
    exp_hmc_ctrl_rt, spec$stan, rt_mouse_lookup$mouse,
    log_pop_mean = spec$centre, seed = spec$seed
  )
  plot_induced_prior_posterior(
    dat_induced,
    paste0("Stage 2 exponential MFULL: mouse-specific ", nm, " prior vs posterior"),
    paste0("stage2_exp_ctrl_rt_induced_prior_vs_posterior_", nm, "_by_RT_mouse.png")
  )
}

individual_growth_plot <- bayesplot::mcmc_intervals(
  exp_hmc_ctrl_rt$draws(variables = c("V0_mouse", "r_mouse")),
  prob = 0.5,
  prob_outer = 0.95
)
ggsave(
  file.path(figures_dir, "stage2_exp_ctrl_rt_individual_V0_r_posteriors.png"),
  individual_growth_plot, width = 12, height = 16, dpi = 300
)

individual_rt_plot <- bayesplot::mcmc_intervals(
  exp_hmc_ctrl_rt$draws(
    variables = c(
      "delta_mouse", "alpha_mouse",
      "beta_a_mouse", "beta_q_mouse", "k_rep_mouse"
    )
  ),
  prob = 0.5,
  prob_outer = 0.95
)
ggsave(
  file.path(figures_dir, "stage2_exp_ctrl_rt_individual_RT_parameter_posteriors.png"),
  individual_rt_plot, width = 14, height = 24, dpi = 300
)

yrep_post <- as_draws_matrix(
  exp_hmc_ctrl_rt$draws(variables = "y_rep")
)
yrep_base <- as_draws_matrix(
  exp_hmc_ctrl_rt$draws(variables = "y_rep_baseline")
)

observed_all <- c(rt_baseline$V_t0, rt_obs_postbaseline$volume)
yrep_all <- cbind(yrep_base, yrep_post)

set.seed(123)
ppc_idx <- sample(seq_len(nrow(yrep_all)), min(100L, nrow(yrep_all)))

ppc_density <- bayesplot::ppc_dens_overlay(
  y = observed_all,
  yrep = yrep_all[ppc_idx, , drop = FALSE]
)
ggsave(
  file.path(figures_dir, "stage2_exp_ctrl_rt_ppc_density_overlay.png"),
  ppc_density, width = 10, height = 7, dpi = 300
)

mu_post <- as_draws_matrix(
  exp_hmc_ctrl_rt$draws(variables = "mu_out")
)
mu_summary <- data.frame(
  obs_index = seq_len(N_obs),
  fitted_median = apply(mu_post, 2, median),
  fitted_q025 = apply(mu_post, 2, quantile, probs = 0.025),
  fitted_q975 = apply(mu_post, 2, quantile, probs = 0.975)
)

trajectory_data <- rt_obs_postbaseline %>%
  left_join(mu_summary, by = "obs_index") %>%
  left_join(
    mouse_lookup %>% select(mouse, is_rt),
    by = "mouse"
  ) %>%
  mutate(group = ifelse(is_rt, "RT-only", "Control"))

baseline_plot_data <- rt_baseline %>%
  left_join(
    mouse_lookup %>% select(mouse, is_rt),
    by = "mouse"
  ) %>%
  mutate(group = ifelse(is_rt, "RT-only", "Control"))

trajectory_plot <- ggplot() +
  geom_ribbon(
    data = trajectory_data,
    aes(
      x = time_days,
      ymin = fitted_q025,
      ymax = fitted_q975,
      group = mouse
    ),
    alpha = 0.20
  ) +
  geom_line(
    data = trajectory_data,
    aes(x = time_days, y = fitted_median, group = mouse),
    linewidth = 0.5
  ) +
  geom_point(
    data = trajectory_data,
    aes(x = time_days, y = volume),
    size = 0.7
  ) +
  geom_point(
    data = baseline_plot_data,
    aes(x = t0, y = V_t0),
    size = 0.7
  ) +
  facet_wrap(~ group + mouse, scales = "free_y") +
  labs(
    title = "Stage 2 exponential MFULL: fitted trajectories",
    x = "Time (days)",
    y = expression("Tumour volume (mm"^3*")")
  ) +
  theme_minimal()

ggsave(
  file.path(figures_dir, "stage2_exp_ctrl_rt_fitted_trajectories_all_mice.png"),
  trajectory_plot, width = 18, height = 24, dpi = 300
)

post_mean_mu_log <- colMeans(
  as_draws_matrix(exp_hmc_ctrl_rt$draws(variables = "mu_log_out"))
)
sigma_mean <- mean(
  as_draws_matrix(exp_hmc_ctrl_rt$draws(variables = "sigma"))[, "sigma"]
)

residual_data <- rt_obs_postbaseline %>%
  mutate(
    fitted_log_volume = post_mean_mu_log,
    log_residual = log(volume) - fitted_log_volume,
    standardised_log_residual = log_residual / sigma_mean
  ) %>%
  left_join(
    mouse_lookup %>% select(mouse, is_rt),
    by = "mouse"
  ) %>%
  mutate(group = ifelse(is_rt, "RT-only", "Control"))

residual_time_plot <- ggplot(
  residual_data,
  aes(x = time_days, y = standardised_log_residual)
) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ group) +
  labs(
    title = "Stage 2 exponential MFULL: standardised log residuals over time",
    x = "Time (days)",
    y = "Standardised log residual"
  ) +
  theme_minimal()

ggsave(
  file.path(figures_dir, "stage2_exp_ctrl_rt_standardised_residuals_over_time.png"),
  residual_time_plot, width = 10, height = 7, dpi = 300
)

residual_hist <- ggplot(
  residual_data,
  aes(x = standardised_log_residual)
) +
  geom_histogram(bins = 30) +
  labs(
    title = "Stage 2 exponential MFULL: standardised residual distribution",
    x = "Standardised log residual",
    y = "Count"
  ) +
  theme_minimal()

ggsave(
  file.path(figures_dir, "stage2_exp_ctrl_rt_standardised_residual_histogram.png"),
  residual_hist, width = 8, height = 6, dpi = 300
)

residual_qq <- ggplot(
  residual_data,
  aes(sample = standardised_log_residual)
) +
  stat_qq() +
  stat_qq_line() +
  labs(
    title = "Stage 2 exponential MFULL: normal Q-Q plot of standardised residuals",
    x = "Theoretical quantiles",
    y = "Observed quantiles"
  ) +
  theme_minimal()

ggsave(
  file.path(figures_dir, "stage2_exp_ctrl_rt_standardised_residual_qq.png"),
  residual_qq, width = 8, height = 6, dpi = 300
)

write.csv(
  residual_data,
  file.path(results_dir, "residual_diagnostics_postbaseline.csv"),
  row.names = FALSE
)

log_lik_draws <- as_draws_matrix(
  exp_hmc_ctrl_rt$draws(variables = "log_lik")
)
loo_exp_ctrl_rt <- loo::loo(log_lik_draws)

capture.output(
  print(loo_exp_ctrl_rt),
  file = file.path(results_dir, "loo_summary.txt")
)

pareto_k <- loo::pareto_k_values(loo_exp_ctrl_rt)
pareto_table <- data.frame(
  observation_index = seq_along(pareto_k),
  pareto_k = as.numeric(pareto_k)
)

write.csv(
  pareto_table,
  file.path(results_dir, "pareto_k_by_observation.csv"),
  row.names = FALSE
)

full_summary <- exp_hmc_ctrl_rt$summary()
full_summary_no_z <- full_summary %>%
  filter(!grepl("^z_", variable))

write.csv(
  as.data.frame(full_summary_no_z),
  file.path(results_dir, "mcmc_summary_all_non_z_parameters.csv"),
  row.names = FALSE
)

mcmc_overview <- data.frame(
  maximum_Rhat = max(full_summary_no_z$rhat, na.rm = TRUE),
  minimum_bulk_ESS = min(full_summary_no_z$ess_bulk, na.rm = TRUE),
  minimum_tail_ESS = min(full_summary_no_z$ess_tail, na.rm = TRUE)
)

write.csv(
  mcmc_overview,
  file.path(results_dir, "mcmc_convergence_overview.csv"),
  row.names = FALSE
)

cat("\nStage-1-style diagnostic figures saved to:\n",
    normalizePath(figures_dir, mustWork = FALSE), "\n")
