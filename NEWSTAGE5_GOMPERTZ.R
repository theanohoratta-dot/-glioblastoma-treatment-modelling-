############################################################
# NEW STAGE-5 TREATMENT-SCHEDULE OPTIMISATION — GOMPERTZ
############################################################

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(deSolve)
})

set.seed(123)


N_POSTERIOR_DRAWS <- 500L
FOLLOWUP_DAY <- 15
V0_REFERENCE <- 30
OUTPUT_DT <- 0.05  # days = 1.2 hours
TIMING_OFFSETS_HOURS <- c(-12, -8, -6, -4, -2, -1, 0, 1, 2, 4, 6, 8, 12)
ODE_RTOL <- 1e-8
ODE_ATOL <- 1e-10


fit_candidates <- c(
  "fit_stage4NEW_gomp_reducesum.rds",
  "results/stage4new_combination_mfull_gompertz_reducesum/fit_stage4NEW_gomp_reducesum.rds",
  "~/Scratch/dissertation/results/stage4new_combination_mfull_gompertz_reducesum/fit_stage4NEW_gomp_reducesum.rds",
  "~/Desktop/STADIO4/FINAL_FITS/stage4new_combination_mfull_gompertz_reducesum/fit_stage4NEW_gomp_reducesum.rds"
)

fit_path <- fit_candidates[file.exists(path.expand(fit_candidates))][1]
if (is.na(fit_path)) {
  stop(
    "Could not find the FINAL NEW mouse-specific-PK Stage-4 Gompertz fit. Tried:\n",
    paste(fit_candidates, collapse = "\n")
  )
}
fit_path <- path.expand(fit_path)


data_candidates <- c(
  "~/Desktop/STADIO4/lemassonParpi.csv",
  "data/lemassonParpi.csv",
  "~/Scratch/dissertation/data/lemassonParpi.csv",
  "lemassonParpi.csv",
  "lemassonParpi (11).csv",
  "/mnt/data/lemassonParpi (11).csv"
)

data_path <- data_candidates[file.exists(path.expand(data_candidates))][1]
if (is.na(data_path)) {
  stop(
    "Could not find lemassonParpi.csv. Tried:\n",
    paste(data_candidates, collapse = "\n")
  )
}
data_path <- path.expand(data_path)

results_dir <- "results/NEWSTAGE5_gompertz"
figures_dir <- file.path(results_dir, "figures")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

cat("Using posterior fit:\n  ", fit_path, "\n", sep = "")
cat("Using dataset:\n  ", data_path, "\n\n", sep = "")

############################################################
#  READ DATA AND VERIFY THE OBSERVED COMBINATION SCHEDULE
############################################################

data <- read.csv(data_path)

required_cols <- c(
  "mouse", "time", "volume", "observation",
  "dose", "modality", "study"
)
stopifnot(all(required_cols %in% names(data)))

data <- data %>%
  mutate(time_days = time / 24)

mouse_groups <- data %>%
  group_by(mouse) %>%
  summarise(
    has_drug = any(modality == 1, na.rm = TRUE),
    has_rt   = any(modality == 2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    group = case_when(
      !has_drug & !has_rt ~ "control",
       has_drug & !has_rt ~ "drug_only",
      !has_drug &  has_rt ~ "rt_only",
       has_drug &  has_rt ~ "combination"
    )
  ) %>%
  arrange(mouse) %>%
  mutate(mouse_id = row_number())

combo_mice <- mouse_groups %>%
  filter(group == "combination") %>%
  pull(mouse)

stopifnot(length(combo_mice) == 9)

# Check that all nine combination mice have exactly the same treatment schedule.
combo_schedule_by_mouse <- data %>%
  filter(mouse %in% combo_mice, modality %in% c(1, 2)) %>%
  arrange(mouse, time, modality) %>%
  select(mouse, time, time_days, modality, dose)

reference_schedule_raw <- combo_schedule_by_mouse %>%
  filter(mouse == combo_mice[1]) %>%
  select(-mouse)

for (m in combo_mice[-1]) {
  tmp <- combo_schedule_by_mouse %>%
    filter(mouse == m) %>%
    select(-mouse)

  if (!isTRUE(all.equal(tmp, reference_schedule_raw, check.attributes = FALSE))) {
    stop("Combination treatment schedules are not identical across all nine mice.")
  }
}

# Convert to the schedule format used by the forward simulator.
# modality: 1 = drug, 2 = RT
reference_schedule <- reference_schedule_raw %>%
  transmute(
    time_days = time_days,
    drug_dose = ifelse(modality == 1, dose, 0),
    rt_dose   = ifelse(modality == 2, dose, 0)
  ) %>%
  group_by(time_days) %>%
  summarise(
    drug_dose = sum(drug_dose),
    rt_dose   = sum(rt_dose),
    .groups = "drop"
  ) %>%
  arrange(time_days)

# Explicit checks for the known experimental combination regimen.
ref_drug <- reference_schedule_raw %>% filter(modality == 1)
ref_rt   <- reference_schedule_raw %>% filter(modality == 2)

stopifnot(
  nrow(ref_drug) == 20,
  nrow(ref_rt) == 10,
  all(ref_drug$dose == 25),
  all(ref_rt$dose == 1)
)

# First daily drug dose is one hour before RT and second dose six hours after RT.
rt_hours <- ref_rt$time
first_drug_hours <- ref_drug$time[seq(1, nrow(ref_drug), by = 2)]
second_drug_hours <- ref_drug$time[seq(2, nrow(ref_drug), by = 2)]
stopifnot(
  all(rt_hours - first_drug_hours == 1),
  all(second_drug_hours - rt_hours == 6)
)

cat("Verified observed combination regimen across all 9 mice:\n")
cat("  10 RT fractions of 1 Gy\n")
cat("  25 mg/kg drug 1 h before each RT fraction\n")
cat("  25 mg/kg second drug dose 6 h after each RT fraction\n\n")

write.csv(
  reference_schedule_raw,
  file.path(results_dir, "observed_combination_schedule.csv"),
  row.names = FALSE
)


############################################################
# LOAD POSTERIOR AND SELECT DRAWS
############################################################

fit <- readRDS(fit_path)

fit_model_params <- fit$metadata()$model_params
if ("logK_extra_pop" %in% fit_model_params) {
  stop("Loaded an obsolete Stage-4 fit containing logK_extra_pop. Use the FINAL analytic Stage-4 fit with logK_pop.")
}
if (!("logK_pop" %in% fit_model_params)) {
  stop("Loaded Stage-4 fit does not contain logK_pop; it is not compatible with this final Stage-5 script.")
}

posterior_variables <- c(
  "logr_pop",
  "sigma_logr",
  "logK_pop",
  "sigma_logK",
  "logdelta_pop",
  "sigma_logdelta",
  "logalpha_pop",
  "sigma_logalpha",
  "logbeta_a_pop",
  "sigma_logbeta_a",
  "logbeta_q_pop",
  "sigma_logbeta_q",
  "logk_rep_pop",
  "sigma_logk_rep",
  "logk_abs_pop",
  "sigma_logk_abs",
  "logk_cl_pop",
  "sigma_logk_cl",
  "logk_repi_pop",
  "sigma_logk_repi"
)

posterior_draws <- fit$draws(
  variables = posterior_variables,
  format = "draws_df"
) %>%
  posterior::as_draws_df() %>%
  as.data.frame()

missing_vars <- setdiff(posterior_variables, names(posterior_draws))
if (length(missing_vars) > 0) {
  stop("Posterior is missing variables: ", paste(missing_vars, collapse = ", "))
}

n_available <- nrow(posterior_draws)
N_POSTERIOR_DRAWS <- min(N_POSTERIOR_DRAWS, n_available)

# Deterministic, approximately even coverage across the full posterior draw set.
selected_rows <- unique(round(seq(1, n_available, length.out = N_POSTERIOR_DRAWS)))
posterior_draws <- posterior_draws[selected_rows, , drop = FALSE]
posterior_draws$posterior_draw_id <- seq_len(nrow(posterior_draws))

set.seed(12345)

posterior_draws$z_r_new       <- rnorm(nrow(posterior_draws))
posterior_draws$z_K_new <- rnorm(nrow(posterior_draws))
posterior_draws$z_delta_new   <- rnorm(nrow(posterior_draws))
posterior_draws$z_alpha_new   <- rnorm(nrow(posterior_draws))
posterior_draws$z_beta_a_new  <- rnorm(nrow(posterior_draws))
posterior_draws$z_beta_q_new  <- rnorm(nrow(posterior_draws))
posterior_draws$z_k_rep_new   <- rnorm(nrow(posterior_draws))
posterior_draws$z_k_abs_new   <- rnorm(nrow(posterior_draws))
posterior_draws$z_k_cl_new    <- rnorm(nrow(posterior_draws))
posterior_draws$z_k_repi_new  <- rnorm(nrow(posterior_draws))

posterior_parameters <- posterior_draws %>%
  transmute(
    posterior_draw_id,
    r       = exp(logr_pop       + sigma_logr       * z_r_new),
    K       = exp(logK_pop + sigma_logK * z_K_new),
    delta   = exp(logdelta_pop   + sigma_logdelta   * z_delta_new),
    alpha   = exp(logalpha_pop   + sigma_logalpha   * z_alpha_new),
    beta_a  = exp(logbeta_a_pop  + sigma_logbeta_a  * z_beta_a_new),
    beta_q  = exp(logbeta_q_pop  + sigma_logbeta_q  * z_beta_q_new),
    k_rep   = exp(logk_rep_pop   + sigma_logk_rep   * z_k_rep_new),
    k_abs   = exp(logk_abs_pop   + sigma_logk_abs   * z_k_abs_new),
    k_cl    = exp(logk_cl_pop    + sigma_logk_cl    * z_k_cl_new),
    k_repi  = exp(logk_repi_pop + sigma_logk_repi  * z_k_repi_new)
  )

stopifnot(
  all(is.finite(as.matrix(posterior_parameters[, -1]))),
  all(as.matrix(posterior_parameters[, -1]) > 0)
)

write.csv(
  posterior_parameters,
  file.path(results_dir, "new_mouse_posterior_parameter_draws.csv"),
  row.names = FALSE
)

cat("Posterior draws available: ", n_available, "\n", sep = "")
cat("Posterior draws used for simulation: ", nrow(posterior_parameters), "\n", sep = "")
cat("Stage-5 individual parameters generated from FINAL Stage-4 mouse-specific hierarchies.\n\n")



pk_c1_exact <- function(dt, C1_start, k_abs) {
  C1_start * exp(-k_abs * dt)
}

pk_c2_exact <- function(dt, C1_start, C2_start, k_abs, k_cl) {
  if (dt <= 0) return(C2_start)

  if (abs(k_abs - k_cl) <= 1e-8) {
    # Same limiting branch as the final Stage-4 Stan implementation.
    exp(-k_cl * dt) * (C2_start + k_abs * C1_start * dt)
  } else {
    C2_start * exp(-k_cl * dt) +
      (k_abs * C1_start / (k_abs - k_cl)) *
      (exp(-k_cl * dt) - exp(-k_abs * dt))
  }
}

pk_int_c2_exact <- function(dt, C1_start, C2_start, k_abs, k_cl) {
  if (dt <= 0) return(0)

  if (abs(k_abs - k_cl) <= 1e-8) {
    # Same limiting branch as the final Stage-4 Stan implementation.
    term0 <- C2_start * (-expm1(-k_cl * dt)) / k_cl
    term0 + C1_start *
      (1 - exp(-k_cl * dt) * (1 + k_cl * dt)) / k_cl
  } else {
    int_cl  <- (-expm1(-k_cl  * dt)) / k_cl
    int_abs <- (-expm1(-k_abs * dt)) / k_abs
    C2_start * int_cl +
      (k_abs * C1_start / (k_abs - k_cl)) * (int_cl - int_abs)
  }
}

damage_exact <- function(dt, D_start, C1_start, C2_start,
                         k_abs, k_cl, k_rep, k_repi) {
  if (dt <= 0 || D_start == 0) return(D_start)

  int_c2 <- pk_int_c2_exact(
    dt = dt,
    C1_start = C1_start,
    C2_start = C2_start,
    k_abs = k_abs,
    k_cl = k_cl
  )

  exponent <- -k_rep * dt + k_rep * k_repi * int_c2
  D_start * exp(exponent)
}

combo_reduced_rhs <- function(t, y, pars) {
  dt <- t - pars[["segment_t0"]]
  D <- damage_exact(
    dt = dt,
    D_start = pars[["D_start"]],
    C1_start = pars[["C1_start"]],
    C2_start = pars[["C2_start"]],
    k_abs = pars[["k_abs"]],
    k_cl = pars[["k_cl"]],
    k_rep = pars[["k_rep"]],
    k_repi = pars[["k_repi"]]
  )

  P <- y[["P"]]
  A <- y[["A"]]
  V <- y[["V"]]
  V_safe <- max(V, 1e-10)

  growth <- pars[["r"]] * P * log(pars[["K"]] / V_safe)
  dP <- growth - (pars[["beta_a"]] + pars[["beta_q"]]) * D * P
  dA <- pars[["beta_a"]] * D * P - pars[["delta"]] * A
  dV <- growth - pars[["delta"]] * A

  list(c(dP, dA, dV))
}


collapse_schedule <- function(schedule) {
  schedule %>%
    filter(drug_dose != 0 | rt_dose != 0) %>%
    group_by(time_days) %>%
    summarise(
      drug_dose = sum(drug_dose),
      rt_dose = sum(rt_dose),
      .groups = "drop"
    ) %>%
    arrange(time_days)
}

make_schedule <- function(drug_times_days, drug_doses,
                          rt_times_days = ref_rt$time_days,
                          rt_doses = ref_rt$dose) {

  stopifnot(length(drug_times_days) == length(drug_doses))
  stopifnot(length(rt_times_days) == length(rt_doses))

  d <- if (length(drug_times_days) > 0) {
    data.frame(
      time_days = drug_times_days,
      drug_dose = drug_doses,
      rt_dose = 0
    )
  } else {
    data.frame(time_days = numeric(0), drug_dose = numeric(0), rt_dose = numeric(0))
  }

  r <- if (length(rt_times_days) > 0) {
    data.frame(
      time_days = rt_times_days,
      drug_dose = 0,
      rt_dose = rt_doses
    )
  } else {
    data.frame(time_days = numeric(0), drug_dose = numeric(0), rt_dose = numeric(0))
  }

  collapse_schedule(bind_rows(d, r))
}

############################################################
# FORWARD SIMULATOR USING FINAL STAGE-4 ANALYTIC REDUCTION
############################################################

simulate_schedule <- function(pars,
                              schedule,
                              followup_day = FOLLOWUP_DAY,
                              V0 = V0_REFERENCE,
                              output_dt = OUTPUT_DT) {

  TIME_TOL <- 1e-10

  schedule <- collapse_schedule(schedule) %>%
    filter(time_days >= 0, time_days <= followup_day) %>%
    arrange(time_days)

  regular_grid <- seq(0, followup_day, by = output_dt)
  output_times <- sort(c(0, regular_grid, schedule$time_days, followup_day))
  if (length(output_times) > 1) {
    output_times <- output_times[c(TRUE, diff(output_times) > TIME_TOL)]
  }

  state <- c(C1 = 0, C2 = 0, D = 0, P = V0, A = 0, V = V0)

  current_time <- 0

  ev0 <- schedule %>% filter(abs(time_days - current_time) <= TIME_TOL)
  if (nrow(ev0) > 0) {
    state[["C1"]] <- state[["C1"]] + sum(ev0$drug_dose)
    state[["D"]]  <- state[["D"]] + pars[["alpha"]] * sum(ev0$rt_dose)
  }

  out <- data.frame(
    time_days = 0,
    C1 = state[["C1"]],
    C2 = state[["C2"]],
    D  = state[["D"]],
    P  = state[["P"]],
    A  = state[["A"]],
    V  = state[["V"]]
  )

  event_times <- sort(unique(schedule$time_days[schedule$time_days > TIME_TOL &
                                                  schedule$time_days < followup_day - TIME_TOL]))
  segment_ends <- sort(unique(c(event_times, followup_day)))

  for (segment_end in segment_ends) {
    if (segment_end <= current_time + TIME_TOL) next

    segment_start_state <- state

    requested <- output_times[
      output_times > current_time + TIME_TOL &
      output_times <= segment_end + TIME_TOL
    ]
    requested <- sort(unique(c(requested, segment_end)))

    segment_pars <- c(
      pars,
      list(
        segment_t0 = current_time,
        C1_start = unname(segment_start_state[["C1"]]),
        C2_start = unname(segment_start_state[["C2"]]),
        D_start = unname(segment_start_state[["D"]])
      )
    )

    y0 <- c(
      P = segment_start_state[["P"]],
      A = segment_start_state[["A"]],
      V = segment_start_state[["V"]]
    )

    sol <- deSolve::ode(
      y = y0,
      times = c(current_time, requested),
      func = combo_reduced_rhs,
      parms = segment_pars,
      method = "lsoda",
      rtol = ODE_RTOL,
      atol = ODE_ATOL
    )

    sol_df <- as.data.frame(sol)
    sol_df <- sol_df[-1, , drop = FALSE]

    if (nrow(sol_df) > 0) {
      analytic_rows <- lapply(sol_df$time, function(tt) {
        dt <- tt - current_time
        c(
          C1 = pk_c1_exact(dt, segment_start_state[["C1"]], pars[["k_abs"]]),
          C2 = pk_c2_exact(dt, segment_start_state[["C1"]], segment_start_state[["C2"]],
                           pars[["k_abs"]], pars[["k_cl"]]),
          D = damage_exact(dt, segment_start_state[["D"]],
                           segment_start_state[["C1"]], segment_start_state[["C2"]],
                           pars[["k_abs"]], pars[["k_cl"]],
                           pars[["k_rep"]], pars[["k_repi"]])
        )
      })
      analytic_mat <- do.call(rbind, analytic_rows)

      new_rows <- data.frame(
        time_days = sol_df$time,
        C1 = analytic_mat[, "C1"],
        C2 = analytic_mat[, "C2"],
        D = analytic_mat[, "D"],
        P = sol_df$P,
        A = sol_df$A,
        V = sol_df$V
      )

      if (any(!is.finite(as.matrix(new_rows[, -1])))) {
        stop("Non-finite state encountered before t = ", segment_end)
      }

      out <- bind_rows(out, new_rows)
      last <- new_rows[nrow(new_rows), ]
      state <- c(
        C1 = last$C1,
        C2 = last$C2,
        D  = last$D,
        P  = last$P,
        A  = last$A,
        V  = last$V
      )
    }

    current_time <- segment_end

    ev <- schedule %>% filter(abs(time_days - current_time) <= TIME_TOL)
    if (nrow(ev) > 0 && current_time < followup_day - TIME_TOL) {
      state[["C1"]] <- state[["C1"]] + sum(ev$drug_dose)
      state[["D"]]  <- state[["D"]] + pars[["alpha"]] * sum(ev$rt_dose)
    }
  }

  out %>%
    arrange(time_days) %>%
    distinct(time_days, .keep_all = TRUE)
}

trapz_auc <- function(time, y) {
  idx <- order(time)
  time <- time[idx]
  y <- y[idx]
  sum(diff(time) * (head(y, -1) + tail(y, -1)) / 2)
}

summarise_trajectory <- function(traj) {
  data.frame(
    final_volume = traj$V[which.max(traj$time_days)],
    tumour_auc = trapz_auc(traj$time_days, traj$V),
    min_volume = min(traj$V),
    max_volume = max(traj$V)
  )
}

############################################################
#  BUILD THE PRIMARY CANDIDATE SCHEDULES
############################################################

rt_times <- ref_rt$time_days
rt_doses <- ref_rt$dose
first_drug_times <- first_drug_hours / 24
second_drug_times <- second_drug_hours / 24

# A. Observed/reference regimen: 25 mg/kg 1 h before + 25 mg/kg 6 h after RT.
schedule_reference <- reference_schedule

# B. Equal total daily dose: 50 mg/kg QD vs 25+25 mg/kg BID.
schedule_qd50 <- make_schedule(
  drug_times_days = first_drug_times,
  drug_doses = rep(50, length(first_drug_times))
)

schedule_bid25x2 <- schedule_reference

# C. Equal dose per administration: 25 mg/kg QD vs 25+25 mg/kg BID.
schedule_qd25 <- make_schedule(
  drug_times_days = first_drug_times,
  drug_doses = rep(25, length(first_drug_times))
)

# D. One single drug administration over the entire RT course.
# Both 25 and 50 mg/kg are included so "single dose" is not silently
# tied to one dose-size interpretation.
schedule_single_course_25 <- make_schedule(
  drug_times_days = first_drug_times[1],
  drug_doses = 25
)

schedule_single_course_50 <- make_schedule(
  drug_times_days = first_drug_times[1],
  drug_doses = 50
)

primary_schedules <- list(
  "Observed BID: 25+25 mg/kg" = schedule_reference,
  "QD: 50 mg/kg" = schedule_qd50,
  "QD: 25 mg/kg" = schedule_qd25,
  "Single course dose: 25 mg/kg" = schedule_single_course_25,
  "Single course dose: 50 mg/kg" = schedule_single_course_50
)

# Write all primary schedules so the exact intervention definitions are auditable.
primary_schedule_table <- imap_dfr(
  primary_schedules,
  ~ mutate(.x, schedule = .y, .before = 1)
)
write.csv(
  primary_schedule_table,
  file.path(results_dir, "primary_candidate_schedules.csv"),
  row.names = FALSE
)

############################################################
#  TIMING-OFFSET SCHEDULES
#
# One 25 mg/kg drug dose is associated with each 1 Gy RT fraction.
# Only drug-RT timing changes; dose amount and number of doses are fixed.
############################################################

timing_schedules <- setNames(
  lapply(TIMING_OFFSETS_HOURS, function(offset_h) {
    make_schedule(
      drug_times_days = rt_times + offset_h / 24,
      drug_doses = rep(25, length(rt_times)),
      rt_times_days = rt_times,
      rt_doses = rt_doses
    )
  }),
  paste0("dt_", TIMING_OFFSETS_HOURS, "h")
)

# Ensure no timing candidate starts before t=0.
min_candidate_time <- min(unlist(lapply(timing_schedules, function(x) x$time_days)))
if (min_candidate_time < 0) {
  stop("At least one timing candidate contains a treatment before t=0.")
}


median_pars <- posterior_parameters %>%
  summarise(across(-posterior_draw_id, median)) %>%
  as.list()

validation_reference <- simulate_schedule(
  pars = median_pars,
  schedule = schedule_reference
)

stopifnot(
  nrow(validation_reference) > 1,
  all(is.finite(validation_reference$V)),
  validation_reference$time_days[1] == 0,
  abs(max(validation_reference$time_days) - FOLLOWUP_DAY) < 1e-10
)

write.csv(
  validation_reference,
  file.path(results_dir, "validation_reference_median_parameter_trajectory.csv"),
  row.names = FALSE
)

cat("Single-draw reference simulation passed.\n")
cat("Median-parameter reference V(T) = ",
    signif(tail(validation_reference$V, 1), 6), " mm^3\n\n", sep = "")


run_one_schedule <- function(schedule_name, schedule, parameter_table) {
  message("Running: ", schedule_name)

  map_dfr(seq_len(nrow(parameter_table)), function(i) {
    pars <- as.list(parameter_table[i, setdiff(names(parameter_table), "posterior_draw_id")])

    ans <- tryCatch({
      traj <- simulate_schedule(pars = pars, schedule = schedule)
      smry <- summarise_trajectory(traj)
      cbind(
        posterior_draw_id = parameter_table$posterior_draw_id[i],
        schedule = schedule_name,
        status = "ok",
        smry
      )
    }, error = function(e) {
      data.frame(
        posterior_draw_id = parameter_table$posterior_draw_id[i],
        schedule = schedule_name,
        status = paste0("ERROR: ", conditionMessage(e)),
        final_volume = NA_real_,
        tumour_auc = NA_real_,
        min_volume = NA_real_,
        max_volume = NA_real_
      )
    })

    ans
  })
}

primary_results <- imap_dfr(
  primary_schedules,
  ~ run_one_schedule(.y, .x, posterior_parameters)
)

write.csv(
  primary_results,
  file.path(results_dir, "primary_schedule_posterior_outcomes.csv"),
  row.names = FALSE
)

error_summary <- primary_results %>%
  count(schedule, status) %>%
  arrange(schedule, desc(n))
print(error_summary)

if (any(primary_results$status != "ok")) {
  warning("Some primary simulations failed. Inspect primary_schedule_posterior_outcomes.csv before interpretation.")
}


timing_results <- imap_dfr(timing_schedules, function(schedule, nm) {
  offset_h <- as.numeric(sub("^dt_", "", sub("h$", "", nm)))
  tmp <- run_one_schedule(
    schedule_name = paste0("Timing offset ", offset_h, " h"),
    schedule = schedule,
    parameter_table = posterior_parameters
  )
  tmp$offset_hours <- offset_h
  tmp
})

write.csv(
  timing_results,
  file.path(results_dir, "timing_sweep_posterior_outcomes.csv"),
  row.names = FALSE
)

if (any(timing_results$status != "ok")) {
  warning("Some timing simulations failed. Inspect timing_sweep_posterior_outcomes.csv before interpretation.")
}

############################################################
#  POSTERIOR SUMMARIES
############################################################

summarise_outcome <- function(df, grouping_cols) {
  df %>%
    filter(status == "ok") %>%
    group_by(across(all_of(grouping_cols))) %>%
    summarise(
      n_draws = n(),
      final_volume_median = median(final_volume),
      final_volume_q025 = quantile(final_volume, 0.025),
      final_volume_q975 = quantile(final_volume, 0.975),
      auc_median = median(tumour_auc),
      auc_q025 = quantile(tumour_auc, 0.025),
      auc_q975 = quantile(tumour_auc, 0.975),
      .groups = "drop"
    )
}

primary_summary <- summarise_outcome(primary_results, "schedule")
timing_summary <- summarise_outcome(timing_results, "offset_hours") %>%
  arrange(offset_hours)

write.csv(
  primary_summary,
  file.path(results_dir, "primary_schedule_summary.csv"),
  row.names = FALSE
)
write.csv(
  timing_summary,
  file.path(results_dir, "timing_sweep_summary.csv"),
  row.names = FALSE
)

print(primary_summary)
print(timing_summary)

############################################################
# PAIRED POSTERIOR COMPARISONS AGAINST OBSERVED BID REGIMEN
############################################################

reference_name <- "Observed BID: 25+25 mg/kg"
reference_outcomes <- primary_results %>%
  filter(schedule == reference_name, status == "ok") %>%
  select(
    posterior_draw_id,
    ref_final_volume = final_volume,
    ref_auc = tumour_auc
  )

paired_vs_reference <- primary_results %>%
  filter(schedule != reference_name, status == "ok") %>%
  inner_join(reference_outcomes, by = "posterior_draw_id") %>%
  mutate(
    delta_final_volume = final_volume - ref_final_volume,
    delta_auc = tumour_auc - ref_auc,
    beats_reference_final = final_volume < ref_final_volume,
    beats_reference_auc = tumour_auc < ref_auc
  )

paired_vs_reference_summary <- paired_vs_reference %>%
  group_by(schedule) %>%
  summarise(
    n_paired_draws = n(),
    delta_final_median = median(delta_final_volume),
    delta_final_q025 = quantile(delta_final_volume, 0.025),
    delta_final_q975 = quantile(delta_final_volume, 0.975),
    prob_lower_final_volume_than_reference = mean(beats_reference_final),
    delta_auc_median = median(delta_auc),
    delta_auc_q025 = quantile(delta_auc, 0.025),
    delta_auc_q975 = quantile(delta_auc, 0.975),
    prob_lower_auc_than_reference = mean(beats_reference_auc),
    .groups = "drop"
  )

write.csv(
  paired_vs_reference_summary,
  file.path(results_dir, "paired_comparisons_vs_observed_bid.csv"),
  row.names = FALSE
)

print(paired_vs_reference_summary)

############################################################
#  DIRECT QD-vs-BID COMPARISONS
############################################################

paired_compare <- function(df, schedule_A, schedule_B, label) {
  A <- df %>%
    filter(schedule == schedule_A, status == "ok") %>%
    select(posterior_draw_id, final_A = final_volume, auc_A = tumour_auc)

  B <- df %>%
    filter(schedule == schedule_B, status == "ok") %>%
    select(posterior_draw_id, final_B = final_volume, auc_B = tumour_auc)

  joined <- inner_join(A, B, by = "posterior_draw_id") %>%
    mutate(
      delta_final_A_minus_B = final_A - final_B,
      delta_auc_A_minus_B = auc_A - auc_B
    )

  data.frame(
    comparison = label,
    schedule_A = schedule_A,
    schedule_B = schedule_B,
    n_paired_draws = nrow(joined),
    median_delta_final_A_minus_B = median(joined$delta_final_A_minus_B),
    q025_delta_final_A_minus_B = quantile(joined$delta_final_A_minus_B, 0.025),
    q975_delta_final_A_minus_B = quantile(joined$delta_final_A_minus_B, 0.975),
    prob_A_lower_final_than_B = mean(joined$final_A < joined$final_B),
    median_delta_auc_A_minus_B = median(joined$delta_auc_A_minus_B),
    q025_delta_auc_A_minus_B = quantile(joined$delta_auc_A_minus_B, 0.025),
    q975_delta_auc_A_minus_B = quantile(joined$delta_auc_A_minus_B, 0.975),
    prob_A_lower_auc_than_B = mean(joined$auc_A < joined$auc_B)
  )
}

qd_bid_comparisons <- bind_rows(
  paired_compare(
    primary_results,
    schedule_A = "QD: 50 mg/kg",
    schedule_B = "Observed BID: 25+25 mg/kg",
    label = "Equal total daily dose: QD 50 vs BID 25+25"
  ),
  paired_compare(
    primary_results,
    schedule_A = "QD: 25 mg/kg",
    schedule_B = "Observed BID: 25+25 mg/kg",
    label = "Equal dose per administration: QD 25 vs BID 25+25"
  )
)

write.csv(
  qd_bid_comparisons,
  file.path(results_dir, "qd_vs_bid_paired_comparisons.csv"),
  row.names = FALSE
)

print(qd_bid_comparisons)

############################################################
# IDENTIFY BEST TIMING OFFSET FOR EACH POSTERIOR DRAW
############################################################

best_timing_by_draw_final <- timing_results %>%
  filter(status == "ok") %>%
  group_by(posterior_draw_id) %>%
  slice_min(final_volume, n = 1, with_ties = TRUE) %>%
  ungroup()

best_timing_by_draw_auc <- timing_results %>%
  filter(status == "ok") %>%
  group_by(posterior_draw_id) %>%
  slice_min(tumour_auc, n = 1, with_ties = TRUE) %>%
  ungroup()

best_timing_probability <- full_join(
  best_timing_by_draw_final %>%
    count(offset_hours, name = "n_best_final") %>%
    mutate(prob_best_final = n_best_final / nrow(posterior_parameters)),
  best_timing_by_draw_auc %>%
    count(offset_hours, name = "n_best_auc") %>%
    mutate(prob_best_auc = n_best_auc / nrow(posterior_parameters)),
  by = "offset_hours"
) %>%
  arrange(offset_hours) %>%
  mutate(
    across(c(n_best_final, n_best_auc), ~ replace_na(.x, 0L)),
    across(c(prob_best_final, prob_best_auc), ~ replace_na(.x, 0))
  )

write.csv(
  best_timing_probability,
  file.path(results_dir, "timing_offset_probability_of_being_best.csv"),
  row.names = FALSE
)

print(best_timing_probability)

############################################################
#  DISSERTATION-STYLE FIGURES
############################################################

p_final <- primary_summary %>%
  ggplot(aes(x = reorder(schedule, final_volume_median), y = final_volume_median)) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = final_volume_q025, ymax = final_volume_q975),
    width = 0.15
  ) +
  coord_flip() +
  labs(
    title = "Treatment-schedule optimisation: dosing-frequency comparison",
    subtitle = paste0("Posterior predictions at day ", FOLLOWUP_DAY, " under the Gompertz model"),
    x = NULL,
    y = expression(paste("Final tumour volume (mm"^3, ")"))
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(figures_dir, "primary_schedules_final_volume.png"),
  p_final,
  width = 9,
  height = 5.5,
  dpi = 300
)

p_auc <- primary_summary %>%
  ggplot(aes(x = reorder(schedule, auc_median), y = auc_median)) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(ymin = auc_q025, ymax = auc_q975),
    width = 0.15
  ) +
  coord_flip() +
  labs(
    title = "Treatment-schedule optimisation: cumulative tumour burden",
    subtitle = paste0("Posterior tumour-volume AUC from day 0 to day ", FOLLOWUP_DAY),
    x = NULL,
    y = expression(paste("Tumour-volume AUC (mm"^3, " day)"))
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(figures_dir, "primary_schedules_tumour_auc.png"),
  p_auc,
  width = 9,
  height = 5.5,
  dpi = 300
)

# 16C. Timing sweep: final tumour volume.
p_timing_final <- timing_summary %>%
  ggplot(aes(x = offset_hours, y = final_volume_median)) +
  geom_ribbon(
    aes(ymin = final_volume_q025, ymax = final_volume_q975),
    alpha = 0.18
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_vline(xintercept = -1, linetype = "dashed") +
  labs(
    title = "Treatment-schedule optimisation: drug-radiotherapy timing",
    subtitle = "One 25 mg/kg drug dose per RT fraction; dashed line marks the observed -1 h timing",
    x = expression(Delta*t == t[drug] - t[RT]~"(hours)"),
    y = expression(paste("Final tumour volume (mm"^3, ")"))
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(figures_dir, "timing_sweep_final_volume.png"),
  p_timing_final,
  width = 8,
  height = 5.5,
  dpi = 300
)

# 16D. Timing sweep: AUC.
p_timing_auc <- timing_summary %>%
  ggplot(aes(x = offset_hours, y = auc_median)) +
  geom_ribbon(
    aes(ymin = auc_q025, ymax = auc_q975),
    alpha = 0.18
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_vline(xintercept = -1, linetype = "dashed") +
  labs(
    title = "Treatment-schedule optimisation: drug-radiotherapy timing",
    subtitle = "Cumulative tumour burden across the posterior timing sweep",
    x = expression(Delta*t == t[drug] - t[RT]~"(hours)"),
    y = expression(paste("Tumour-volume AUC (mm"^3, " day)"))
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(figures_dir, "timing_sweep_tumour_auc.png"),
  p_timing_auc,
  width = 8,
  height = 5.5,
  dpi = 300
)

metadata <- data.frame(
  quantity = c(
    "model",
    "posterior_fit",
    "dataset",
    "n_posterior_draws",
    "followup_day",
    "initial_volume_mm3",
    "K_population_definition",
    "timing_offset_definition",
    "ode_rtol",
    "ode_atol"
  ),
  value = c(
    "Full joint Gompertz treatment model",
    fit_path,
    data_path,
    nrow(posterior_parameters),
    FOLLOWUP_DAY,
    V0_REFERENCE,
    "K = exp(logK_pop), propagated draw-by-draw from the final Stage-4 posterior",
    "delta_t = t_drug - t_RT; negative means drug before RT",
    ODE_RTOL,
    ODE_ATOL
  )
)

write.csv(
  metadata,
  file.path(results_dir, "analysis_metadata.csv"),
  row.names = FALSE
)

############################################################
# FINAL CONSOLE REPORT
############################################################

cat("\n============================================================\n")
cat("TREATMENT-SCHEDULE POSTERIOR SIMULATION COMPLETE\n")
cat("============================================================\n")
cat("Model: Gompertz full joint treatment model\n")
cat("Posterior draws: ", nrow(posterior_parameters), "\n", sep = "")
cat("Follow-up: day ", FOLLOWUP_DAY, "\n", sep = "")
cat("V0 reference: ", V0_REFERENCE, " mm^3\n", sep = "")
cat("K reference: posterior draw-specific exp(logK_pop); no observed-volume anchor used.\n")
cat("Results folder: ", results_dir, "\n\n", sep = "")

cat("QD vs BID comparisons:\n")
print(qd_bid_comparisons)

cat("\nTiming-sweep summary:\n")
print(timing_summary)

cat("\nProbability each timing offset is best:\n")
print(best_timing_probability)

cat("\nNEXT ANALYSIS AFTER REVIEWING THESE RESULTS:\n")
cat("  1. Refine the timing grid locally if a clear optimum appears.\n")
cat("  2. Repeat key comparisons over alternative V0 reference values.\n")
cat("  3. Repeat the key schedules under Exponential and Gompertz models.\n")
cat("============================================================\n")


combo_observations_stage5 <- data %>%
  filter(mouse %in% combo_mice, observation == 1) %>%
  arrange(mouse, time_days)

combo_baselines_stage5 <- combo_observations_stage5 %>%
  group_by(mouse) %>%
  slice_min(time_days, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(mouse, time_days, baseline_volume = volume)

stopifnot(
  nrow(combo_baselines_stage5) == 9L,
  all(is.finite(combo_baselines_stage5$baseline_volume)),
  all(combo_baselines_stage5$baseline_volume > 0)
)

V0_SENSITIVITY_VALUES <- as.numeric(quantile(
  combo_baselines_stage5$baseline_volume,
  probs = c(0, 0.5, 1),
  names = FALSE,
  type = 7
))
V0_SENSITIVITY_VALUES <- sort(unique(V0_SENSITIVITY_VALUES))

# Helper: exactly one result per draw/schedule is expected.
assert_unique_draw_schedule <- function(x, schedule_cols) {
  check <- x %>%
    count(across(all_of(c("posterior_draw_id", schedule_cols)))) %>%
    filter(n != 1L)
  if (nrow(check) > 0) stop("Duplicate/missing draw-schedule result rows detected.")
  invisible(TRUE)
}

complete_success_draws <- function(x, candidate_col, expected_n) {
  x %>%
    filter(status == "ok") %>%
    group_by(posterior_draw_id) %>%
    summarise(n_candidates = n_distinct(.data[[candidate_col]]), .groups = "drop") %>%
    filter(n_candidates == expected_n) %>%
    pull(posterior_draw_id)
}

############################################################
# REFINED DRUG-RT TIMING SEARCH
############################################################

REFINED_OFFSETS_HOURS <- c(
  -2.00, -1.50, -1.25, -1.00, -0.75, -0.50,
  -0.25,  0.00,  0.25,  0.50,  0.75,  1.00
)

refined_dir <- file.path(results_dir, "refined_timing")
dir.create(refined_dir, recursive = TRUE, showWarnings = FALSE)

refined_schedules <- setNames(
  lapply(REFINED_OFFSETS_HOURS, function(offset_h) {
    make_schedule(
      drug_times_days = rt_times + offset_h / 24,
      drug_doses = rep(25, length(rt_times)),
      rt_times_days = rt_times,
      rt_doses = ref_rt$dose
    )
  }),
  paste0("dt_", REFINED_OFFSETS_HOURS, "h")
)

if (min(unlist(lapply(refined_schedules, function(z) z$time_days))) < 0) {
  stop("A refined timing candidate contains treatment before t=0.")
}

refined_schedule_table <- imap_dfr(refined_schedules, function(schedule, nm) {
  offset_h <- as.numeric(sub("^dt_", "", sub("h$", "", nm)))
  mutate(schedule, offset_hours = offset_h, .before = 1)
})
write.csv(refined_schedule_table,
          file.path(refined_dir, "refined_timing_candidate_schedules.csv"),
          row.names = FALSE)

refined_results <- imap_dfr(refined_schedules, function(schedule, nm) {
  offset_h <- as.numeric(sub("^dt_", "", sub("h$", "", nm)))
  tmp <- run_one_schedule(
    schedule_name = paste0("Refined timing ", offset_h, " h"),
    schedule = schedule,
    parameter_table = posterior_parameters
  )
  tmp$offset_hours <- offset_h
  tmp
})

assert_unique_draw_schedule(refined_results, "offset_hours")
write.csv(refined_results,
          file.path(refined_dir, "refined_timing_posterior_outcomes.csv"),
          row.names = FALSE)

refined_summary <- refined_results %>%
  filter(status == "ok") %>%
  group_by(offset_hours) %>%
  summarise(
    n_draws = n(),
    final_volume_median = median(final_volume),
    final_volume_q025 = quantile(final_volume, 0.025),
    final_volume_q975 = quantile(final_volume, 0.975),
    auc_median = median(tumour_auc),
    auc_q025 = quantile(tumour_auc, 0.025),
    auc_q975 = quantile(tumour_auc, 0.975),
    .groups = "drop"
  ) %>%
  arrange(offset_hours)
write.csv(refined_summary,
          file.path(refined_dir, "refined_timing_summary.csv"),
          row.names = FALSE)

complete_refined_draws <- complete_success_draws(
  refined_results, "offset_hours", length(REFINED_OFFSETS_HOURS)
)
if (length(complete_refined_draws) == 0L) {
  stop("No posterior draw completed every refined timing candidate.")
}

best_refined_final <- refined_results %>%
  filter(status == "ok", posterior_draw_id %in% complete_refined_draws) %>%
  group_by(posterior_draw_id) %>%
  slice_min(final_volume, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  count(offset_hours, name = "n_best_final") %>%
  mutate(prob_best_final = n_best_final / length(complete_refined_draws))

best_refined_auc <- refined_results %>%
  filter(status == "ok", posterior_draw_id %in% complete_refined_draws) %>%
  group_by(posterior_draw_id) %>%
  slice_min(tumour_auc, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  count(offset_hours, name = "n_best_auc") %>%
  mutate(prob_best_auc = n_best_auc / length(complete_refined_draws))

refined_best_probability <- full_join(best_refined_final, best_refined_auc,
                                      by = "offset_hours") %>%
  complete(
    offset_hours = REFINED_OFFSETS_HOURS,
    fill = list(n_best_final = 0L, prob_best_final = 0,
                n_best_auc = 0L, prob_best_auc = 0)
  ) %>%
  arrange(offset_hours)
write.csv(refined_best_probability,
          file.path(refined_dir, "refined_timing_probability_of_being_best.csv"),
          row.names = FALSE)

reference_minus1 <- refined_results %>%
  filter(offset_hours == -1, status == "ok") %>%
  select(posterior_draw_id,
         ref_final_volume = final_volume,
         ref_auc = tumour_auc)

refined_paired_minus1 <- refined_results %>%
  filter(offset_hours != -1, status == "ok") %>%
  inner_join(reference_minus1, by = "posterior_draw_id") %>%
  mutate(delta_final = final_volume - ref_final_volume,
         delta_auc = tumour_auc - ref_auc) %>%
  group_by(offset_hours) %>%
  summarise(
    n_paired_draws = n(),
    delta_final_median = median(delta_final),
    delta_final_q025 = quantile(delta_final, 0.025),
    delta_final_q975 = quantile(delta_final, 0.975),
    prob_lower_final_than_minus1 = mean(delta_final < 0),
    delta_auc_median = median(delta_auc),
    delta_auc_q025 = quantile(delta_auc, 0.025),
    delta_auc_q975 = quantile(delta_auc, 0.975),
    prob_lower_auc_than_minus1 = mean(delta_auc < 0),
    .groups = "drop"
  ) %>% arrange(offset_hours)
write.csv(refined_paired_minus1,
          file.path(refined_dir, "refined_timing_paired_vs_minus1h.csv"),
          row.names = FALSE)

median_best_refined_final <- refined_summary %>%
  slice_min(final_volume_median, n = 1, with_ties = TRUE)
median_best_refined_auc <- refined_summary %>%
  slice_min(auc_median, n = 1, with_ties = TRUE)

p_refined_final <- ggplot(refined_summary,
                          aes(offset_hours, final_volume_median)) +
  geom_ribbon(aes(ymin = final_volume_q025, ymax = final_volume_q975), alpha = 0.18) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  geom_vline(xintercept = -1, linetype = "dashed") +
  labs(title = "Refined drug-radiotherapy timing search",
       subtitle = "One 25 mg/kg dose per 1 Gy RT fraction; dashed line = observed -1 h timing",
       x = expression(Delta*t == t[drug] - t[RT]~"(hours)"),
       y = expression(paste("Final tumour volume (mm"^3, ")"))) +
  theme_bw(base_size = 12)
ggsave(file.path(refined_dir, "refined_timing_final_volume.png"),
       p_refined_final, width = 8, height = 5.5, dpi = 300)

############################################################
# PK / REPAIR-FACTOR DIAGNOSTIC
############################################################

repair_dir <- file.path(results_dir, "repair_factor_diagnostic")
dir.create(repair_dir, recursive = TRUE, showWarnings = FALSE)

repair_schedules <- c(primary_schedules, refined_schedules)

run_repair_diagnostic <- function(schedule, schedule_name, parameter_table) {
  map_dfr(seq_len(nrow(parameter_table)), function(i) {
    pars <- as.list(parameter_table[i, setdiff(names(parameter_table), "posterior_draw_id")])
    tryCatch({
      traj <- simulate_schedule(pars = pars, schedule = schedule)
      ratio <- pars[["k_repi"]] * traj$C2
      data.frame(
        posterior_draw_id = parameter_table$posterior_draw_id[i],
        schedule = schedule_name,
        status = "ok",
        max_k_repi_C2 = max(ratio),
        min_repair_factor = min(1 - ratio),
        ever_k_repi_C2_gt_1 = any(ratio > 1),
        final_volume = tail(traj$V, 1)
      )
    }, error = function(e) {
      data.frame(
        posterior_draw_id = parameter_table$posterior_draw_id[i],
        schedule = schedule_name,
        status = paste0("ERROR: ", conditionMessage(e)),
        max_k_repi_C2 = NA_real_,
        min_repair_factor = NA_real_,
        ever_k_repi_C2_gt_1 = NA,
        final_volume = NA_real_
      )
    })
  })
}

repair_results <- imap_dfr(repair_schedules, run_repair_diagnostic,
                           parameter_table = posterior_parameters)
assert_unique_draw_schedule(repair_results, "schedule")
write.csv(repair_results,
          file.path(repair_dir, "repair_factor_diagnostic_by_draw.csv"),
          row.names = FALSE)

repair_summary <- repair_results %>%
  filter(status == "ok") %>%
  group_by(schedule) %>%
  summarise(
    n_draws = n(),
    prob_k_repi_C2_gt_1 = mean(ever_k_repi_C2_gt_1),
    max_ratio_median = median(max_k_repi_C2),
    max_ratio_q025 = quantile(max_k_repi_C2, 0.025),
    max_ratio_q975 = quantile(max_k_repi_C2, 0.975),
    min_repair_factor_median = median(min_repair_factor),
    .groups = "drop"
  ) %>% arrange(desc(prob_k_repi_C2_gt_1), schedule)
write.csv(repair_summary,
          file.path(repair_dir, "repair_factor_diagnostic_summary.csv"),
          row.names = FALSE)

############################################################
#  REFERENCE-STATE SENSITIVITY ANALYSIS
############################################################

SENSITIVITY_OFFSETS_HOURS <- c(-1, -0.5, -0.25, 0)

sensitivity_schedules <- setNames(
  lapply(SENSITIVITY_OFFSETS_HOURS, function(offset_h) {
    make_schedule(
      drug_times_days = rt_times + offset_h / 24,
      drug_doses = rep(25, length(rt_times)),
      rt_times_days = rt_times,
      rt_doses = ref_rt$dose
    )
  }),
  paste0("dt_", SENSITIVITY_OFFSETS_HOURS, "h")
)


sensitivity_dir <- file.path(results_dir, "v0_sensitivity")
dir.create(sensitivity_dir, recursive = TRUE, showWarnings = FALSE)

run_sensitivity_case <- function(V0_value, offset_h, schedule) {
  pars_table <- posterior_parameters

  map_dfr(seq_len(nrow(pars_table)), function(i) {
    pars <- as.list(pars_table[i, setdiff(names(pars_table), "posterior_draw_id")])
    tryCatch({
      traj <- simulate_schedule(pars = pars, schedule = schedule, V0 = V0_value)
      smry <- summarise_trajectory(traj)
      data.frame(posterior_draw_id = pars_table$posterior_draw_id[i],
                 V0 = V0_value, offset_hours = offset_h, status = "ok",
                 final_volume = smry$final_volume, tumour_auc = smry$tumour_auc)
    }, error = function(e) {
      data.frame(posterior_draw_id = pars_table$posterior_draw_id[i],
                 V0 = V0_value, offset_hours = offset_h,
                 status = paste0("ERROR: ", conditionMessage(e)),
                 final_volume = NA_real_, tumour_auc = NA_real_)
    })
  })
}

sensitivity_grid <- expand.grid(
  V0 = V0_SENSITIVITY_VALUES,
  offset_hours = SENSITIVITY_OFFSETS_HOURS,
  KEEP.OUT.ATTRS = FALSE
)

sensitivity_results <- map_dfr(seq_len(nrow(sensitivity_grid)), function(j) {
  cc <- sensitivity_grid[j, ]
  nm <- paste0("dt_", cc$offset_hours, "h")
  run_sensitivity_case(cc$V0, cc$offset_hours, sensitivity_schedules[[nm]])
})
write.csv(sensitivity_results,
          file.path(sensitivity_dir, "v0_sensitivity_posterior_outcomes.csv"),
          row.names = FALSE)

sensitivity_summary <- sensitivity_results %>%
  filter(status == "ok") %>%
  group_by(V0, offset_hours) %>%
  summarise(n_draws = n(),
            final_volume_median = median(final_volume),
            final_volume_q025 = quantile(final_volume, .025),
            final_volume_q975 = quantile(final_volume, .975),
            auc_median = median(tumour_auc),
            auc_q025 = quantile(tumour_auc, .025),
            auc_q975 = quantile(tumour_auc, .975),
            .groups = "drop")
write.csv(sensitivity_summary,
          file.path(sensitivity_dir, "v0_sensitivity_summary.csv"), row.names = FALSE)

robustness_table <- sensitivity_summary %>%
  group_by(V0) %>%
  summarise(best_offset_final = offset_hours[which.min(final_volume_median)],
            best_offset_auc = offset_hours[which.min(auc_median)],
            .groups = "drop")
write.csv(robustness_table,
          file.path(sensitivity_dir, "v0_robustness_summary.csv"), row.names = FALSE)

write.csv(combo_baselines_stage5,
          file.path(sensitivity_dir, "combination_mouse_baseline_volumes.csv"),
          row.names = FALSE)

all_status_objects <- list(
  primary = primary_results$status,
  coarse_timing = timing_results$status,
  refined_timing = refined_results$status,
  repair_factor = repair_results$status,
  sensitivity = sensitivity_results$status
)

integrity_report <- bind_rows(lapply(names(all_status_objects), function(nm) {
  x <- all_status_objects[[nm]]
  data.frame(
    analysis = nm,
    n_rows = length(x),
    n_ok = sum(x == "ok"),
    n_failed = sum(x != "ok")
  )
}))
write.csv(integrity_report,
          file.path(results_dir, "STAGE5_complete_integrity_report.csv"),
          row.names = FALSE)

if (any(integrity_report$n_failed > 0)) {
  warning("At least one Stage-5 simulation failed. Inspect the integrity report and corresponding outcome CSVs before interpretation.")
}

stage5_metadata <- data.frame(
  quantity = c(
    "posterior_fit", "dataset", "posterior_draws_used", "followup_day",
    "coarse_timing_candidates", "refined_timing_candidates",
    "sensitivity_timing_candidates", "V0_sensitivity_values",
    "ode_rtol", "ode_atol"
  ),
  value = c(
    fit_path, data_path, nrow(posterior_parameters), FOLLOWUP_DAY,
    paste(TIMING_OFFSETS_HOURS, collapse = ","),
    paste(REFINED_OFFSETS_HOURS, collapse = ","),
    paste(SENSITIVITY_OFFSETS_HOURS, collapse = ","),
    paste(signif(V0_SENSITIVITY_VALUES, 8), collapse = ","),
    ODE_RTOL, ODE_ATOL
  )
)
write.csv(stage5_metadata,
          file.path(results_dir, "STAGE5_complete_metadata.csv"),
          row.names = FALSE)

cat("\n============================================================\n")
cat("COMPLETE STAGE-5 ANALYSIS FINISHED\n")
cat("============================================================\n")
print(integrity_report, row.names = FALSE)
cat("\nRefined median-optimal timing (final volume):\n")
print(median_best_refined_final, row.names = FALSE)
cat("\nRefined median-optimal timing (AUC):\n")
print(median_best_refined_auc, row.names = FALSE)
cat("\nReference-state sensitivity best timings:\n")
print(as.data.frame(robustness_table), row.names = FALSE)
cat("\nRepair-factor diagnostic summary:\n")
print(as.data.frame(repair_summary), row.names = FALSE)
cat("\nAll results saved under:\n", results_dir, "\n", sep = "")
cat("============================================================\n")
