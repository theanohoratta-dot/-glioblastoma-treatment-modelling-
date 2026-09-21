##################################################
# STAGE 1: UNTREATED TUMOUR GROWTH
##################################################
library(cmdstanr)
library(dplyr)
library(posterior)
library(loo)
library(bayesplot)
library(ggplot2)

setwd("~/Desktop/dissertation")

data <- read.csv("~/Desktop/dissertation/data/lemassonParpi.csv")

#output folders

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

control_mice <- data %>%
  group_by(mouse) %>%
  summarise(
    received_treatment = any(
      modality %in% c(1, 2),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  filter(!received_treatment) %>%
  pull(mouse)

control_data <- data %>%
  filter(mouse %in% control_mice, observation == 1) %>%
  mutate(
    time_days = time / 24,
    mouse_id = as.integer(factor(mouse))
  )

# Validate untreated-control data 


stopifnot(
  nrow(control_data) > 0,
  all(!is.na(control_data$mouse_id)),
  all(is.finite(control_data$time_days)),
  all(is.finite(control_data$volume)),
  all(control_data$time_days >= 0),
  all(control_data$volume > 0)
)

duplicate_rows <- control_data %>%
  group_by(mouse, time_days) %>%
  filter(n() > 1) %>%
  arrange(mouse, time_days, volume) %>%
  ungroup()

print(duplicate_rows)

stopifnot(nrow(duplicate_rows) == 0)

# Check baseline observation times
control_baseline_check <-
  control_data %>%
  group_by(
    mouse,
    mouse_id
  ) %>%
  summarise(
    first_time_days =
      min(time_days),
    
    first_volume =
      volume[
        which.min(time_days)
      ],
    
    number_of_observations =
      n(),
    
    .groups = "drop"
  ) %>%
  arrange(
    mouse_id
  )

print(
  control_baseline_check
)

write.csv(
  control_baseline_check,
  "results/control_baseline_check.csv",
  row.names = FALSE
)

# Construct Stan data


stan_data_exp <- list(
  N = nrow(control_data),
  M = n_distinct(control_data$mouse_id),
  mouse = control_data$mouse_id,
  time = control_data$time_days,
  volume = control_data$volume
)

# Stan data for Logistic and Gompertz models
stan_data_bounded <- stan_data_exp
##################################################
# PRIOR PREDICTIVE CHECKS-
# simulations performed before fitting and
# verify that the final weakly informative
# priors generate plausible tumour-volume paths
###################################################


dir.create(
  file.path("figures", "prior_predictive"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path("results", "prior_predictive"),
  recursive = TRUE,
  showWarnings = FALSE
)

set.seed(123)

n_prior_draws <- 1000
n_prior_plot_draws <- 100

prior_time_grid <- seq(
  from = 0,
  to = max(control_data$time_days),
  length.out = 100
)

rhalfnorm <- function(n, sd) {
  abs(rnorm(n, mean = 0, sd = sd))
}

simulate_prior_model <- function(
    model_name,
    logr_centre,
    use_K = FALSE
) {
  
  logV0_pop <-
    rnorm(
      n_prior_draws,
      mean = log(30),
      sd = 0.5
    )
  
  logr_pop <-
    rnorm(
      n_prior_draws,
      mean = log(logr_centre),
      sd = 0.5
    )
  
  sigma_logV0 <-
    rhalfnorm(
      n_prior_draws,
      sd = 0.5
    )
  
  sigma_logr <-
    rhalfnorm(
      n_prior_draws,
      sd = 0.5
    )
  
  sigma_obs <-
    rhalfnorm(
      n_prior_draws,
      sd = 0.5
    )
  
  z_logV0 <-
    rnorm(n_prior_draws)
  
  z_logr <-
    rnorm(n_prior_draws)
  
  V0 <-
    exp(
      logV0_pop +
        sigma_logV0 * z_logV0
    )
  
  r <-
    exp(
      logr_pop +
        sigma_logr * z_logr
    )
  
  if (use_K) {
    
    logK_pop <-
      rnorm(
        n_prior_draws,
        mean = log(500),
        sd = 0.7
      )
    
    sigma_logK <-
      rhalfnorm(
        n_prior_draws,
        sd = 0.5
      )
    
    z_logK <-
      rnorm(
        n_prior_draws
      )
    
    K <-
      exp(
        logK_pop +
          sigma_logK * z_logK
      )
    
  } else {
    
    K <-
      rep(
        NA_real_,
        n_prior_draws
      )
  }
  
  trajectory_list <-
    vector(
      "list",
      n_prior_draws
    )
  
  for (i in seq_len(n_prior_draws)) {
    
    t <-
      prior_time_grid
    
    if (model_name == "Exponential") {
      
      mu <-
        V0[i] *
        exp(
          r[i] * t
        )
      
    } else if (model_name == "Logistic") {
      
      denominator <-
        1 +
        ((K[i] - V0[i]) / V0[i]) *
        exp( -r[i] * t )
      mu <-
        K[i] / denominator
      
    } else if (model_name == "Gompertz") {
      
      mu <-
        K[i] *
        exp(
          log(
            V0[i] /
              K[i]
          ) *
            exp(
              -r[i] * t
            )
        )
      
    } else {
      
      stop("Unknown model_name")
    }
    
    valid_mu <-
      is.finite(mu) &
      mu > 0
    
    y_rep <-
      rep(
        NA_real_,
        length(mu)
      )
    
    if (all(valid_mu)) {
      
      y_rep <-
        rlnorm(
          length(mu),
          meanlog = log(mu),
          sdlog = sigma_obs[i]
        )
    }
    
    trajectory_list[[i]] <-
      data.frame(
        model = model_name,
        draw = i,
        time_days = t,
        mu = mu,
        y_rep = y_rep,
        valid =
          all(valid_mu)
      )
  }
  
  trajectories <-
    bind_rows(
      trajectory_list
    )
  
  parameter_draws <-
    data.frame(
      model = model_name,
      draw = seq_len(n_prior_draws),
      
      logV0_pop =
        logV0_pop,
      
      logr_pop =
        logr_pop,
      
      sigma_logV0 =
        sigma_logV0,
      
      sigma_logr =
        sigma_logr,
      
      V0 =
        V0,
      
      r =
        r,
      
      K =
        K,
      
      sigma =
        sigma_obs
    )
  list(
    trajectories = trajectories,
    parameters = parameter_draws
  )
}

##################################################
# Simulating all three Stage 1 models
##################################################

prior_exp <-
  simulate_prior_model(
    model_name = "Exponential",
    logr_centre = 0.3,
    use_K = FALSE
  )

prior_logistic <-
  simulate_prior_model(
    model_name = "Logistic",
    logr_centre = 0.4,
    use_K = TRUE
  )

prior_gompertz <-
  simulate_prior_model(
    model_name = "Gompertz",
    logr_centre = 0.2,
    use_K = TRUE
  )

prior_trajectories <-
  bind_rows(
    prior_exp$trajectories,
    prior_logistic$trajectories,
    prior_gompertz$trajectories
  )

prior_parameter_draws <-
  bind_rows(
    prior_exp$parameters,
    prior_logistic$parameters,
    prior_gompertz$parameters
  )

# Numerical prior-predictive diagnostics
prior_validity_summary <-
  prior_trajectories %>%
  group_by(
    model,
    draw
  ) %>%
  summarise(
    valid =
      first(valid),
    max_mu =
      if (all(is.finite(mu))) {
        max(mu)
      } else {
        NA_real_
      },
    max_y_rep =
      if (all(is.finite(y_rep))) {
        max(y_rep)
      } else {
        NA_real_
      },
    .groups = "drop"
  ) %>%
  group_by(model) %>%
  summarise(
    n_draws = n(),
    n_valid = sum(valid),
    proportion_valid = mean(valid),
    median_max_mu =
      median(
        max_mu,
        na.rm = TRUE
      ),
    q95_max_mu =
      quantile(
        max_mu,
        0.95,
        na.rm = TRUE
      ),
    median_max_y_rep =
      median(
        max_y_rep,
        na.rm = TRUE
      ),
    q95_max_y_rep =
      quantile(
        max_y_rep,
        0.95,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

print(
  prior_validity_summary
)

write.csv(
  prior_validity_summary,
  file.path(
    "results",
    "prior_predictive",
    "prior_predictive_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  prior_parameter_draws,
  file.path(
    "results",
    "prior_predictive",
    "prior_parameter_draws.csv"
  ),
  row.names = FALSE
)

# Prior parameter summaries


prior_parameter_summary <-
  prior_parameter_draws %>%
  group_by(model) %>%
  summarise(
    V0_median =
      median(V0),
    V0_q025 =
      quantile(V0, 0.025),
    V0_q975 =
      quantile(V0, 0.975),
    r_median =
      median(r),
    r_q025 =
      quantile(r, 0.025),
    r_q975 =
      quantile(r, 0.975),
    K_median =
      if (all(is.na(K))) {
        NA_real_
      } else {
        median(K, na.rm = TRUE)
      },
    K_q025 =
      if (all(is.na(K))) {
        NA_real_
      } else {
        quantile(K, 0.025, na.rm = TRUE)
      },
    K_q975 =
      if (all(is.na(K))) {
        NA_real_
      } else {
        quantile(K, 0.975, na.rm = TRUE)
      },
    sigma_median =
      median(sigma),
    sigma_q025 =
      quantile(sigma, 0.025),
    sigma_q975 =
      quantile(sigma, 0.975),
    .groups = "drop"
  )

print(
  prior_parameter_summary
)

write.csv(
  prior_parameter_summary,
  file.path(
    "results",
    "prior_predictive",
    "prior_parameter_summary.csv"
  ),
  row.names = FALSE
)

# Prior diagnostic:
# probability of K < = V0


prior_K_V0_diagnostic <-
  prior_parameter_draws %>%
  filter(
    model %in%
      c(
        "Logistic",
        "Gompertz"
      )
  ) %>%
  group_by(model) %>%
  summarise(
    n_prior_draws = n(),
    
    n_K_le_V0 =
      sum(
        K <= V0,
        na.rm = TRUE
      ),
    
    proportion_K_le_V0 =
      mean(
        K <= V0,
        na.rm = TRUE
      ),
    
    proportion_K_gt_V0 =
      mean(
        K > V0,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

print(
  prior_K_V0_diagnostic
)

write.csv(
  prior_K_V0_diagnostic,
  file.path(
    "results",
    "prior_predictive",
    "prior_K_vs_V0_diagnostic.csv"
  ),
  row.names = FALSE
)

# Prior predictive structural trajectories
prior_plot_data <-
  prior_trajectories %>%
  filter(
    draw <= n_prior_plot_draws,
    valid,
    is.finite(mu),
    mu > 0
  )

prior_predictive_plot <-
  ggplot(
    prior_plot_data,
    aes(
      x = time_days,
      y = mu,
      group =
        interaction(
          model,
          draw
        )
    )
  ) +
  geom_line(
    alpha = 0.12,
    linewidth = 0.35
  ) +
  geom_point(
    data = control_data,
    aes(
      x = time_days,
      y = volume
    ),
    inherit.aes = FALSE,
    alpha = 0.35,
    size = 1
  ) +
  facet_wrap(
    ~ model,
    scales = "free_y"
  ) +
  scale_y_log10() +
  labs(
    title =
      "Stage 1 prior predictive structural trajectories",
    subtitle =
      "100 prior trajectories per model; observed control data shown only as a scale reference",
    x =
      "Time (days)",
    y =
      "Tumour volume (mm^3, log scale)"
  ) +
  theme_bw()

print(
  prior_predictive_plot
)

ggsave(
  file.path(
    "figures",
    "prior_predictive",
    "prior_predictive_structural_trajectories.png"
  ),
  prior_predictive_plot,
  width = 11,
  height = 5.5,
  dpi = 300
)


# Prior predictive noisy observationss
prior_noisy_plot_data <-
  prior_trajectories %>%
  filter(
    draw <= n_prior_plot_draws,
    valid,
    is.finite(y_rep),
    y_rep > 0
  )

prior_predictive_noisy_plot <-
  ggplot(
    prior_noisy_plot_data,
    aes(
      x = time_days,
      y = y_rep,
      group =
        interaction(
          model,
          draw
        )
    )
  ) +
  geom_line(
    alpha = 0.08,
    linewidth = 0.3
  ) +
  geom_point(
    data = control_data,
    aes(
      x = time_days,
      y = volume
    ),
    inherit.aes = FALSE,
    alpha = 0.35,
    size = 1
  ) +
  facet_wrap(
    ~ model,
    scales = "free_y"
  ) +
  scale_y_log10() +
  labs(
    title =
      "Stage 1 prior predictive replicated observations",
    subtitle =
      "Includes residual lognormal observation variability",
    x =
      "Time (days)",
    y =
      "Replicated tumour volume (mm^3, log scale)"
  ) +
  theme_bw()

print(
  prior_predictive_noisy_plot
)

ggsave(
  file.path(
    "figures",
    "prior_predictive",
    "prior_predictive_replicated_observations.png"
  ),
  prior_predictive_noisy_plot,
  width = 11,
  height = 5.5,
  dpi = 300
)

cat(
  "\nPRIOR PREDICTIVE CHECKS COMPLETED.\n",
  "See results/prior_predictive and figures/prior_predictive.\n"
)


##################################################
#1. Exponential model
##################################################

writeLines('
data {
  int<lower=1> N;
  int<lower=1> M;
  array[N] int<lower=1, upper=M> mouse;
  vector<lower=0>[N] time;
  vector<lower=0>[N] volume;
}

parameters {
  real logV0_pop;
  real logr_pop;
  real<lower=0> sigma_logV0;
  real<lower=0> sigma_logr;
  vector[M] z_logV0;
  vector[M] z_logr;
  real<lower=0> sigma;
}

transformed parameters {
  vector[M] logV0_mouse = logV0_pop + sigma_logV0 * z_logV0;
  vector[M] logr_mouse = logr_pop + sigma_logr * z_logr;
}

model {
  vector[N] mu_log;

  logV0_pop ~ normal(log(30), 0.5);
  logr_pop ~ normal(log(0.3), 0.5);
  sigma_logV0 ~ normal(0, 0.5);
  sigma_logr ~ normal(0, 0.5);
  z_logV0 ~ normal(0, 1);
  z_logr ~ normal(0, 1);
  sigma ~ normal(0, 0.5);

  for (n in 1:N) {
    mu_log[n] = logV0_mouse[mouse[n]] + exp(logr_mouse[mouse[n]]) * time[n];
  }

  log(volume) ~ normal(mu_log, sigma);
}

generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  vector[N] mu_log_out;
  vector[N] mu_out;

  for (n in 1:N) {
    real mu_log_n =
      logV0_mouse[mouse[n]] +
      exp(logr_mouse[mouse[n]]) * time[n];

    mu_log_out[n] = mu_log_n;
    mu_out[n] = exp(mu_log_n);

log_lik[n] =
  lognormal_lpdf(volume[n] | mu_log_n, sigma);

    y_rep[n] =
      lognormal_rng(mu_log_n, sigma);
  }
}
', "exp_control_mixed_loo.stan")

##################################################
# 2. Logistic model
##################################################

writeLines('
data {
  int<lower=1> N;
  int<lower=1> M;

  array[N] int<lower=1, upper=M> mouse;

  vector[N] time;

  vector<lower=0>[N] volume;


}

parameters {
  real logV0_pop;
  real logr_pop;
  real logK_pop;
  real<lower=0> sigma_logV0;
  real<lower=0> sigma_logr;
  real<lower=0> sigma_logK;
  vector[M] z_logV0;
  vector[M] z_logr;
  vector[M] z_logK;

  real<lower=0> sigma;
}

transformed parameters {

  vector[M] V0_mouse;
  vector[M] r_mouse;
  vector[M] K_mouse;
  for (m in 1:M) {
    V0_mouse[m] =
      exp(
        logV0_pop +
        sigma_logV0 *
        z_logV0[m]
      );
    r_mouse[m] =
      exp(
        logr_pop +
        sigma_logr *
        z_logr[m]
      );

K_mouse[m] =
  exp(
    logK_pop +
    sigma_logK *
    z_logK[m]
  );
  }
}

model {

  vector[N] mu_log;

  logV0_pop ~
    normal(
      log(30),
      0.5
    );

  logr_pop ~
    normal(
      log(0.4),
      0.5
    );

  logK_pop ~
    normal(
      log(500),
      0.7
    );

  sigma_logV0 ~
    normal(
      0,
      0.5
    );

  sigma_logr ~
    normal(
      0,
      0.5
    );

  sigma_logK ~
    normal(
      0,
      0.5
    );

  z_logV0 ~
    std_normal();

  z_logr ~
    std_normal();

  z_logK ~
    std_normal();
    
  sigma ~
    normal(
      0,
      0.5
    );

  for (n in 1:N) {

    real V0 =
      V0_mouse[
        mouse[n]
      ];

    real r =
      r_mouse[
        mouse[n]
      ];

    real K =
      K_mouse[
        mouse[n]
      ];

    real mu =
      K /
      (1 +((K - V0) / V0) *exp( -r * time[n])
      );

    mu_log[n] =
      log(mu);
  }

  log(volume) ~
    normal(
      mu_log,
      sigma
    );
}

generated quantities {

  vector[N] log_lik;
  vector[N] y_rep;
  vector[N] mu_log_out;
  vector[N] mu_out;

  for (n in 1:N) {

    real V0 =
      V0_mouse[
        mouse[n]
      ];

    real r =
      r_mouse[
        mouse[n]
      ];

    real K =
      K_mouse[
        mouse[n]
      ];

    real mu =
      K /
      (1 +((K - V0) / V0) *exp(-r * time[n])
      );

    mu_log_out[n] =
      log(mu);

    mu_out[n] =
      mu;

log_lik[n] =
  lognormal_lpdf(
    volume[n] |
    mu_log_out[n],
    sigma
  );

    y_rep[n] =
      lognormal_rng(
        mu_log_out[n],
        sigma
      );
  }
}
',
"logistic_control_mixed_loo.stan"
)
##################################################
# 3. Gompertz model
##################################################

writeLines('
data {
  int<lower=1> N;
  int<lower=1> M;

  array[N] int<lower=1, upper=M> mouse;

  vector[N] time;

  vector<lower=0>[N] volume;
}

parameters {
  real logV0_pop;
  real logr_pop;
  real logK_pop;
  real<lower=0> sigma_logV0;
  real<lower=0> sigma_logr;
  real<lower=0> sigma_logK;
  vector[M] z_logV0;
  vector[M] z_logr;
  vector[M] z_logK;
  
  real<lower=0> sigma;
}

transformed parameters {

  vector[M] V0_mouse;
  vector[M] r_mouse;
  vector[M] K_mouse;
  for (m in 1:M) {
    V0_mouse[m] =
      exp(
        logV0_pop +
        sigma_logV0 * z_logV0[m]
      );

    r_mouse[m] =
      exp(
        logr_pop +
        sigma_logr * z_logr[m]
      );

K_mouse[m] =
  exp(
    logK_pop +
    sigma_logK * z_logK[m]
  );
  }
}
model {
  vector[N] mu_log;
 logV0_pop ~ normal(log(30), 0.5);
logr_pop ~ normal(log(0.2), 0.5);
logK_pop ~ normal(log(500), 0.7);

sigma_logV0 ~ normal(0, 0.5);
sigma_logr ~ normal(0, 0.5);
sigma_logK ~ normal(0, 0.5);

z_logV0 ~ std_normal();
z_logr ~ std_normal();
z_logK ~ std_normal();

sigma ~ normal(0, 0.5);

  for (n in 1:N) {
  real V0 = V0_mouse[mouse[n]];
  real r  = r_mouse[mouse[n]];
  real K  = K_mouse[mouse[n]];

  real mu =
    K *
    exp(
      log(V0 / K) *
      exp(-r * time[n])
    );

  mu_log[n] = log(mu);
}

log(volume) ~ normal(mu_log, sigma);
}

generated quantities {

  vector[N] log_lik;
  vector[N] y_rep;
  vector[N] mu_log_out;
  vector[N] mu_out;

  for (n in 1:N) {

  real V0 = V0_mouse[mouse[n]];
  real r  = r_mouse[mouse[n]];
  real K  = K_mouse[mouse[n]];

  real mu =
    K *
    exp(
      log(V0 / K) *
      exp(-r * time[n])
    );

  mu_log_out[n] = log(mu);
  mu_out[n] = mu;

log_lik[n] =
  lognormal_lpdf(
    volume[n] |
    mu_log_out[n],
    sigma
  );

  y_rep[n] =
    lognormal_rng(
      mu_log_out[n],
      sigma
    );
}
}
',
"gompertz_control_mixed_loo.stan"
)

# Fit models


fit_one <- function(stan_file, stan_data) {
  mod <- cmdstan_model(stan_file)
  mod$sample(
    data = stan_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 1000,
    adapt_delta = 0.99,
    seed = 123,
    refresh = 100
  )
}

fit_exp <- fit_one(
  "exp_control_mixed_loo.stan",
  stan_data_exp
)

fit_logistic <- fit_one(
  "logistic_control_mixed_loo.stan",
  stan_data_bounded
)

fit_gompertz <- fit_one(
  "gompertz_control_mixed_loo.stan",
  stan_data_bounded
)

# Posterior diagnostic:
calculate_K_V0_posterior_diagnostic <-
  function(
    fit,
    model_name,
    control_data
  ) {
    
    K_draws <-
      as_draws_matrix(
        fit$draws(
          "K_mouse"
        )
      )
    
    V0_draws <-
      as_draws_matrix(
        fit$draws(
          "V0_mouse"
        )
      )
    
    stopifnot(
      nrow(K_draws) ==
        nrow(V0_draws),
      
      ncol(K_draws) ==
        ncol(V0_draws),
      
      all(
        is.finite(
          K_draws
        )
      ),
      
      all(
        is.finite(
          V0_draws
        )
      ),
      
      all(
        K_draws > 0
      ),
      
      all(
        V0_draws > 0
      )
    )
    
    mouse_lookup <-
      control_data %>%
      distinct(
        mouse,
        mouse_id
      ) %>%
      arrange(
        mouse_id
      )
    
    stopifnot(
      nrow(mouse_lookup) ==
        ncol(K_draws)
    )
    
    probability_K_gt_V0 <-
      colMeans(
        K_draws >
          V0_draws
      )
    
    probability_K_le_V0 <-
      colMeans(
        K_draws <=
          V0_draws
      )
    
    diagnostic <-
      mouse_lookup %>%
      mutate(
        model =
          model_name,
        
        posterior_probability_K_gt_V0 =
          probability_K_gt_V0,
        
        posterior_probability_K_le_V0 =
          probability_K_le_V0
      ) %>%
      select(
        model,
        mouse,
        mouse_id,
        posterior_probability_K_gt_V0,
        posterior_probability_K_le_V0
      )
    
    return(
      diagnostic
    )
  }


logistic_K_V0_posterior_diagnostic <-
  calculate_K_V0_posterior_diagnostic(
    fit =
      fit_logistic,
    model_name =
      "Logistic",
    control_data =
      control_data
  )


gompertz_K_V0_posterior_diagnostic <-
  calculate_K_V0_posterior_diagnostic(
    fit =
      fit_gompertz,
    model_name =
      "Gompertz",
    control_data =
      control_data
  )


posterior_K_V0_diagnostic <-
  bind_rows(
    logistic_K_V0_posterior_diagnostic,
    gompertz_K_V0_posterior_diagnostic
  )


print(
  as.data.frame(
    posterior_K_V0_diagnostic
  ),
  row.names = FALSE
)

write.csv(
  posterior_K_V0_diagnostic,
  file.path(
    "results",
    "posterior_K_vs_V0_diagnostic.csv"
  ),
  row.names = FALSE
)

fit_exp$save_object(
  "results/fit_exp_control.rds"
)

fit_logistic$save_object(
  "results/fit_logistic_control.rds"
)

fit_gompertz$save_object(
  "results/fit_gompertz_control.rds"
)

# Posterior summaries - 95% credible intervals

summarise_model_posterior <- function(fit, variables) {
  
  posterior::summarise_draws(
    fit$draws(variables),
    
    mean,
    median,
    sd,
    
    q2.5 = function(x) {
      quantile(x, probs = 0.025)
    },
    
    q97.5 = function(x) {
      quantile(x, probs = 0.975)
    },
    
    posterior::rhat,
    posterior::ess_bulk,
    posterior::ess_tail
  )
}

# Exponential model

summary_exp <- summarise_model_posterior(
  fit_exp,
  variables = c(
    "logV0_pop",
    "logr_pop",
    "sigma_logV0",
    "sigma_logr",
    "sigma"
  )
)

print(summary_exp)

# Logistic model

summary_logistic <- summarise_model_posterior(
  fit_logistic,
  variables = c(
    "logV0_pop",
    "logr_pop",
    "logK_pop",
    "sigma_logV0",
    "sigma_logr",
    "sigma_logK",
    "sigma"
  )
)

print(summary_logistic)

# Gompertz model

summary_gompertz <- summarise_model_posterior(
  fit_gompertz,
  variables = c(
    "logV0_pop",
    "logr_pop",
    "logK_pop",
    "sigma_logV0",
    "sigma_logr",
    "sigma_logK",
    "sigma"
  )
)

print(summary_gompertz)

# Residual diagnostics for all models


create_residual_diagnostics <- function(
    fit,
    data,
    model_name
) {
  
  mu_log_draws <- as_draws_matrix(
    fit$draws("mu_log_out")
  )
  
  mu_draws <- as_draws_matrix(
    fit$draws("mu_out")
  )
  
  fitted_log <- colMeans(mu_log_draws)
  fitted_original <- colMeans(mu_draws)
  
  sigma_hat <- fit$summary("sigma")$mean
  
  stopifnot(
    length(fitted_log) == nrow(data),
    length(fitted_original) == nrow(data)
  )
  
  diagnostic_data <- data %>%
    mutate(
      model = model_name,
      fitted_log = fitted_log,
      fitted_volume = fitted_original,
      log_residual =
        log(volume) - fitted_log,
      standardised_log_residual =
        log_residual / sigma_hat
    )  
  
  return(diagnostic_data)
}
residuals_exp <- create_residual_diagnostics(
  fit = fit_exp,
  data = control_data,
  model_name = "Exponential"
)

residuals_logistic <- create_residual_diagnostics(
  fit = fit_logistic,
  data = control_data,
  model_name = "Logistic"
)

residuals_gompertz <- create_residual_diagnostics(
  fit = fit_gompertz,
  data = control_data,
  model_name = "Gompertz"
)
residuals_all_models <- bind_rows(
  residuals_exp,
  residuals_logistic,
  residuals_gompertz
)
residual_fitted_plot <- ggplot(
  residuals_all_models,
  aes(
    x = fitted_log,
    y = standardised_log_residual,
    colour = factor(mouse_id)
  )
) +
  geom_point(alpha = 0.7) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = c(-2, 2),
    linetype = "dotted"
  ) +
  facet_wrap(
    ~ model,
    scales = "free_x"
  ) +
  labs(
    x = "Posterior mean log fitted tumour volume",
    y = "Standardised log residual",
    colour = "Mouse",
    title = "Residual diagnostics for untreated tumour-growth models"
  ) +
  theme_minimal()

residual_fitted_plot

# Bayesian residual diagnostics


create_bayesian_residual_diagnostics <- function(
    fit,
    data,
    model_name
) {
  
  # Posterior draws of fitted log volumes -
  # rows = posterior draws, columns = observations
  mu_log_draws <- as_draws_matrix(
    fit$draws("mu_log_out")
  )
  
  # Posterior draws of residual SD 
  # we convert to an ordinary numeric vector
  sigma_draws <- as.numeric(
    as_draws_matrix(fit$draws("sigma"))[, 1]
  )
  
  stopifnot(
    ncol(mu_log_draws) == nrow(data),
    nrow(mu_log_draws) == length(sigma_draws)
  )
  
  # Observed log volumes
  observed_log <- log(data$volume)
  
  # Residual for every posterior draw and observation
  log_residual_draws <- sweep(
    mu_log_draws,
    MARGIN = 2,
    STATS = observed_log,
    FUN = function(mu, observed) observed - mu
  )
  
  standardised_residual_draws <- sweep(
    log_residual_draws,
    MARGIN = 1,
    STATS = sigma_draws,
    FUN = "/"
  )
  
  # summarise each observation's posterior residual distribution
  residual_summary <- data %>%
    mutate(
      model = model_name,
      
      fitted_log = colMeans(mu_log_draws),
      
      residual_median = apply(
        standardised_residual_draws,
        2,
        median
      ),
      
      residual_lower = apply(
        standardised_residual_draws,
        2,
        quantile,
        probs = 0.025
      ),
      
      residual_upper = apply(
        standardised_residual_draws,
        2,
        quantile,
        probs = 0.975
      ),
      
      probability_positive = colMeans(
        standardised_residual_draws > 0
      )
    )
  
  return(residual_summary)
}
bayesian_residuals_exp <-
  create_bayesian_residual_diagnostics(
    fit = fit_exp,
    data = control_data,
    model_name = "Exponential"
  )

bayesian_residuals_logistic <-
  create_bayesian_residual_diagnostics(
    fit = fit_logistic,
    data = control_data,
    model_name = "Logistic"
  )

bayesian_residuals_gompertz <-
  create_bayesian_residual_diagnostics(
    fit = fit_gompertz,
    data = control_data,
    model_name = "Gompertz"
  )

bayesian_residuals_all <- bind_rows(
  bayesian_residuals_exp,
  bayesian_residuals_logistic,
  bayesian_residuals_gompertz
)
bayesian_residual_plot <- ggplot(
  bayesian_residuals_all,
  aes(
    x = fitted_log,
    y = residual_median,
    colour = factor(mouse_id)
  )
) +
  geom_errorbar(
    aes(
      ymin = residual_lower,
      ymax = residual_upper
    ),
    width = 0,
    alpha = 0.35
  ) +
  geom_point(alpha = 0.8) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = c(-2, 2),
    linetype = "dotted"
  ) +
  facet_wrap(~ model) +
  labs(
    x = "Posterior mean log fitted tumour volume",
    y = "Posterior median standardised residual",
    colour = "Mouse",
    title = "Bayesian residual diagnostics for untreated growth models"
  ) +
  theme_minimal()

bayesian_residual_plot

ggsave(
  filename = "figures/bayesian_residual_diagnostics.png",
  plot = bayesian_residual_plot,
  width = 10,
  height = 7,
  dpi = 300
)

# Standardised residuals over time

residual_time_plot <- ggplot(
  residuals_all_models,
  aes(
    x = time_days,
    y = standardised_log_residual,
    colour = factor(mouse_id)
  )
) +
  geom_point(alpha = 0.7) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = c(-2, 2),
    linetype = "dotted"
  ) +
  facet_wrap(~ model) +
  labs(
    x = "Time (days)",
    y = "Standardised log residual",
    colour = "Mouse",
    title = "Standardised residuals over time"
  ) +
  theme_minimal()

residual_time_plot

ggsave(
  filename = "figures/residuals_vs_fitted.png",
  plot = residual_fitted_plot,
  width = 10,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "figures/residuals_vs_time.png",
  plot = residual_time_plot,
  width = 10,
  height = 6,
  dpi = 300
)
# Observed versus fitted tumour volumes

observed_fitted_plot <- ggplot(
  residuals_all_models,
  aes(
    x = fitted_volume,
    y = volume,
    colour = factor(mouse_id)
  )
) +
  geom_point(alpha = 0.7) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  facet_wrap(~ model) +
  coord_equal() +
  labs(
    x = "Posterior mean fitted tumour volume",
    y = "Observed tumour volume",
    colour = "Mouse",
    title = "Observed versus fitted tumour volumes"
  ) +
  theme_minimal()

observed_fitted_plot

ggsave(
  filename = "figures/observed_vs_fitted.png",
  plot = observed_fitted_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# Individual observed and fitted trajectories

create_individual_fit_plot <- function(
    data,
    selected_model
) {
  
  data %>%
    filter(model == selected_model) %>%
    ggplot(aes(x = time_days)) +
    geom_point(
      aes(y = volume),
      alpha = 0.7,
      size = 1.5
    ) +
    geom_line(
      aes(
        y = fitted_volume,
        group = mouse_id
      ),
      linewidth = 0.8,
      linetype = "dashed"
    ) +
    facet_wrap(
      ~ mouse_id,
      ncol = 5,
      scales = "free_y"
    ) +
    labs(
      x = "Time (days)",
      y = expression("Tumour volume (mm"^3*")"),
      title = paste(
        "Observed and fitted untreated tumour trajectories:",
        selected_model,
        "model"
      )
    ) +
    theme_minimal() +
    theme(
      panel.spacing = grid::unit(0.8, "lines"),
      strip.text = element_text(size = 9)
    )
}

individual_fit_exp <- create_individual_fit_plot(
  data = residuals_all_models,
  selected_model = "Exponential"
)

individual_fit_gompertz <- create_individual_fit_plot(
  data = residuals_all_models,
  selected_model = "Gompertz"
)

individual_fit_logistic <- create_individual_fit_plot(
  data = residuals_all_models,
  selected_model = "Logistic"
)

individual_fit_exp
individual_fit_gompertz
individual_fit_logistic

# Numerical residual summaries

residual_summary <- residuals_all_models %>%
  group_by(model) %>%
  summarise(
    mean_log_residual =
      mean(log_residual),
    
    sd_log_residual =
      sd(log_residual),
    
    mean_absolute_log_residual =
      mean(abs(log_residual)),
    
    root_mean_squared_log_residual =
      sqrt(mean(log_residual^2)),
    
    maximum_absolute_log_residual =
      max(abs(log_residual)),
    
    .groups = "drop"
  )

print(
  residual_summary,
  width = Inf
)

# Prior-posterior density plots on fitted model scale

create_prior_posterior_plot <- function(
    fit,
    model_name,
    prior_logr_mean,
    seed = 123
) {
  
  set.seed(seed)
  
  # Exponential

  if (model_name == "Exponential") {
    
    posterior_draws <- as_draws_matrix(
      fit$draws(
        c(
          "logV0_pop",
          "logr_pop",
          "sigma_logV0",
          "sigma_logr",
          "sigma"
        )
      )
    )
    
    n_draws <- nrow(posterior_draws)
    
    prior_data <- data.frame(
      value = c(
        rnorm(
          n_draws,
          log(30),
          0.5
        ),
        
        rnorm(
          n_draws,
          prior_logr_mean,
          0.5
        ),
        
        abs(
          rnorm(
            n_draws,
            0,
            0.5
          )
        ),
        
        abs(
          rnorm(
            n_draws,
            0,
            0.5
          )
        ),
        
        abs(
          rnorm(
            n_draws,
            0,
            0.5
          )
        )
      ),
      
      parameter = rep(
        c(
          "log population initial volume",
          "log population growth rate",
          "Between-mouse variability in V0",
          "Between-mouse variability in growth rate",
          "Residual SD"
        ),
        each = n_draws
      ),
      
      distribution = "Prior"
    )
    
    posterior_data <- data.frame(
      value = c(
        posterior_draws[, "logV0_pop"],
        posterior_draws[, "logr_pop"],
        posterior_draws[, "sigma_logV0"],
        posterior_draws[, "sigma_logr"],
        posterior_draws[, "sigma"]
      ),
      
      parameter = rep(
        c(
          "log population initial volume",
          "log population growth rate",
          "Between-mouse variability in V0",
          "Between-mouse variability in growth rate",
          "Residual SD"
        ),
        each = n_draws
      ),
      
      distribution = "Posterior"
    )
  }
  
  # Logistic
  
  if (model_name == "Logistic") {
    
    posterior_draws <- as_draws_matrix(
      fit$draws(
        c(
          "logV0_pop",
          "logr_pop",
          "logK_pop",
          "sigma_logV0",
          "sigma_logr",
          "sigma_logK",
          "sigma"
        )
      )
    )
    
    n_draws <- nrow(posterior_draws)
    
    prior_data <- data.frame(
      value = c(
        rnorm(
          n_draws,
          log(30),
          0.5
        ),
        
        rnorm(
          n_draws,
          prior_logr_mean,
          0.5
        ),
        
        rnorm(
          n_draws,
          log(500),
          0.7
        ),
        
        abs(
          rnorm(
            n_draws,
            0,
            0.5
          )
        ),
        
        abs(
          rnorm(
            n_draws,
            0,
            0.5
          )
        ),
        
        abs(
          rnorm(
            n_draws,
            0,
            0.5
          )
        ),
        
        abs(
          rnorm(
            n_draws,
            0,
            0.5
          )
        )
      ),
      
      parameter = rep(
        c(
          "log population initial volume",
          "log population growth rate",
          "log population carrying capacity",
          "Between-mouse variability in V0",
          "Between-mouse variability in growth rate",
          "Between-mouse variability in K",
          "Residual SD"
        ),
        each = n_draws
      ),
      
      distribution = "Prior"
    )
    
    posterior_data <- data.frame(
      value = c(
        posterior_draws[, "logV0_pop"],
        posterior_draws[, "logr_pop"],
        posterior_draws[, "logK_pop"],
        posterior_draws[, "sigma_logV0"],
        posterior_draws[, "sigma_logr"],
        posterior_draws[, "sigma_logK"],
        posterior_draws[, "sigma"]
      ),
      
      parameter = rep(
        c(
          "log population initial volume",
          "log population growth rate",
          "log population carrying capacity",
          "Between-mouse variability in V0",
          "Between-mouse variability in growth rate",
          "Between-mouse variability in K",
          "Residual SD"
        ),
        each = n_draws
      ),
      
      distribution = "Posterior"
    )
  }
  
  # Gompertz
  
  if (model_name == "Gompertz") {
    
    posterior_draws <- as_draws_matrix(
      fit$draws(
        c(
          "logV0_pop",
          "logr_pop",
          "logK_pop",
          "sigma_logV0",
          "sigma_logr",
          "sigma_logK",
          "sigma"
        )
      )
    )
    
    n_draws <- nrow(posterior_draws)
    
    # Prior draws
    
    prior_data <- data.frame(
      value = c(
        
        rnorm(
          n_draws,
          mean = log(30),
          sd = 0.5
        ),
        
        rnorm(
          n_draws,
          mean = prior_logr_mean,
          sd = 0.5
        ),
        
        rnorm(
          n_draws,
          mean = log(500),
          sd = 0.7
        ),
        
        abs(
          rnorm(
            n_draws,
            mean = 0,
            sd = 0.5
          )
        ),
        
        abs(
          rnorm(
            n_draws,
            mean = 0,
            sd = 0.5
          )
        ),
        
        abs(
          rnorm(
            n_draws,
            mean = 0,
            sd = 0.5
          )
        ),
        
        abs(
          rnorm(
            n_draws,
            mean = 0,
            sd = 0.5
          )
        )
      ),
      
      parameter = rep(
        c(
          "log population initial volume",
          "log population growth rate",
          "log population carrying capacity",
          "Between-mouse variability in V0",
          "Between-mouse variability in growth rate",
          "Between-mouse variability in K",
          "Residual SD"
        ),
        each = n_draws
      ),
      
      distribution = "Prior"
    )

    # Posterior draws
    posterior_data <- data.frame(
      value = c(
        posterior_draws[, "logV0_pop"],
        posterior_draws[, "logr_pop"],
        posterior_draws[, "logK_pop"],
        posterior_draws[, "sigma_logV0"],
        posterior_draws[, "sigma_logr"],
        posterior_draws[, "sigma_logK"],
        posterior_draws[, "sigma"]
      ),
      
      parameter = rep(
        c(
          "log population initial volume",
          "log population growth rate",
          "log population carrying capacity",
          "Between-mouse variability in V0",
          "Between-mouse variability in growth rate",
          "Between-mouse variability in K",
          "Residual SD"
        ),
        each = n_draws
      ),
      
      distribution = "Posterior"
    )
  }

  # Combine them

  
  combined_data <- bind_rows(
    prior_data,
    posterior_data
  )
  
  combined_data$distribution <- factor(
    combined_data$distribution,
    levels = c(
      "Prior",
      "Posterior"
    )
  )
  
  ggplot(
    combined_data,
    aes(
      x = value,
      fill = distribution,
      colour = distribution
    )
  ) +
    geom_density(
      alpha = 0.30,
      adjust = 1
    ) +
    facet_wrap(
      ~ parameter,
      scales = "free",
      ncol = 2
    ) +
    labs(
      x = "Parameter value on fitted model scale",
      y = "Density",
      fill = "Distribution",
      colour = "Distribution",
      title = paste(
        "Prior and posterior distributions:",
        model_name,
        "model"
      )
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom"
    )
}

prior_posterior_exp_plot <- create_prior_posterior_plot(
  fit = fit_exp,
  model_name = "Exponential",
  prior_logr_mean = log(0.3)
)

prior_posterior_logistic_plot <- create_prior_posterior_plot(
  fit = fit_logistic,
  model_name = "Logistic",
  prior_logr_mean = log(0.4)
)

prior_posterior_gompertz_plot <- create_prior_posterior_plot(
  fit = fit_gompertz,
  model_name = "Gompertz",
  prior_logr_mean = log(0.2)
)

prior_posterior_exp_plot
prior_posterior_logistic_plot
prior_posterior_gompertz_plot

# Posterior distributions

# Exponential model

mcmc_dens(
  
  fit_exp$draws(
    
    c(
      
      "logV0_pop",
      
      "logr_pop",
      
      "sigma"
      
    )
    
  )
  
)

# Logistic model

mcmc_dens(
  fit_logistic$draws(
    c(
      "logV0_pop",
      "logr_pop",
      "logK_pop",
      "sigma_logK",
      "sigma"
    )
  )
)

# Gompertz model

mcmc_dens(
  fit_gompertz$draws(
    c(
      "logV0_pop",
      "logr_pop",
      "logK_pop",
      "sigma_logK",
      "sigma"
    )
  )
)

# Trace plots
# Exponential model

mcmc_trace(
  fit_exp$draws(
    c(
      "logV0_pop",
      "logr_pop",
      "sigma"
    )
  )
)

# Logistic model

mcmc_trace(
  fit_logistic$draws(
    c(
      "logV0_pop",
      "logr_pop",
      "logK_pop",
      "sigma_logK",
      "sigma"
    )
  )
)

#Gompertz model
mcmc_trace(
  fit_gompertz$draws(
    c(
      "logV0_pop",
      "logr_pop",
      "logK_pop",
      "sigma_logK",
      "sigma"
    )
  )
) +
  labs(
    title = "Trace plots for the Gompertz model",
    x = "Post-warm-up iteration",
    y = "Parameter value",
    colour = "Chain"
  ) +
  theme_bw(
    base_size = 14
  )

##################################################
# Observation-level PSIS-LOO cross-validation
##################################################

loo_exp <- loo(
  fit_exp$draws("log_lik")
)

loo_logistic <- loo(
  fit_logistic$draws("log_lik")
)

loo_gompertz <- loo(
  fit_gompertz$draws("log_lik")
)

loo_comparison_observation <- loo_compare(
  list(
    Exponential = loo_exp,
    Logistic = loo_logistic,
    Gompertz = loo_gompertz
  )
)

loo_comparison_observation

# Pareto-k diagnostics

pareto_k_table(loo_exp)
pareto_k_table(loo_logistic)
pareto_k_table(loo_gompertz)
print(loo_exp)
print(loo_logistic)
print(loo_gompertz)
print(loo_comparison_observation)
pareto_k_exp <- pareto_k_values(loo_exp)
pareto_k_logistic <- pareto_k_values(loo_logistic)
pareto_k_gompertz <- pareto_k_values(loo_gompertz)

# Posterior predictive checks
yrep_exp <- as_draws_matrix(
  fit_exp$draws("y_rep")
)

yrep_logistic <- as_draws_matrix(
  fit_logistic$draws("y_rep")
)

yrep_gompertz <- as_draws_matrix(
  fit_gompertz$draws("y_rep")
)

set.seed(123)

ppc_draw_ids <- sample(
  seq_len(nrow(yrep_exp)),
  size = min(100, nrow(yrep_exp))
)

# Log-scale density 

ppc_exp_log <- ppc_dens_overlay(
  y = log(control_data$volume),
  yrep = log(yrep_exp[ppc_draw_ids, , drop = FALSE])
) +
  ggtitle("Exponential") +
  xlab("Log tumour volume")

ppc_logistic_log <- ppc_dens_overlay(
  y = log(control_data$volume),
  yrep = log(yrep_logistic[ppc_draw_ids, , drop = FALSE])
) +
  ggtitle("Logistic") +
  xlab("Log tumour volume")

ppc_gompertz_log <- ppc_dens_overlay(
  y = log(control_data$volume),
  yrep = log(yrep_gompertz[ppc_draw_ids, , drop = FALSE])
) +
  ggtitle("Gompertz") +
  xlab("Log tumour volume")

ppc_exp_log
ppc_logistic_log
ppc_gompertz_log

# Combined posterior predictive checks
library(patchwork)

combined_ppc_plot <-
  ppc_exp_log +
  ppc_logistic_log +
  ppc_gompertz_log +
  plot_layout(
    ncol = 3,
    guides = "collect"
  ) +
  plot_annotation(
    title = "Posterior predictive checks for untreated tumour-growth models"
  ) &
  theme(
    legend.position = "bottom"
  )

combined_ppc_plot
dev.new(
  width = 13,
  height = 5
)

print(combined_ppc_plot)
ggsave(
  filename = "figures/combined_posterior_predictive_checks.pdf",
  plot = combined_ppc_plot,
  width = 13,
  height = 5
)

#raw-scale PPC - maybe for appendix
ppc_gompertz_raw <- ppc_dens_overlay(
  y = control_data$volume,
  yrep = yrep_gompertz[ppc_draw_ids, , drop = FALSE]
) +
  ggtitle("Posterior predictive check: Gompertz model, original scale") +
  xlab(expression("Tumour volume (mm"^3*")"))

ppc_gompertz_raw

# Observed and fitted trajectories by mouse

trajectory_all_models_plot <- ggplot(
  residuals_all_models,
  aes(x = time_days)
) +
  geom_point(
    aes(y = volume),
    alpha = 0.65
  ) +
  geom_line(
    aes(
      y = fitted_volume,
      group = interaction(model, mouse_id)
    ),
    linewidth = 0.8
  ) +
  facet_grid(
    model ~ mouse_id,
    scales = "free_y"
  ) +
  labs(
    x = "Time (days)",
    y = expression("Tumour volume (mm"^3*")"),
    title = "Observed and fitted untreated tumour trajectories"
  ) +
  theme_minimal()

trajectory_all_models_plot

# Time-binned posterior predictive check
# -maybe for appendix

control_data_ppc <- control_data %>%
  mutate(
    time_bin = cut(
      time_days,
      breaks = 5,
      include.lowest = TRUE
    )
  )

time_bin_id <- as.integer(control_data_ppc$time_bin)

ppc_time_gompertz <- ppc_stat_grouped(
  y = log(control_data_ppc$volume),
  yrep = log(yrep_gompertz[ppc_draw_ids, , drop = FALSE]),
  group = time_bin_id,
  stat = "median"
) +
  ggtitle("Time-binned posterior predictive check: Gompertz model") +
  xlab("Time bin") +
  ylab("Median log tumour volume")

ppc_time_gompertz


trajectory_gompertz_plot <- residuals_gompertz %>%
  ggplot(aes(x = time_days)) +
  geom_point(
    aes(y = volume),
    alpha = 0.65,
    size = 1.5
  ) +
  geom_line(
    aes(
      y = fitted_volume,
      group = mouse_id
    ),
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  facet_wrap(
    ~ mouse_id,
    ncol = 5,
    scales = "free_y"
  ) +
  labs(
    x = "Time (days)",
    y = expression("Tumour volume (mm"^3*")"),
    title = "Observed and fitted Gompertz trajectories by mouse"
  ) +
  theme_minimal()

trajectory_gompertz_plot

# Observation-level Pareto-k diagnostics

pareto_diagnostics <- control_data %>%
  mutate(
    observation_index = row_number(),
    k_exponential = pareto_k_values(loo_exp),
    k_logistic = pareto_k_values(loo_logistic),
    k_gompertz = pareto_k_values(loo_gompertz)
  ) %>%
  mutate(
    maximum_k = pmax(
      k_exponential,
      k_logistic,
      k_gompertz
    )
  )

pareto_diagnostics %>%
  filter(maximum_k > 0.7) %>%
  arrange(desc(maximum_k)) %>%
  select(
    observation_index,
    mouse,
    mouse_id,
    time_days,
    volume,
    k_exponential,
    k_logistic,
    k_gompertz,
    maximum_k
  )

# PSIS-LOO reliability check

loo_reliability_summary <- data.frame(
  model = c(
    "Exponential",
    "Logistic",
    "Gompertz"
  ),
  n_k_above_0_7 = c(
    sum(pareto_k_values(loo_exp) > 0.7),
    sum(pareto_k_values(loo_logistic) > 0.7),
    sum(pareto_k_values(loo_gompertz) > 0.7)
  ),
  n_k_above_1 = c(
    sum(pareto_k_values(loo_exp) >= 1),
    sum(pareto_k_values(loo_logistic) >= 1),
    sum(pareto_k_values(loo_gompertz) >= 1)
  ),
  maximum_k = c(
    max(pareto_k_values(loo_exp)),
    max(pareto_k_values(loo_logistic)),
    max(pareto_k_values(loo_gompertz))
  )
)

loo_reliability_summary

# Population-typical and new-mouse prediction:

population_gompertz_draws <-
  as_draws_matrix(
    fit_gompertz$draws(
      c(
        "logV0_pop",
        "logr_pop",
        "logK_pop",
        "sigma_logV0",
        "sigma_logr",
        "sigma_logK",
        "sigma"
      )
    )
  )

population_time_grid <-
  seq(
    from = 0,
    to = max(control_data$time_days),
    length.out = 200
  )

logV0_population_draws <-
  population_gompertz_draws[, "logV0_pop"]

logr_population_draws <-
  population_gompertz_draws[, "logr_pop"]

logK_population_draws <-
  population_gompertz_draws[
    ,
    "logK_pop"
  ]

sigma_logV0_draws <-
  population_gompertz_draws[, "sigma_logV0"]

sigma_logr_draws <-
  population_gompertz_draws[, "sigma_logr"]

sigma_logK_draws <-
  population_gompertz_draws[, "sigma_logK"]

sigma_residual_draws <-
  population_gompertz_draws[, "sigma"]

n_draws <-
  nrow(population_gompertz_draws)

V0_typical_draws <-
  exp(logV0_population_draws)

r_typical_draws <-
  exp(logr_population_draws)

K_typical_draws <-
  exp(
    logK_population_draws
  )

# 5. Population-typical trajectory


population_typical_matrix <-
  sapply(
    population_time_grid,
    function(t) {
      
      K_typical_draws *
        exp(
          log(V0_typical_draws / K_typical_draws) *
            exp(-r_typical_draws * t)
        )
    }
  )

population_typical_matrix <-
  as.matrix(population_typical_matrix)

# Simulate one new mouse for each posterior draw

set.seed(123)

z_V0_new <-
  rnorm(n_draws)

z_r_new <-
  rnorm(n_draws)

z_K_new <-
  rnorm(n_draws)

V0_new_mouse_draws <-
  exp(
    logV0_population_draws +
      sigma_logV0_draws * z_V0_new
  )

r_new_mouse_draws <-
  exp(
    logr_population_draws +
      sigma_logr_draws * z_r_new
  )

K_new_mouse_draws <-
  exp(
    logK_population_draws +
      sigma_logK_draws * z_K_new
  )

# Mean trajectory for the new mouse

new_mouse_mean_matrix <-
  sapply(
    population_time_grid,
    function(t) {
      
      K_new_mouse_draws *
        exp(
          log(V0_new_mouse_draws / K_new_mouse_draws) *
            exp(-r_new_mouse_draws * t)
        )
    }
  )

new_mouse_mean_matrix <-
  as.matrix(new_mouse_mean_matrix)

# observations for the new mouse
new_mouse_predictive_matrix <-
  matrix(
    NA_real_,
    nrow = n_draws,
    ncol = length(population_time_grid)
  )

for (j in seq_along(population_time_grid)) {
  
  new_mouse_predictive_matrix[, j] <-
    rlnorm(
      n = n_draws,
      meanlog = log(new_mouse_mean_matrix[, j]),
      sdlog = sigma_residual_draws
    )
}


#  uncertainty

population_trajectory_summary <-
  data.frame(
    
    time_days = population_time_grid,
    
    typical_median =
      apply(
        population_typical_matrix,
        2,
        median
      ),
    
    typical_lower_95 =
      apply(
        population_typical_matrix,
        2,
        quantile,
        probs = 0.025
      ),
    
    typical_upper_95 =
      apply(
        population_typical_matrix,
        2,
        quantile,
        probs = 0.975
      ),
    
    new_mouse_median =
      apply(
        new_mouse_predictive_matrix,
        2,
        median
      ),
    
    new_mouse_lower_95 =
      apply(
        new_mouse_predictive_matrix,
        2,
        quantile,
        probs = 0.025
      ),
    
    new_mouse_upper_95 =
      apply(
        new_mouse_predictive_matrix,
        2,
        quantile,
        probs = 0.975
      )
  )


# population-typical and new-mouse predictive uncertainty
# plots

gompertz_population_plot <-
  ggplot(
    population_trajectory_summary,
    aes(x = time_days)
  ) +
  
  geom_ribbon(
    aes(
      ymin = new_mouse_lower_95,
      ymax = new_mouse_upper_95
    ),
    alpha = 0.18
  ) +
  
  geom_ribbon(
    aes(
      ymin = typical_lower_95,
      ymax = typical_upper_95
    ),
    alpha = 0.45
  ) +
  
  geom_line(
    aes(y = typical_median),
    linewidth = 1
  ) +
  
  geom_point(
    data = control_data,
    aes(
      x = time_days,
      y = volume
    ),
    inherit.aes = FALSE,
    alpha = 0.35
  ) +
  
  labs(
    x = "Time (days)",
    y = expression("Tumour volume (mm"^3*")"),
    title = "Population-level Gompertz predictions",
    subtitle =
      "Dark ribbon: population-typical uncertainty; light ribbon: prediction for a new mouse"
  ) +
  
  theme_minimal()

gompertz_population_plot

ggsave(
  filename =
    "figures/gompertz_population_predictions.pdf",
  plot =
    gompertz_population_plot,
  width =
    7.5,
  height =
    5.5
)

# Population-typical trajectories for all models

exp_population_draws <-
  as_draws_matrix(
    fit_exp$draws(
      c(
        "logV0_pop",
        "logr_pop"
      )
    )
  )

logistic_population_draws <-
  as_draws_matrix(
    fit_logistic$draws(
      c(
        "logV0_pop",
        "logr_pop",
        "logK_pop"
      )
    )
  )

gompertz_comparison_draws <-
  as_draws_matrix(
    fit_gompertz$draws(
      c(
        "logV0_pop",
        "logr_pop",
        "logK_pop"
      )
    )
  )

comparison_time_grid <-
  seq(
    from = 0,
    to = max(control_data$time_days),
    length.out = 200
  )

# Exponential population-typical trajectory

exp_V0_draws <-
  exp(
    exp_population_draws[
      ,
      "logV0_pop"
    ]
  )

exp_r_draws <-
  exp(
    exp_population_draws[
      ,
      "logr_pop"
    ]
  )

exp_trajectory_matrix <-
  sapply(
    comparison_time_grid,
    function(t) {
      
      exp_V0_draws *
        exp(
          exp_r_draws *
            t
        )
    }
  )

exp_trajectory_matrix <-
  as.matrix(
    exp_trajectory_matrix
  )

# Logistic 

logistic_V0_typical_draws <-
  exp(
    logistic_population_draws[
      ,
      "logV0_pop"
    ]
  )

logistic_r_typical_draws <-
  exp(
    logistic_population_draws[
      ,
      "logr_pop"
    ]
  )

logistic_K_typical_comparison_draws <-
  exp(
    logistic_population_draws[
      ,
      "logK_pop"
    ]
  )

logistic_trajectory_matrix <-
  sapply(
    comparison_time_grid,
    function(t) {
      
      logistic_K_typical_comparison_draws /
        (
          1 +
            (
              (
                logistic_K_typical_comparison_draws -
                  logistic_V0_typical_draws
              ) /
                logistic_V0_typical_draws
            ) *
            exp(
              -logistic_r_typical_draws *
                t
            )
        )
    }
  )

logistic_trajectory_matrix <-
  as.matrix(
    logistic_trajectory_matrix
  )


# Gompertz


gompertz_V0_typical_draws <-
  exp(
    gompertz_comparison_draws[, "logV0_pop"]
  )

gompertz_r_typical_comparison_draws <-
  exp(
    gompertz_comparison_draws[, "logr_pop"]
  )

gompertz_K_typical_comparison_draws <-
  exp(
    gompertz_comparison_draws[
      ,
      "logK_pop"
    ]
  )

gompertz_trajectory_matrix <-
  sapply(
    comparison_time_grid,
    function(t) {
      
      gompertz_K_typical_comparison_draws *
        exp(
          log(
            gompertz_V0_typical_draws /
              gompertz_K_typical_comparison_draws
          ) *
            exp(
              -gompertz_r_typical_comparison_draws * t
            )
        )
    }
  )

gompertz_trajectory_matrix <-
  as.matrix(gompertz_trajectory_matrix)


summarise_trajectory_matrix <- function(
    trajectory_matrix,
    model_name
) {
  
  data.frame(
    
    time_days =
      comparison_time_grid,
    
    model =
      model_name,
    
    median =
      apply(
        trajectory_matrix,
        2,
        median
      ),
    
    lower_95 =
      apply(
        trajectory_matrix,
        2,
        quantile,
        probs = 0.025
      ),
    
    upper_95 =
      apply(
        trajectory_matrix,
        2,
        quantile,
        probs = 0.975
      )
  )
}

trajectory_comparison_data <-
  bind_rows(
    
    summarise_trajectory_matrix(
      exp_trajectory_matrix,
      "Exponential"
    ),
    
    summarise_trajectory_matrix(
      logistic_trajectory_matrix,
      "Logistic"
    ),
    
    summarise_trajectory_matrix(
      gompertz_trajectory_matrix,
      "Gompertz"
    )
  )

trajectory_comparison_data$model <-
  factor(
    trajectory_comparison_data$model,
    levels = c(
      "Exponential",
      "Logistic",
      "Gompertz"
    )
  )

population_model_comparison_plot <-
  ggplot(
    trajectory_comparison_data,
    aes(
      x = time_days,
      y = median,
      colour = model,
      fill = model
    )
  ) +
  
  geom_ribbon(
    aes(
      ymin = lower_95,
      ymax = upper_95
    ),
    alpha = 0.12,
    colour = NA
  ) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_point(
    data = control_data,
    aes(
      x = time_days,
      y = volume
    ),
    inherit.aes = FALSE,
    alpha = 0.25
  ) +
  
  labs(
    x =
      "Time (days)",
    
    y =
      expression(
        "Tumour volume (mm"^3*")"
      ),
    
    colour =
      "Model",
    
    fill =
      "Model",
    
    title =
      "Population-typical untreated tumour-growth trajectories",
    
    subtitle =
      "Lines show posterior medians; ribbons show 95% credible intervals"
  ) +
  
  theme_minimal()

population_model_comparison_plot

ggsave(
  filename =
    "figures/population_model_comparison.pdf",
  
  plot =
    population_model_comparison_plot,
  
  width =
    7.5,
  
  height =
    5.5
)

print(
  residual_summary,
  width = Inf
)

write.csv(
  residual_summary,
  "results/residual_summary.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(
    loo_comparison_observation
  ),
  "results/loo_comparison_observation.csv",
  row.names = TRUE
)

write.csv(
  loo_reliability_summary,
  "results/loo_reliability_summary.csv",
  row.names = FALSE
)

# Untreated-control dataset summary

control_dataset_summary <- data.frame(
  number_of_observations = nrow(control_data),
  number_of_mice = n_distinct(control_data$mouse),
  minimum_time_days = min(control_data$time_days),
  maximum_time_days = max(control_data$time_days),
  minimum_volume = min(control_data$volume),
  maximum_volume = max(control_data$volume)
)

print(control_dataset_summary)

write.csv(
  control_dataset_summary,
  "results/control_dataset_summary.csv",
  row.names = FALSE
)



sink(
  "results/lomo_console_output.txt",
  split = TRUE
)

# Leave-One-Mouse-Out Cross-Validation


mod_exp_lomo <- cmdstan_model(
  "exp_control_mixed_loo.stan"
)

mod_logistic_lomo <- cmdstan_model(
  "logistic_control_mixed_loo.stan"
)

mod_gompertz_lomo <- cmdstan_model(
  "gompertz_control_mixed_loo.stan"
)

log_mean_exp <- function(x) {
  
  x_max <- max(x)
  
  x_max + log(
    mean(
      exp(x - x_max)
    )
  )
}

make_lomo_stan_data <- function(
    training_data,
    model_name
) {
  
  training_data <-
    training_data %>%
    mutate(
      mouse_id_lomo =
        as.integer(
          factor(mouse)
        )
    )
  
  stan_data <-
    list(
      
      N =
        nrow(
          training_data
        ),
      
      M =
        n_distinct(
          training_data$mouse_id_lomo
        ),
      
      mouse =
        training_data$mouse_id_lomo,
      
      time =
        training_data$time_days,
      
      volume =
        training_data$volume
    )
  
  
  return(
    list(
      stan_data =
        stan_data
    )
  )
}

score_held_out_mouse <- function(
    fit,
    held_out_data,
    model_name,
    random_seed
) {

  number_of_test_observations <-
    nrow(
      held_out_data
    )
  
  
  stopifnot(
    number_of_test_observations > 0,
    all(
      is.finite(
        held_out_data$time_days
      )
    ),
    all(
      is.finite(
        held_out_data$volume
      )
    ),
    all(
      held_out_data$volume > 0
    )
  )
  
  set.seed(
    random_seed
  )
  
  # EXPONENTIAL
  
  if (
    model_name ==
    "Exponential"
  ) {
    
    draws <-
      as_draws_matrix(
        fit$draws(
          c(
            "logV0_pop",
            "logr_pop",
            "sigma_logV0",
            "sigma_logr",
            "sigma"
          )
        )
      )
    
    
    number_of_draws <-
      nrow(
        draws
      )
    
    z_V0_new <-
      rnorm(
        number_of_draws
      )
    
    z_r_new <-
      rnorm(
        number_of_draws
      )
    
    V0_new <-
      exp(
        draws[
          ,
          "logV0_pop"
        ] +
          draws[
            ,
            "sigma_logV0"
          ] *
          z_V0_new
      )
    
    r_new <-
      exp(
        draws[
          ,
          "logr_pop"
        ] +
          draws[
            ,
            "sigma_logr"
          ] *
          z_r_new
      )
    
    sigma_new <-
      draws[
        ,
        "sigma"
      ]

    
    expected_volume <-
      matrix(
        NA_real_,
        nrow =
          number_of_draws,
        ncol =
          number_of_test_observations
      )
    
    
    for (
      j in seq_len(
        number_of_test_observations
      )
    ) {
      
      t_j <-
        held_out_data$time_days[
          j
        ]
      
      expected_volume[
        ,
        j
      ] <-
        V0_new *
        exp(
          r_new *
            t_j
        )
    }
  }

  # LOGISTIC 
  
  else if (
    model_name ==
    "Logistic"
  ) {
    
    draws <-
      as_draws_matrix(
        fit$draws(
          c(
            "logV0_pop",
            "logr_pop",
            "logK_pop",
            "sigma_logV0",
            "sigma_logr",
            "sigma_logK",
            "sigma"
          )
        )
      )
    
    
    number_of_draws <-
      nrow(
        draws
      )
    
    z_V0_new <-
      rnorm(
        number_of_draws
      )
    
    z_r_new <-
      rnorm(
        number_of_draws
      )
    
    z_K_new <-
      rnorm(
        number_of_draws
      )
  
    V0_new <-
      exp(
        draws[
          ,
          "logV0_pop"
        ] +
          draws[
            ,
            "sigma_logV0"
          ] *
          z_V0_new
      )
    
    r_new <-
      exp(
        draws[
          ,
          "logr_pop"
        ] +
          draws[
            ,
            "sigma_logr"
          ] *
          z_r_new
      )
    
    K_new <-
      exp(
        draws[
          ,
          "logK_pop"
        ] +
          draws[
            ,
            "sigma_logK"
          ] *
          z_K_new
      )
    
    sigma_new <-
      draws[
        ,
        "sigma"
      ]
    
    stopifnot(
      all(
        is.finite(
          V0_new
        )
      ),
      all(
        is.finite(
          r_new
        )
      ),
      all(
        is.finite(
          K_new
        )
      ),
      all(
        V0_new > 0
      ),
      all(
        r_new > 0
      ),
      all(
        K_new > 0
      )
    )
    
    expected_volume <-
      matrix(
        NA_real_,
        nrow =
          number_of_draws,
        ncol =
          number_of_test_observations
      )
    
    
    for (
      j in seq_len(
        number_of_test_observations
      )
    ) {
      
      t_j <-
        held_out_data$time_days[
          j
        ]
      
      expected_volume[
        ,
        j
      ] <-
        K_new /
        (
          1 +
            (
              (
                K_new -
                  V0_new
              ) /
                V0_new
            ) *
            exp(
              -r_new *
                t_j
            )
        )
    }
  }
  
  # GOMPERTZ 
  else if (
    model_name ==
    "Gompertz"
  ) {
    
    draws <-
      as_draws_matrix(
        fit$draws(
          c(
            "logV0_pop",
            "logr_pop",
            "logK_pop",
            "sigma_logV0",
            "sigma_logr",
            "sigma_logK",
            "sigma"
          )
        )
      )
    
    number_of_draws <-
      nrow(draws)
    
    z_V0_new <-
      rnorm(number_of_draws)
    
    z_r_new <-
      rnorm(number_of_draws)
    
    z_K_new <-
      rnorm(number_of_draws)
    
    V0_new <-
      exp(
        draws[, "logV0_pop"] +
          draws[, "sigma_logV0"] *
          z_V0_new
      )
    
    r_new <-
      exp(
        draws[, "logr_pop"] +
          draws[, "sigma_logr"] *
          z_r_new
      )
    
    K_new <-
      exp(
        draws[, "logK_pop"] +
          draws[, "sigma_logK"] *
          z_K_new
      )
    
    sigma_new <-
      draws[, "sigma"]
    
    stopifnot(
      all(is.finite(V0_new)),
      all(is.finite(r_new)),
      all(is.finite(K_new)),
      all(V0_new > 0),
      all(r_new > 0),
      all(K_new > 0)
    )
    
    expected_volume <-
      matrix(
        NA_real_,
        nrow = number_of_draws,
        ncol = number_of_test_observations
      )
    
    for (
      j in seq_len(
        number_of_test_observations
      )
    ) {
      
      t_j <-
        held_out_data$time_days[j]
      
      expected_volume[, j] <-
        K_new *
        exp(
          log(V0_new / K_new) *
            exp(-r_new * t_j)
        )
    }
  }

  stopifnot(
    all(
      is.finite(
        expected_volume
      )
    ),
    all(
      expected_volume > 0
    ),
    all(
      is.finite(
        sigma_new
      )
    ),
    all(
      sigma_new > 0
    )
  )
  
  
  joint_log_likelihood <-
    rep(
      0,
      number_of_draws
    )
  
  
  for (
    j in seq_len(
      number_of_test_observations
    )
  ) {
    
    joint_log_likelihood <-
      joint_log_likelihood +
      dlnorm(
        held_out_data$volume[
          j
        ],
        meanlog =
          log(
            expected_volume[
              ,
              j
            ]
          ),
        sdlog =
          sigma_new,
        log =
          TRUE
      )
  }

  
  mouse_elpd <-
    log_mean_exp(
      joint_log_likelihood
    )

  
  data.frame(
    
    mouse =
      unique(
        held_out_data$mouse
      ),
    
    number_of_observations =
      number_of_test_observations,
    
    elpd_mouse =
      mouse_elpd
  ) 
  
}

run_lomo_cv <- function(
    model_object,
    model_name,
    full_data,
    output_file,
    seed = 123
) {
  
  mouse_values <-
    sort(
      unique(full_data$mouse)
    )
  
  ################################################
  # Resume from an existing partial result to prevent completed folds 
  # from being lost if the R stops.
  ################################################
  
  if (file.exists(output_file)) {
    
    lomo_results <- read.csv(
      output_file
    )
    
    completed_mice <-
      lomo_results$mouse
    
  } else {
    
    lomo_results <- data.frame()
    
    completed_mice <- numeric(0)
  }
  
  for (fold_index in seq_along(mouse_values)) {
    
    held_out_mouse <-
      mouse_values[fold_index]
    
    if (held_out_mouse %in% completed_mice) {
      
      message(
        model_name,
        ": mouse ",
        held_out_mouse,
        " already completed; skipping."
      )
      
      next
    }
    
    message(
      "\n",
      model_name,
      " LOMO fold ",
      fold_index,
      " of ",
      length(mouse_values),
      ": holding out mouse ",
      held_out_mouse
    )
    
    training_data <- full_data %>%
      filter(
        mouse != held_out_mouse
      )
    
    held_out_data <- full_data %>%
      filter(
        mouse == held_out_mouse
      ) %>%
      arrange(time_days)
    
    fold_data <-
      make_lomo_stan_data(
        training_data =
          training_data,
        
        model_name =
          model_name
      )
    
    
    fold_fit <- model_object$sample(
      data = fold_data$stan_data,
      chains = 4,
      parallel_chains = 4,
      iter_warmup = 500,
      iter_sampling = 1000,
      adapt_delta = 0.99,
      seed = seed + fold_index,
      refresh = 250
    )
    
    fold_diagnostics <-
      fold_fit$diagnostic_summary()
    
    if (
      any(fold_diagnostics$num_divergent > 0) ||
      any(fold_diagnostics$num_max_treedepth > 0)
    ) {
      warning(
        model_name,
        ", held-out mouse ",
        held_out_mouse,
        ": sampling diagnostics require inspection."
      )
    }
    
    fold_score <- score_held_out_mouse(
      fit =
        fold_fit,
      
      held_out_data =
        held_out_data,
      
      model_name =
        model_name,
      
      random_seed =
        seed +
        10000 +
        fold_index
    )
    
    fold_score$model <- model_name
    fold_score$fold <- fold_index
    
    lomo_results <- bind_rows(
      lomo_results,
      fold_score
    )
    
    write.csv(
      lomo_results,
      output_file,
      row.names = FALSE
    )
    
    saveRDS(
      fold_fit,
      file.path(
        "results/lomo_fits",
        paste0(
          tolower(model_name),
          "_held_out_mouse_",
          held_out_mouse,
          ".rds"
        )
      )
    )
    
    rm(fold_fit)
    gc()
  }
  
  lomo_results %>%
    arrange(mouse)
}
# LOMO smoke test

smoke_test_mouse <-
  sort(
    unique(
      control_data$mouse
    )
  )[1]

smoke_training_data <-
  control_data %>%
  filter(
    mouse !=
      smoke_test_mouse
  )

smoke_held_out_data <-
  control_data %>%
  filter(
    mouse ==
      smoke_test_mouse
  ) %>%
  arrange(
    time_days
  )


run_lomo_smoke_test <- function(
    model_object,
    model_name,
    seed = 999
) {
  

  
  smoke_fold_data <-
    make_lomo_stan_data(
      training_data =
        smoke_training_data,
      
      model_name =
        model_name
    )

  smoke_fit <-
    model_object$sample(
      
      data =
        smoke_fold_data$stan_data,
      
      chains =
        2,
      
      parallel_chains =
        2,
      
      iter_warmup =
        100,
      
      iter_sampling =
        100,
      
      adapt_delta =
        0.99,
      
      seed =
        seed,
      
      refresh =
        50
    )
  
  print(
    smoke_fit$diagnostic_summary()
  )

  smoke_score <-
    score_held_out_mouse(
      
      fit =
        smoke_fit,
      
      held_out_data =
        smoke_held_out_data,
      
      model_name =
        model_name,
      
      random_seed =
        seed +
        10000
    )

  
  print(
    smoke_score
  )
  
  stopifnot(
    nrow(
      smoke_score
    ) ==
      1,
    
    is.finite(
      smoke_score$elpd_mouse
    )
  )
  
  
  cat(
    "\n",
    model_name,
    " smoke test PASSED.\n",
    sep = ""
  )
  
  
  return(
    invisible(
      list(
        fit =
          smoke_fit,
        
        score =
          smoke_score
      )
    )
  )
}


##################################################
# Run smoke tests
##################################################

smoke_exp <-
  run_lomo_smoke_test(
    model_object =
      mod_exp_lomo,
    
    model_name =
      "Exponential",
    
    seed =
      901
  )


smoke_logistic <-
  run_lomo_smoke_test(
    model_object =
      mod_logistic_lomo,
    
    model_name =
      "Logistic",
    
    seed =
      902
  )


smoke_gompertz <-
  run_lomo_smoke_test(
    model_object =
      mod_gompertz_lomo,
    
    model_name =
      "Gompertz",
    
    seed =
      903
  )

# Gompertz LOMO-CV

lomo_gompertz <- run_lomo_cv(
  model_object = mod_gompertz_lomo,
  model_name = "Gompertz",
  full_data = control_data,
  output_file =
    "results/lomo_gompertz_individualK.csv",
  seed = 123
)

print(
  lomo_gompertz
)

# Logistic LOMO-CV


lomo_logistic <- run_lomo_cv(
  model_object = mod_logistic_lomo,
  model_name = "Logistic",
  full_data = control_data,
  output_file =
    "results/lomo_logistic_individualK.csv",
  seed = 123
)

# Exponential LOMO-CV

lomo_exponential <- run_lomo_cv(
  model_object = mod_exp_lomo,
  model_name = "Exponential",
  full_data = control_data,
  output_file =
    "results/lomo_exponential_individualK.csv",
  seed = 123
)


lomo_all_models <- bind_rows(
  lomo_exponential,
  lomo_logistic,
  lomo_gompertz
) %>%
  select(
    model,
    mouse,
    number_of_observations,
    elpd_mouse
  )

write.csv(
  lomo_all_models,
  "results/lomo_all_models.csv",
  row.names = FALSE
)

lomo_model_totals <- lomo_all_models %>%
  group_by(model) %>%
  summarise(
    number_of_mice = n(),
    total_elpd = sum(elpd_mouse),
    mean_elpd_per_mouse = mean(elpd_mouse),
    .groups = "drop"
  ) %>%
  arrange(desc(total_elpd))

print(
  lomo_model_totals
)

write.csv(
  lomo_model_totals,
  "results/lomo_model_totals.csv",
  row.names = FALSE
)


lomo_wide <- lomo_all_models %>%
  select(
    mouse,
    model,
    elpd_mouse
  ) %>%
  tidyr::pivot_wider(
    names_from = model,
    values_from = elpd_mouse
  )

best_model <-
  lomo_model_totals$model[1]

best_mouse_scores <-
  lomo_wide[[best_model]]

lomo_comparison <- lomo_model_totals %>%
  transmute(
    model = model,
    total_elpd = total_elpd,
    elpd_difference =
      total_elpd -
      max(total_elpd)
  )

lomo_comparison$se_difference <- sapply(
  lomo_comparison$model,
  function(current_model) {
    
    mouse_difference <-
      lomo_wide[[current_model]] -
      best_mouse_scores
    
    sqrt(
      nrow(lomo_wide) *
        var(mouse_difference)
    )
  }
)

print(
  lomo_comparison
)

write.csv(
  lomo_comparison,
  "results/lomo_model_comparison.csv",
  row.names = FALSE
)

sink()

report_file <- "results/lomo_results_report.txt"

sink(report_file)

sink()


# Raw untreated-control tumour trajectories


raw_control_trajectories <- ggplot(
  control_data,
  aes(
    x = time_days,
    y = volume,
    group = mouse_id
  )
) +
  geom_line(
    colour = "grey60",
    linewidth = 0.6
  ) +
  geom_point(
    size = 1.8
  ) +
  facet_wrap(
    ~ mouse_id,
    scales = "free_y"
  ) +
  labs(
    x = "Time (days)",
    y = expression("Tumour volume (mm"^3*")"),
    title = "Observed tumour trajectories for untreated control mice"
  ) +
  theme_minimal()

raw_control_trajectories

ggsave(
  "figures/raw_control_trajectories.pdf",
  raw_control_trajectories,
  width = 8,
  height = 7
)

# Complete posterior-density figures

make_density_plot <- function(
    draw_data,
    model_name
) {
  
  if (model_name == "Exponential") {
    
    parameter_levels <- c(
      "V0_pop",
      "r_pop",
      "sigma_logV0",
      "sigma_logr",
      "sigma"
    )
    
    parameter_labels <- c(
      "Initial tumour volume (V0)",
      "Growth rate (r)",
      "Between-mouse variability in V0",
      "Between-mouse variability in r",
      "Residual standard deviation"
    )
    
  } else if (model_name == "Logistic") {
    
    parameter_levels <- c(
      "V0_pop",
      "r_pop",
      "K_typical",
      "sigma_logV0",
      "sigma_logr",
      "sigma_logK",
      "sigma"
    )
    
    parameter_labels <- c(
      "Initial tumour volume (V0)",
      "Growth rate (r)",
      "Typical mouse carrying capacity (K)",
      "Between-mouse variability in V0",
      "Between-mouse variability in r",
      "Between-mouse variability in K",
      "Residual standard deviation"
    )
    
  } else if (model_name == "Gompertz") {
    
    parameter_levels <- c(
      "V0_pop",
      "r_pop",
      "K_typical",
      "sigma_logV0",
      "sigma_logr",
      "sigma_logK",
      "sigma"
    )
    
    parameter_labels <- c(
      "Initial tumour volume (V0)",
      "Growth rate (r)",
      "Typical mouse asymptotic volume (K)",
      "Between-mouse variability in V0",
      "Between-mouse variability in r",
      "Between-mouse variability in K",
      "Residual standard deviation"
    )
    
  } else {
    
    stop(
      "model_name must be Exponential, Logistic or Gompertz."
    )
  }
  
  long_data <- draw_data %>%
    select(
      all_of(
        parameter_levels
      )
    ) %>%
    tidyr::pivot_longer(
      cols = everything(),
      names_to = "parameter",
      values_to = "value"
    ) %>%
    mutate(
      parameter = factor(
        parameter,
        levels = parameter_levels,
        labels = parameter_labels
      )
    )
  
  ggplot(
    long_data,
    aes(
      x = value
    )
  ) +
    geom_density(
      fill = "grey70",
      alpha = 0.7,
      linewidth = 0.8
    ) +
    facet_wrap(
      ~ parameter,
      scales = "free",
      ncol = 4
    ) +
    labs(
      x = "Parameter value",
      y = "Posterior density",
      title = paste0(
        "Posterior distributions: ",
        tolower(model_name),
        " model"
      )
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      strip.text = element_text(
        face = "bold",
        size = 10
      ),
      axis.title = element_text(
        size = 11
      ),
      panel.grid.minor = element_blank()
    )
}

# Exponential 

exp_density_data <- as_draws_df(
  fit_exp$draws(
    c(
      "logV0_pop",
      "logr_pop",
      "sigma_logV0",
      "sigma_logr",
      "sigma"
    )
  )
) %>%
  transmute(
    V0_pop = exp(logV0_pop),
    r_pop = exp(logr_pop),
    sigma_logV0 = sigma_logV0,
    sigma_logr = sigma_logr,
    sigma = sigma
  )

posterior_density_exp <- make_density_plot(
  exp_density_data,
  "Exponential"
)

# Logistic 

logistic_density_draws <- as_draws_df(
  fit_logistic$draws(
    c(
      "logV0_pop",
      "logr_pop",
      "sigma_logV0",
      "sigma_logr",
      "sigma_logK",
      "sigma"
    )
  )
)

logistic_K_mouse_draws <- as_draws_matrix(
  fit_logistic$draws(
    "K_mouse"
  )
)

logistic_density_data <- logistic_density_draws %>%
  transmute(
    V0_pop = exp(logV0_pop),
    r_pop = exp(logr_pop),
    
    K_typical = apply(
      logistic_K_mouse_draws,
      1,
      median
    ),
    
    sigma_logV0 = sigma_logV0,
    sigma_logr = sigma_logr,
    sigma_logK = sigma_logK,
    sigma = sigma
  )

posterior_density_logistic <- make_density_plot(
  logistic_density_data,
  "Logistic"
)


# Gompertz 

gompertz_density_draws <- as_draws_df(
  fit_gompertz$draws(
    c(
      "logV0_pop",
      "logr_pop",
      "sigma_logV0",
      "sigma_logr",
      "sigma_logK",
      "sigma"
    )
  )
)

gompertz_K_mouse_draws <- as_draws_matrix(
  fit_gompertz$draws(
    "K_mouse"
  )
)

gompertz_density_data <- gompertz_density_draws %>%
  transmute(
    V0_pop = exp(logV0_pop),
    r_pop = exp(logr_pop),
    
    K_typical = apply(
      gompertz_K_mouse_draws,
      1,
      median
    ),
    
    sigma_logV0 = sigma_logV0,
    sigma_logr = sigma_logr,
    sigma_logK = sigma_logK,
    sigma = sigma
  )

posterior_density_gompertz <- make_density_plot(
  gompertz_density_data,
  "Gompertz"
)


posterior_density_exp
posterior_density_logistic
posterior_density_gompertz

dir.create(
  "figures",
  showWarnings = FALSE
)

ggsave(
  "figures/posterior_density_exponential.pdf",
  posterior_density_exp,
  width = 8,
  height = 5
)

ggsave(
  "figures/posterior_density_logistic.pdf",
  posterior_density_logistic,
  width = 10,
  height = 5.5
)

ggsave(
  "figures/posterior_density_gompertz.pdf",
  posterior_density_gompertz,
  width = 10,
  height = 5.5
)

# Original-scale posterior summaries - medians and 95% credible intervals


summarise_posterior_95 <- function(draw_data) {
  
  draw_data %>%
    as.data.frame() %>%
    summarise(
      across(
        everything(),
        list(
          median =
            ~ median(.x),
          
          lower95 =
            ~ quantile(
              .x,
              probs = 0.025
            ),
          
          upper95 =
            ~ quantile(
              .x,
              probs = 0.975
            )
        )
      )
    )
}

# Exponential 

exp_original_scale_draws <-
  fit_exp$draws(
    c(
      "logV0_pop",
      "logr_pop",
      "sigma_logV0",
      "sigma_logr",
      "sigma"
    )
  ) %>%
  as_draws_df() %>%
  transmute(
    
    V0_pop =
      exp(logV0_pop),
    
    r_pop =
      exp(logr_pop),
    
    sigma_logV0 =
      sigma_logV0,
    
    sigma_logr =
      sigma_logr,
    
    sigma =
      sigma
  )

exp_original_scale_95 <-
  summarise_posterior_95(
    exp_original_scale_draws
  )

# Logistic 

logistic_original_scale_draws <-
  fit_logistic$draws(
    c(
      "logV0_pop",
      "logr_pop",
      "sigma_logV0",
      "sigma_logr",
      "sigma_logK",
      "sigma"
    )
  ) %>%
  as_draws_df()

logistic_K_mouse_draws_summary <-
  as_draws_matrix(
    fit_logistic$draws(
      "K_mouse"
    )
  )


logistic_K_typical_draws <-
  apply(
    logistic_K_mouse_draws_summary,
    1,
    median
  )

logistic_original_scale_draws <-
  logistic_original_scale_draws %>%
  transmute(
    
    V0_pop =
      exp(logV0_pop),
    
    r_pop =
      exp(logr_pop),
    
    K_typical =
      logistic_K_typical_draws,
    
    sigma_logV0 =
      sigma_logV0,
    
    sigma_logr =
      sigma_logr,
    
    sigma_logK =
      sigma_logK,
    
    sigma =
      sigma
  )

logistic_original_scale_95 <-
  summarise_posterior_95(
    logistic_original_scale_draws
  )

# Gompertz 


gompertz_original_scale_draws <-
  fit_gompertz$draws(
    c(
      "logV0_pop",
      "logr_pop",
      "sigma_logV0",
      "sigma_logr",
      "sigma_logK",
      "sigma"
    )
  ) %>%
  as_draws_df()

gompertz_K_mouse_draws_summary <-
  as_draws_matrix(
    fit_gompertz$draws(
      "K_mouse"
    )
  )

gompertz_K_typical_draws <-
  apply(
    gompertz_K_mouse_draws_summary,
    1,
    median
  )


gompertz_original_scale_draws <-
  gompertz_original_scale_draws %>%
  transmute(
    
    V0_pop =
      exp(logV0_pop),
    
    r_pop =
      exp(logr_pop),
    
    K_typical =
      gompertz_K_typical_draws,
    
    sigma_logV0 =
      sigma_logV0,
    
    sigma_logr =
      sigma_logr,
    
    sigma_logK =
      sigma_logK,
    
    sigma =
      sigma
  )

gompertz_original_scale_95 <-
  summarise_posterior_95(
    gompertz_original_scale_draws
  )

print(
  exp_original_scale_95
)

print(
  logistic_original_scale_95
)

print(
  gompertz_original_scale_95
)

dir.create(
  "results",
  showWarnings = FALSE
)

write.csv(
  exp_original_scale_95,
  "results/exp_original_scale_95_summary.csv",
  row.names = FALSE
)

write.csv(
  logistic_original_scale_95,
  "results/logistic_original_scale_95_summary.csv",
  row.names = FALSE
)

write.csv(
  gompertz_original_scale_95,
  "results/gompertz_original_scale_95_summary.csv",
  row.names = FALSE
)


exp_original_scale_95_long <-
  exp_original_scale_95 %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "quantity",
    values_to = "value"
  )

logistic_original_scale_95_long <-
  logistic_original_scale_95 %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "quantity",
    values_to = "value"
  )

gompertz_original_scale_95_long <-
  gompertz_original_scale_95 %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "quantity",
    values_to = "value"
  )

print(
  exp_original_scale_95_long
)

print(
  logistic_original_scale_95_long
)

print(
  gompertz_original_scale_95_long
)
##################################################
# MCMC convergence diagnostics
##################################################

library(posterior)

extract_mcmc_diagnostics <- function(fit, model_name){
  
  summary <- fit$summary()
  
  parameters <- summary
  
  diagnostics <- data.frame(
    
    Model = model_name,
    
    Maximum_Rhat =
      max(parameters$rhat, na.rm = TRUE),
    
    Minimum_Bulk_ESS =
      min(parameters$ess_bulk, na.rm = TRUE),
    
    Minimum_Tail_ESS =
      min(parameters$ess_tail, na.rm = TRUE)
    
  )
  
  diagnostics
}

mcmc_diagnostics <- dplyr::bind_rows(
  
  extract_mcmc_diagnostics(
    fit_exp,
    "Exponential"
  ),
  
  extract_mcmc_diagnostics(
    fit_logistic,
    "Logistic"
  ),
  
  extract_mcmc_diagnostics(
    fit_gompertz,
    "Gompertz"
  )
  
)

print(mcmc_diagnostics)

write.csv(
  mcmc_diagnostics,
  "results/mcmc_diagnostics.csv",
  row.names = FALSE
)


fit_exp$diagnostic_summary()

fit_logistic$diagnostic_summary()

fit_gompertz$diagnostic_summary()

#E-BFMI
fit_exp$diagnostic_summary()$ebfmi

fit_logistic$diagnostic_summary()$ebfmi

fit_gompertz$diagnostic_summary()$ebfmi

##################################################
# Residual-distribution diagnostics
##################################################

# 1. Histogram of standardised log residuals


residual_histogram_plot <- ggplot(
  residuals_all_models,
  aes(x = standardised_log_residual)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 15,
    fill = "grey80",
    colour = "black",
    linewidth = 0.3
  ) +
  stat_function(
    fun = dnorm,
    args = list(
      mean = 0,
      sd = 1
    ),
    linewidth = 0.9,
    linetype = "dashed"
  ) +
  facet_wrap(
    ~ model,
    nrow = 1
  ) +
  labs(
    x = "Standardised log residual",
    y = "Density",
    title = "Distribution of standardised log residuals",
    subtitle = "Dashed curves show the standard Normal density"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    panel.grid.minor = element_blank()
  )

residual_histogram_plot

# 2. Normal Q-Q plots of standardised log residuals


residual_qq_plot <- ggplot(
  residuals_all_models,
  aes(sample = standardised_log_residual)
) +
  stat_qq(
    alpha = 0.7,
    size = 1.5
  ) +
  stat_qq_line(
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  facet_wrap(
    ~ model,
    nrow = 1
  ) +
  labs(
    x = "Theoretical Normal quantiles",
    y = "Observed residual quantiles",
    title = "Normal Q-Q plots of standardised log residuals"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    panel.grid.minor = element_blank()
  )

residual_qq_plot

ggsave(
  filename = "figures/standardised_residual_histograms.pdf",
  plot = residual_histogram_plot,
  width = 8,
  height = 4
)

ggsave(
  filename = "figures/standardised_residual_qq_plots.pdf",
  plot = residual_qq_plot,
  width = 8,
  height = 4
)


residual_summary_complete <- residuals_all_models %>%
  group_by(model) %>%
  summarise(
    mean_log_residual =
      mean(log_residual),
    
    sd_log_residual =
      sd(log_residual),
    
    mean_absolute_log_residual =
      mean(abs(log_residual)),
    
    root_mean_squared_log_residual =
      sqrt(mean(log_residual^2)),
    
    maximum_absolute_log_residual =
      max(abs(log_residual)),
    
    number_outside_2 =
      sum(abs(standardised_log_residual) > 2),
    
    proportion_outside_2 =
      mean(abs(standardised_log_residual) > 2),
    
    .groups = "drop"
  )

print(
  as.data.frame(residual_summary_complete),
  row.names = FALSE
)


write.csv(
  residual_summary_complete,
  "results/residual_summary_complete.csv",
  row.names = FALSE
)
##################################################
# Investigation of potentially unusual observations
##################################################

outlier_investigation <- residuals_all_models %>%
  mutate(
    absolute_standardised_residual =
      abs(standardised_log_residual)
  ) %>%
  arrange(
    model,
    desc(absolute_standardised_residual)
  )


# Observations outside +/- 2 residual SD

residual_outliers <- outlier_investigation %>%
  filter(
    absolute_standardised_residual > 2
  ) %>%
  select(
    model,
    mouse,
    mouse_id,
    time_days,
    volume,
    fitted_volume,
    log_residual,
    standardised_log_residual
  )

print(
  as.data.frame(residual_outliers),
  row.names = FALSE
)


# Largest residuals in each model


largest_residuals <- outlier_investigation %>%
  group_by(model) %>%
  slice_max(
    order_by = absolute_standardised_residual,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  select(
    model,
    mouse,
    mouse_id,
    time_days,
    volume,
    fitted_volume,
    standardised_log_residual
  )

print(
  as.data.frame(largest_residuals),
  row.names = FALSE
)

# Detailed investigation of mouse 5


mouse5_diagnostics <- residuals_all_models %>%
  filter(
    mouse_id == 5
  ) %>%
  arrange(
    model,
    time_days
  ) %>%
  select(
    model,
    mouse,
    mouse_id,
    time_days,
    volume,
    fitted_volume,
    log_residual,
    standardised_log_residual
  )

print(
  as.data.frame(mouse5_diagnostics),
  row.names = FALSE
)

# Mouse 5 observed versus fitted trajectories

mouse5_plot <- residuals_all_models %>%
  filter(
    mouse_id == 5
  ) %>%
  ggplot(
    aes(
      x = time_days
    )
  ) +
  geom_point(
    aes(
      y = volume
    ),
    size = 2.5
  ) +
  geom_line(
    aes(
      y = fitted_volume
    ),
    linewidth = 1
  ) +
  facet_wrap(
    ~ model
  ) +
  labs(
    x = "Time (days)",
    y = expression("Tumour volume (mm"^3*")"),
    title = "Observed and fitted trajectory for mouse 5"
  ) +
  theme_minimal()

mouse5_plot

# Mouse 5: fitted trajectory with posterior credible band


get_mu_quantiles <- function(fit, data, model_name) {
  
  mu_draws <- as_draws_matrix(fit$draws("mu_out"))
  
  stopifnot(ncol(mu_draws) == nrow(data))
  
  data %>%
    mutate(
      model = model_name,
      fitted_volume = colMeans(mu_draws),
      fitted_lower  = apply(mu_draws, 2, quantile, probs = 0.025),
      fitted_upper  = apply(mu_draws, 2, quantile, probs = 0.975)
    )
}

mu_band_exp       <- get_mu_quantiles(fit_exp,       control_data, "Exponential")
mu_band_logistic  <- get_mu_quantiles(fit_logistic,  control_data, "Logistic")
mu_band_gompertz  <- get_mu_quantiles(fit_gompertz,  control_data, "Gompertz")

mu_band_all <- bind_rows(mu_band_exp, mu_band_logistic, mu_band_gompertz) %>%
  mutate(model = factor(model, levels = c("Exponential", "Logistic", "Gompertz")))

mouse5_band <- mu_band_all %>%
  filter(mouse_id == 5) %>%
  mutate(is_day2 = time_days == 2)

mouse5_plot <- ggplot(mouse5_band, aes(x = time_days)) +
  geom_ribbon(
    aes(ymin = fitted_lower, ymax = fitted_upper, fill = model),
    alpha = 0.18
  ) +
  geom_line(
    aes(y = fitted_volume, colour = model),
    linewidth = 1
  ) +
  geom_point(
    data = filter(mouse5_band, !is_day2),
    aes(y = volume),
    size = 2.2, colour = "black"
  ) +
  geom_point(
    data = filter(mouse5_band, is_day2),
    aes(y = volume),
    size = 3.2, colour = "firebrick", shape = 17
  ) +
  facet_wrap(~ model) +
  scale_colour_manual(values = c("Exponential" = "#E69F00",
                                 "Logistic"    = "#009E73",
                                 "Gompertz"    = "#0072B2")) +
  scale_fill_manual(values = c("Exponential" = "#E69F00",
                               "Logistic"    = "#009E73",
                               "Gompertz"    = "#0072B2")) +
  labs(
    x = "Time (days)",
    y = expression("Tumour volume (mm"^3*")"),
    title = "Observed and fitted trajectory for mouse 5",
    subtitle = "Shaded band: 95% posterior credible interval; red triangle: day-2 observation"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

mouse5_plot

ggsave(
  filename = "figures/mouse5_observed_fitted.png",
  plot = mouse5_plot,
  width = 10, height = 4.5, dpi = 300
)
ggsave(
  filename = "figures/mouse5_observed_fitted.png",
  plot = mouse5_plot,
  width = 10, height = 4.5, dpi = 300
)

# Residual behaviour by mouse


mouse_residual_summary <- residuals_all_models %>%
  group_by(
    model,
    mouse,
    mouse_id
  ) %>%
  summarise(
    n_observations = n(),
    
    mean_absolute_standardised_residual =
      mean(
        abs(standardised_log_residual)
      ),
    
    maximum_absolute_standardised_residual =
      max(
        abs(standardised_log_residual)
      ),
    
    n_outside_2 =
      sum(
        abs(standardised_log_residual) > 2
      ),
    
    .groups = "drop"
  ) %>%
  arrange(
    model,
    desc(maximum_absolute_standardised_residual)
  )

print(
  as.data.frame(mouse_residual_summary),
  row.names = FALSE
)

# Investigate shape of pooled observed PPC density


control_data_ppc_shape <- control_data %>%
  mutate(
    time_group = cut(
      time_days,
      breaks = 4,
      include.lowest = TRUE
    )
  )

observed_density_by_time <- ggplot(
  control_data_ppc_shape,
  aes(
    x = log(volume),
    colour = time_group
  )
) +
  geom_density(
    linewidth = 1
  ) +
  labs(
    x = "Log tumour volume",
    y = "Density",
    colour = "Time interval",
    title = "Observed log tumour-volume distributions by time"
  ) +
  theme_minimal()

observed_density_by_time



# Log tumour volume versus time


volume_time_shape_plot <- ggplot(
  control_data,
  aes(
    x = time_days,
    y = log(volume),
    colour = factor(mouse_id)
  )
) +
  geom_point(
    alpha = 0.7
  ) +
  geom_line(
    aes(group = mouse_id),
    alpha = 0.4
  ) +
  labs(
    x = "Time (days)",
    y = "Log tumour volume",
    colour = "Mouse",
    title = "Observed log tumour volume over time"
  ) +
  theme_minimal()

volume_time_shape_plot

#histogram-ECDF
ggplot(
  control_data,
  aes(x = log(volume))
) +
  geom_histogram(
    bins = 15,
    boundary = 0
  ) +
  labs(
    x = "Log tumour volume",
    y = "Count",
    title = "Distribution of observed log tumour volumes"
  ) +
  theme_minimal()


#trace_plot
mcmc_trace(
  fit_gompertz$draws(
    c(
      "logV0_pop",
      "logr_pop",
      "logK_pop",
      "sigma_logK",
      "sigma"
    )
  )
) +
  labs(
    title = "Trace plots for the Gompertz model",
    x = "Post-warm-up iteration",
    y = "Parameter value",
    colour = "Chain"
  ) +
  theme_bw(
    base_size = 14
  )
trace_gompertz <- mcmc_trace(
  fit_gompertz$draws(
    c(
      "logV0_pop",
      "logr_pop",
      "logK_pop",
      "sigma_logK",
      "sigma"
    )
  )
) +
  labs(
    title = "Trace plots for the Gompertz model",
    x = "Post-warm-up iteration",
    y = "Parameter value",
    colour = "Chain"
  ) +
  theme_bw(base_size = 14)

trace_gompertz

ggsave(
  "figures/traceplot_gompertz.png",
  trace_gompertz,
  width = 9,
  height = 7,
  dpi = 300
)


# Residual-distribution diagnostics


##################################################
# 1. Histogram of standardised log residuals
# with standard Normal density overlay
##################################################

residual_histogram_plot <- ggplot(
  residuals_all_models,
  aes(
    x = standardised_log_residual
  )
) +
  geom_histogram(
    aes(
      y = after_stat(density)
    ),
    bins = 15,
    fill = "grey80",
    colour = "black",
    linewidth = 0.3
  ) +
  stat_function(
    fun = dnorm,
    args = list(
      mean = 0,
      sd = 1
    ),
    linewidth = 0.9,
    linetype = "dashed"
  ) +
  facet_wrap(
    ~ model,
    nrow = 1
  ) +
  labs(
    x = "Standardised log residual",
    y = "Density",
    title = "Distribution of standardised log residuals",
    subtitle = "Dashed curves show the standard Normal density"
  ) +
  theme_minimal(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    panel.grid.minor = element_blank()
  )

residual_histogram_plot


##################################################
# 2. Normal Q-Q plots of standardised log residuals
##################################################

residual_qq_plot <- ggplot(
  residuals_all_models,
  aes(
    sample = standardised_log_residual
  )
) +
  stat_qq(
    alpha = 0.7,
    size = 1.5
  ) +
  stat_qq_line(
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  facet_wrap(
    ~ model,
    nrow = 1
  ) +
  labs(
    x = "Theoretical Normal quantiles",
    y = "Observed residual quantiles",
    title = "Normal Q-Q plots of standardised log residuals"
  ) +
  theme_minimal(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    panel.grid.minor = element_blank()
  )

residual_qq_plot


ggsave(
  filename =
    "figures/standardised_residual_histograms.pdf",
  plot =
    residual_histogram_plot,
  width =
    8,
  height =
    4
)

ggsave(
  filename =
    "figures/standardised_residual_qq_plots.pdf",
  plot =
    residual_qq_plot,
  width =
    8,
  height =
    4
)


# Time-binned posterior predictive checks
# Clear for dissertation

control_data_ppc <- control_data %>%
  mutate(
    time_bin = cut(
      time_days,
      breaks = 5,
      include.lowest = TRUE
    ),
    time_bin_id = as.integer(time_bin)
  )

observed_time_summary <- control_data_ppc %>%
  group_by(
    time_bin_id,
    time_bin
  ) %>%
  summarise(
    median_time_days =
      median(time_days),
    
    observed_median_log_volume =
      median(log(volume)),
    
    n_observations =
      n(),
    
    .groups = "drop"
  )

print(observed_time_summary)

create_time_binned_ppc_summary <- function(
    yrep,
    model_name,
    data,
    draw_ids
) {
  
  yrep_selected <-
    yrep[
      draw_ids,
      ,
      drop = FALSE
    ]
  
  log_yrep_selected <-
    log(
      yrep_selected
    )
  
  predictive_bin_medians <-
    lapply(
      sort(
        unique(
          data$time_bin_id
        )
      ),
      function(current_bin) {
        
        observation_indices <-
          which(
            data$time_bin_id ==
              current_bin
          )
        
        apply(
          log_yrep_selected[
            ,
            observation_indices,
            drop = FALSE
          ],
          1,
          median
        )
      }
    )
  
  predictive_bin_medians <-
    do.call(
      cbind,
      predictive_bin_medians
    )
  
  predictive_summary <-
    data.frame(
      
      time_bin_id =
        sort(
          unique(
            data$time_bin_id
          )
        ),
      
      predictive_median =
        apply(
          predictive_bin_medians,
          2,
          median
        ),
      
      predictive_lower_95 =
        apply(
          predictive_bin_medians,
          2,
          quantile,
          probs = 0.025
        ),
      
      predictive_upper_95 =
        apply(
          predictive_bin_medians,
          2,
          quantile,
          probs = 0.975
        ),
      
      predictive_lower_50 =
        apply(
          predictive_bin_medians,
          2,
          quantile,
          probs = 0.25
        ),
      
      predictive_upper_50 =
        apply(
          predictive_bin_medians,
          2,
          quantile,
          probs = 0.75
        ),
      
      model =
        model_name
    )
  
  
  predictive_summary %>%
    left_join(
      observed_time_summary,
      by = "time_bin_id"
    )
}


time_ppc_exp <-
  create_time_binned_ppc_summary(
    yrep =
      yrep_exp,
    
    model_name =
      "Exponential",
    
    data =
      control_data_ppc,
    
    draw_ids =
      ppc_draw_ids
  )


time_ppc_logistic <-
  create_time_binned_ppc_summary(
    yrep =
      yrep_logistic,
    
    model_name =
      "Logistic",
    
    data =
      control_data_ppc,
    
    draw_ids =
      ppc_draw_ids
  )


time_ppc_gompertz <-
  create_time_binned_ppc_summary(
    yrep =
      yrep_gompertz,
    
    model_name =
      "Gompertz",
    
    data =
      control_data_ppc,
    
    draw_ids =
      ppc_draw_ids
  )


time_binned_ppc_data <-
  bind_rows(
    time_ppc_exp,
    time_ppc_logistic,
    time_ppc_gompertz
  )

time_binned_ppc_data$model <-
  factor(
    time_binned_ppc_data$model,
    levels = c(
      "Exponential",
      "Logistic",
      "Gompertz"
    )
  )

print(
  time_binned_ppc_data
)


combined_time_binned_ppc <-
  ggplot(
    time_binned_ppc_data,
    aes(
      x = median_time_days
    )
  ) +
  

geom_ribbon(
  aes(
    ymin =
      predictive_lower_95,
    
    ymax =
      predictive_upper_95
  ),
  alpha = 0.18
) +


geom_ribbon(
  aes(
    ymin =
      predictive_lower_50,
    
    ymax =
      predictive_upper_50
  ),
  alpha = 0.30
) +

geom_line(
  aes(
    y =
      predictive_median
  ),
  linewidth = 0.9
) +
  
  geom_point(
    aes(
      y =
        predictive_median
    ),
    size = 2
  ) +

geom_point(
  aes(
    y =
      observed_median_log_volume
  ),
  shape = 4,
  size = 3,
  stroke = 1.2
) +
  
  facet_wrap(
    ~ model,
    nrow = 1
  ) +
  
  labs(
    x =
      "Time (days)",
    
    y =
      "Median log tumour volume",
    
    title =
      "Time-binned posterior predictive checks",
    
    subtitle =
      "Crosses show observed medians; lines show posterior predictive medians"
  ) +
  
  theme_minimal(
    base_size = 13
  ) +
  
  theme(
    plot.title =
      element_text(
        face = "bold"
      ),
    
    panel.grid.minor =
      element_blank()
  )


combined_time_binned_ppc

ggsave(
  filename =
    "figures/combined_time_binned_posterior_predictive_checks.pdf",
  
  plot =
    combined_time_binned_ppc,
  
  width =
    10,
  
  height =
    4.5
)


combined_time_binned_ppc <-
  ggplot(
    time_binned_ppc_data,
    aes(
      x = median_time_days,
      colour = model,
      fill = model
    )
  ) +
  
  geom_ribbon(
    aes(
      ymin = predictive_lower_95,
      ymax = predictive_upper_95
    ),
    alpha = 0.15,
    colour = NA
  ) +
  
  geom_ribbon(
    aes(
      ymin = predictive_lower_50,
      ymax = predictive_upper_50
    ),
    alpha = 0.30,
    colour = NA
  ) +
  
  geom_line(
    aes(
      y = predictive_median
    ),
    linewidth = 1
  ) +
  
  geom_point(
    aes(
      y = predictive_median
    ),
    size = 2.3
  ) +
  
  geom_point(
    data = time_binned_ppc_data,
    mapping = aes(
      x = median_time_days,
      y = observed_median_log_volume
    ),
    inherit.aes = FALSE,
    shape = 4,
    size = 3.2,
    stroke = 1.3,
    colour = "black"
  ) +
  
  facet_wrap(
    ~ model,
    nrow = 1
  ) +
  
  labs(
    x = "Time (days)",
    y = "Median log tumour volume"
  ) +
  
  theme_minimal(
    base_size = 13
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

combined_time_binned_ppc

dir.create(
  "figures",
  showWarnings = FALSE,
  recursive = TRUE
)

ggsave(
  filename =
    "figures/combined_time_binned_posterior_predictive_checks.pdf",
  plot =
    combined_time_binned_ppc,
  width =
    10,
  height =
    4.5,
  units =
    "in"
)


# Mouse 5: numerical credible-interval diagnostics


mouse5_interval_diagnostics <- mouse5_band %>%
  mutate(
    residual_volume = volume - fitted_volume,
    
    absolute_residual_volume = abs(residual_volume),
    
    relative_residual_percent =
      100 * residual_volume / volume,
    
    inside_95_CrI =
      volume >= fitted_lower &
      volume <= fitted_upper,
    
    below_95_CrI =
      volume < fitted_lower,
    
    above_95_CrI =
      volume > fitted_upper,
    
    distance_outside_CrI = case_when(
      volume < fitted_lower ~ volume - fitted_lower,
      volume > fitted_upper ~ volume - fitted_upper,
      TRUE ~ 0
    )
  ) %>%
  select(
    model,
    mouse,
    mouse_id,
    time_days,
    volume,
    fitted_volume,
    fitted_lower,
    fitted_upper,
    residual_volume,
    absolute_residual_volume,
    relative_residual_percent,
    inside_95_CrI,
    below_95_CrI,
    above_95_CrI,
    distance_outside_CrI
  ) %>%
  arrange(
    model,
    time_days
  )

print(
  as.data.frame(mouse5_interval_diagnostics),
  row.names = FALSE
)


# Mouse 5: interval diagnostic summary by model

mouse5_interval_summary <- mouse5_interval_diagnostics %>%
  group_by(model) %>%
  summarise(
    n_observations = n(),
    
    n_inside_95_CrI =
      sum(inside_95_CrI),
    
    n_outside_95_CrI =
      sum(!inside_95_CrI),
    
    proportion_inside_95_CrI =
      mean(inside_95_CrI),
    
    maximum_absolute_residual =
      max(absolute_residual_volume),
    
    time_of_maximum_absolute_residual =
      time_days[
        which.max(absolute_residual_volume)
      ],
    
    .groups = "drop"
  )

print(
  as.data.frame(mouse5_interval_summary),
  row.names = FALSE
)


# Mouse 5: numerical fitted-trajectory
# credible-interval diagnostics

mouse5_credible_diagnostics <- mouse5_band %>%
  mutate(
    inside_95_credible_interval =
      volume >= fitted_lower &
      volume <= fitted_upper,
    
    below_95_credible_interval =
      volume < fitted_lower,
    
    above_95_credible_interval =
      volume > fitted_upper,
    
    distance_from_credible_interval = case_when(
      volume < fitted_lower ~ volume - fitted_lower,
      volume > fitted_upper ~ volume - fitted_upper,
      TRUE ~ 0
    )
  ) %>%
  select(
    model,
    mouse,
    mouse_id,
    time_days,
    volume,
    fitted_volume,
    fitted_lower,
    fitted_upper,
    inside_95_credible_interval,
    below_95_credible_interval,
    above_95_credible_interval,
    distance_from_credible_interval
  ) %>%
  arrange(
    model,
    time_days
  )

print(
  as.data.frame(mouse5_credible_diagnostics),
  row.names = FALSE
)

# Mouse 5: posterior predictive interval

get_predictive_quantiles <- function(
    fit,
    data,
    model_name
) {
  
  yrep_draws <- as_draws_matrix(
    fit$draws("y_rep")
  )
  
  stopifnot(
    ncol(yrep_draws) == nrow(data)
  )
  
  data %>%
    mutate(
      model = model_name,
      
      predictive_median =
        apply(
          yrep_draws,
          2,
          median
        ),
      
      predictive_lower =
        apply(
          yrep_draws,
          2,
          quantile,
          probs = 0.025
        ),
      
      predictive_upper =
        apply(
          yrep_draws,
          2,
          quantile,
          probs = 0.975
        )
    )
}


predictive_band_exp <-
  get_predictive_quantiles(
    fit_exp,
    control_data,
    "Exponential"
  )

predictive_band_logistic <-
  get_predictive_quantiles(
    fit_logistic,
    control_data,
    "Logistic"
  )

predictive_band_gompertz <-
  get_predictive_quantiles(
    fit_gompertz,
    control_data,
    "Gompertz"
  )


predictive_band_all <-
  bind_rows(
    predictive_band_exp,
    predictive_band_logistic,
    predictive_band_gompertz
  ) %>%
  mutate(
    model = factor(
      model,
      levels = c(
        "Exponential",
        "Logistic",
        "Gompertz"
      )
    )
  )

# Mouse 5: posterior predictive diagnostics

mouse5_predictive_diagnostics <-
  predictive_band_all %>%
  filter(
    mouse_id == 5
  ) %>%
  mutate(
    inside_95_predictive_interval =
      volume >= predictive_lower &
      volume <= predictive_upper,
    
    below_95_predictive_interval =
      volume < predictive_lower,
    
    above_95_predictive_interval =
      volume > predictive_upper,
    
    distance_from_predictive_interval =
      case_when(
        volume < predictive_lower ~
          volume - predictive_lower,
        
        volume > predictive_upper ~
          volume - predictive_upper,
        
        TRUE ~ 0
      )
  ) %>%
  select(
    model,
    mouse,
    mouse_id,
    time_days,
    volume,
    predictive_median,
    predictive_lower,
    predictive_upper,
    inside_95_predictive_interval,
    below_95_predictive_interval,
    above_95_predictive_interval,
    distance_from_predictive_interval
  ) %>%
  arrange(
    model,
    time_days
  )

print(
  as.data.frame(
    mouse5_predictive_diagnostics
  ),
  row.names = FALSE
)

# Mouse 5: combine credible and predictive
# interval diagnostics

mouse5_interval_comparison <-
  mouse5_credible_diagnostics %>%
  select(
    model,
    mouse,
    mouse_id,
    time_days,
    volume,
    fitted_volume,
    fitted_lower,
    fitted_upper,
    inside_95_credible_interval
  ) %>%
  left_join(
    mouse5_predictive_diagnostics %>%
      select(
        model,
        mouse,
        mouse_id,
        time_days,
        predictive_median,
        predictive_lower,
        predictive_upper,
        inside_95_predictive_interval
      ),
    by = c(
      "model",
      "mouse",
      "mouse_id",
      "time_days"
    )
  ) %>%
  arrange(
    model,
    time_days
  )

print(
  as.data.frame(
    mouse5_interval_comparison
  ),
  row.names = FALSE
)


mouse5_interval_summary <-
  mouse5_interval_comparison %>%
  group_by(model) %>%
  summarise(
    n_observations = n(),
    
    n_outside_95_credible_interval =
      sum(
        !inside_95_credible_interval
      ),
    
    n_outside_95_predictive_interval =
      sum(
        !inside_95_predictive_interval
      ),
    
    .groups = "drop"
  )

print(
  as.data.frame(
    mouse5_interval_summary
  ),
  row.names = FALSE
)