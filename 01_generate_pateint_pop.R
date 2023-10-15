library(tidyverse)

# function to simulate liver biopsy staging for a patient with given fibrosis area
simulate_biopsy <- function(simulated_data) {
  
  left_join(simulated_data, prob_data) %>%
    mutate(
      est_stage = 
        case_when(
          throw_est < prob_f0 ~ "F0",
          throw_est > prob_f0 & throw_est < prob_f0 + prob_f1 ~ "F1",
          throw_est > prob_f0 + prob_f1 & throw_est < prob_f0 + prob_f1 + prob_f2 ~ "F2",
          throw_est > prob_f0 + prob_f1 + prob_f2 & throw_est < prob_f0 + prob_f1 + prob_f2 + prob_f3 ~ "F3",
          throw_est > prob_f0 + prob_f1 + prob_f2 + prob_f3 ~ "F4",
          TRUE ~ "error"
        ),
      est_stage = as_factor(est_stage),
      est_stage = fct_relevel(est_stage, c("F0", "F1", "F2", "F3", "F4")),
      true_stage = 
        case_when(
          throw_true < prob_f0 ~ "F0",
          throw_true > prob_f0 & throw_true < prob_f0 + prob_f1 ~ "F1",
          throw_true > prob_f0 + prob_f1 & throw_true < prob_f0 + prob_f1 + prob_f2 ~ "F2",
          throw_true > prob_f0 + prob_f1 + prob_f2 & throw_true < prob_f0 + prob_f1 + prob_f2 + prob_f3 ~ "F3",
          throw_true > prob_f0 + prob_f1 + prob_f2 + prob_f3 ~ "F4",
          TRUE ~ "error"
        )
    ) %>%
    mutate(
      advanced = 
        if_else(true_stage == "F3" | true_stage == "F4", "advanced", "early"),
      significant =
        if_else(true_stage == "F2" | true_stage == "F3" | true_stage == "F4", "significant", "early"),
      cpa = as.numeric(cpa)
    ) 
  
}

# generate distribution for test population, n = 100,000
n_sim <- 200000

## plot distribution of fibrosis
dist_plot <-
  ggplot() +
  geom_density(aes(x = rbeta(n_sim, 1.5, 4) * 30)) +
  theme_classic() + 
  xlab("Collagen proportionate area (%)") +
  ylab("Probability density") +
  ggtitle("Estimated distribution of fibrosis in cohorts enrolled in biomarkler studies")


## generate 100,000 simulated patients across that distribution
sim_biopsy_grey <-
  tibble(
    cpa = round(rbeta(n_sim, 1.5, 4) * 30, 1),
    throw_true = runif(n_sim),
    throw_est = runif(n_sim)
  ) %>%
  mutate(cpa = as.character(cpa)) %>%
  filter(cpa > 1) %>%
  slice_head(n = 100000)

grey <- simulate_biopsy(sim_biopsy_grey)

# write_csv(grey, "grey_patients.csv") # not run, included in repo
