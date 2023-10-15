library(tidyverse)
library(pROC)
library(patchwork)


#### functions ####
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

# get biopsy performance characteristics from population
## auc for advanced fibrosis
roc_advanced <- function(x) {
  
  x %>% 
    mutate(
      est_stage = 
        case_when(
          throw_est < prob_f0 ~ 0,
          throw_est > prob_f0 & throw_est < prob_f0 + prob_f1 ~ 1,
          throw_est > prob_f0 + prob_f1 & throw_est < prob_f0 + prob_f1 + prob_f2 ~ 2,
          throw_est > prob_f0 + prob_f1 + prob_f2 & throw_est < prob_f0 + prob_f1 + prob_f2 + prob_f3 ~ 3,
          throw_est > prob_f0 + prob_f1 + prob_f2 + prob_f3 ~ 4,
          TRUE ~ 5
        )
    ) %>% 
    filter(est_stage != 5) %>%
    roc(advanced, est_stage, levels = c("early", "advanced"))
  
}

## sensitivity and specificity for advanced fibrosis
roc_cuts_advanced <- function(x) {
  
  x %>% 
    mutate(
      est_stage = 
        case_when(
          throw_est < prob_f0 ~ 0,
          throw_est > prob_f0 & throw_est < prob_f0 + prob_f1 ~ 1,
          throw_est > prob_f0 + prob_f1 & throw_est < prob_f0 + prob_f1 + prob_f2 ~ 2,
          throw_est > prob_f0 + prob_f1 + prob_f2 & throw_est < prob_f0 + prob_f1 + prob_f2 + prob_f3 ~ 3,
          throw_est > prob_f0 + prob_f1 + prob_f2 + prob_f3 ~ 4,
          TRUE ~ 5
        )
    ) %>% 
    filter(est_stage != 5) %>%
    roc(advanced, est_stage, levels = c("early", "advanced")) %>%
    coords(
      transpose = FALSE,
      c(1, 2, 3, 4)
    )
  
}

## auc for significant fibrosis
roc_significant <- function(x) {
  
  x %>% 
    mutate(
      est_stage = 
        case_when(
          throw_est < prob_f0 ~ 0,
          throw_est > prob_f0 & throw_est < prob_f0 + prob_f1 ~ 1,
          throw_est > prob_f0 + prob_f1 & throw_est < prob_f0 + prob_f1 + prob_f2 ~ 2,
          throw_est > prob_f0 + prob_f1 + prob_f2 & throw_est < prob_f0 + prob_f1 + prob_f2 + prob_f3 ~ 3,
          throw_est > prob_f0 + prob_f1 + prob_f2 + prob_f3 ~ 4,
          TRUE ~ 5
        )
    ) %>% 
    filter(est_stage != 5) %>%
    roc(significant, est_stage, levels = c("early", "significant"))
  
}

## sensitivity and specificity for significant fibrosis
roc_cuts_significant <- function(x) {
  
  x %>% 
    mutate(
      est_stage = 
        case_when(
          throw_est < prob_f0 ~ 0,
          throw_est > prob_f0 & throw_est < prob_f0 + prob_f1 ~ 1,
          throw_est > prob_f0 + prob_f1 & throw_est < prob_f0 + prob_f1 + prob_f2 ~ 2,
          throw_est > prob_f0 + prob_f1 + prob_f2 & throw_est < prob_f0 + prob_f1 + prob_f2 + prob_f3 ~ 3,
          throw_est > prob_f0 + prob_f1 + prob_f2 + prob_f3 ~ 4,
          TRUE ~ 5
        )
    ) %>% 
    filter(est_stage != 5) %>%
    roc(significant, est_stage, levels = c("early", "significant")) %>%
    coords(
      transpose = FALSE,
      c(1, 2, 3, 4)
    )
  
}

# select patients with appropriate fibrosis stages to match observational data
## identifies patients from overall 100,000 patient cohort 
## generates performance characteristics for advanced and significant fibrosis
bx_performance_fn <- function(x, n_pat, prob_f0, prob_f1, prob_f2, prob_f3, prob_f4){
  
  data_f0 <-
    x %>%
    filter(true_stage == "F0") %>%
    slice_sample(n = prob_f0 * n_pat, replace = TRUE)
  
  data_f1 <-
    x %>%
    filter(true_stage == "F1") %>%
    slice_sample(n = prob_f1 * n_pat, replace = TRUE)
  
  data_f2 <-
    x %>%
    filter(true_stage == "F2") %>%
    slice_sample(n = prob_f2 * n_pat, replace = TRUE)
  
  data_f3 <-
    x %>%
    filter(true_stage == "F3") %>%
    slice_sample(n = prob_f3 * n_pat, replace = TRUE)
  
  data_f4 <-
    x %>%
    filter(true_stage == "F4") %>%
    slice_sample(n = prob_f4 * n_pat, replace = TRUE)
  
  data <-
    bind_rows(data_f0, data_f1, data_f2, data_f3, data_f4)
  
  auc_adv <- roc_advanced(data)
  
  auc_adv_tib <-
    as_tibble(auc_adv$auc) %>%
    mutate(auc_adv = as.double(value)) %>%
    select(auc_adv)
  
  f3_cuts <-
    roc_cuts_advanced(data) %>%
    filter(threshold == 3) %>%
    rename(
      spec_f3_adv = specificity,
      sens_f3_adv = sensitivity) %>%
    select(-threshold)
  
  auc_sig <- roc_significant(data)
  
  auc_sig_tib <-
    as_tibble(auc_sig$auc) %>%
    mutate(auc_sig = as.double(value)) %>%
    select(auc_sig)
  
  f2_cuts <-
    roc_cuts_significant(data) %>%
    filter(threshold == 2) %>%
    rename(
      spec_f2_sig = specificity,
      sens_f2_sig = sensitivity) %>%
    select(-threshold)
  
  out <- 
    bind_cols(auc_adv_tib, f3_cuts, auc_sig_tib, f2_cuts)
  
}

#### analysis ####
# import population data
grey <- read_csv("grey_patients.csv")

# estimate biopsy performance
## biopsy performance of single study of n_pat
n_pat = 500
prob_f0 = 0.25
prob_f1 = 0.25
prob_f2 = 0.2
prob_f3 = 0.2
prob_f4 = 0.1

data_f0 <-
  grey %>%
  filter(true_stage == "F0") %>%
  slice_sample(n = prob_f0 * n_pat, replace = TRUE)

data_f1 <-
  grey %>%
  filter(true_stage == "F1") %>%
  slice_sample(n = prob_f1 * n_pat, replace = TRUE)

data_f2 <-
  grey %>%
  filter(true_stage == "F2") %>%
  slice_sample(n = prob_f2 * n_pat, replace = TRUE)

data_f3 <-
  grey %>%
  filter(true_stage == "F3") %>%
  slice_sample(n = prob_f3 * n_pat, replace = TRUE)

data_f4 <-
  grey %>%
  filter(true_stage == "F4") %>%
  slice_sample(n = prob_f4 * n_pat, replace = TRUE)

single_bx_study <-
  bind_rows(data_f0, data_f1, data_f2, data_f3, data_f4)

auc_adv <- roc_advanced(single_bx_study)
cuts_adv <- roc_cuts_advanced(single_bx_study)

auc_sig <- roc_significant(single_bx_study)
cuts_sig <- roc_cuts_significant(single_bx_study)

auc_adv
cuts_adv

auc_sig
cuts_sig

## simulate n_sim biopsy cohorts of n_pat
n_pat = 500
n_sim = 1000

bx_sim <-
  replicate(
    n = n_sim,
    bx_performance_fn(grey, 500, 0.25, 0.25, 0.2, 0.2, 0.1),
    simplify = FALSE
  )

bx_sim_tib <- bind_rows(bx_sim)


## descibe biopsy performance characteristics
summary(bx_sim_tib)

## plots to show distribution of 
sum_vals <-
  bx_sim_tib %>%
  pivot_longer(
    everything(), 
    names_to = "characteristic", values_to = "value") %>%
  group_by(characteristic) %>%
  summarise(
    mean = mean(value),
    sd = sd(value),
    min = min(value),
    max = max(value)
  )

cut_f3 <-
  bx_sim_tib %>%
  select(sens_f3_adv, spec_f3_adv) %>%
  rename(
    sensitivity = sens_f3_adv,
    specificity = spec_f3_adv
  ) %>%
  pivot_longer(everything(), names_to = "characteristic", values_to = "value")

plot_f3 <-
  ggplot(cut_f3) +
  geom_boxplot(aes(x = characteristic, y = value)) +
  geom_jitter(aes(x = characteristic, y = value), alpha = 0.1) +
  geom_hline(yintercept = 0.8, linetype = 2) +
  ylim(0.5, 1) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(
    title = "Identification of advanced fibrosis (F3+)",
    subtitle = glue::glue("Mean AUROC {round(sum_vals[[1,2]], 2)}")
  )


cut_f2 <-
  bx_sim_tib %>%
  select(sens_f2_sig, spec_f2_sig) %>%
  rename(
    sensitivity = sens_f2_sig,
    specificity = spec_f2_sig
  ) %>%
  pivot_longer(everything(), names_to = "characteristic", values_to = "value")

plot_f2 <-
  ggplot(cut_f2) +
  geom_boxplot(aes(x = characteristic, y = value)) +
  geom_jitter(aes(x = characteristic, y = value), alpha = 0.1) +
  geom_hline(yintercept = 0.8, linetype = 2) +
  ylim(0.5, 1) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(
    title = "Identification of significant fibrosis (F2+)",
    subtitle = glue::glue("Mean AUROC {round(sum_vals[[2,2]], 2)}")
  )
  
## manuscript figure 2
plot_f3 / plot_f2 + plot_annotation(tag_levels = "A")

ggsave("Figure 2.jpeg", device = "jpeg", width = 4.5, height = 9, units = "in")


