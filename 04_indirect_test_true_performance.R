library(tidyverse)
library(pROC)
library(patchwork)


#### functions ####
# function to simulate liver biopsy staging for a patient with given fibrosis area
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


# estimate true biomarker performance
## solution derived from Waikar et al, doi: 10.1681/ASN.2010111124 
make_data_est_true_nit <- function(x, prevalence, app_sens, app_spec){
  
  data <- 
    x %>%
    mutate(
      Prev = prevalence,
      AppSens = app_sens,
      AppSpec = app_spec
    )
  
  data_analyse <-
    tibble(
      a = data$Prev * data$SensG,
      b = (1-data$Prev) * (1-data$SpecG),
      c = data$AppSens,
      d = data$AppSpec,
      k = (1-data$Prev) * data$SpecG,
      l = data$Prev * (1-data$SensG) 
    )
  
}

est_true_nit <- function(c, a, b, d, k, l){
  
  sens_lhs = c * (a + b) - b
  spec_lhs = d * (k + l) - l
  
  rhs_matrix <-
    matrix(c(a, -b, -l, k), nrow = 2, byrow = TRUE)
  
  lhs_matrix <-
    matrix(c(sens_lhs, spec_lhs))
  
  j <- solve(rhs_matrix, lhs_matrix)
  
  as_tibble(j) %>%
    mutate(characteristic = c("sensitivity", "specificity")) %>%
    pivot_wider(id_cols = everything(), names_from = characteristic, values_from = V1)
  
}

## read in patient population (from 01_generate_patient_pop.R)
grey <- read_csv("grey_patients.csv")


## simulate n_sim biopsy cohorts of n_pat
### fibrosis distribution from mozes et al, dx.doi.org/10.1136/gutjnl-2021-324243
n_pat = 500
n_sim = 100

## f3
### fib-4
prob_f0 = 0.25
prob_f1 = 0.25
prob_f2 = 0.20
prob_f3 = 0.20
prob_f4 = 0.10


bx_sim_fib4 <-
  replicate(
    n = n_sim,
    bx_performance_fn(
      grey, n_pat, prob_f0, prob_f1, prob_f2, prob_f3, prob_f4),
    simplify = FALSE
  )

bx_sim_fib4_tib <- bind_rows(bx_sim_fib4)


f3_bx_perform_fib4 <-
  bx_sim_fib4_tib %>%
  select(spec_f3_adv, sens_f3_adv) %>%
  rename(
    SpecG = spec_f3_adv,
    SensG = sens_f3_adv)


est_nit_data_fib4 <-
  make_data_est_true_nit(
    f3_bx_perform_fib4,
    prevalence = 0.3, # 0.3
    app_sens = 0.69,
    app_spec = 0.70
  ) %>%
  select(c, a, b, d, k, l)

est_nit_f3_fib4 <-
  pmap(est_nit_data_fib4, est_true_nit) %>%
  bind_rows() %>%
  pivot_longer(
    everything(), names_to = "characteristic", values_to = "value"
  ) %>%
  mutate(test = "fib4")


### apri
prob_f0 = 0.25
prob_f1 = 0.25
prob_f2 = 0.20
prob_f3 = 0.20
prob_f4 = 0.10


bx_sim_apri <-
  replicate(
    n = n_sim,
    bx_performance_fn(
      grey, n_pat, prob_f0, prob_f1, prob_f2, prob_f3, prob_f4),
    simplify = FALSE
  )

bx_sim_apri_tib <- bind_rows(bx_sim_apri)


f3_bx_perform_apri <-
  bx_sim_apri_tib %>%
  select(spec_f3_adv, sens_f3_adv) %>%
  rename(
    SpecG = spec_f3_adv,
    SensG = sens_f3_adv)


est_nit_data_apri <-
  make_data_est_true_nit(
    f3_bx_perform_apri,
    prevalence = 0.3, # 0.30
    app_sens = 0.67,
    app_spec = 0.63
  ) %>%
  select(c, a, b, d, k, l)

est_nit_f3_apri <-
  pmap(est_nit_data_apri, est_true_nit) %>%
  bind_rows() %>%
  pivot_longer(
    everything(), names_to = "characteristic", values_to = "value"
  ) %>%
  mutate(test = "apri")


### proc3
prob_f0 = 0.23
prob_f1 = 0.27
prob_f2 = 0.21
prob_f3 = 0.24
prob_f4 = 0.05


bx_sim_proc3 <-
  replicate(
    n = n_sim,
    bx_performance_fn(
      grey, n_pat, prob_f0, prob_f1, prob_f2, prob_f3, prob_f4),
    simplify = FALSE
  )

bx_sim_proc3_tib <- bind_rows(bx_sim_proc3)


f3_bx_perform_proc3 <-
  bx_sim_proc3_tib %>%
  select(spec_f3_adv, sens_f3_adv) %>%
  rename(
    SpecG = spec_f3_adv,
    SensG = sens_f3_adv)


est_nit_data_proc3 <-
  make_data_est_true_nit(
    f3_bx_perform_proc3,
    prevalence = 0.3, # 0.30
    app_sens = 0.77,
    app_spec = 0.59
  ) %>%
  select(c, a, b, d, k, l)

est_nit_f3_proc3 <-
  pmap(est_nit_data_proc3, est_true_nit) %>%
  bind_rows() %>%
  pivot_longer(
    everything(), names_to = "characteristic", values_to = "value"
  ) %>%
  mutate(test = "pro-c3")


### adapt
prob_f0 = 0.23
prob_f1 = 0.27
prob_f2 = 0.21
prob_f3 = 0.24
prob_f4 = 0.05


bx_sim_adapt <-
  replicate(
    n = n_sim,
    bx_performance_fn(
      grey, n_pat, prob_f0, prob_f1, prob_f2, prob_f3, prob_f4),
    simplify = FALSE
  )

bx_sim_adapt_tib <- bind_rows(bx_sim_adapt)


f3_bx_perform_adapt <-
  bx_sim_adapt_tib %>%
  select(spec_f3_adv, sens_f3_adv) %>%
  rename(
    SpecG = spec_f3_adv,
    SensG = sens_f3_adv)


est_nit_data_adapt <-
  make_data_est_true_nit(
    f3_bx_perform_adapt,
    prevalence = 0.3, # 0.30
    app_sens = 0.78,
    app_spec = 0.69
  ) %>%
  select(c, a, b, d, k, l)

est_nit_f3_adapt <-
  pmap(est_nit_data_adapt, est_true_nit) %>%
  bind_rows() %>%
  pivot_longer(
    everything(), names_to = "characteristic", values_to = "value"
  ) %>%
  mutate(test = "adapt")


### elf
prob_f0 = 0.23
prob_f1 = 0.27
prob_f2 = 0.21
prob_f3 = 0.24
prob_f4 = 0.05


bx_sim_elf <-
  replicate(
    n = n_sim,
    bx_performance_fn(
      grey, n_pat, prob_f0, prob_f1, prob_f2, prob_f3, prob_f4),
    simplify = FALSE
  )

bx_sim_elf_tib <- bind_rows(bx_sim_elf)


f3_bx_perform_elf <-
  bx_sim_elf_tib %>%
  select(spec_f3_adv, sens_f3_adv) %>%
  rename(
    SpecG = spec_f3_adv,
    SensG = sens_f3_adv)


est_nit_data_elf <-
  make_data_est_true_nit(
    f3_bx_perform_elf,
    prevalence = 0.3, # 0.30
    app_sens = 0.71,
    app_spec = 0.76
  ) %>%
  select(c, a, b, d, k, l)

est_nit_f3_elf <-
  pmap(est_nit_data_elf, est_true_nit) %>%
  bind_rows() %>%
  pivot_longer(
    everything(), names_to = "characteristic", values_to = "value"
  ) %>%
  mutate(test = "elf")


## collate and tabulate advanced fibrosis
### make manuscript table 1
true_f3_summary_indirect_performance <-
  bind_rows(
    est_nit_f3_fib4, est_nit_f3_apri, est_nit_f3_proc3, est_nit_f3_adapt, est_nit_f3_elf
  ) %>%
  mutate(
    n_sim = rep(1:500, each = 2)
  ) %>%
  pivot_wider(
    id_cols = c(test, n_sim), names_from = characteristic, values_from = value
  ) %>%
  mutate(
    test = as_factor(test),
    test = fct_relevel(test, "fib4", "apri", "pro-c3", "adapt", "elf"),
    sensitivity = if_else(sensitivity >1, 1, sensitivity),
    specificity = if_else(specificity >1, 1, specificity),
    auc = ((sensitivity * specificity) + (sensitivity + (1 - specificity) * specificity)) / 2
  ) %>%
  select(-n_sim) %>%
  pivot_longer(
    -test, names_to = "characteristic", values_to = "value"
  ) %>%
  group_by(test, characteristic) %>%
  summarise(
    mean = round(mean(value), 2),
    sd = round(sd(value), 2),
    min = round(min(value), 2),
    max = round(max(value), 2)
  ) %>%
  mutate(
    max = if_else(max > 1, 1, max)
  ) %>%
  gt::gt()


gt::gtsave(true_f3_summary_indirect_performance, "blood_biomarker_performance.docx")
