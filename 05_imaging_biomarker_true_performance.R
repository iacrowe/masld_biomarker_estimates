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
### fibrosis distribution from selveraj, j hepatol
n_pat = 500
n_sim = 100

## f3
### vcte

prob_f0 = 0.27
prob_f1 = 0.28
prob_f2 = 0.20
prob_f3 = 0.16
prob_f4 = 0.09


bx_sim_lsm <-
  replicate(
    n = n_sim,
    bx_performance_fn(
      grey, n_pat, prob_f0, prob_f1, prob_f2, prob_f3, prob_f4),
    simplify = FALSE
  )

bx_sim_lsm_tib <- bind_rows(bx_sim_lsm)


f3_bx_perform_lsm <-
  bx_sim_lsm_tib %>%
  select(spec_f3_adv, sens_f3_adv) %>%
  rename(
    SpecG = spec_f3_adv,
    SensG = sens_f3_adv)


est_nit_data_lsm <-
  make_data_est_true_nit(
    f3_bx_perform_lsm,
    prevalence = 0.275,
    app_sens = 0.80,
    app_spec = 0.77
  ) %>%
  select(c, a, b, d, k, l)

est_nit_f3_lsm <-
  pmap(est_nit_data_lsm, est_true_nit) %>%
  bind_rows() %>%
  pivot_longer(
    everything(), names_to = "characteristic", values_to = "value"
  ) %>%
  mutate(test = "vcte")


### mre
prob_f0 = 0.30
prob_f1 = 0.31
prob_f2 = 0.20
prob_f3 = 0.11
prob_f4 = 0.08


bx_sim_mre <-
  replicate(
    n = n_sim,
    bx_performance_fn(
      grey, n_pat, prob_f0, prob_f1, prob_f2, prob_f3, prob_f4),
    simplify = FALSE
  )

bx_sim_mre_tib <- bind_rows(bx_sim_mre)


f3_bx_perform_mre <-
  bx_sim_mre_tib %>%
  select(spec_f3_adv, sens_f3_adv) %>%
  rename(
    SpecG = spec_f3_adv,
    SensG = sens_f3_adv)


est_nit_data_mre <-
  make_data_est_true_nit(
    f3_bx_perform_mre,
    prevalence = 0.275,
    app_sens = 0.83,
    app_spec = 0.89
  ) %>%
  select(c, a, b, d, k, l)

est_nit_f3_mre <-
  pmap(est_nit_data_mre, est_true_nit) %>%
  bind_rows() %>%
  pivot_longer(
    everything(), names_to = "characteristic", values_to = "value"
  ) %>%
  mutate(test = "mre")


### pswe
prob_f0 = 0.27
prob_f1 = 0.27
prob_f2 = 0.16
prob_f3 = 0.13
prob_f4 = 0.17


bx_sim_pswe <-
  replicate(
    n = n_sim,
    bx_performance_fn(
      grey, n_pat, prob_f0, prob_f1, prob_f2, prob_f3, prob_f4),
    simplify = FALSE
  )

bx_sim_pswe_tib <- bind_rows(bx_sim_pswe)


f3_bx_perform_pswe <-
  bx_sim_pswe_tib %>%
  select(spec_f3_adv, sens_f3_adv) %>%
  rename(
    SpecG = spec_f3_adv,
    SensG = sens_f3_adv)


est_nit_data_pswe <-
  make_data_est_true_nit(
    f3_bx_perform_pswe,
    prevalence = 0.275,
    app_sens = 0.80,
    app_spec = 0.86
  ) %>%
  select(c, a, b, d, k, l)

est_nit_f3_pswe <-
  pmap(est_nit_data_pswe, est_true_nit) %>%
  bind_rows() %>%
  pivot_longer(
    everything(), names_to = "characteristic", values_to = "value"
  ) %>%
  mutate(test = "pswe")


### 2dswe
prob_f0 = 0.22
prob_f1 = 0.21
prob_f2 = 0.21
prob_f3 = 0.21
prob_f4 = 0.15


bx_sim_2dswe <-
  replicate(
    n = n_sim,
    bx_performance_fn(
      grey, n_pat, prob_f0, prob_f1, prob_f2, prob_f3, prob_f4),
    simplify = FALSE
  )

bx_sim_2dswe_tib <- bind_rows(bx_sim_2dswe)


f3_bx_perform_2dswe <-
  bx_sim_2dswe_tib %>%
  select(spec_f3_adv, sens_f3_adv) %>%
  rename(
    SpecG = spec_f3_adv,
    SensG = sens_f3_adv)


est_nit_data_2dswe <-
  make_data_est_true_nit(
    f3_bx_perform_2dswe,
    prevalence = 0.275,
    app_sens = 0.72,
    app_spec = 0.72
  ) %>%
  select(c, a, b, d, k, l)

est_nit_f3_2dswe <-
  pmap(est_nit_data_2dswe, est_true_nit) %>%
  bind_rows() %>%
  pivot_longer(
    everything(), names_to = "characteristic", values_to = "value"
  ) %>%
  mutate(test = "2dswe")


# collate and plot
## manuscript table 2
true_f3_summary_performance <-
  bind_rows(
    est_nit_f3_lsm, est_nit_f3_mre, est_nit_f3_pswe, est_nit_f3_2dswe
  ) %>%
  mutate(
    n_sim = rep(1:400, each = 2)
  ) %>%
  pivot_wider(
    id_cols = c(test, n_sim), names_from = characteristic, values_from = value
  ) %>%
  mutate(
    test = as_factor(test),
    test = fct_relevel(test, "vcte", "mre", "pswe", "2dswe"),
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
  ) 

true_f3_summary_table <-
  true_f3_summary_performance %>%
  gt::gt()

gt::gtsave(true_f3_summary_table, "imaging_biomarker_performance.docx")

## manuscript figure 3
true_f3_imaging_plot_sum <-
  true_f3_summary_performance %>%
  select(test, characteristic, mean) %>%
  pivot_wider(
    id_cols = test, names_from = characteristic, values_from = mean
  ) %>%
  mutate(
    performance = "true",
    test = as_factor(test),
    test = fct_relevel(test, "mre", "pswe", "vcte", "2dswe")
  )

app_f3_imaging <-
  tibble(
    test = true_f3_imaging_plot_sum$test, # vcte, mre, pswe, 2dswe
    sensitivity = c(0.798, 0.83, 0.802, 0.72),
    specificity = c(0.77, 0.89, 0.86, 0.72),
    auc = c(0.85, 0.9, 0.89, 0.72),
    performance = "apparent"
  )

plot_app_v_true_perform_f3 <-
  true_f3_imaging_plot_sum %>%
  bind_rows(app_f3_imaging) %>%
  mutate(
    sensitivity = if_else(sensitivity >1, 1, sensitivity),
    specificity = if_else(specificity >1, 1, specificity)
  )


color_fill_imaging <-
  c("#528fad", "#aadce0", "#ffd06f", "#ef8a47")


f3_dot_auc <-
  ggplot(
    plot_app_v_true_perform_f3,
    aes(x = performance, y = auc, colour = test, group = test)) +
  geom_point(
    size = 4) +
  geom_line(
    linewidth = 1,
    show.legend = FALSE) +
  geom_hline(yintercept = 0.92, linetype = 2) +
  scale_colour_manual(
    values = color_fill_imaging, 
    guide_legend(title = "Imaging\nbiomarker",
                 title.position = "top")) +
  scale_y_continuous(limits = c(0.5, 1), expand = c(0,0)) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(
    title = "AUC"
  )

f3_dot_sens <-
  ggplot(
    plot_app_v_true_perform_f3,
    aes(x = performance, y = sensitivity, colour = test, group = test)) +
  geom_point(
    size = 4) +
  geom_line(
    linewidth = 1,
    show.legend = FALSE) +
  geom_hline(yintercept = 0.83, linetype = 2) +
  scale_colour_manual(
    values = color_fill_imaging, 
    guide_legend(title = "Imaging\nbiomarker",
                 title.position = "top")) +
  scale_y_continuous(limits = c(0.5, 1), expand = c(0,0)) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(
    title = "Sensitivity"
  )

f3_dot_spec <-
  ggplot(
    plot_app_v_true_perform_f3,
    aes(x = performance, y = specificity, colour = test, group = test)) +
  geom_point(
    size = 4) +
  geom_line(
    linewidth = 1,
    show.legend = FALSE) +
  geom_hline(yintercept = 0.93, linetype = 2) +
  scale_colour_manual(
    values = color_fill_imaging, 
    guide_legend(title = "Imaging\nbiomarker",
                 title.position = "top")) +
  scale_y_continuous(limits = c(0.5, 1), expand = c(0,0)) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(
    title = "Specificity"
  )

f3_dot_auc + f3_dot_sens + f3_dot_spec + 
  plot_layout(guide = "collect") + plot_annotation(tag_levels = "A")

ggsave("Figure 3.jpeg", device = "jpeg", width = 9, height = 4, units = "in")

