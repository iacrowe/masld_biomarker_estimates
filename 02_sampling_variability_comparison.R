library(tidyverse)
# functions
## select patients with appropriate fibrosis stages to match observational data
### identifies patients from overall 100,000 patient cohort
bx_select_fn <- function(x, n_pat, prob_f0, prob_f1, prob_f2, prob_f3, prob_f4){
  
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
  
}

## estimates variability due to sampling for 500 patients across all fibrosis stages
bx_sim_ratziu_all_fn <- function(x){
  
  a <-
    bx_select_fn(x, 500, 0.27, 0.27, 0.2, 0.2, 0.06) %>% # F stage distribution from Ratziu
    select(cpa, est_stage, true_stage) %>%
    rename(
      biopsy_1 = est_stage, 
      biopsy_2 = true_stage)
  
  
  a %>%
    mutate(
      f_change = 
        case_when(
          biopsy_1 == "F0" & biopsy_2 == "F0" ~ "no change",
          biopsy_1 == "F0" & biopsy_2 == "F1" ~ "worsened",
          biopsy_1 == "F0" & biopsy_2 == "F2" ~ "worsened",
          biopsy_1 == "F0" & biopsy_2 == "F3" ~ "worsened",
          biopsy_1 == "F0" & biopsy_2 == "F4" ~ "worsened",
          biopsy_1 == "F1" & biopsy_2 == "F0" ~ "improved",
          biopsy_1 == "F1" & biopsy_2 == "F1" ~ "no change",
          biopsy_1 == "F1" & biopsy_2 == "F2" ~ "worsened",
          biopsy_1 == "F1" & biopsy_2 == "F3" ~ "worsened",
          biopsy_1 == "F1" & biopsy_2 == "F4" ~ "worsened",
          biopsy_1 == "F2" & biopsy_2 == "F0" ~ "improved",
          biopsy_1 == "F2" & biopsy_2 == "F1" ~ "improved",
          biopsy_1 == "F2" & biopsy_2 == "F2" ~ "no change",
          biopsy_1 == "F2" & biopsy_2 == "F3" ~ "worsened",
          biopsy_1 == "F2" & biopsy_2 == "F4" ~ "worsened",
          biopsy_1 == "F3" & biopsy_2 == "F0" ~ "improved",
          biopsy_1 == "F3" & biopsy_2 == "F1" ~ "improved",
          biopsy_1 == "F3" & biopsy_2 == "F2" ~ "improved",
          biopsy_1 == "F3" & biopsy_2 == "F3" ~ "no change",
          biopsy_1 == "F3" & biopsy_2 == "F4" ~ "worsened",
          biopsy_1 == "F4" & biopsy_2 == "F0" ~ "improved",
          biopsy_1 == "F4" & biopsy_2 == "F1" ~ "improved",
          biopsy_1 == "F4" & biopsy_2 == "F2" ~ "improved",
          biopsy_1 == "F4" & biopsy_2 == "F3" ~ "improved",
          biopsy_1 == "F4" & biopsy_2 == "F4" ~ "improved"
        )
    ) %>%
    count(f_change)
  
}

## estimates variability due to sampling for 500 patients with f2 or f3
bx_sim_ratziu_f2f3_fn <- function(x){
  
  a <-
    bx_select_fn(x, 500, 0.27, 0.27, 0.2, 0.2, 0.06) %>% # F stage distribution from Ratziu
    select(cpa, est_stage, true_stage) %>%
    rename(
      biopsy_1 = est_stage, 
      biopsy_2 = true_stage) %>%
    filter(biopsy_1 == "F2" | biopsy_1 == "F3")
  
  
  a %>%
    mutate(
      f_change = 
        case_when(
          biopsy_1 == "F0" & biopsy_2 == "F0" ~ "no change",
          biopsy_1 == "F0" & biopsy_2 == "F1" ~ "worsened",
          biopsy_1 == "F0" & biopsy_2 == "F2" ~ "worsened",
          biopsy_1 == "F0" & biopsy_2 == "F3" ~ "worsened",
          biopsy_1 == "F0" & biopsy_2 == "F4" ~ "worsened",
          biopsy_1 == "F1" & biopsy_2 == "F0" ~ "improved",
          biopsy_1 == "F1" & biopsy_2 == "F1" ~ "no change",
          biopsy_1 == "F1" & biopsy_2 == "F2" ~ "worsened",
          biopsy_1 == "F1" & biopsy_2 == "F3" ~ "worsened",
          biopsy_1 == "F1" & biopsy_2 == "F4" ~ "worsened",
          biopsy_1 == "F2" & biopsy_2 == "F0" ~ "improved",
          biopsy_1 == "F2" & biopsy_2 == "F1" ~ "improved",
          biopsy_1 == "F2" & biopsy_2 == "F2" ~ "no change",
          biopsy_1 == "F2" & biopsy_2 == "F3" ~ "worsened",
          biopsy_1 == "F2" & biopsy_2 == "F4" ~ "worsened",
          biopsy_1 == "F3" & biopsy_2 == "F0" ~ "improved",
          biopsy_1 == "F3" & biopsy_2 == "F1" ~ "improved",
          biopsy_1 == "F3" & biopsy_2 == "F2" ~ "improved",
          biopsy_1 == "F3" & biopsy_2 == "F3" ~ "no change",
          biopsy_1 == "F3" & biopsy_2 == "F4" ~ "worsened",
          biopsy_1 == "F4" & biopsy_2 == "F0" ~ "improved",
          biopsy_1 == "F4" & biopsy_2 == "F1" ~ "improved",
          biopsy_1 == "F4" & biopsy_2 == "F2" ~ "improved",
          biopsy_1 == "F4" & biopsy_2 == "F3" ~ "improved",
          biopsy_1 == "F4" & biopsy_2 == "F4" ~ "improved"
        )
    ) %>%
    count(f_change)
  
}

# analysis
## import patient population
grey <- read_csv("grey_patients.csv")

n_sim = 1

## estimate sampling across all fibrosis stages
multi_cpa_sampl_all <-
  replicate(
    n_sim,
    bx_sim_ratziu_all_fn(grey), 
    simplify = FALSE)

cpa_est_sampl_error_all <-
  bind_rows(multi_cpa_sampl_all) %>%
  mutate(sim = rep(1:n_sim, each = 3))

all <- cpa_est_sampl_error_all %>% select(n) %>% mutate(n2 = c(21, 59, 21))
all
chisq.test(all)


## estimate sampling for fibrosis stages 2 and 3
multi_cpa_sampl_f2f3 <-
  replicate(
    n_sim,
    bx_sim_ratziu_f2f3_fn(grey), 
    simplify = FALSE)

cpa_est_sampl_error_f2f3 <-
  bind_rows(multi_cpa_sampl_f2f3) %>%
  mutate(sim = rep(1:n_sim, each = 3))

f2f3 <- cpa_est_sampl_error_f2f3 %>% select(n) %>% mutate(n2 = c(15, 26, 7))
f2f3
chisq.test(f2f3)

## manipulate to plot
f2f3_prop <-
  f2f3 %>%
  mutate(
    cpa_prop = n / sum(f2f3$n),
    ratziu_prop = n2 / sum(f2f3$n2)
  ) %>%
  select(ends_with("prop")) %>%
  mutate(f_change = c("improved", "no change", "worsened")) %>%
  pivot_longer(
    cols = -f_change, names_to = "dataset", values_to = "proportion"
  ) %>%
  mutate(
    dataset = if_else(dataset == "cpa_prop", "Simulated", "Observed")
  )


## plot
color_fill <-
  c("#1e466e", "#aadce0", "#ef8a47")

ggplot(f2f3_prop) +
  geom_col(aes(x = dataset, y = proportion, fill = f_change)) +
  scale_fill_manual(
    values = color_fill, 
    guide_legend(title = "Fibrosis change", title.position = "left", title.hjust = 0)) +
  ylab("Proportion") +
  theme_classic() +
  theme(axis.title.x = element_blank())

## supplemental figure 2
ggsave("Supplementary Figure 2.jpeg", device = "jpeg", width = 4, height = 6, units = "in")




