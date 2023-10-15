library(tidyverse)
library(patchwork)

nash_data <-
  read_csv("virtual_slide_manual_threshold_1.25x.csv") %>%
  rename(
    length = `length (mm)`,
    f_stage = score,
    cpa = `cpa (%)`
  ) %>%
  filter(length > 17 & length < 23)

sum_original <-
  nash_data %>%
  group_by(f_stage) %>%
  summarise(
    mean_cpa = mean(cpa),
    sd_cpa = sd(cpa)
  )

original <-
  nash_data %>%
  group_by(f_stage) %>%
  count() %>%
  left_join(sum_original)


color_fill <-
  c("#1e466e", "#528fad", "#aadce0", "#ffd06f", "#ef8a47")

ggplot(nash_data) +
  ggdist::stat_slab(aes(y = cpa, x = as.character(f_stage)), fill = "skyblue", alpha = 0.5) +
  ggdist::stat_dotsinterval(aes(y = cpa, x = as.character(f_stage)), fill = "deepskyblue", side = "topright") +
  theme_classic()

ggplot(nash_data, aes(x = cpa, y = as.character(f_stage), fill = as.character(f_stage))) +
  ggdist::stat_slab(alpha = 0.2) +
  ggdist::stat_dotsinterval(side = "topright",  slab_linewidth = NA) +
  scale_fill_manual(values = color_fill) +
  theme_classic() +
  theme(legend.position = "none")


groups_cpa <-
  tibble(
    cpa = seq(0, 40, 0.1)
  ) %>%
  mutate(cpa = as.character(cpa))


## f0
f0_data <-
  nash_data %>% 
  filter(f_stage == 0)

MASS::fitdistr(f0_data$cpa, "lognormal")

f0_plot <-
  ggplot() +
  geom_density(aes(x = rlnorm(10000, 1.599, 0.355)), fill = "blue") +
  geom_density(data = f0_data, aes(x = cpa)) +
  xlim(0, 35) +
  ylim(0, 0.25)

F0 <- 
  tibble(
    cpa = round(rlnorm(100000, 1.599, 0.355), 1),
    stage = "F0"
  ) %>%
  count(cpa, name = "count_f0") %>%
  mutate(cpa = as.character(cpa))

## f1
f1_data <-
  nash_data %>% 
  filter(f_stage == 1)

MASS::fitdistr(f1_data$cpa, "lognormal")

f1_plot <-
  ggplot() +
  geom_density(aes(x = rlnorm(10000, 1.875, 0.288)), fill = "blue") +
  geom_density(data = f1_data, aes(x = cpa)) +
  xlim(0, 35) +
  ylim(0, 0.25)

F1 <- 
  tibble(
    cpa = round(rlnorm(100000, 1.875, 0.288), 1),
    stage = "F1"
  ) %>%
  count(cpa, name = "count_f1") %>%
  mutate(cpa = as.character(cpa))

## f2
f2_data <-
  nash_data %>% 
  filter(f_stage == 2)

MASS::fitdistr(f2_data$cpa, "lognormal")

f2_plot <-
  ggplot() +
  geom_density(aes(x = rlnorm(10000, 2.154, 0.25)), fill = "blue") +
  geom_density(data = f2_data, aes(x = cpa)) +
  xlim(0, 35) +
  ylim(0, 0.25)

F2 <- 
  tibble(
    cpa = round(rlnorm(100000, 2.154, 0.25), 1),
    stage = "F2"
  ) %>%
  count(cpa, name = "count_f2") %>%
  mutate(cpa = as.character(cpa))

## f3
f3_data <-
  nash_data %>% 
  filter(f_stage == 3)

MASS::fitdistr(f3_data$cpa, "lognormal")

f3_plot <-
  ggplot() +
  geom_density(aes(x = rlnorm(10000, 2.79, 0.188)), fill = "blue") +
  geom_density(data = f3_data, aes(x = cpa)) +
  xlim(0, 35) +
  ylim(0, 0.25)

F3 <- 
  tibble(
    cpa = round(rlnorm(100000, 2.79, 0.188), 1),
    stage = "F3"
  ) %>%
  count(cpa, name = "count_f3") %>%
  mutate(cpa = as.character(cpa))

## f4
f4_data <-
  nash_data %>% 
  filter(f_stage == 4)

MASS::fitdistr(f4_data$cpa, "normal")

f4_plot <-
  ggplot() +
  geom_density(aes(x = rnorm(10000, 18.2, 3.87)), fill = "blue") +
  geom_density(data = f4_data, aes(x = cpa)) +
  xlim(0, 35) +
  ylim(0, 0.25)

F4 <- 
  tibble(
    cpa = round(rnorm(100000, 18.2, 3.87), 1),
    stage = "F4"
  ) %>%
  count(cpa, name = "count_f4") %>%
  mutate(cpa = as.character(cpa))

f0_plot / f1_plot / f2_plot / f3_plot / f4_plot


sim_data_long <-
  tibble(
    f0 = rlnorm(10000, 1.599, 0.355),
    f1 = rlnorm(10000, 1.875, 0.288),
    f2 = rlnorm(10000, 2.154, 0.25),
    f3 = rlnorm(10000, 2.79, 0.188),
    f4 = rnorm(10000, 18.2, 3.87)
  ) %>%
  pivot_longer(
    everything(),
    names_to = "stage",
    values_to = "cpa"
  ) %>%
  mutate(
    f_stage = rep(seq(0, 4, 1), times = 10000)
  )

measured_plus_sim_cpa <-
  ggplot() +
  ggdist::stat_slab(
    data =sim_data_long,
    aes(x = cpa, y = as.character(f_stage), fill = as.character(f_stage)),
    alpha = 0.2) +
  ggdist::stat_dotsinterval(
    data = nash_data, 
    aes(x = cpa, y = as.character(f_stage), fill = as.character(f_stage)),
    side = "topright",  slab_linewidth = NA) +
  scale_fill_manual(values = color_fill) +
  xlim(1, 28) +
  xlab("Collagen proportionate area (%)") +
  ylab("True fibrosis stage") +
  theme_classic() +
  theme(legend.position = "none")



## generate prob data
data_wide <-
  groups_cpa %>%
  left_join(F0) %>%
  left_join(F1) %>%
  left_join(F2) %>%
  left_join(F3) %>%
  left_join(F4)

prob_data <-
  data_wide %>%
  replace_na(
    list(count_f0 = 0, count_f1 = 0, count_f2 = 0, count_f3 = 0, count_f4 = 0)) %>%
  mutate(
    prob_f0 = count_f0 / (count_f0 + count_f1 + count_f2 + count_f3 + count_f4),
    prob_f1 = count_f1 / (count_f0 + count_f1 + count_f2 + count_f3 + count_f4),
    prob_f2 = count_f2 / (count_f0 + count_f1 + count_f2 + count_f3 + count_f4),
    prob_f3 = count_f3 / (count_f0 + count_f1 + count_f2 + count_f3 + count_f4),
    prob_f4 = count_f4 / (count_f0 + count_f1 + count_f2 + count_f3 + count_f4)
  ) %>%
  select(
    cpa,
    starts_with("prob")
  )

# write_csv(prob_data, "prob_data.csv")

probs <- read_csv("prob_data.csv")

prob_data_plot <-
  probs %>%
  pivot_longer(
    -cpa,
    names_to = "f_stage",
    values_to = "probability"
  ) %>%
  mutate(
    cpa = as.numeric(cpa),
    f_stage =
      case_when(
        f_stage == "prob_f0" ~ "0",
        f_stage == "prob_f1" ~ "1",
        f_stage == "prob_f2" ~ "2",
        f_stage == "prob_f3" ~ "3",
        f_stage == "prob_f4" ~ "4",
      )
  ) %>%
  filter(
    cpa > 1.8 & cpa <30)


prob_plot <-
  ggplot(prob_data_plot %>% filter(cpa <28)) +
  geom_col(aes(x = cpa, y = probability, fill = f_stage), width = .1) +
  scale_fill_manual(
    values = color_fill, 
    guide_legend(title = "Fibrosis\nstage", title.position = "left")) +
  xlab("Collagen proportionate area (%)") +
  ylab("Probability") +
  xlim(1,28) +
  theme_classic() +
  theme(
    legend.position = "right")


measured_plus_sim_cpa / prob_plot + 
  plot_annotation(tag_levels = 'A')

ggsave("Figure 1.jpeg", device = "jpeg", width = 9, height = 6, units = "in")
