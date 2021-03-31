# load packages you need
library(tidyverse)


# read in data
si <- read.csv("data/champagne_creekchub_isotope_data.csv")



# subset only creek chub data to a file 
si_cc <-
  si %>% 
  filter(!is.na(sitecode)) %>%
  filter(role == "cc") %>%
  mutate(delta15n = round(delta15n, 2),
         delta2h = round(delta2h, 2))




# organize and calculate baselines

# calculate algal baseline min for each site
baseline_d2h_algae <-
  si %>% 
  filter(!is.na(sitecode)) %>%
  filter(role == "algae") %>%
  select(sitecode, role, delta2h) %>%
  mutate(delta2h = as.numeric(delta2h)) %>%
  group_by(sitecode, role) %>%
  summarise(d2h_algae = mean(delta2h, na.rm = T),
            d2h_algae_sd = sd(delta2h, na.rm = T),
            n_d2h_algae = n())

baseline_d15n_algae <-
  si %>% 
  filter(!is.na(sitecode)) %>%
  filter(role == "algae") %>%
  select(sitecode, role, delta15n) %>%
  mutate(delta15n = as.numeric(delta15n)) %>%
  group_by(sitecode, role) %>%
  summarise(d15n_algae = mean(delta15n, na.rm = T),
            d15n_algae_sd = sd(delta15n, na.rm = T),
            n_d15n_algae = n())

baseline_d2h_det <-
  si %>% 
  filter(!is.na(sitecode)) %>%
  filter(role == "det") %>%
  select(sitecode, role, delta2h) %>%
  mutate(delta2h = as.numeric(delta2h)) %>%
  group_by(sitecode, role) %>%
  summarise(d2h_det = mean(delta2h, na.rm = T),
            d2h_det_sd = sd(delta2h, na.rm = T),
            n_d2h_det = n())

baseline_d15n_det <-
  si %>% 
  filter(!is.na(sitecode)) %>%
  filter(role == "det") %>%
  select(sitecode, role, delta15n) %>%
  mutate(delta15n = as.numeric(delta15n)) %>%
  group_by(sitecode, role) %>%
  summarise(d15n_det = mean(delta15n, na.rm = T),
            d15n_det_sd = sd(delta15n, na.rm = T),
            n_d15n_det = n())

# combine baseline files
baselines_d15n <- full_join(baseline_d15n_algae, baseline_d2h_algae, by = c("sitecode", "role"))

baselines_d2h <- full_join(baseline_d15n_det, baseline_d2h_det, by = c("sitecode", "role"))

baselines <- merge(baselines_d15n, baselines_d2h, by = c("sitecode")) %>%
  select(-role.x, -role.y) %>%
  arrange(sitecode) 



# data missing for 2018 terrestrial detritus for ep2 and ep3 so was assumed to be the same as ep1
baselines[12, 11] <- baselines[11, 11] #H for P2 det from P1
baselines[13, 11] <- baselines[11, 11] #H for P3 det from P1



# merge cc data to baseline data by site
si_est <- left_join(si_cc, baselines) %>%
  arrange(sitecode) %>%
  mutate(delta2h = as.numeric(delta2h)) %>%
  filter(! is.na(delta2h))



# adjust for water for d2h in creek chub
d2h_water = -62.7
w1 = 0.17

# adjust creek chub for water for d2h
w2 = 1 - (1 - w1)^2

#correct d2h for CC
si_est$delta2h_cor <- (si_est$delta2h - (w2 * d2h_water))/(1 - w2)

# check correction
plot(si_est$delta2h_cor ~ si_est$delta2h, xlim = c(-180, -80), ylim = c(-180, -80))
abline(0, 1)



# read in environmental data for regressions
p_ag <- read.csv("data/watershed_local_prop_ag.csv")


# remove extra files
rm(list=setdiff(ls(), c("si_est", "p_ag")))

