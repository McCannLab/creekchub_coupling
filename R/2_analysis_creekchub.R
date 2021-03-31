# run code to prepare data for analysis
source("R/1_prep_data_creekchub.R")

###############################################################################################################
# proportion terrestrial and trophic position estimates based on primary producer values

# prop terrestrial energy estimation
si_est$ter_energy_pp <- (si_est$delta2h_cor - si_est$d2h_algae)/(si_est$d2h_det - si_est$d2h_algae)

# constrain between 0-1
si_est <-
  si_est %>%
  mutate(ter_energy_pp = case_when(
    ter_energy_pp <= 0 ~ 0.001,
    ter_energy_pp >=1 ~ 0.999,
    ter_energy_pp > 0 & ter_energy_pp  < 1 ~ ter_energy_pp,
  ))


# trophic position estimation
si_est$tp <- 1 + ((si_est$delta15n - si_est$d15n_det)/3.4)

######################################################################################################################
#data prep

# join fish meta to full isotope dataset
si_est <- 
  si_est %>%
  left_join(fish_meta, by = c("sitecode", "sampleid"))

# quick summary to check
si_est %>%
  filter(!is.na(totalmm)) %>%
  summarise(mean_tl = mean(totalmm),
            sd_tl = sd(totalmm),
            n = n())

# calculate means for each site of tp and prop terrestrial
si_x <- 
  si_est %>% group_by(sitecode) %>%
  summarise(prop_ter_mean = mean(ter_energy_pp, na.rm = T),
            tp = mean(tp, na.rm = T),
            n = n(),
            tl_mean = mean(totalmm, na.rm = T),
            tl_sd = sd(totalmm, na.rm = T))
            

# merge si data to p ag data 
si_est <- right_join(si_est, p_ag)
si_x<-  right_join(si_x, p_ag)



##########################################################################################

# statistics

# linear regressions mean creek chub length vs local and regional ag
summary(lm(tl_mean ~ p_ag_250, data = si_x))
summary(lm(tl_mean ~ p_ag_watershed, data = si_x))

# proportion terrestrial energy vs local ag
mod_local_ter <- lm(qlogis(prop_ter_mean) ~ p_ag_250 + tl_mean, data = si_x)
summary(mod_local_ter) # mean length not sig
mod_local_ter <- update(mod_local_ter, .~. -tl_mean)
summary(mod_local_ter) # local ag sig

# proportion terrestrial energy vs regional ag
mod_reg_ter <- lm(qlogis(prop_ter_mean) ~ p_ag_watershed + tl_mean, data = si_x)
summary(mod_reg_ter)
mod_reg_ter <- update(mod_reg_ter, .~. -tl_mean)
summary(mod_reg_ter)

# trophic position vs local ag
mod_local_tp <- lm(tp ~ p_ag_250 + tl_mean, data = si_x)
summary(mod_local_tp) # mean length not sig
mod_local_tp <- update(mod_local_tp, .~. -tl_mean)
summary(mod_local_tp) # local ag sig

# trophic position vs regional ag
mod_reg_tp <- lm(tp ~ p_ag_watershed + tl_mean, data = si_x)
summary(mod_reg_tp)
mod_reg_tp <- update(mod_reg_tp, .~. -tl_mean)
summary(mod_reg_tp)

# try some basic predictions based on local ag
new <- data.frame(p_ag_250 = c(0, 1))
plogis(predict(mod_local_ter, new))
predict(mod_local_tp, new)

# make trophic response panel plot

a=
  ggplot(si_x, aes(x = p_ag_250, y = qlogis(prop_ter_mean))) +
  geom_point(data = si_est, aes(x = p_ag_250, y = qlogis(ter_energy_pp)), colour = "gray", pch = 1) +
  geom_point(size = 3) +
  geom_smooth(se = F, method = "lm", colour = "black") +
  scale_y_continuous(labels=function(x)round(plogis(x),digits=2)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .2), labels = seq(0, 1, .2)) +
  theme_classic() +
  xlab("proportion local agriculture (250 m buffer)") +
  ylab("proportion terrestrial energy")

b=
  ggplot(si_x, aes(x = p_ag_watershed, y = qlogis(prop_ter_mean))) +
  geom_point(data = si_est, aes(x = p_ag_watershed, y = qlogis(ter_energy_pp)), colour = "gray", pch = 1) +
  geom_point(size = 3) +
  scale_y_continuous(labels=function(x)round(plogis(x),digits=2)) +
  scale_x_continuous(limits = c(0.38, 1), breaks = seq(0.4, 1, .2), labels = seq(0.4, 1, .2)) +
  theme_classic() +
  xlab("proportion watershed agriculture") +
  ylab("proportion terrestrial energy")

c=
  ggplot(si_x, aes(x = p_ag_250, y = tp)) +
  geom_point(data = si_est, aes(x = p_ag_250, y = tp), colour = "gray", pch = 1) +
  geom_point(size = 3) +
  geom_smooth(se = F, method = "lm", colour = "black") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .2), labels = seq(0, 1, .2)) +
  theme_classic() +
  xlab("proportion local agriculture (250 m buffer)") +
  ylab("trophic position")

d=
  ggplot(si_x, aes(x = p_ag_watershed, y = tp)) +
  geom_point(data = si_est, aes(x = p_ag_watershed, y = tp), colour = "gray", pch = 1) +
  geom_point(size = 3) +
  theme_classic() +
  scale_x_continuous(limits = c(0.38, 1), breaks = seq(0.4, 1, .2), labels = seq(0.4, 1, .2)) +
  xlab("proportion watershed agriculture") +
  ylab("trophic position")


# combine plots
cowplot::plot_grid(a, b, c, d, ncol = 2, labels = c("A", "B", "C", "D", label_size = 18))

# save plot at tiff
#ggsave("trophic_response.tiff", units = "in", width = 8, height = 6, dpi = 600, compression = "lzw")


###########################################################################################
# alternative isotope regressions

# without the two sites around 40% regional ag

# remove the two sites with regional % agriculture at around 40%
si_xx <- si_x %>%
  filter(!p_ag_watershed < 0.6) 

# proportion terrestrial energy vs regional ag without ~40% values
mod_reg_terx <- lm(qlogis(prop_ter_mean) ~ p_ag_watershed + tl_mean, data = si_xx)
summary(mod_reg_terx)
mod_reg_terx <- update(mod_reg_terx, .~. -tl_mean)
summary(mod_reg_terx)

# trophic position vs regional ag without ~40% values
mod_reg_tpx <- lm(tp ~ p_ag_watershed + tl_mean, data = si_xx)
summary(mod_reg_tpx)
mod_reg_tpx <- update(mod_reg_tpx, .~. -tl_mean)
summary(mod_reg_tpx)



# regressions for terrestrial energy without sites P2, P3 
# for which we assumed terrestrial baseline d2h was same as P1

si_xxx <- 
  si_x %>%
  filter(!sitecode %in% c("P2", "P3"))

# proportion terrestrial energy vs local ag
mod_local_terxx <- lm(qlogis(prop_ter_mean) ~ p_ag_250 + tl_mean, data = si_xxx)
summary(mod_local_terxx) # mean length not sig
mod_local_terxx <- update(mod_local_terxx, .~. -tl_mean)
summary(mod_local_terxx) # local ag sig

# proportion terrestrial energy vs regional ag
mod_reg_terxx <- lm(qlogis(prop_ter_mean) ~ p_ag_watershed + tl_mean, data = si_xxx)
summary(mod_reg_terxx)
mod_reg_terxx <- update(mod_reg_terxx, .~. -tl_mean)
summary(mod_reg_terxx)





##########################################################################################

# pca of env data from study sites included in paper

library(vegan)
library(corrgram)
library(Hmisc)
library(DMwR)
library(ggfortify)


# read in site data
envdata <- read.csv("data/2018PCAaggrad.csv") %>%
  filter(sitecode %in% unique(si_x$sitecode)) %>%
  right_join(p_ag)

# log transform variables to be used in pca
envdata <- 
  envdata %>% 
  mutate(
    p_ag_watershed = p_ag_watershed,
    log_sin = log10(sinuosity),
    log_turbidity = log10(turbidity),
    log_p = log10(july.phosphorus),
    log_n = log10(july.nitrogen),
    p_ag_250 = p_ag_250,
    log_water_t = log10(temp)
  )


hist(envdata$temp)

# create a dataframe with data for pca
pcadata <- 
  envdata %>% 
  select(p_ag_watershed, buf_width, p_ag_250, log_turbidity, log_n, log_p, log_sin, log_water_t)
rownames(pcadata) <- envdata$sitecode 

# no P data for EP1 so we use K-nearest neighbor imputation to fill in this missing value so we can include the rest of the
# data from this site in the PCA
pcadata1 <- knnImputation(pcadata, k = 10, scale = T, meth = "weighAvg", distData = NULL)

#####
# stop here unless running PCA
# run pca
pca <- prcomp(pcadata1, scale = TRUE, center = TRUE)
#biplot(pca, xlim  = c(-0.5, 0.6))
pca$x
summary(pca)
pca$rotation

autoplot(pca, data = pcadata1,
         loadings = T, loadings.colour = 'blue',
         loadings.label = T, loadings.label.size = 5, label = T, shape = F) +
  theme_classic(base_size = 16)

#ggsave("output/pca_env_data.pdf", units = "in", width = 8, height = 6, dpi = 600)

# correlogram of variables in the pca
corrgram(pcadata, lower.panel=panel.cor, upper.panel=panel.pts, cex.labels = 1.2)

