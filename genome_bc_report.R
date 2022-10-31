
# Import packages ---------------------------------------------------------


library(tidyverse)
library(nlme)
library(janitor)
library(ggeffects)
library(mgcv)
library(glue)



# Import Data -------------------------------------------------------------


enzyme <- read.csv("../../enzyme_cal_combined.csv") %>% 
  clean_names()

df_enz <- enzyme %>% 
  mutate(biosolid = str_extract(details, pattern = "[0-9]+$"),
         age = case_when(
           str_detect(season, pattern = "2021") ~ (2021 - as.numeric(biosolid)),
           str_detect(season, pattern = "2022") ~ (2022 - as.numeric(biosolid))
         ),
         age = replace_na(as.character(age), "ref"),
         age = factor(age,
                      levels = c("ref", "3" ,  "4" ,  "6" ,  "7" ,  "8", "15",  "16",  "21",  "22",  "26",  "27")),
         year = str_extract(season, pattern = "[0-9]+$"),
         cn_ratio = carbon/nitrogen)
df_enz_noref <- df_enz %>% 
  filter(age != "ref") %>% 
  mutate(age = as.numeric(paste(age)))



# NLME models with age ----------------------------------------------------


theme_plot <- theme_bw() +
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        panel.grid.minor = element_blank())

mod1_lap <- lme(sqrt(lap) ~ poly(age,2)  ,
                random = ~ 1|season/quadrate_no,
                na.action = na.omit,
                data = df_enz_noref)
summary(mod1_lap)
anova.lme(mod1_lap)
plot(ggpredict(mod1_lap))
plot(mod1_lap)


##plot LAP##
lap <- plot(ggpredict(mod1_lap, terms = "age")) +
  theme_plot +
  labs(title = NULL,
       x = "Reclamation age (years)",
       y = expression(paste("Leucine-amino-peptidase activity ",
                            "(",mu, "M/g/hr)"))) +
  scale_x_continuous(breaks = seq(0, 30, by = 5)) +
  annotate(geom = "text", label = glue("Rsq = {round(MuMIn::r.squaredGLMM(mod1_lap)[1], 2)},
                                       pval < 0.0001"),
           x = 5, y = 25, size = 5)
lap

ggsave(lap, filename = "../../plot/lap.jpeg", width = 8, height = 6, units = "in", dpi =1000)





mod1_cb <- lme(cb ~ poly(age,2) ,
                random = ~ 1|season/quadrate_no,
                na.action = na.omit,
                data = df_enz_noref)
anova.lme(mod1_cb)
plot(ggpredict(mod1_cb, terms="age [all]"))

##plot cb##
cb <- plot(ggpredict(mod1_cb, terms = "age")) +
  theme_plot +
  labs(title = NULL,
       x = "Reclamation age (years)",
       y = expression(paste("Cellobiohydrolase activity ",
                            "(",mu, "M/g/hr)"))) +
  scale_x_continuous(breaks = seq(0, 30, by = 5)) +
  annotate(geom = "text", label = glue("Rsq = {round(MuMIn::r.squaredGLMM(mod1_cb)[1], 2)},
                                       pval = 0.33"),
           x = 5, y = 10, size = 5)
cb

ggsave(cb, filename = "../../plot/cb.jpeg", width = 8, height = 6, units = "in", dpi =1000)


mod1_phos <- lme(phos ~ poly(age,2) ,
                random = ~ 1|season/quadrate_no,
                na.action = na.omit,
                data = df_enz_noref)
anova.lme(mod1_phos)
plot(ggpredict(mod1_phos))
MuMIn::r.squaredGLMM(mod1_phos)

##plot phos##
phos <- plot(ggpredict(mod1_phos, terms = "age")) +
  theme_plot +
  labs(title = NULL,
       x = "Reclamation age (years)",
       y = expression(paste("Phosphatase activity ",
                            "(",mu, "M/g/hr)"))) +
  scale_x_continuous(breaks = seq(0, 30, by = 5)) +
  annotate(geom = "text", label = glue("Rsq = {round(MuMIn::r.squaredGLMM(mod1_phos)[1], 2)},
                                      pval < 0.0001"),
           x = 5, y = 10, size = 5)
phos

ggsave(phos, filename = "../../plot/phos.jpeg", width = 8, height = 6, units = "in", dpi =1000)


##plot NAG

mod1_nag <- lme(nag ~ age,
                random = ~ 1|season/quadrate_no,
                na.action = na.omit,
                data = df_enz_noref)
anova.lme(mod1_nag)
plot(ggpredict(mod1_nag))
MuMIn::r.squaredGLMM(mod1_nag)

nag <- plot(ggpredict(mod1_nag, terms = "age")) +
  theme_plot +
  labs(title = NULL,
       x = "Reclamation age (years)",
       y = expression(paste("N-acetyl-",β,"-Glucosaminidase ",
                            "(",mu, "M/g/hr)"))) +
  scale_x_continuous(breaks = seq(0, 30, by = 5)) +
  annotate(geom = "text", label = glue("Rsq = {round(MuMIn::r.squaredGLMM(mod1_nag)[1], 2)},
                                      pval < 0.0001"),
           x = 5, y = 5, size = 5)
nag

ggsave(nag, filename = "../../plot/nag.jpeg", width = 8, height = 6, units = "in", dpi =1000)


##AG
mod1_ag <- lme(ag ~ age,
                random = ~ 1|season/quadrate_no,
                na.action = na.omit,
                data = df_enz_noref)
##plot ag##
ag <- plot(ggpredict(mod1_ag, terms = "age")) +
  theme_plot +
  labs(title = NULL,
       x = "Reclamation age (years)",
       y = expression(paste(alpha,"-Glucosidase activity ",
                            "(",mu, "M/g/hr)"))) +
  scale_x_continuous(breaks = seq(0, 30, by = 5)) +
  annotate(geom = "text", label = glue("Rsq = {round(MuMIn::r.squaredGLMM(mod1_ag)[1], 2)},
                                      pval < 0.0001"),
           x = 5, y = 5, size = 5)
ag

ggsave(ag, filename = "../../plot/ag.jpeg", width = 8, height = 6, units = "in", dpi =1000)



# NLME model with cn_ratio ------------------------------------------------

##LAP##
mod2_lap <- lme(sqrt(lap) ~ poly(cn_ratio,2) ,
                random = ~ 1|season/quadrate_no,
                na.action = na.omit,
                data = df_enz_noref)
summary(mod2_lap)
anova.lme(mod2_lap)$`p-value`[2]
plot(ggpredict(mod2_lap))
plot(mod2_lap)


##plot LAP##
lap_cn_ratio <- plot(ggpredict(mod2_lap, terms = "cn_ratio")) +
  theme_plot +
  labs(title = NULL,
       x = "C/N ratio",
       y = expression(paste("Leucine-amino-peptidase activity ",
                            "(",mu, "M/g/hr)"))) +
  scale_x_continuous(breaks = seq(7, 30, by = 5)) +
  annotate(geom = "text", label = glue("Rsq = {round(MuMIn::r.squaredGLMM(mod2_lap)[1], 2)},
                                       pval < {round(anova.lme(mod2_lap)$`p-value`[2], 2)}"),
           x = 9, y = 25, size = 5)
lap_cn_ratio
ggsave(lap_cn_ratio, filename = "../../plot/lap_cn_ratio.jpeg", width = 8, height = 6, units = "in", dpi =1000)


##CB##
mod2_cb <- lme(sqrt(cb) ~ poly(cn_ratio,2) ,
                random = ~ 1|season/quadrate_no,
                na.action = na.omit,
                data = df_enz_noref)
summary(mod2_cb)
anova.lme(mod2_cb)$`p-value`[2]
plot(ggpredict(mod2_cb))
plot(mod2_cb)


##plot cb##
cb_cn_ratio <- plot(ggpredict(mod2_cb, terms = "cn_ratio")) +
  theme_plot +
  labs(title = NULL,
       x = "C/N ratio",
       y = expression(paste("Cellobiohydrolase activity ",
                            "(",mu, "M/g/hr)"))) +
  scale_x_continuous(breaks = seq(7, 30, by = 5)) +
  annotate(geom = "text", label = glue("Rsq = {round(MuMIn::r.squaredGLMM(mod2_cb)[1], 2)},
                                       pval < {round(anova.lme(mod2_cb)$`p-value`[2], 2)}"),
           x = 9, y = 25, size = 5)
cb_cn_ratio
ggsave(cb_cn_ratio, filename = "../../plot/cb_cn_ratio.jpeg", width = 8, height = 6, units = "in", dpi =1000)


##Phos##
mod2_phos <- lme(sqrt(phos) ~ poly(cn_ratio,2) ,
               random = ~ 1|season/quadrate_no,
               na.action = na.omit,
               data = df_enz_noref)
summary(mod2_phos)
anova.lme(mod2_phos)$`p-value`[2]
plot(ggpredict(mod2_phos))
plot(mod2_phos)


##plot phos##
phos_cn_ratio <- plot(ggpredict(mod2_phos, terms = "cn_ratio")) +
  theme_plot +
  labs(title = NULL,
       x = "C/N ratio",
       y = expression(paste("Phosphatase activity ",
                            "(",mu, "M/g/hr)"))) +
  scale_x_continuous(breaks = seq(7, 30, by = 5)) +
  annotate(geom = "text", label = glue("Rsq = {round(MuMIn::r.squaredGLMM(mod2_phos)[1], 2)},
                                       pval < {round(anova.lme(mod2_phos)$`p-value`[2], 2)}"),
           x = 9, y = 60, size = 5)
phos_cn_ratio
ggsave(phos_cn_ratio, filename = "../../plot/phos_cn_ratio.jpeg", width = 8, height = 6, units = "in", dpi =1000)



##NAG##
mod2_nag <- lme(sqrt(nag) ~ poly(cn_ratio,2) ,
                 random = ~ 1|season/quadrate_no,
                 na.action = na.omit,
                 data = df_enz_noref)
summary(mod2_nag)
anova.lme(mod2_nag)$`p-value`[2]
plot(ggpredict(mod2_nag))
plot(mod2_nag)


##plot nag##
nag_cn_ratio <- plot(ggpredict(mod2_nag, terms = "cn_ratio")) +
  theme_plot +
  labs(title = NULL,
       x = "C/N ratio",
       y = expression(paste("N-acetyl-",β,"-Glucosaminidase ",
                            "(",mu, "M/g/hr)"))) +
  scale_x_continuous(breaks = seq(7, 30, by = 5)) +
  annotate(geom = "text", label = glue("Rsq = {round(MuMIn::r.squaredGLMM(mod2_nag)[1], 2)},
                                       pval < {round(anova.lme(mod2_nag)$`p-value`[2], 6)}"),
           x = 9, y = 60, size = 5)
nag_cn_ratio
ggsave(nag_cn_ratio, filename = "../../plot/nag_cn_ratio.jpeg", width = 8, height = 6, units = "in", dpi =1000)



##AG##
mod2_ag <- lme(sqrt(ag) ~ poly(cn_ratio,2, raw=TRUE) ,
                random = ~ 1|season/quadrate_no,
                na.action = na.omit,
                data = df_enz_noref)
summary(mod2_ag)
anova.lme(mod2_ag)$`p-value`[2]
plot(ggpredict(mod2_ag))
plot(mod2_ag)


##plot ag##
ag_cn_ratio <- plot(ggpredict(mod2_ag, terms = "cn_ratio")) +
 # geom_jitter(data = df_enz_noref, aes(x = cn_ratio, y = ag, fill = age, inherit.aes = F))+
  theme_plot +
  labs(title = NULL,
       x = "C/N ratio",
       y = expression(paste("N-acetyl-",β,"-Glucosaminidase ",
                            "(",mu, "M/g/hr)"))) +
  scale_x_continuous(breaks = seq(7, 30, by = 5)) +
  annotate(geom = "text", label = glue("Rsq = {round(MuMIn::r.squaredGLMM(mod2_ag)[1], 2)},
                                       pval < {round(anova.lme(mod2_ag)$`p-value`[2], 3)}"),
           x = 9, y = 10, size = 5)
ag_cn_ratio
ggsave(ag_cn_ratio, filename = "../../plot/ag_cn_ratio.jpeg", width = 8, height = 6, units = "in", dpi =1000)




# CN_ratio ----------------------------------------------------------------

mod_cn <- gam(cn_ratio ~ s(age, k=4),
              data = df_enz_noref)
p <-  broom.mixed::tidy(mod_cn)
plot(mod_cn)

cn <- ggplot(df_enz_noref, aes(x = age, y = cn_ratio)) +
  stat_smooth(method="gam", formula = y ~ s(x, bs = "cs", k=4),
              color = "black",
              size = 0.75)+
  scale_x_continuous(breaks = seq(7, 30, by = 5))+
   theme_plot +
  annotate(geom = "text", label = glue("Rsq = {round(MuMIn::r.squaredGLMM(mod_cn)[1], 2)},
                                       pval < 0.001"),
           x =25, y = 16, size = 5)+
  labs(x = "Reclamation age (years)",
       y = "C/N ratio")
ggsave(cn, filename = "../../plot/cn_ratio.jpeg", width = 8, height = 6, units = "in", dpi =1000)

