
# Import packages ---------------------------------------------------------------------------------------
library(tidyverse)
library(nlme)


# Import and wrangle data -------------------------------------------------------------------------------
soil_data_combined <- read.csv("soil_data_combined.csv")

soil_df <- soil_data_combined %>% 
  janitor::clean_names() %>% 
  mutate(site_code = str_remove(code, patter = "[0-9]+"),
         cn_ratio = carbon/nitrogen,
         age = as.numeric(str_extract(season, pattern = "[0-9]+")) - as.numeric(str_extract(details, pattern = "[0-9]+"))) %>% 
  filter(site_code != "LL")%>% 
  group_by(season) %>% 
  mutate(zscore_cn_ratio = (cn_ratio - mean(cn_ratio))/sd(cn_ratio)) %>% 
  ungroup()

df <- soil_df %>% 
  janitor::clean_names() %>% 
  select(moisture:carbon, p_h,age) %>% 
  correlation::correlation() %>% 
  summary() %>% plot()
df  

ggplot(soil_df, aes(x = site_code, y = p_h)) +
  geom_boxplot() +
  facet_grid( ~ season)

names(soil_df)
plot(log(ag) ~ nitrogen, soil_df)
abline(lm(log(ag) ~ nitrogen, soil_df))

soil_df <- as.data.frame(soil_df)
m1 <- lme(sqrt(phos) ~ carbon + nitrogen + p_h + age ,
          random = ~ 1|season,
          soil_df, na.action = na.omit)
anova.lme(m1, verbose = T)
MuMIn::r.squaredGLMM(m1)
hist(resid(m1), breaks = 10)
shapiro.test(resid(m1))
ggeffects::ggpredict(m1) %>% plot()


m2 <- lme(log1p(lap) ~ carbon + nitrogen + p_h + age,
          random = ~ 1|season,
          soil_df, na.action = na.omit)
anova.lme(m2)
MuMIn::r.squaredGLMM(m2)
hist(resid(m2))
shapiro.test(resid(m2))
ggeffects::ggpredict(m2) %>% plot()

m3 <- lme(log2(cb) ~ carbon + nitrogen + p_h + age,
          random = ~ 1|season,
          soil_df, na.action = na.omit)
anova.lme(m3)
MuMIn::r.squaredGLMM(m3)
hist(resid(m3))
shapiro.test(resid(m3))
ggeffects::ggpredict(m3) %>% plot()

m4 <- lme(log1p(nag) ~ carbon + nitrogen + p_h + age,
          random = ~ 1|season,
          soil_df, na.action = na.omit)
anova.lme(m4)
MuMIn::r.squaredGLMM(m4)
hist(resid(m4))
shapiro.test(resid(m4))
ggeffects::ggpredict(m4) %>% plot()


m5 <- lme(sqrt(ag) ~ carbon + nitrogen + p_h + age,
          random = ~ 1|season,
          soil_df, na.action = na.omit)
anova.lme(m5)
MuMIn::r.squaredGLMM(m5)
hist(resid(m5))
shapiro.test(resid(m5))
ggeffects::ggpredict(m5) %>% plot()


ggplot(soil_df, aes(x = age, y  = carbon/nitrogen)) +
  geom_point()+
  geom_smooth(method = "lm",
              formula=y ~ poly(x, 2, raw=TRUE))


