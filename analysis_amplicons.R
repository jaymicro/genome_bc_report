
library(tidyverse)
library(vegan)
library(ggsci)
library(emmeans)
library(ggpubr)


# Import data --------------------------------------------------------


count_protist <- read.csv("protist/ASVs_counts.tsv", sep = "\t")
samp_name <- colnames(count_protist)


count_bac <- read.csv("bac/ASVs_counts.tsv", sep = "\t")

bac_count <- count_bac %>% 
  column_to_rownames(., var = "X" ) %>% 
  select(!any_of(samp_name)) %>% 
  select(!starts_with("con")) %>% 
  filter(rowSums(.)> 2000) %>% 
  t() %>% 
  as.data.frame() %>% 
  filter(rowSums(.)> 5000)

resampe_count <- min(rowSums(bac_count))

dist_bray <- bac_count %>% 
  avgdist(., sample = resampe_count, meanfun = median, transf = "sqrt", iterations = 100)

long_dist <- dist_bray %>% 
  as.matrix(.) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "rw_name") %>% 
  pivot_longer(-rw_name,
              names_to = "col_name")

set.seed(1234)
mds <- metaMDS(dist_bray, trymax = 999)
plot(mds)  
text(mds)

bac_mds_score <- scores(mds) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "sample_id") %>% 
  separate(sample_id, into = c("season", "site_id", "treatment")) %>% 
  mutate(site = str_remove(site_id, pattern = "[0-9]+$"),
         season = factor(season, level = c("S", "F", "Sp"), ordered = T),
         site = factor(site, level = c( "GE", "GC", "GD", "GA", "GB", "GF", "LL"),  ordered = T))

bac_mds_summer <- bac_mds_score %>% 
  filter(season == "S")

trt <- c("PMA untreated", "PMA treated")
season <- c("Summer", "Fall", "Spring")
names(trt) <- c("NT", "T")
names(season) <- c("S", "F", "Sp")

ggplot(bac_mds_score, aes(x = NMDS1, y = NMDS2, fill = site)) +
  geom_point(pch = 21, size = 5) +
  facet_grid(season ~ treatment,
             labeller = labeller(treatment = trt, season = season)) +
  scale_fill_viridis_d() +
  theme_bw()

bac_beta_div <- ggplot(bac_mds_summer, aes(x = NMDS1, y = NMDS2, fill = site)) +
  geom_point(pch = 21, size = 5) +
  facet_grid(~ treatment,
             labeller = labeller(treatment = trt, season = season)) +
  scale_fill_viridis_d(name = "Reclamation\nage",
                       labels = c("3", "6", "7", "15","21", "26", "reference")
                       ) +
  theme_bw()+
  theme(axis.title.y = element_text(size = 16, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        strip.text = element_text(size = 14, color = "black")) +
  labs(title = "Bacteria")

bac_beta_div

# Alphadiv ----------------------------------------------------------------

alphadiv <- bac_count %>% 
  mutate(shann = diversity(., index = "shannon"),
         invsimp = diversity(., index = "invsimpson"),
         richness = specnumber(.)) %>% 
  rownames_to_column(., var = "sample_id") %>% 
  separate(sample_id, into = c("season", "site_id", "treatment")) %>% 
  mutate(site = str_remove(site_id, pattern = "[0-9]+$"),
         season = factor(season, level = c("S", "F", "Sp"), ordered = T))

alphadiv_summer <- alphadiv %>% 
  filter(season == "S")

library(lme4)
library(lmerTest)

m1 <- glmer(richness ~ treatment + (1|site), family = "poisson",
          data = alphadiv_summer)         
car::Anova(m1)
DHARMa::simulateResiduals(m1, plot = TRUE)
hist(resid(m1))


bac_div <- ggplot(alphadiv_summer, aes(x = treatment, y = richness)) +
  geom_boxplot(aes(fill = treatment)) +
  ggpubr::geom_bracket(xmin = "NT", xmax = "T",
                       y.position = 375,
                       label = "***") +
  theme_bw() +
  scale_x_discrete(label = c("PMA untreated", "PMA treated")) +
  ggsci::scale_fill_jama() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text = element_text(size = 14, color = "black")) +
  labs(y = "Observed species\nrichness",
       title = "Bacteria")


# fungal diversity --------------------------------------------------------

count_fungi <- read.csv("fungi/ASVs_counts.tsv", sep = "\t")

fungi_count <- count_fungi %>% 
  column_to_rownames(., var = "X" ) %>% 
  select(!any_of(samp_name)) %>% 
  select(!starts_with("con")) %>%
  filter(rowSums(.)> 2000) %>% 
  t() %>% 
  as.data.frame() %>% 
  filter(rowSums(.)> 5000)


alphadiv_fungi <- fungi_count %>% 
  mutate(shann = diversity(., index = "shannon"),
         invsimp = diversity(., index = "invsimpson"),
         richness = specnumber(.)) %>% 
  rownames_to_column(., var = "sample_id") %>% 
  separate(sample_id, into = c("season", "site_id", "treatment")) %>% 
  mutate(site = str_remove(site_id, pattern = "[0-9]+$"),
         season = factor(season, level = c("S", "F", "Sp"), ordered = T)) %>% 
  drop_na()

alphadiv_summer_fungi <- alphadiv_fungi %>% 
  filter(season == "S")


m1 <- glmer(richness ~ treatment + (1|site), family = "poisson",
            data = alphadiv_summer_fungi)         
car::Anova(m1)
DHARMa::simulateResiduals(m1, plot = TRUE)
hist(resid(m1))

fungi_div <- ggplot(alphadiv_summer_fungi, aes(x = treatment, y = richness)) +
  geom_boxplot(aes(fill = treatment)) +
  ggpubr::geom_bracket(xmin = "NT", xmax = "T",
                       y.position = 110,
                       label = "***") +
  theme_bw() +
  scale_x_discrete(label = c("PMA untreated", "PMA treated")) +
  ggsci::scale_fill_jama() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text = element_text(size = 14, color = "black")) +
  labs(y = "Observed species\nrichness",
       title = "Fungi")


library(patchwork)

micro_div <- bac_div / fungi_div

ggsave(micro_div, file = "alpha_div.jpg", width = 6, height = 8, units = "in", dpi = 800)


# fungal betadiversity ----------------------------------------------------

resampe_count_fungi <- min(rowSums(fungi_count))

dist_bray_fungi <- fungi_count %>% 
  avgdist(., sample = resampe_count_fungi, meanfun = median, transf = "sqrt", iterations = 100)


set.seed(12345)
mds_fungi <- metaMDS(dist_bray_fungi,k=3, trymax = 999)
plot(mds_fungi)  
text(mds_fungi)



fungi_mds_score <- scores(mds_fungi) %>% 
  as.data.frame() %>% 
  rownames_to_column(., var = "sample_id") %>% 
  separate(sample_id, into = c("season", "site_id", "treatment")) %>% 
  mutate(site = str_remove(site_id, pattern = "[0-9]+$"),
         season = factor(season, level = c("S", "F", "Sp"), ordered = T),
         site = factor(site, level = c( "GE", "GC", "GD", "GA", "GB", "GF", "LL"),  ordered = T))

fungi_mds_summer <- fungi_mds_score %>% 
  filter(season == "S")

trt <- c("PMA untreated", "PMA treated")
season <- c("Summer", "Fall", "Spring")
names(trt) <- c("NT", "T")
names(season) <- c("S", "F", "Sp")

ggplot(fungi_mds_score, aes(x = NMDS1, y = NMDS2, fill = site)) +
  geom_point(pch = 21, size = 5) +
  facet_grid(season ~ treatment,
             labeller = labeller(treatment = trt, season = season)) +
  scale_fill_viridis_d() +
  theme_bw() +
  theme()

fungi_beta_div <- ggplot(fungi_mds_summer, aes(x = NMDS1, y = NMDS2, fill = site)) +
  geom_point(pch = 21, size = 5) +
  facet_grid(~ treatment,
             labeller = labeller(treatment = trt, season = season)) +
  scale_fill_viridis_d(name = "Reclamation\nage",
                       labels = c("3", "6", "7", "15","21", "26", "reference")
  ) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 16, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        strip.text = element_text(size = 14, color = "black")) +
  labs(title = "Fungi")
fungi_beta_div

beta_div_micro <- bac_beta_div/fungi_beta_div +  plot_layout(guides = "collect")

ggsave(beta_div_micro, file = "beta_div_micro.jpg", width = 8, height = 8, units = "in", dpi = 800)
