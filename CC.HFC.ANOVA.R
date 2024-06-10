library(here) # here makes a project transportable
library(janitor) # clean_names
library(readxl) # read excel, duh!
library(data.table) # magical data frames
library(magrittr) # pipes
library(stringr) # string functions
library(forcats) # factor functions

# analysis packages
library(emmeans) # the workhorse for inference
library(nlme) # gls and some lmm
library(lme4) # linear mixed models
library(lmerTest) # linear mixed model inference
library(afex) # ANOVA linear models
library(glmmTMB) # generalized linear models
library(MASS) # negative binomial and some other functions
library(car) # model checking and ANOVA
library(DHARMa) # model checking

# graphing packages
library(ggsci) # color palettes
library(ggpubr) # publication quality plots
library(ggforce) # better jitter

library(knitr) # kable tables
library(kableExtra) # kable_styling tables

library(insight)
library(lazyWeave)

library(cowplot) # combine plots
library(patchwork)# combine plots
# Data Read In

library(reshape2)
library(matrixStats)
library(dplyr)
library(AICcmodavg)
library(rstatix)
library(ggdark)
library(multcompView)
library(tidyverse)
#STRIPED BARS
library(ggpattern)
##Functions for normalization
source_path <- here("raw-R work", "ggplot_the_model.R")
source(source_path)

## Read in MRI Data

c_dt <- read.csv("./Data/CC&HFC MRI.csv")

# Get summary of MRI data
c.summary <- c_dt %>% group_by(group, m.start) %>%
  summarise(mean_mass = mean(weight),
            se_mass = sd(weight)/sqrt(length(weight)),
            mean_fat = mean(fat),
            se_fat = sd(fat)/sqrt(length(fat)),
            mean_lean = mean(lean),
            se_lean = sd(lean)/sqrt(length(lean)))

# All fixed-effects models for total mass
c6.lm <- lm(weight ~ fatc * sex, data = c_dt[which(c_dt$m.start == 6),])
c12.lm <- lm(weight ~ fatc * sex, data = c_dt[which(c_dt$m.start == 12),])
c18.lm <- lm(weight ~ fatc * sex, data = c_dt[which(c_dt$m.start == 18),])

# 6 months first
anova(c6.lm)
summary(c6.lm)

# Estimated marginal means (emmeans)
c6.em <- emmeans(c6.lm, specs = c("fatc", "sex"))
c6.em

c6.em_dt <- data.table(summary(c6.em))
c6.em_dt

# Contrast "pairwise comparisions"
c6.pairs <- contrast(c6.em,
                     method = "revpairwise",
                     simple = "each",
                     combine = T,
                     adjust = "sidak") %>%
  summary(infer = T, )

c6.pairs_dt <- data.table(c6.pairs)

# Coordinates for significance indicators
group1 <- c6.pairs_dt$contrast

group1 <- gsub(" - fatc0", "", as.character(group1))
group1 <- gsub(" - sex0", "", as.character(group1))
group1 <- as.factor(group1)

group2 <- c6.pairs_dt$contrast
group2 <- gsub("fatc1 - ", "", as.character(group2))
group2 <- gsub("sex1 - ", "", as.character(group2))
group2 <- as.factor(group2)
group2

c6.pairs_dt$group1 <- group1
c6.pairs_dt$group2 <- group2
c6.pairs_dt

# Nice p-values for graphs
c6.pairs_dt[, p_rounded := p_round(p.value,
                                   digits = 2)]
c6.pairs_dt[, p_pretty := p_format(p_rounded,
                                   digits = 2,
                                   accuracy = 1e-04,
                                   add.p = TRUE)]
c6.pairs_dt

# Assumptions and Residuals
c6.lm_plot <- ggplot_the_model(
  fit = c6.lm,
  fit_emm = c6.em,
  fit_pairs = c6.pairs,
  palette = pal_okabe_ito_blue,
  y_label = "Weight",
  g_label = "sex",
  dots = "jitter"
  
)
c6.lm_plot

# Graphing labels (group 2 = CC | group 6 = HFC, f = Female | m = Male)
c6.em_dt$group <- c("group 2 f", "group 6 f", "group 2 m", "group 6 m")
c6.pairs_dt$group1 <- c("group 6 f", "group 6 m", "group 2 m", "group 6 m")
c6.pairs_dt$group2 <- c("group 2 f", "group 2 m", "group 2 f", "group 6 f")

#Plot first time-point total mass

ggbarplot(data = c6.em_dt, 
          x = "group", 
          y = "emmean", 
          fill = "group") +
  theme_pubr() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  geom_errorbar(data = c6.em_dt, aes(y = emmean, ymin = lower.CL, ymax = upper.CL), width = 0.1, color = "black") +
  stat_pvalue_manual(data = c6.pairs_dt,
                     label = "p_pretty",
                     xmin = "group1",
                     xmax = "group2",
                     y.position = c(seq(from = 50, to = 65, by = 5)),
                     tip.length = 0.01) +
  labs(title = "Month 6 Mass", y ="Weight")

# Repeat for 12 and 18
anova(c12.lm)

# Estimated marginal means (emmeans)
c12.em <- emmeans(c12.lm, specs = c("fatc", "sex"))
c12.em

c12.em_dt <- data.table(summary(c12.em))
c12.em_dt

# Contrast "pairwise comparisions"
c12.pairs <- contrast(c12.em,
                      method = "revpairwise",
                      simple = "each",
                      combine = T,
                      adjust = "sidak") %>%
  summary(infer = T, )

c12.pairs_dt <- data.table(c12.pairs)

# Assumptions and Residuals
c12.lm_plot <- ggplot_the_model(
  fit = c12.lm,
  fit_emm = c12.em,
  fit_pairs = c12.pairs,
  palette = pal_okabe_ito_blue,
  y_label = "Weight",
  g_label = "sex",
  dots = "jitter"
  
)
c12.lm_plot

###
anova(c18.lm)

c18.em <- emmeans(c18.lm, specs = c("fatc", "sex"))
c18.em

c18.em_dt <- data.table(summary(c18.em))
c18.em_dt

c18.pairs <- contrast(c18.em,
                      method = "revpairwise",
                      simple = "each",
                      combine = T,
                      adjust = "sidak") %>%
  summary(infer = T, )

c18.pairs_dt <- data.table(c18.pairs)

c18.lm_plot <- ggplot_the_model(
  fit = c18.lm,
  fit_emm = c18.em,
  fit_pairs = c18.pairs,
  palette = pal_okabe_ito_blue,
  y_label = "Weight",
  g_label = "sex",
  dots = "jitter"
  
)
c18.lm_plot

## ALL MASS ##
c.em_dt <- rbind(c6.em_dt, c12.em_dt, c18.em_dt, fill = T)
# Add missing graphing labels
c.em_dt$group <- rep(c("group 2 f", "group 6 f", "group 2 m", "group 6 m"), times = 3)
# Add age
c.em_dt$age <- rep(c(6, 12, 18), each = 4)

# All mass significance bars
c.pairs_dt <- rbind(c6.pairs_dt, c12.pairs_dt, c18.pairs_dt, fill = T)
c.pairs_dt$age <- rep(c(6, 12, 18), each = 4)
c.pairs_dt$group1 <- rep(c("group 6 f", "group 6 m", "group 2 m", "group 6 m"), times = 3)
c.pairs_dt$group2 <- rep(c("group 2 f", "group 2 m", "group 2 f", "group 6 f"), times = 3)
c.pairs_dt[, p_rounded := p_round(p.value,
                                   digits = 2)]
c.pairs_dt[, p_pretty := p_format(p_rounded,
                                   digits = 2,
                                   accuracy = 1e-04,
                                   add.p = TRUE)]
# X-axis label
c.em_dt$xlab <- paste(c.em_dt$group," ",c.em_dt$age,"", sep = "")

# Significance bar coordinates
c.pairs_dt$xlab1 <- paste(c.pairs_dt$group1, c.pairs_dt$age)
c.pairs_dt$xlab2 <- paste(c.pairs_dt$group2, c.pairs_dt$age)

# Use only those that ar significant
c.pairs.sig <- c.pairs_dt[which(c.pairs_dt$p.value <= 0.05),]

# Graph by diet and sex
c.em_dt <- c.em_dt[order(c.em_dt$group),]

#All timepoints total mass
ggbarplot(data = c.em_dt, 
          x = "xlab", 
          y = "emmean", 
          fill = "group") +
  theme_pubr() +
  geom_errorbar(data = c.em_dt, aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE), width = 0.1, color = "black") +
  stat_pvalue_manual(data = c.pairs.sig,
                     label = "p_pretty",
                     xmin = "xlab1",
                     xmax = "xlab2",
                     y.position = c(seq(from = 50, to = 85, by = 5)),
                     tip.length = 0.01) +
  scale_x_discrete(labels = rep(c(6, 12, 18), times = 4)) +
  scale_fill_discrete(name = "Diets:", labels = c("CC F", "CC M", "HFC F", "HFC M")) +
  labs(title = "C3H/HeJ Total Mass Over Time", y ="Mass (g)", x = "Month")

### Repeat for lean Mass ###
# Begin at fixed-effects model

cl6.lm <- lm(lean ~ fatc * sex, data = c_dt[which(c_dt$m.start == 6),])
cl12.lm <- lm(lean ~ fatc * sex, data = c_dt[which(c_dt$m.start == 12),])
cl18.lm <- lm(lean ~ fatc * sex, data = c_dt[which(c_dt$m.start == 18),])

# 6 months first
anova(cl6.lm)
summary(cl6.lm)

# Estimated marginal means (emmeans)
cl6.em <- emmeans(cl6.lm, specs = c("fatc", "sex"))
cl6.em

cl6.em_dt <- data.table(summary(cl6.em))
cl6.em_dt

# Contrast "pairwise comparisions"
cl6.pairs <- contrast(cl6.em,
                     method = "revpairwise",
                     simple = "each",
                     combine = T,
                     adjust = "sidak") %>%
  summary(infer = T, )

cl6.pairs_dt <- data.table(cl6.pairs)

# Coordinates for significance indicators
cl6.pairs_dt$group1 <- group1
cl6.pairs_dt$group2 <- group2
cl6.pairs_dt

# Nice p-values for graphs
cl6.pairs_dt[, p_rounded := p_round(p.value,
                                   digits = 2)]
cl6.pairs_dt[, p_pretty := p_format(p_rounded,
                                   digits = 2,
                                   accuracy = 1e-04,
                                   add.p = TRUE)]
cl6.pairs_dt

# Assumptions and Residuals
cl6.lm_plot <- ggplot_the_model(
  fit = cl6.lm,
  fit_emm = cl6.em,
  fit_pairs = cl6.pairs,
  palette = pal_okabe_ito_blue,
  y_label = "Weight",
  g_label = "sex",
  dots = "jitter"
  
)
cl6.lm_plot

# Graphing labels (group 2 = CC | group 6 = HFC, f = Female | m = Male)
cl6.em_dt$group <- c("group 2 f", "group 6 f", "group 2 m", "group 6 m")
cl6.pairs_dt$group1 <- c("group 6 f", "group 6 m", "group 2 m", "group 6 m")
cl6.pairs_dt$group2 <- c("group 2 f", "group 2 m", "group 2 f", "group 6 f")

#Plot first time-point lean mass

ggbarplot(data = cl6.em_dt, 
          x = "group", 
          y = "emmean", 
          fill = "group") +
  theme_pubr() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  geom_errorbar(data = cl6.em_dt, aes(y = emmean, ymin = lower.CL, ymax = upper.CL), width = 0.1, color = "black") +
  stat_pvalue_manual(data = cl6.pairs_dt,
                     label = "p_pretty",
                     xmin = "group1",
                     xmax = "group2",
                     y.position = c(seq(from = 50, to = 65, by = 5)),
                     tip.length = 0.01) +
  labs(title = "Month 6 Lean Mass", y ="Mass")

# Repeat for 12 and 18
anova(cl12.lm)

# Estimated marginal means (emmeans)
cl12.em <- emmeans(cl12.lm, specs = c("fatc", "sex"))
cl12.em

cl12.em_dt <- data.table(summary(cl12.em))
cl12.em_dt

# Contrast "pairwise comparisions"
cl12.pairs <- contrast(c12.em,
                      method = "revpairwise",
                      simple = "each",
                      combine = T,
                      adjust = "sidak") %>%
  summary(infer = T, )

cl12.pairs_dt <- data.table(cl12.pairs)

# Assumptions and Residuals
cl12.lm_plot <- ggplot_the_model(
  fit = cl12.lm,
  fit_emm = cl12.em,
  fit_pairs = cl12.pairs,
  palette = pal_okabe_ito_blue,
  y_label = "Weight",
  g_label = "sex",
  dots = "jitter"
  
)
cl12.lm_plot

###
anova(cl18.lm)

cl18.em <- emmeans(cl18.lm, specs = c("fatc", "sex"))
cl18.em

cl18.em_dt <- data.table(summary(cl18.em))
c18.em_dt

cl18.pairs <- contrast(cl18.em,
                      method = "revpairwise",
                      simple = "each",
                      combine = T,
                      adjust = "sidak") %>%
  summary(infer = T, )

cl18.pairs_dt <- data.table(cl18.pairs)

cl18.lm_plot <- ggplot_the_model(
  fit = cl18.lm,
  fit_emm = cl18.em,
  fit_pairs = cl18.pairs,
  palette = pal_okabe_ito_blue,
  y_label = "Weight",
  g_label = "sex",
  dots = "jitter"
  
)
cl18.lm_plot

## ALL LEAN MASS ##
cl.em_dt <- rbind(cl6.em_dt, cl12.em_dt, cl18.em_dt, fill = T)
# Add missing graphing labels
cl.em_dt$group <- rep(c("group 2 f", "group 6 f", "group 2 m", "group 6 m"), times = 3)
# Add age
cl.em_dt$age <- rep(c(6, 12, 18), each = 4)

# All mass significance bars
cl.pairs_dt <- rbind(cl6.pairs_dt, cl12.pairs_dt, cl18.pairs_dt, fill = T)
cl.pairs_dt$age <- rep(c(6, 12, 18), each = 4)
cl.pairs_dt$group1 <- rep(c("group 6 f", "group 6 m", "group 2 m", "group 6 m"), times = 3)
cl.pairs_dt$group2 <- rep(c("group 2 f", "group 2 m", "group 2 f", "group 6 f"), times = 3)
cl.pairs_dt[, p_rounded := p_round(p.value,
                                  digits = 2)]
cl.pairs_dt[, p_pretty := p_format(p_rounded,
                                  digits = 2,
                                  accuracy = 1e-04,
                                  add.p = TRUE)]
# X-axis label
cl.em_dt$xlab <- paste(cl.em_dt$group," ",cl.em_dt$age,"", sep = "")

# Significance bar coordinates
cl.pairs_dt$xlab1 <- paste(cl.pairs_dt$group1, cl.pairs_dt$age)
cl.pairs_dt$xlab2 <- paste(cl.pairs_dt$group2, cl.pairs_dt$age)

# Use only those that ar significant
cl.pairs.sig <- cl.pairs_dt[which(cl.pairs_dt$p.value <= 0.05),]

# Graph by diet and sex
cl.em_dt <- cl.em_dt[order(cl.em_dt$group),]

#All timepoints lean mass
ggbarplot(data = cl.em_dt, 
          x = "xlab", 
          y = "emmean", 
          fill = "group") +
  theme_pubr() +
  geom_errorbar(data = cl.em_dt, aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE), width = 0.1, color = "black") +
  stat_pvalue_manual(data = cl.pairs.sig,
                     label = "p_pretty",
                     xmin = "xlab1",
                     xmax = "xlab2",
                     y.position = c(seq(from = 50, to = 95, by = 5)), #Adjusted based on # of sig
                     tip.length = 0.01) +
  scale_x_discrete(labels = rep(c(6, 12, 18), times = 4)) +
  scale_fill_discrete(name = "Diets:", labels = c("CC F", "CC M", "HFC F", "HFC M")) +
  labs(title = "C3H/HeJ Lean Mass Over Time", y ="Mass (g)", x = "Month")

### Repeat for fat mass ###

cf6.lm <- lm(fat ~ fatc * sex, data = c_dt[which(c_dt$m.start == 6),])
cf12.lm <- lm(fat ~ fatc * sex, data = c_dt[which(c_dt$m.start == 12),])
cf18.lm <- lm(fat ~ fatc * sex, data = c_dt[which(c_dt$m.start == 18),])

# 6 months first
anova(cf6.lm)
summary(cf6.lm)

# Estimated marginal means (emmeans)
cf6.em <- emmeans(cf6.lm, specs = c("fatc", "sex"))
cf6.em

cf6.em_dt <- data.table(summary(cf6.em))
cf6.em_dt

# Contrast "pairwise comparisions"
cf6.pairs <- contrast(cf6.em,
                      method = "revpairwise",
                      simple = "each",
                      combine = T,
                      adjust = "sidak") %>%
  summary(infer = T, )

cf6.pairs_dt <- data.table(cf6.pairs)

# Coordinates for significance indicators
cf6.pairs_dt$group1 <- group1
cf6.pairs_dt$group2 <- group2
cf6.pairs_dt

# Nice p-values for graphs
cf6.pairs_dt[, p_rounded := p_round(p.value,
                                    digits = 2)]
cf6.pairs_dt[, p_pretty := p_format(p_rounded,
                                    digits = 2,
                                    accuracy = 1e-04,
                                    add.p = TRUE)]
cf6.pairs_dt

# Assumptions and Residuals
cf6.lm_plot <- ggplot_the_model(
  fit = cf6.lm,
  fit_emm = cf6.em,
  fit_pairs = cf6.pairs,
  palette = pal_okabe_ito_blue,
  y_label = "Weight",
  g_label = "sex",
  dots = "jitter"
  
)
cf6.lm_plot

# Graphing labels (group 2 = CC | group 6 = HFC, f = Female | m = Male)
cf6.em_dt$group <- c("group 2 f", "group 6 f", "group 2 m", "group 6 m")
cf6.pairs_dt$group1 <- c("group 6 f", "group 6 m", "group 2 m", "group 6 m")
cf6.pairs_dt$group2 <- c("group 2 f", "group 2 m", "group 2 f", "group 6 f")

#Plot first time-point total mass

ggbarplot(data = cf6.em_dt, 
          x = "group", 
          y = "emmean", 
          fill = "group") +
  theme_pubr() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  geom_errorbar(data = cf6.em_dt, aes(y = emmean, ymin = lower.CL, ymax = upper.CL), width = 0.1, color = "black") +
  stat_pvalue_manual(data = cf6.pairs_dt,
                     label = "p_pretty",
                     xmin = "group1",
                     xmax = "group2",
                     y.position = c(seq(from = 50, to = 65, by = 5)),
                     tip.length = 0.01) +
  labs(title = "Month 6 Lean Mass", y ="Mass")

# Repeat for 12 and 18
anova(cf12.lm)

# Estimated marginal means (emmeans)
cf12.em <- emmeans(cf12.lm, specs = c("fatc", "sex"))
cf12.em

cf12.em_dt <- data.table(summary(cf12.em))
cf12.em_dt

# Contrast "pairwise comparisions"
cf12.pairs <- contrast(c12.em,
                       method = "revpairwise",
                       simple = "each",
                       combine = T,
                       adjust = "sidak") %>%
  summary(infer = T, )

cf12.pairs_dt <- data.table(cf12.pairs)

# Assumptions and Residuals
cf12.lm_plot <- ggplot_the_model(
  fit = cf12.lm,
  fit_emm = cf12.em,
  fit_pairs = cf12.pairs,
  palette = pal_okabe_ito_blue,
  y_label = "Weight",
  g_label = "sex",
  dots = "jitter"
  
)
cf12.lm_plot

###
anova(cf18.lm)

cf18.em <- emmeans(cf18.lm, specs = c("fatc", "sex"))
cf18.em

cf18.em_dt <- data.table(summary(cf18.em))
c18.em_dt

cf18.pairs <- contrast(cf18.em,
                       method = "revpairwise",
                       simple = "each",
                       combine = T,
                       adjust = "sidak") %>%
  summary(infer = T, )

cf18.pairs_dt <- data.table(cf18.pairs)

cf18.lm_plot <- ggplot_the_model(
  fit = cf18.lm,
  fit_emm = cf18.em,
  fit_pairs = cf18.pairs,
  palette = pal_okabe_ito_blue,
  y_label = "Weight",
  g_label = "sex",
  dots = "jitter"
  
)
cf18.lm_plot

## ALL MASS ##
cf.em_dt <- rbind(cf6.em_dt, cf12.em_dt, cf18.em_dt, fill = T)
# Add missing graphing labels
cf.em_dt$group <- rep(c("group 2 f", "group 6 f", "group 2 m", "group 6 m"), times = 3)
# Add age
cf.em_dt$age <- rep(c(6, 12, 18), each = 4)

# All mass significance bars
cf.pairs_dt <- rbind(cf6.pairs_dt, cf12.pairs_dt, cf18.pairs_dt, fill = T)
cf.pairs_dt$age <- rep(c(6, 12, 18), each = 4)
cf.pairs_dt$group1 <- rep(c("group 6 f", "group 6 m", "group 2 m", "group 6 m"), times = 3)
cf.pairs_dt$group2 <- rep(c("group 2 f", "group 2 m", "group 2 f", "group 6 f"), times = 3)
cf.pairs_dt[, p_rounded := p_round(p.value,
                                   digits = 2)]
cf.pairs_dt[, p_pretty := p_format(p_rounded,
                                   digits = 2,
                                   accuracy = 1e-04,
                                   add.p = TRUE)]
# X-axis label
cf.em_dt$xlab <- paste(cf.em_dt$group," ",cf.em_dt$age,"", sep = "")

# Significance bar coordinates
cf.pairs_dt$xlab1 <- paste(cf.pairs_dt$group1, cf.pairs_dt$age)
cf.pairs_dt$xlab2 <- paste(cf.pairs_dt$group2, cf.pairs_dt$age)

# Use only those that ar significant
cf.pairs.sig <- cf.pairs_dt[which(cf.pairs_dt$p.value <= 0.05),]

# Graph by diet and sex
cf.em_dt <- cf.em_dt[order(cf.em_dt$group),]

#All timepoints fat mass
ggbarplot(data = cf.em_dt, 
          x = "xlab", 
          y = "emmean", 
          fill = "group") +
  theme_pubr() +
  geom_errorbar(data = cf.em_dt, aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE), width = 0.1, color = "black") +
  stat_pvalue_manual(data = cf.pairs.sig,
                     label = "p_pretty",
                     xmin = "xlab1",
                     xmax = "xlab2",
                     y.position = c(seq(from = 50, to = 80, by = 5)), #Adjusted based on # of sig
                     tip.length = 0.01) +
  scale_x_discrete(labels = rep(c(6, 12, 18), times = 4)) +
  scale_fill_discrete(name = "Diets:", labels = c("CC F", "CC M", "HFC F", "HFC M")) +
  labs(title = "C3H/HeJ Fat Mass Over Time", y ="Mass (g)", x = "Month")

### Combine all mass types
c.em_dt$cat <- "Total Mass (g)"
cl.em_dt$cat <- "Lean Mass (g)"
cf.em_dt$cat <- "Fat Mass (g)"
all.em_dt <- rbind(c.em_dt, cl.em_dt, cf.em_dt)

all.em_dt$group <- as.factor(all.em_dt$group)

all.em_dt <- all.em_dt[order(all.em_dt$age),]

all.em_dt$xlab <- paste(all.em_dt$age, rep(c("T", "L", "F"), each = 4, times = 3))

all.em_dt$xlab <- as.factor(all.em_dt$xlab)

is.numeric(all.em_dt$emmean)

glab <- c(`group 2 f` = "CC Females", `group 2 m` = "CC Males", 
          `group 6 f` = "HFC Females", `group 6 m` = "HFC Males")

names(glab) <-levels(all.em_dt$group)
names(glab) <- as.factor(names(glab))

names(glab)
is.factor(all.em_dt$group)

all.em_dt$group <- as.factor(all.em_dt$group)

gcolor <- c("brown1", "deepskyblue3", "chartreuse4", "darkgoldenrod1")

#### Plot Combining ####

all.em_dt$cat <- as.factor(all.em_dt$cat)
## Legend order ##
all.em_dt$cat <- factor(all.em_dt$cat, levels = rev(levels(all.em_dt$cat)))
all.em_dt$group <- factor(all.em_dt$group, levels = c("group 2 f", "group 6 f", 
                                                      "group 2 m", "group 6 m"))
levels(all.em_dt$group)

#geom_col_pattern(aes(pattern = cat, fill = group), 
#color = 'black', pattern_fill = "white",
#pattern_angle = 45,
#pattern_density = rep(c(.35, .35, .75), times = 12)

m.types <- ggplot(data = all.em_dt, aes(x = factor(xlab, levels = c("6 T", "6 L", "6 F",
                                                         "12 T", "12 L", "12 F",
                                                         "18 T", "18 L", "18 F")),
                             fill = group,
                             y = emmean)) +
  geom_col_pattern(aes(pattern = cat, pattern_angle = cat, fill = group), 
                   colour = 'black',
                   fill = 'brown1', # Flip fill and pattern_fill for then cut/paste from 
                   pattern_fill = 'white', #total finals for lean
                   pattern_color = "black", pattern_angle = 45,
                   pattern_spacing = 0.05, # Set key density
                   pattern_density = .75, # Graph density
                   pattern_key_scale_factor = 0.3) +
  scale_pattern_discrete(choices = c("none", "stripe", "circle"), name = "Diets:") +
  scale_pattern_spacing_discrete(range = c(0.01, 0.1)) +
  facet_wrap(~factor(group, levels = c("group 2 f", "group 6 f", 
                                       "group 2 m", "group 6 m")),
             labeller = as_labeller(glab)) +
  scale_x_discrete(labels = rep(c(6, 12, 18), times = 4)) +
  geom_errorbar(data = all.em_dt, aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE, x = xlab), width = 0.1, color = "black") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 52)) +
  theme_pubr() +
  theme(legend.box = "vertical") +
  labs(title = "C3H/HeJ Echo MRI Analysis", y ="Mass (g)", x = "Age (Months)")

m.types

types = get_plot_component(last_plot(), "guide-box-top") #extract legend

#No legend
plt <- ggplot(data = all.em_dt, aes(x = factor(xlab, levels = c("6 T", "6 L", "6 F",
                                                                "12 T", "12 L", "12 F",
                                                                "18 T", "18 L", "18 F")), 
                                    y = emmean)) +
  geom_col_pattern(aes(pattern = cat, fill = group), 
                   color = 'black', pattern_fill = "white",
                   pattern_angle = 45,
                   pattern_density = rep(c(.35, .35, .75), times = 12)) +
  scale_pattern_discrete(choices = c("none", "stripe", "circle"), name = "Mass Type") +
  scale_pattern_spacing_discrete(range = c(0.01, 0.1)) +
  facet_wrap(~factor(group, levels = c("group 2 f", "group 6 f", 
                                       "group 2 m", "group 6 m")),
             labeller = as_labeller(glab)) +
  scale_fill_manual(values = gcolor, name = "Diets", labels = c("CC Females", "CC Males", 
                                                                "HFC Females", "HFC Males")) +
  geom_errorbar(data = all.em_dt, aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE, x = xlab), width = 0.1, color = "black") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 52)) +
  scale_x_discrete(labels = rep(c(6, 12, 18), each = 3)) +
  guides(pattern_angle = "none") +
  theme_pubr() +
  theme(legend.position = "none",
        title = element_text(size = 16),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  labs(title = "Echo MRI Analysis", y ="Mass (g)", x = "Age (Months)")

plot_grid(plt, types, nrow = 2, ncol = 1, rel_widths = c(.8, .3), rel_heights = c(.8, .1))
plot_grid(plt, types, plot_spacer(), nrow = 2, ncol = 1, rel_heights = c(4, 1))
