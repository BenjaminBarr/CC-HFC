install.packages("ggdark")
install.packages("ggplot2")
install.packages("readxl")
install.packages("reshape2")
install.packages("matrixStats")
install.packages("dplyr")
install.packages("AICcmodavg")
install.packages("tidyverse")
install.packages("ggpubr")
install.packages("rstatix")
install.packages("emmeans")
install.packages("lme4")



#install.packages(c("lubridate", "ggsurvfit", "gtsummary", "tidycmprsk"))
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)

#devtools::install_github("zabore/condsurv")
#library(condsurv)

# Data Read In

library(ggplot2)
library(readxl)
library(reshape2)
library(matrixStats)
library(dplyr)
library(AICcmodavg)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(ggdark)
library(emmeans)
library(lme4)
library(survival)
library(knitr)
library(tibble)
######### From Scratch ###############

## G26 KM

g26.km <- read.csv("./g26.live.csv")

g26.model <- survfit(data = g26.km, Surv(week, death) ~ label)
g26.km$label <- as.factor(g26.km$label)

survfit2(Surv(week, death) ~ label, data = g26.km) %>%
  ggsurvfit(linewidth = 2, linetype_aes = T) +
  scale_color_manual(values = c("brown1", "chartreuse4", "deepskyblue3", "darkgoldenrod1"), 
                     labels = c("CC F", "CC M", "HFC F", "HFC M")) +
  labs(x = "Weeks",
       y = "Overall Survival Probability",
       title = "Survival Assessment") +
  scale_ggsurvfit(y_scales = list(breaks = seq(.2, 1, by = .2))) +
  guides(linetype = "none") +
  add_pvalue(location = "annotation", y = .25, x = 40, caption = "Log-rank {p.value}", size = 5) +
  theme_ggsurvfit_KMunicate()+
  theme(axis.title = element_text(size = 16),
        title = element_text(size = 20),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14))

  
survdiff(formula = Surv(week, death) ~ label, data = g26.km)
coxph(formula = Surv(week, death) ~ label, data = g26.km) %>%
  tbl_regression(exp = T)

survfit(Surv(week, death) ~ label, data = g26.km) %>% 
  tbl_survfit(
    times = 73,
    label_header = "**18-month survival (95% CI)**"
  )
summary(g26.km)


