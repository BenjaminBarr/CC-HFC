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

## Read in Data
g26 <- read.csv("./Data/CC&HFC Weekly Mass.csv")

# Calculate weekly mean
g26.avg <- g26 %>%
  group_by(sex, group, week) %>%
  summarise("mean" = mean(weight), "sd" = sd(weight))
# Generate color label for weekly means
g26.avg$series <- paste(g26.avg$group," ",g26.avg$sex,"", sep = "")

# Plot the data using LOESS
ggplot(data = g26, aes(week, weight, group = series, col = series)) +
  geom_smooth() +
  geom_point(data = g26.avg, aes(week, mean, group = series, col = series)) +
  scale_x_continuous(limits = c(0, 73),
                     breaks = c(4, 24, 48, 72)) +
  scale_color_manual(label = c("CC Females", "CC Males", "HFC Females", "HFC Males"),
                     values = c("brown1", "chartreuse4", "deepskyblue3", "darkgoldenrod1")) +
  labs(title = "Weekly Total Mass",
       y = "Median Mass (g)", x = "Week", color = "Diet") +
  theme(plot.title = element_text(family = "Fira Sans Condensed", hjust = 0.5, size = 20),
        plot.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text = element_text(size = 15, color = "black"),
        axis.ticks = element_blank(),
        legend.key = element_blank(),
        legend.position = "right")
