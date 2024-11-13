## Author: David Kottelenberg
## Institute: Wageningen University & Research
## Last edited: 13/11/2024
## Summary: Analyses of 2022-2024 field experiments, weed biomass

rm(list = ls())

library(tidyverse)
library(readxl)
library(RColorBrewer)
library(lme4)
library(car)
library(scales)
library(caret)
library(gridExtra)
library(segmented)
library(lmtest)
library(stringr)
library(ggbeeswarm)
library(xtable)
library(ggpubr)
library(patchwork)

#setwd() # Set working directory to data location

sowingDate2022 <- as.Date("2022-04-19")
sowingDate2023 <- as.Date("2023-03-02")
sowingDate2024 <- as.Date("2023-04-12")

source("WSFWeedFunctions.R")

colourVector <- c(brewer.pal(n=8,"Set1"), "#FDC086")

modelColors <- c("Cereal" = colourVector[1],
                 "Legume" = colourVector[2],
                 "Weed" = colourVector[4],
                 "T" = colourVector[1],
                 "Triticale" = colourVector[1],
                 "F" = colourVector[2],
                 "Faba" = colourVector[2],
                 "TF" = colourVector[3],
                 "C" = colourVector[1],
                 "L" = colourVector[2],
                 "CL" = colourVector[3],
                 "TA" = colourVector[5],
                 "TM" = colourVector[9],
                 "FA" = colourVector[7],
                 "FM" =  colourVector[8],
                 "Intercrop" = colourVector[3],
                 "Sole cereal" = colourVector[1],
                 "Sole legume" = colourVector[2],
                 "Sole triticale" = colourVector[1],
                 "Sole faba" = colourVector[2],
                 "Sole weeds" = colourVector[4],
                 "Weed infested" = colourVector[8],
                 "Herbicide treated" = colourVector[7],
                 "W" = colourVector[4],
                 "Y" = colourVector[8],
                 "N" = colourVector[7])


######################
##                  ##
## Experiment 2022  ##
##                  ##
######################

##
## First harvest
##

data2022_1 <- read_xlsx("WSF_2022_data.xlsx", sheet = "Biomass1", range = "A1:F77", col_names = TRUE) %>% 
  rename("TreatmentN" = "TreatmentN",
         "W" = "BiomassWeed") %>%
  arrange(Plot)

data2022_1Weeds <- data2022_1 %>% 
  dplyr::select(Plot, Treatment, Block, W) %>% 
  na.omit()

weedCLD2022_1 <- calcWeeds(data2022_1Weeds)
data2022_1WeedsGroup <- sapply(data2022_1Weeds$Treatment, function(t)ifelse(t %in% c("Barley", "Rye", "Triticale", "Wheat"), "Sole triticale", ifelse(t %in% c("Pea", "Lupine", "Faba"), "Sole faba", "Intercrop")))

data2022_1WeedsGroup <- data2022_1Weeds %>% 
  mutate(Group = data2022_1WeedsGroup)

# Create plot
(plotWeeds2022_1 <- ggplot(data = data2022_1WeedsGroup) +
    geom_boxplot(aes(x = reorder(Treatment, -W, function(x)mean(x, na.rm = TRUE)), y = W, fill = Group)) +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")"),
         fill = "Crop type") +
    theme_classic(base_size = 30) +
    labs(title = "Weed biomass, 2022",
         subtitle = "First harvest") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    scale_fill_manual(name = "Crop", values = modelColors) +
    annotate("text", x = 9.5, y = 4.2, label = "n.s.", size = 4))

cldWeedsGrouped <- calcWeedsGrouped(data2022_1WeedsGroup)

(plotWeedsBox2022_1 <- ggplot(data = data2022_1WeedsGroup, aes(x = reorder(Group, -W, mean), y = W)) +
    geom_boxplot() +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")")) +
    theme_classic(base_size = 30) +
    labs(title = "Weed biomass, 2022",
         subtitle = "First harvest") +
    annotate("text", x = 1:3, y = 4.5, label = sapply(cldWeedsGrouped$.group, function(x)gsub(" ", "", x)), size = 4))

# Weed prediction
predWeeds2022_1 <- predictWeeds(data2022_1Weeds)

(plotWeedPred2022_1 <- predWeeds2022_1[[2]] +
    labs(subtitle = "2022, first harvest"))

predWeeds2022Long_1 <- predWeeds2022_1[[1]] %>% 
  pivot_longer(cols = c("W", "Arithmetic", "Harmonic"), names_to = "Type", values_to = "WeedBiomass")

newType <- sapply(predWeeds2022Long_1$Type, function(t)ifelse(t == "W", "Observed", t))

predWeeds2022Long_1 <- mutate(predWeeds2022Long_1, Type = newType) %>% 
  mutate(Type = factor(Type, levels = c("Observed", "Arithmetic", "Harmonic")))

(plotWeedPredBox2022_1 <- ggplot(data = predWeeds2022Long_1, aes(x = Type, y = WeedBiomass)) +
    geom_boxplot() +
    labs(x = "", y = bquote("Weed dry weight (g "~m^-2~")")) +
    theme_classic(base_size = 30) +
    labs(title = "Predicted and observed weed biomass, 2022",
         subtitle = "First harvest") +
    annotate("text", x = 2, y = 4.0, label = "n.s.", size = 4))

fit0 <- lm(W ~ 0 + offset(Harmonic), data = predWeeds2022_1[[1]])
fit1 <- lm(W ~ Harmonic, data = predWeeds2022_1[[1]])
fit2 <- lm(W ~ 0 + offset(Arithmetic), data = predWeeds2022_1[[1]])
fit3 <- lm(W ~ Arithmetic, data = predWeeds2022_1[[1]])
lrtest(fit0, fit1) # Harmonic line is not significantly different from 0:1 line, p = 0.5725
lrtest(fit2, fit3) # Arithmetic line is significantly different from 0:1 line, p = 0.01326
lrtest(fit1, fit3) # Harmonic and arithmetic lines are significantly different from each other, p < 2.2e-16

##
## Second harvest
##

# Read data
data2022_2 <- read_xlsx("WSF_2022_data.xlsx", sheet = "Biomass2", range = "A1:F77", col_names = TRUE) %>% 
  dplyr::select(Plot, Treatment, "TreatmentN", "Block", "Date", "BiomassWeed") %>% 
  rename("TreatmentN" = "TreatmentN",
         "WeedWW" = "BiomassWeed",
         "Time" = "Date") %>% 
  mutate(WeedWW = as.double(WeedWW)) %>% 
  dplyr::select(!Time) %>% 
  arrange(Plot)

data2022_2Weeds <- data2022_2 %>% 
  rename(W = WeedWW) %>% 
  dplyr::select(Plot, Treatment, Block, W) %>% 
  na.omit()

weedCLD2022_2 <- calcWeeds(data2022_2Weeds)
data2022_2WeedsGroup <- sapply(data2022_2Weeds$Treatment, function(t)ifelse(t %in% c("Barley", "Rye", "Triticale", "Wheat"), "Sole cereal", ifelse(t %in% c("Pea", "Lupine", "Faba"), "Sole legume", "Intercrop")))

data2022_2WeedsGroup <- data2022_2Weeds %>% 
  mutate(Group = data2022_2WeedsGroup)

# Create plot
(plotWeeds2022_2 <- ggplot(data = data2022_2WeedsGroup) +
    geom_boxplot(aes(x = reorder(Treatment, -W, function(x)mean(x, na.rm = TRUE)), y = W, fill = Group)) +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")"),
         fill = "Crop type") +
    theme_classic(base_size = 30) +
    labs(title = "Weed biomass, 2022",
         subtitle = "Second harvest") +
    scale_fill_manual(name = "Crop", values = modelColors) +
    annotate("text", x = 1:19, y = c(rep(570, 4), seq(540, 210, -30), rep(180, 3)), label = sapply(weedCLD2022_2$.group, function(x)gsub(" ", "", x)), size = 2.5) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))) 

# Weed prediction
predWeeds2022_2 <- predictWeeds(data2022_2Weeds)
(plotWeedPred2022_2 <- predWeeds2022_2[[2]] +
    labs(subtitle = "2022, second harvest"))

cldWeedsGrouped <- calcWeedsGrouped(data2022_2WeedsGroup)

(plotWeedsBox2022_2 <- ggplot(data = data2022_2WeedsGroup, aes(x = reorder(Group, -W, mean), y = W)) +
    geom_boxplot() +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")")) +
    theme_classic(base_size = 30) +
    labs(title = "Weed biomass, 2022",
         subtitle = "Second harvest") +
    annotate("text", x = 1:3, y = 570, label = sapply(cldWeedsGrouped$.group, function(x)gsub(" ", "", x)), size = 4))

fit0 <- lm(W ~ 0 + offset(Harmonic), data = predWeeds2022_2[[1]])
fit1 <- lm(W ~ Harmonic, data = predWeeds2022_2[[1]])
fit2 <- lm(W ~ 0 + offset(Arithmetic), data = predWeeds2022_2[[1]])
fit3 <- lm(W ~ Arithmetic, data = predWeeds2022_2[[1]])
lrtest(fit0, fit1) # Harmonic line is not significantly different from 0:1 line, p = 0.2746
lrtest(fit2, fit3) # Arithmetic line is significantly different from 0:1 line, p = 2.877e-8
lrtest(fit1, fit3) # Harmonic and arithmetic lines are significantly different from each other, p < 2.2e-16


predWeeds2022Long_2 <- predWeeds2022_2[[1]] %>% 
  pivot_longer(cols = c("W", "Arithmetic", "Harmonic"), names_to = "Type", values_to = "WeedBiomass")

newType <- sapply(predWeeds2022Long_2$Type, function(t)ifelse(t == "W", "Observed", t))

predWeeds2022Long_2 <- mutate(predWeeds2022Long_2, Type = newType) %>% 
  mutate(Type = factor(Type, levels = c("Observed", "Arithmetic", "Harmonic")))

(plotWeedPredBox2022_2 <- ggplot(data = predWeeds2022Long_2, aes(x = factor(Type, levels = c("Arithmetic", "Harmonic", "Observed")), y = WeedBiomass)) +
    geom_boxplot() +
    labs(x = "", y = bquote("Weed dry weight (g "~m^-2~")")) +
    theme_classic(base_size = 30) +
    labs(title = "Predicted and observed weed biomass, 2022",
         subtitle = "Second harvest") +
    annotate("text", x = 1:3, y = 485, label = c("a", "b", "b"), size = 4))


# Split per cereal
predWeeds2022_2C <- sapply(predWeeds2022_2[[1]]$Treatment, function(Treatment)strsplit(Treatment, split = "_")[[1]][1])

predWeeds2022_2Cereal <- predWeeds2022_2[[1]] %>% 
  mutate(Cereal = predWeeds2022_2C) %>% 
  rename(Observed = W) %>% 
  pivot_longer(cols = c(Observed, Arithmetic, Harmonic), names_to = "Model", values_to = "WeedBiomass") %>% 
  mutate(Model = factor(Model, levels = c("Arithmetic", "Harmonic", "Observed")))

pAOVBarley <- summary(aov(WeedBiomass ~ Model + (1|Block), data = filter(predWeeds2022_2Cereal, Cereal == "Barley")))[[1]]$`Pr(>F)`[1]
modWeedBarley <- lmer(WeedBiomass ~ Model + (1|Block), data = filter(predWeeds2022_2Cereal, Cereal == "Barley"))
PHWeedsBarley <- emmeans(modWeedBarley, list(pairwise ~ Model), adjust = "tukey")
CLDWeedsBarley <- cld(PHWeedsBarley$emmeans,
                Letters = letters,
                decreasing = TRUE)

(plotWeedsBox2022_2Barley <- ggplot(data = filter(predWeeds2022_2Cereal, Cereal == "Barley"), aes(x = Model, y = WeedBiomass)) +
    geom_boxplot() +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")")) +
    theme_classic(base_size = 8) +
    labs(title = "Predicted and observed weed biomass",
         subtitle = "2022, second harvest, Barley",
         x = "") +
    annotate("text", x = 1:3, y = 340, label = sapply(CLDWeedsBarley$.group, function(x)gsub(" ", "", x)), size = 3))

pAOVRye <- summary(aov(WeedBiomass ~ Model + (1|Block), data = filter(predWeeds2022_2Cereal, Cereal == "Rye")))[[1]]$`Pr(>F)`[1]
modWeedRye <- lmer(WeedBiomass ~ Model + (1|Block), data = filter(predWeeds2022_2Cereal, Cereal == "Rye"))
PHWeedsRye <- emmeans(modWeedRye, list(pairwise ~ Model), adjust = "tukey")
CLDWeedsRye <- cld(PHWeedsRye$emmeans,
                      Letters = letters,
                      decreasing = TRUE)

(plotWeedsBox2022_2Rye <- ggplot(data = filter(predWeeds2022_2Cereal, Cereal == "Rye"), aes(x = Model, y = WeedBiomass)) +
    geom_boxplot() +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")")) +
    theme_classic(base_size = 8) +
    labs(title = "Predicted and observed weed biomass",
         subtitle = "2022, second harvest, Rye",
         x = "") +
    annotate("text", x = 1:3, y = 340, label = sapply(CLDWeedsRye$.group, function(x)gsub(" ", "", x)), size = 3))

pAOVTriticale <- summary(aov(WeedBiomass ~ Model + (1|Block), data = filter(predWeeds2022_2Cereal, Cereal == "Triticale")))[[1]]$`Pr(>F)`[1]
modWeedTriticale <- lmer(WeedBiomass ~ Model + (1|Block), data = filter(predWeeds2022_2Cereal, Cereal == "Triticale"))
PHWeedsTriticale <- emmeans(modWeedTriticale, list(pairwise ~ Model), adjust = "tukey")
CLDWeedsTriticale <- cld(PHWeedsTriticale$emmeans,
                      Letters = letters,
                      decreasing = TRUE)

(plotWeedsBox2022_2Triticale <- ggplot(data = filter(predWeeds2022_2Cereal, Cereal == "Triticale"), aes(x = Model, y = WeedBiomass)) +
    geom_boxplot() +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")")) +
    theme_classic(base_size = 8) +
    labs(title = "Predicted and observed weed biomass",
         subtitle = "2022, second harvest, Triticale",
         x = "") +
    annotate("text", x = 1:3, y = 330, label = sapply(CLDWeedsTriticale$.group, function(x)gsub(" ", "", x)), size = 3))

pAOVWheat <- summary(aov(WeedBiomass ~ Model + (1|Block), data = filter(predWeeds2022_2Cereal, Cereal == "Wheat")))[[1]]$`Pr(>F)`[1]

(plotWeedsBox2022_2Wheat <- ggplot(data = filter(predWeeds2022_2Cereal, Cereal == "Wheat"), aes(x = Model, y = WeedBiomass)) +
    geom_boxplot() +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")")) +
    theme_classic(base_size = 8) +
    labs(title = "Predicted and observed weed biomass",
         subtitle = "2022, second harvest, Wheat",
         x = "") +
    annotate("text", x = 2, y = 450, label = "n.s.", size = 3))

######################
##                  ##
## Experiment 2023A ##
##                  ##
######################

##
## First harvest
##

# Read data
data2023_A1 <- read_xlsx("WSF_2023_data.xlsx", 
    sheet = "BiomassA1", 
    range = "A1:F41", 
    col_names = TRUE) %>% 
  rename(Date = HarvestDate) %>% 
  dplyr::select(Plot, Treatment, Block, Date, WeedBiomass) %>% 
  mutate(Time = as.integer(as.Date(Date) - sowingDate2023)) %>%
  rename(W = WeedBiomass) %>% 
  mutate(W = as.numeric(W)) %>% 
  dplyr::select(-Date) %>% 
  arrange(Time, Plot)

data2023_A1Weeds <- data2023_A1 %>% 
  dplyr::select(Plot, Treatment, Block, W) #%>% 
  # na.omit()

weedCLD2023_A1 <- calcWeeds(na.omit(data2023_A1Weeds))
data2023_A1WeedsGroup <- sapply(data2023_A1Weeds$Treatment, function(t)ifelse(t %in% c("T", "T+"), "Sole triticale", ifelse(t %in% c("F", "F+"), "Sole faba", "Intercrop")))

data2023_A1WeedsGroup <- data2023_A1Weeds %>% 
  mutate(Group = data2023_A1WeedsGroup)

# Create plot
(plotWeeds2023_A1 <- ggplot(data = data2023_A1WeedsGroup) +
    geom_boxplot(aes(x = reorder(Treatment, -W, function(x)mean(x, na.rm = TRUE)), y = W, fill = Group)) +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")"),
         fill = "Crop type") +
    theme_classic(base_size = 30) +
    scale_fill_manual(name = "Crop", values = modelColors) +
    labs(title = "Weed biomass, 2023 A",
         subtitle = "First harvest") +
    annotate("text", x = 4.5, y = 0.9, label = "n.s.", size = 4))

# Weed prediction
predWeeds2023_A1 <- predictWeeds(data2023_A1Weeds)
(plotWeedPred2023_A1 <- predWeeds2023_A1[[2]] +
    labs(subtitle = "2023 A, first harvest"))


# Read data
data2023_A2 <- read_xlsx("WSF_2023_data.xlsx", sheet = "BiomassA2", range = "A1:F41", col_names = TRUE) %>% 
  rename(Date = HarvestDate) %>% 
  dplyr::select(Plot, Treatment, Block, Date, WeedBiomass) %>% 
  mutate(Time = as.integer(as.Date(Date) - sowingDate2023)) %>%
  rename(W = WeedBiomass) %>% 
  mutate(W = as.numeric(W)) %>% 
  dplyr::select(-c(Date)) %>% 
  arrange(Time, Plot)

data2023_A2Weeds <- data2023_A2 %>% 
  dplyr::select(Plot, Treatment, Block, W)

weedCLD2023_A2 <- calcWeeds(na.omit(data2023_A2Weeds))
data2023_A2WeedsGroup <- sapply(data2023_A2Weeds$Treatment, function(t)ifelse(t %in% c("T", "T+"), "Sole triticale", ifelse(t %in% c("F", "F+"), "Sole faba", "Intercrop")))

data2023_A2WeedsGroup <- data2023_A2Weeds %>% 
  mutate(Group = data2023_A2WeedsGroup)

# Create plot
(plotWeeds2023_A2 <- ggplot(data = data2023_A2WeedsGroup) +
    geom_boxplot(aes(x = reorder(Treatment, -W, function(x)mean(x, na.rm = TRUE)), y = W, fill = Group)) +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")"),
         fill = "Crop type") +
    theme_classic(base_size = 30) +
    scale_fill_manual(name = "Crop", values = modelColors) +
    labs(title = "Weed biomass, 2023 A",
         subtitle = "Final harvest") +
    annotate("text", x = 1:8, y = 120, label = weedCLD2023_A2$.group, size = 2.5))

# Weed prediction
predWeeds2023_A2 <- predictWeeds(data2023_A2Weeds)
(plotWeedPred2023_A2 <- predWeeds2023_A2[[2]] +
    labs(subtitle = "2023 A, second harvest"))


predWeeds2023ALong_2 <- predWeeds2023_A2[[1]] %>% 
  pivot_longer(cols = c("W", "Arithmetic", "Harmonic"), names_to = "Type", values_to = "WeedBiomass")

newType <- sapply(predWeeds2023ALong_2$Type, function(t)ifelse(t == "W", "Observed", t))

predWeeds2023ALong_2 <- mutate(predWeeds2023ALong_2, Type = newType) %>% 
  mutate(Type = factor(Type, levels = c("Observed", "Arithmetic", "Harmonic")))

(plotWeedPredBox2023A_2 <- ggplot(data = predWeeds2023ALong_2, aes(x = factor(Type, levels = c("Arithmetic", "Harmonic", "Observed")), y = WeedBiomass)) +
    geom_boxplot() +
    labs(x = "", y = bquote("Weed dry weight (g "~m^-2~")")) +
    theme_classic(base_size = 30) +
    labs(title = "Predicted and observed weed biomass, 2023A",
         subtitle = "Second harvest") +
    annotate("text", x = 1:3, y = 98, label = c("a", "b", "b"), size = 4))

# Separate per treatment
predWeeds2023_A2Treatment <- predWeeds2023_A2[[1]] %>% 
  na.omit() %>% 
  rename(Observed = W) %>% 
  pivot_longer(cols = c(Observed, Arithmetic, Harmonic), names_to = "Model", values_to = "WeedBiomass") %>% 
  mutate(Model = factor(Model, levels = c("Arithmetic", "Harmonic", "Observed")))

pAOV1T1F <- summary(aov(WeedBiomass ~ Model + (1|Block), data = filter(predWeeds2023_A2Treatment, Treatment == "1T:1F")))[[1]]$`Pr(>F)`[1]
modWeed1T1F <- lmer(WeedBiomass ~ Model + (1|Block), data = filter(predWeeds2023_A2Treatment, Treatment == "1T:1F"))
PHWeeds1T1F <- emmeans(modWeedBarley, list(pairwise ~ Model), adjust = "tukey")
CLDWeeds1T1F <- cld(PHWeeds1T1F$emmeans,
                Letters = letters,
                decreasing = TRUE)

(plotWeedsBox2023A_21T1F <- ggplot(data = filter(predWeeds2023_A2Treatment, Treatment == "1T:1F"), aes(x = Model, y = WeedBiomass)) +
    geom_boxplot() +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")")) +
    theme_classic(base_size = 30) +
    labs(title = "Predicted and observed weed biomass",
         subtitle = "2023A, second harvest, 1T:1F",
         x = "") +
    annotate("text", x = 1:3, y = 67, label = sapply(CLDWeeds1T1F$.group, function(x)gsub(" ", "", x)), size = 4))


pAOVTFM <- summary(aov(WeedBiomass ~ Model + (1|Block), data = filter(predWeeds2023_A2Treatment, Treatment == "TF-M")))[[1]]$`Pr(>F)`[1]
modWeedTFM <- lmer(WeedBiomass ~ Model + (1|Block), data = filter(predWeeds2023_A2Treatment, Treatment == "TF-M"))
PHWeedsTFM <- emmeans(modWeedBarley, list(pairwise ~ Model), adjust = "tukey")
CLDWeedsTFM <- cld(PHWeedsTFM$emmeans,
                Letters = letters,
                decreasing = TRUE)

(plotWeedsBox2023A_2TFM <- ggplot(data = filter(predWeeds2023_A2Treatment, Treatment == "TF-M"), aes(x = Model, y = WeedBiomass)) +
    geom_boxplot() +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")")) +
    theme_classic(base_size = 30) +
    labs(title = "Predicted and observed weed biomass",
         subtitle = "2023A, second harvest, TF-M",
         x = "") +
    annotate("text", x = 1:3, y = 67, label = sapply(CLDWeedsTFM$.group, function(x)gsub(" ", "", x)), size = 4))



pAOV1T3F <- summary(aov(WeedBiomass ~ Model + (1|Block), data = filter(predWeeds2023_A2Treatment, Treatment == "1T:3F")))[[1]]$`Pr(>F)`[1]
modWeed1T3F <- lmer(WeedBiomass ~ Model + (1|Block), data = filter(predWeeds2023_A2Treatment, Treatment == "1T:3F"))
PHWeeds1T3F <- emmeans(modWeedBarley, list(pairwise ~ Model), adjust = "tukey")
CLDWeeds1T3F <- cld(PHWeeds1T3F$emmeans,
                Letters = letters,
                decreasing = TRUE)

(plotWeedsBox2023A_21T3F <- ggplot(data = filter(predWeeds2023_A2Treatment, Treatment == "1T:3F"), aes(x = Model, y = WeedBiomass)) +
    geom_boxplot() +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")")) +
    theme_classic(base_size = 30) +
    labs(title = "Predicted and observed weed biomass",
         subtitle = "2023A, second harvest, 1T:3F",
         x = "") +
    annotate("text", x = 1:3, y = 92, label = sapply(CLDWeeds1T3F$.group, function(x)gsub(" ", "", x)), size = 4))



pAOV3T1F <- summary(aov(WeedBiomass ~ Model + (1|Block), data = filter(predWeeds2023_A2Treatment, Treatment == "3T:1F")))[[1]]$`Pr(>F)`[1]
modWeed3T1F <- lmer(WeedBiomass ~ Model + (1|Block), data = filter(predWeeds2023_A2Treatment, Treatment == "3T:1F"))
PHWeeds3T1F <- emmeans(modWeedBarley, list(pairwise ~ Model), adjust = "tukey")
CLDWeeds3T1F <- cld(PHWeeds3T1F$emmeans,
                Letters = letters,
                decreasing = TRUE)

(plotWeedsBox2023A_23T1F <- ggplot(data = filter(predWeeds2023_A2Treatment, Treatment == "3T:1F"), aes(x = Model, y = WeedBiomass)) +
    geom_boxplot() +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")")) +
    theme_classic(base_size = 30) +
    labs(title = "Predicted and observed weed biomass",
         subtitle = "2023A, second harvest, 3T:1F",
         x = "") +
    annotate("text", x = 1:3, y = 42, label = sapply(CLDWeeds3T1F$.group, function(x)gsub(" ", "", x)), size = 4))

predWeeds2023_A2[[1]]
fit0 <- lm(W ~ 0 + offset(Harmonic), data = predWeeds2023_A2[[1]])
fit1 <- lm(W ~ Harmonic, data = predWeeds2023_A2[[1]])
fit2 <- lm(W ~ 0 + offset(Arithmetic), data = predWeeds2023_A2[[1]])
fit3 <- lm(W ~ Arithmetic, data = predWeeds2023_A2[[1]])
lrtest(fit0, fit1) # Harmonic is significantly different from 0:1 line, p = 0.0003345
lrtest(fit2, fit3) # Arithmetic is significantly different from 0:1 line, p = 7.994e-13
lrtest(fit1, fit3) # Harmonic is significantly different from Arithmetic line, p < 2.2e-16


##
## Third harvest
##

# Read data
data2023_A3 <- read_xlsx("WSF_2023_data.xlsx", sheet = "FinalHarvestA", range = "A1:F41", col_names = TRUE) %>%
  dplyr::select(-TreatmentN) %>% 
 arrange(Plot) %>% 
  dplyr::select(Plot, Treatment, Block, WeedBiomass) %>% 
  rename(W = WeedBiomass)

data2023_A3Weeds <- data2023_A3 %>% 
  dplyr::select(Plot, Treatment, Block, W)

weedCLD2023_A3 <- calcWeeds(na.omit(data2023_A3Weeds))
data2023_A3WeedsGroup <- sapply(data2023_A3Weeds$Treatment, function(t)ifelse(t %in% c("T", "T+"), "Sole triticale", ifelse(t %in% c("F", "F+"), "Sole faba", "Intercrop")))

data2023_A3WeedsGroup <- data2023_A3Weeds %>% 
  mutate(Group = data2023_A3WeedsGroup)

# Create plot
(plotWeeds2023_A3 <- ggplot(data = data2023_A3WeedsGroup) +
    geom_boxplot(aes(x = reorder(Treatment, -W, function(x)mean(x, na.rm = TRUE)), y = W, fill = Group)) +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")"),
         fill = "Crop type") +
    theme_classic(base_size = 30) +
    scale_fill_manual(name = "Crop", values = modelColors) +
    labs(title = "Weed biomass, 2023 A",
         subtitle = "Final harvest") +
    annotate("text", x = 1:8, y = 340, label = weedCLD2023_A3$.group, size = 2.5))

# Weed prediction
predWeeds2023_A3 <- predictWeeds(data2023_A3Weeds)
(plotWeedPred2023_A3 <- predWeeds2023_A3[[2]] +
    labs(subtitle = "2023 A, final harvest"))

fit0 <- lm(W ~ 0 + offset(Harmonic), data = predWeeds2023_A3[[1]])
fit1 <- lm(W ~ Harmonic, data = predWeeds2023_A3[[1]])
fit2 <- lm(W ~ 0 + offset(Arithmetic), data = predWeeds2023_A3[[1]])
fit3 <- lm(W ~ Arithmetic, data = predWeeds2023_A3[[1]])
lrtest(fit0, fit1) # Harmonic is not significantly different from 0:1 line, p = 4041
lrtest(fit2, fit3) # Arithmetic is not significantly different from 0:1 line, p = 0.06758
lrtest(fit1, fit3) # Harmonic is significantly different from Arithmetic line, p < 2.2e-16

######################
##                  ##
## Experiment 2023B ##
##                  ##
######################

##
## First harvest
##

# Read data
data2023_B1 <- read_xlsx("WSF_2023_data.xlsx", sheet = "BiomassB1", range = "A1:F45", col_names = TRUE) %>% 
  dplyr::select(Plot, Treatment, Block, WeedBiomass) %>% 
  rename(W = WeedBiomass) %>%  
  arrange(Plot) %>% 
  filter(!Treatment %in% c("T-25", "1T:1F-25", "TF-M-25", "F-25"))


data2023_B1Weeds <- data2023_B1 %>% 
  dplyr::select(Plot, Treatment, Block, W)

weedCLD2023_B1 <- calcWeeds(na.omit(data2023_B1Weeds))
data2023_B1WeedsGroup <- sapply(data2023_B1Weeds$Treatment, function(t)ifelse(t %in% c("T", "T-25", "T-375"), "Sole triticale", ifelse(t %in% c("F", "F-25", "F-375"), "Sole faba", "Intercrop")))

data2023_B1WeedsGroup <- data2023_B1Weeds %>% 
  mutate(Group = data2023_B1WeedsGroup)

# Create plot
(plotWeeds2023_B1 <- ggplot(data = data2023_B1WeedsGroup) +
    geom_boxplot(aes(x = reorder(Treatment, -W, function(x)mean(x, na.rm = TRUE)), y = W, fill = Group)) +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")"),
         fill = "Crop type") +
    theme_classic(base_size = 30) +
    scale_fill_manual(name = "Crop", values = modelColors) +
    labs(title = "Weed biomass, 2023 B",
         subtitle = "First harvest") +
    annotate("text", x = 5.5, y = 1.0, label = "n.s.", size = 4))

# Weed prediction
predWeeds2023_B1 <- predictWeeds(data2023_B1Weeds)
(plotWeedPred2023_B1 <- predWeeds2023_B1[[2]] +
    labs(subtitle = "2023 B, first harvest"))

##
## Second harvest
##

# Read data
data2023_B2 <- read_xlsx("WSF_2023_data.xlsx", sheet = "BiomassB2", range = "A1:F45", col_names = TRUE) %>% 
  dplyr::select(Plot, Treatment, Block, WeedBiomass) %>% 
  rename(W = WeedBiomass) %>%  
  arrange(Plot) %>% 
  filter(!Treatment %in% c("T-25", "1T:1F-25", "TF-M-25", "F-25"))

data2023_B2Weeds <- data2023_B2 %>% 
  dplyr::select(Plot, Treatment, Block, W) %>% 
  na.omit()

weedCLD2023_B2 <- calcWeeds(data2023_B2Weeds)
data2023_B2WeedsGroup <- sapply(data2023_B2Weeds$Treatment, function(t)ifelse(t %in% c("T", "T-25", "T-375"), "Sole triticale", ifelse(t %in% c("F", "F-25", "F-375"), "Sole faba", "Intercrop")))

data2023_B2WeedsGroup <- data2023_B2Weeds %>% 
  mutate(Group = data2023_B2WeedsGroup)

# Create plot
(plotWeeds2023_B2 <- ggplot(data = data2023_B2WeedsGroup) +
    geom_boxplot(aes(x = reorder(Treatment, -W, function(x)mean(x, na.rm = TRUE)), y = W, fill = Group)) +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")"),
         fill = "Crop type") +
    theme_classic(base_size = 30) +
    scale_fill_manual(name = "Crop", values = modelColors) +
    labs(title = "Weed biomass, 2023 A",
         subtitle = "Final harvest") +
    annotate("text", x = 1:8, y = 120, label = weedCLD2023_B2$.group, size = 2.5))

# Weed prediction
predWeeds2023_B2 <- predictWeeds(data2023_B2Weeds)
(plotWeedPred2023_B2 <- predWeeds2023_B2[[2]] +
    labs(subtitle = "2023 B, second harvest"))

predWeeds2023BLong_2 <- predWeeds2023_B2[[1]] %>% 
  pivot_longer(cols = c("W", "Arithmetic", "Harmonic"), names_to = "Type", values_to = "WeedBiomass")

newType <- sapply(predWeeds2023BLong_2$Type, function(t)ifelse(t == "W", "Observed", t))

predWeeds2023BLong_2 <- mutate(predWeeds2023BLong_2, Type = newType) %>% 
  mutate(Type = factor(Type, levels = c("Observed", "Arithmetic", "Harmonic")))

(plotWeedPredBox2023B_2 <- ggplot(data = predWeeds2023BLong_2, aes(x = factor(Type, levels = c("Arithmetic", "Harmonic", "Observed")), y = WeedBiomass)) +
    geom_boxplot() +
    labs(x = "", y = bquote("Weed dry weight (g "~m^-2~")")) +
    theme_classic(base_size = 30) +
    labs(title = "Predicted and observed weed biomass, 2023B",
         subtitle = "Second harvest") +
    annotate("text", x = 1:3, y = 78, label = c("a", "b", "b"), size = 4))

fit0 <- lm(W ~ 0 + offset(Harmonic), data = predWeeds2023_B2[[1]])
fit1 <- lm(W ~ Harmonic, data = predWeeds2023_B2[[1]])
fit2 <- lm(W ~ 0 + offset(Arithmetic), data = predWeeds2023_B2[[1]])
fit3 <- lm(W ~ Arithmetic, data = predWeeds2023_B2[[1]])
lrtest(fit0, fit1) # Harmonic line is not significantly different from the 0:1 line, p = 0.771
lrtest(fit2, fit3) # Arithmetic line is significantly different from the 0:1 line, p = 0.0001649
lrtest(fit1, fit3)# Harmonic and arithmetic lines are significantly different from each other, p < 2.2e-16

# Split per row distance
predWeeds2023_B2ICT <- sapply(predWeeds2023_B2[[1]]$Treatment, function(Treatment)strsplit(Treatment, split = "_")[[1]][1])

predWeeds2023_B2IntType <- predWeeds2023_B2[[1]] %>% 
  mutate(IntType = predWeeds2023_B2ICT) %>% 
  rename(Observed = W) %>% 
  pivot_longer(cols = c(Observed, Arithmetic, Harmonic), names_to = "Model", values_to = "WeedBiomass") %>% 
  mutate(Model = factor(Model, levels = c("Arithmetic", "Harmonic", "Observed")))

pAOV1T1F <- summary(aov(WeedBiomass ~ Model + (1|Block), data = filter(predWeeds2023_B2IntType, IntType == "1T:1F")))[[1]]$`Pr(>F)`[1]
modWeed1T1F <- lmer(WeedBiomass ~ Model + (1|Block), data = filter(predWeeds2023_B2IntType, IntType == "1T:1F"))
PHWeeds1T1F <- emmeans(modWeed1T1F, list(pairwise ~ Model), adjust = "tukey")
CLDWeeds1T1F <- cld(PHWeeds1T1F$emmeans,
                      Letters = letters,
                      decreasing = TRUE)

(plotWeedsBox2023_B21T1F <- ggplot(data = filter(predWeeds2023_B2IntType, IntType == "1T:1F"), aes(x = Model, y = WeedBiomass)) +
    geom_boxplot() +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")")) +
    theme_classic(base_size = 8) +
    labs(title = "Predicted and observed weed biomass",
         subtitle = "2023 B2, second harvest, 1T:1F",
         x = "") +
    annotate("text", x = 1:3, y = 40, label = sapply(CLDWeeds1T1F$.group, function(x)gsub(" ", "", x)), size = 3))

pAOVTFM <- summary(aov(WeedBiomass ~ Model + (1|Block), data = filter(predWeeds2023_B2IntType, IntType == "TF-M")))[[1]]$`Pr(>F)`[1]
modWeedTFM <- lmer(WeedBiomass ~ Model + (1|Block), data = filter(predWeeds2023_B2IntType, IntType == "TF-M"))
PHWeedsTFM <- emmeans(modWeedTFM, list(pairwise ~ Model), adjust = "tukey")
CLDWeedsTFM <- cld(PHWeedsTFM$emmeans,
                     Letters = letters,
                     decreasing = TRUE)

(plotWeedsBox2023_B2TFM <- ggplot(data = filter(predWeeds2023_B2IntType, IntType == "TF-M"), aes(x = Model, y = WeedBiomass)) +
    geom_boxplot() +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")")) +
    theme_classic(base_size = 8) +
    labs(title = "Predicted and observed weed biomass",
         subtitle = "2023 B2, second harvest, TF-M",
         x = "") +
    annotate("text", x = 1:3, y = 40, label = sapply(CLDWeedsTFM$.group, function(x)gsub(" ", "", x)), size = 3))

pAOV1T1F375 <- summary(aov(WeedBiomass ~ Model + (1|Block), data = filter(predWeeds2023_B2IntType, IntType == "1T:1F-375")))[[1]]$`Pr(>F)`[1]
# p = 0.489

(plotWeedsBox2023_B21T1F375 <- ggplot(data = filter(predWeeds2023_B2IntType, IntType == "1T:1F-375"), aes(x = Model, y = WeedBiomass)) +
    geom_boxplot() +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")")) +
    theme_classic(base_size = 8) +
    labs(title = "Predicted and observed weed biomass",
         subtitle = "2023 B2, second harvest, 1T:1F-375",
         x = "") +
    annotate("text", x = 2, y = 65, label = "n.s.", size = 3))

pAOVTFM375 <- summary(aov(WeedBiomass ~ Model + (1|Block), data = filter(predWeeds2023_B2IntType, IntType == "TF-M-375")))[[1]]$`Pr(>F)`[1]
# p = 0.3278

(plotWeedsBox2023_B2TFM375 <- ggplot(data = filter(predWeeds2023_B2IntType, IntType == "TF-M-375"), aes(x = Model, y = WeedBiomass)) +
    geom_boxplot() +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")")) +
    theme_classic(base_size = 8) +
    labs(title = "Predicted and observed weed biomass",
         subtitle = "2023 B2, second harvest, TF-M-375",
         x = "") +
    annotate("text", x = 2, y = 65, label = "n.s.", size = 3))

##
## Third harvest
##

# Read data
data2023_B3 <- read_xlsx("WSF_2023_data.xlsx", sheet = "FinalHarvestB", range = "A1:F45", col_names = TRUE) %>%
  dplyr::select(-TreatmentN) %>% 
  arrange(Plot) %>% 
  dplyr::select(Plot, Treatment, Block, WeedBiomass) %>% 
  rename(W = WeedBiomass) %>% 
  filter(!Treatment %in% c("T-25", "1T:1F-25", "TF-M-25", "F-25"))

data2023_B3Weeds <- data2023_B3 %>% 
  dplyr::select(Plot, Treatment, Block, W) %>% 
  na.omit()

weedCLD2023_B3 <- calcWeeds(data2023_B3Weeds)
data2023_B3WeedsGroup <- sapply(data2023_B3Weeds$Treatment, function(t)ifelse(t %in% c("T", "T-375"), "Sole triticale", ifelse(t %in% c("F", "F-375"), "Sole faba", "Intercrop")))

data2023_B3WeedsGroup <- data2023_B3Weeds %>% 
  mutate(Group = data2023_B3WeedsGroup)

# Create plot
(plotWeeds2023_B3 <- ggplot(data = data2023_B3WeedsGroup) +
    geom_boxplot(aes(x = reorder(Treatment, -W, function(x)mean(x, na.rm = TRUE)), y = W, fill = Group)) +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")"),
         fill = "Crop type") +
    theme_classic(base_size = 30) +
    scale_fill_manual(name = "Crop", values = modelColors) +
    labs(title = "Weed biomass, 2023 B",
         subtitle = "Final harvest") +
    annotate("text", x = 1:8, y = 340, label = weedCLD2023_B3$.group, size = 2.5))

# Weed prediction
predWeeds2023_B3 <- predictWeeds(data2023_B3Weeds)
(plotWeedPred2023_B3 <- predWeeds2023_B3[[2]] +
    labs(subtitle = "2023 B, final harvest"))


######################
##                  ##
## Experiment 2024 ##
##                  ##
######################

##
## First harvest
##

# Read data
data2024_1 <- read_xlsx("WSF_2024_data.xlsx", sheet = "Biomass1", range = "A1:E21", col_names = TRUE) %>% 
  dplyr::select(Plot, Treatment, Block, WeedBiomass) %>%
  rename(W = WeedBiomass) %>% 
  arrange(Plot)

data2024_1Weeds <- data2024_1 %>% 
  dplyr::select(Plot, Treatment, Block, W)

weedCLD2024_1 <- calcWeeds(na.omit(data2024_1Weeds))
data2024_1WeedsGroup <- sapply(data2024_1Weeds$Treatment, function(t)ifelse(t %in% c("T"), "Sole triticale", ifelse(t %in% c("F"), "Sole faba", ifelse(t %in% c("W"), "Sole weeds", "Intercrop"))))

data2024_1WeedsGroup <- data2024_1Weeds %>% 
  mutate(Group = data2024_1WeedsGroup)

# Create plot
(plotWeeds2024_1 <- ggplot(data = data2024_1WeedsGroup) +
    geom_boxplot(aes(x = reorder(Treatment, -W, function(x)mean(x, na.rm = TRUE)), y = W, fill = Group)) +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")"),
         fill = "Crop type") +
    theme_classic(base_size = 30) +
    scale_fill_manual(name = "Crop", values = modelColors) +
    labs(title = "Weed biomass, 2024",
         subtitle = "First harvest") +
    annotate("text", x = 2.5, y = 5.5, label = "n.s.", size = 4))

# Weed prediction
data2024_1WeedsPred <- data2024_1Weeds %>% 
  filter(Treatment != "W")
predWeeds2024_1 <- predictWeeds(data2024_1WeedsPred)
(plotWeedPred2024_1 <- predWeeds2024_1[[2]] +
    labs(subtitle = "2024, first harvest"))

##
## Second harvest
##

# Read data
data2024_2 <- read_xlsx("WSF_2024_data.xlsx", sheet = "Biomass2", range = "A1:E21", col_names = TRUE) %>% 
  dplyr::select(Plot, Treatment, Block, WeedBiomass) %>%
  rename(W = WeedBiomass) %>% 
  arrange(Plot)

data2024_2Weeds <- data2024_2 %>% 
  dplyr::select(Plot, Treatment, Block, W) %>% 
  na.omit()

weedCLD2024_2 <- calcWeeds(data2024_2Weeds)
data2024_2WeedsGroup <- sapply(data2024_2Weeds$Treatment, function(t)ifelse(t %in% c("T"), "Sole triticale", ifelse(t %in% c("F"), "Sole faba", ifelse(t %in% c("W"), "Sole weeds", "Intercrop"))))

data2024_2WeedsGroup <- data2024_2Weeds %>% 
  mutate(Group = data2024_2WeedsGroup)

# Create plot
(plotWeeds2024_2 <- ggplot(data = data2024_2WeedsGroup) +
    geom_boxplot(aes(x = reorder(Treatment, -W, function(x)mean(x, na.rm = TRUE)), y = W, fill = Group)) +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")"),
         fill = "Crop type") +
    theme_classic(base_size = 30) +
    scale_fill_manual(name = "Crop", values = modelColors) +
    labs(title = "Weed biomass, 2024",
         subtitle = "Final harvest") +
    annotate("text", x = 1:4, y = 450, label = weedCLD2024_2$.group, size = 2.5))

# Weed prediction
data2024_2WeedsPred <- data2024_2Weeds %>% 
  filter(Treatment != "W")
predWeeds2024_2 <- predictWeeds(data2024_2WeedsPred)
(plotWeedPred2024_2 <- predWeeds2024_2[[2]] +
    labs(subtitle = "2024, second harvest"))

predWeeds2024Long_2 <- predWeeds2024_2[[1]] %>% 
  pivot_longer(cols = c("W", "Arithmetic", "Harmonic"), names_to = "Type", values_to = "WeedBiomass")

newType <- sapply(predWeeds2024Long_2$Type, function(t)ifelse(t == "W", "Observed", t))

predWeeds2024Long_2 <- mutate(predWeeds2024Long_2, Type = newType) %>% 
  mutate(Type = factor(Type, levels = c("Observed", "Arithmetic", "Harmonic")))

(plotWeedPredBox2024_2 <- ggplot(data = predWeeds2024Long_2, aes(x = factor(Type, levels = c("Arithmetic", "Harmonic", "Observed")), y = WeedBiomass)) +
    geom_boxplot() +
    labs(x = "", y = bquote("Weed dry weight (g "~m^-2~")")) +
    theme_classic(base_size = 30) +
    labs(title = "Predicted and observed weed biomass, 2024",
         subtitle = "Second harvest") +
    annotate("text", x = 1:3, y = 125, label = c("a", "ab", "b"), size = 4))


fit0 <- lm(W ~ 0 + offset(Harmonic), data = predWeeds2024_2[[1]])
fit1 <- lm(W ~ Harmonic, data = predWeeds2024_2[[1]])
fit2 <- lm(W ~ 0 + offset(Arithmetic), data = predWeeds2024_2[[1]])
fit3 <- lm(W ~ Arithmetic, data = predWeeds2024_2[[1]])
lrtest(fit0, fit1) # Harmonic line is significantly different from the 0:1 line, p = 0.008456
lrtest(fit2, fit3) # Arithmetic line is significantly different from the 0:1 line, p = 0.00493
lrtest(fit1, fit3)# Harmonic and arithmetic lines are significantly different from each other, p < 2.2e-16

##
## Third harvest
##

# Read data
data2024_3 <- read_xlsx("WSF_2024_data.xlsx", sheet = "FinalHarvest", range = "A1:E21", col_names = TRUE) %>%
  arrange(Plot) %>% 
  dplyr::select(Plot, Treatment, Block, WeedBiomass) %>% 
  rename(W = WeedBiomass)

data2024_3Weeds <- data2024_3 %>% 
  dplyr::select(Plot, Treatment, Block, W) %>% 
  na.omit()

weedCLD2024_3 <- calcWeeds(data2024_3Weeds)
data2024_3WeedsGroup <- sapply(data2024_3Weeds$Treatment, function(t)ifelse(t %in% c("T"), "Sole triticale", ifelse(t %in% c("F"), "Sole faba", ifelse(t %in% c("W"), "Sole weeds", "Intercrop"))))

data2024_3WeedsGroup <- data2024_3Weeds %>% 
  mutate(Group = data2024_3WeedsGroup)

# Create plot
(plotWeeds2024_3 <- ggplot(data = data2024_3WeedsGroup) +
    geom_boxplot(aes(x = reorder(Treatment, -W, function(x)mean(x, na.rm = TRUE)), y = W, fill = Group)) +
    labs(x = "Treatment", y = bquote("Weed dry weight (g "~m^-2~")"),
         fill = "Crop type") +
    theme_classic(base_size = 30) +
    scale_fill_manual(name = "Crop", values = modelColors) +
    labs(title = "Weed biomass, 2024",
         subtitle = "Final harvest") +
    annotate("text", x = 1:4, y = 820, label = sapply(weedCLD2024_3$.group, function(x)gsub(" ", "", x)), size = 2.5))

# Weed prediction
data2024_3WeedsPred <- data2024_3Weeds %>% 
  filter(Treatment != "W")
predWeeds2024_3 <- predictWeeds(data2024_3WeedsPred)
(plotWeedPred2024_3 <- predWeeds2024_3[[2]] +
    labs(subtitle = "2024, second harvest"))

fit0 <- lm(W ~ 0 + offset(Harmonic), data = predWeeds2024_3[[1]])
fit1 <- lm(W ~ Harmonic, data = predWeeds2024_3[[1]])
fit2 <- lm(W ~ 0 + offset(Arithmetic), data = predWeeds2024_3[[1]])
fit3 <- lm(W ~ Arithmetic, data = predWeeds2024_3[[1]])
lrtest(fit0, fit1) # Harmonic line is significantly different from the 0:1 line, p = 0.008456
lrtest(fit2, fit3) # Arithmetic line is significantly different from the 0:1 line, p = 0.00493
lrtest(fit1, fit3)# Harmonic and arithmetic lines are significantly different from each other, p < 2.2e-16


### 
### Create figures
###

dir.create("figures")
dir.create("figures/weed")

tiff("figures/weed/plotWeedPred2Combined.tiff", units = "mm", width = 174, height = 174, res = 400)
plotWeedPred2022_2 + theme_classic(base_size = 15) + 
  theme(legend.position = "none", plot.subtitle = element_text(hjust = -0.2)) + 
  labs(y = "", x = bquote(Pred[w]), title = "", subtitle = bquote(Obs[w]~~~~~~~~~~~~"2022")) +
  geom_segment(aes(x = 300, y = 400, xend = 350, yend = 350),
               arrow = arrow(length = unit(0.3, "cm")),
               size = 1,
               col = "black") +
  geom_segment(aes(x = 190, y = 300, xend = 275, yend = 250),
               arrow = arrow(length = unit(0.3, "cm")),
               size = 1,
               col = "black") +
  geom_segment(aes(x = 325, y = 100, xend = 285, yend = 170),
               arrow = arrow(length = unit(0.3, "cm")),
               size = 1,plotWeeds2023_B2
               col = "black") +
  annotate("text", x = c(275, 65, 300), y = c(430, 300, 90), label = c("1:1", "Harmonic", "Arithmetic"), size = 4, hjust = 0.0, col = c("black", "#619CFF", "#F8766D")) +
  
plotWeedPred2023_A2 + labs(title = "", subtitle = "2023A, second harvest") + theme_classic(base_size = 15) + 
  theme(legend.position = "none", plot.subtitle = element_text(hjust = -0.2)) + 
  labs(y = "", x = bquote(Pred[w]), title = "", subtitle = bquote(Obs[w]~~~~~~~~~~~~"2023A")) +
plotWeedPred2023_B2 + labs(title = "", subtitle = "2023B, second harvest") + theme_classic(base_size = 15) + 
  theme(legend.position = "none", plot.subtitle = element_text(hjust = -0.2)) + 
  labs(y = "", x = bquote(Pred[w]), title = "", subtitle = bquote(Obs[w]~~~~~~~~~~~~"2023B")) +
plot_layout(design = "
            12
            3#
            ") +
plot_annotation(tag_levels = "a")
dev.off()

a <- plotWeeds2022_1 + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~~~"2022"), x = "") + theme_classic(base_size = 10) + theme(legend.position = "none", plot.subtitle = element_text(hjust = -0.02), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
b <- plotWeeds2023_A1 + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~~~"2023A"), x = "") + theme_classic(base_size = 10) + theme(legend.position = "none", plot.subtitle = element_text(hjust = -0.2), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
c <- plotWeeds2023_B1 + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~~~"2023B"), y = "") + theme_classic(base_size = 10) + theme(legend.position = "none", plot.subtitle = element_text(hjust = -0.2), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
d <- plotWeeds2024_1 + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~~~"2024"), y = "", x = "") + theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.subtitle = element_text(hjust = -0.2))

tiff("figures/weed/plotWeed1Combined.tiff", units = "mm", width = 174, height = 152, res = 400)
a / (b | c | d) +
plot_annotation(tag_levels = "a")
dev.off()

a <- plotWeeds2022_2 + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~~~"2022"), x = "") + theme_classic(base_size = 10) + theme(legend.position = "none", plot.subtitle = element_text(hjust = -0.02), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
b <- plotWeeds2023_A2 + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~~~"2023A"), x = "") + theme_classic(base_size = 10) + theme(legend.position = "none", plot.subtitle = element_text(hjust = -0.2), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
c <- plotWeeds2023_B2 + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~~~"2023B"), y = "") + theme_classic(base_size = 10) + theme(legend.position = "none", plot.subtitle = element_text(hjust = -0.2), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
d <- plotWeeds2024_2 + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~~~"2024"), x = "", y = "") + theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.subtitle = element_text(hjust = -0.2))

tiff("figures/weed/plotWeed2Combined.tiff", units = "mm", width = 174, height = 152, res = 400)
a / (b | c | d) +
plot_annotation(tag_levels = "a")
dev.off()

b <- plotWeeds2023_A3 + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~~~"2023A"), x = "") + theme_classic(base_size = 10) + theme(legend.position = "none", plot.subtitle = element_text(hjust = -0.02), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
c <- plotWeeds2023_B3 + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~~~"2023B"), y = "") + theme_classic(base_size = 10) + theme(legend.position = "none", plot.subtitle = element_text(hjust = -0.2), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
d <- plotWeeds2024_3 + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~~~"2024"), x = "", y = "") + theme_classic(base_size = 10) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.subtitle = element_text(hjust = -0.2))

tiff("figures/weed/plotWeed3Combined.tiff", units = "mm", width = 174, height = 76, res = 400)
(b | c | d) +
plot_annotation(tag_levels = "a",
                theme = theme(plot.title = element_text(size = 35)))
dev.off()

tiff("figures/weed/plotWeedBoxCombined.tiff", units = "mm", width = 129, height = 129, res = 400)
plotWeedPredBox2022_2 + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~~~"2022")) + 
  theme_classic(base_size = 10) + 
  theme(plot.subtitle = element_text(hjust = -0.1), axis.title.x = element_text(face = "bold")) +
plotWeedPredBox2023A_2 + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~~~"2023A")) +  
  theme_classic(base_size = 10) + 
  theme(plot.subtitle = element_text(hjust = -0.1), axis.title.x = element_text(face = "bold")) +
plotWeedPredBox2023B_2 + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~~~"2023B")) + 
  theme_classic(base_size = 10) + 
  theme(plot.subtitle = element_text(hjust = -0.1), axis.title.x = element_text(face = "bold"))  +
plotWeedPredBox2024_2 + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~~~"2024")) +  
  theme_classic(base_size = 10) + 
  theme(plot.subtitle = element_text(hjust = -0.1), axis.title.x = element_text(face = "bold")) +
plot_layout(nrow = 2, ncol = 2) +
plot_annotation(tag_levels = "a")
dev.off()

tiff("figures/weed/plotWeedBox20222023BSep.tiff", units = "mm", width = 84, height = 174, res = 400)
plotWeedsBox2022_2Barley + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~"Barley")) + 
  theme_classic(base_size = 7) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.subtitle = element_text(hjust = -0.4)) + 
  theme_classic(base_size = 7) + 
plotWeedsBox2022_2Rye + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~"Rye")) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.subtitle = element_text(hjust = -0.4)) + 
  theme_classic(base_size = 7) + 
plotWeedsBox2022_2Triticale + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~"Triticale")) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.subtitle = element_text(hjust = -0.4)) + 
  theme_classic(base_size = 7) + 
plotWeedsBox2022_2Wheat + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~"Wheat")) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.subtitle = element_text(hjust = -0.4)) + 
  theme_classic(base_size = 7) + 
plotWeedsBox2023_B21T1F + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~"1T:1F")) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.subtitle = element_text(hjust = -0.4)) + 
  theme_classic(base_size = 7) + 
plotWeedsBox2023_B21T1F375 + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~"1T:1F-375")) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.subtitle = element_text(hjust = -0.4)) + 
  theme_classic(base_size = 7) + 
plotWeedsBox2023_B2TFM + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~~~"TF-M")) + 
  theme_classic(base_size = 7) + 
  theme(plot.subtitle = element_text(hjust = -0.2), axis.title.x = element_text(face = "bold")) + 
plotWeedsBox2023_B2TFM375 + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~~~"TF-M-375")) +  
  theme_classic(base_size = 7) + 
  theme(plot.subtitle = element_text(hjust = -0.2), axis.title.x = element_text(face = "bold")) + 
plot_layout(nrow = 4, ncol = 2, axis_titles = "collect") +
plot_annotation(tag_levels = "a")
dev.off()

tiff("figures/weed/plotWeedBox2023ASep.tiff", units = "mm", width = 129, height = 129, res = 400)
plotWeedsBox2023A_21T1F + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~~~"1T:1F")) + theme(plot.subtitle = element_text(hjust = -0.1), axis.title.x = element_text(face = "bold")) + 
  theme_classic(base_size = 10) + 
plotWeedsBox2023A_2TFM + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~~~"TF-M")) + theme(plot.subtitle = element_text(hjust = -0.1), axis.title.x = element_text(face = "bold")) + 
  theme_classic(base_size = 10) + 
plotWeedsBox2023A_21T3F + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~~~"1T:3F")) + theme(plot.subtitle = element_text(hjust = -0.1), axis.title.x = element_text(face = "bold")) + 
  theme_classic(base_size = 10) + 
plotWeedsBox2023A_23T1F + labs(y = "", title = "", subtitle = bquote(Y[w]~~~~~~~~~~~~"3T:1F")) + theme(plot.subtitle = element_text(hjust = -0.1), axis.title.x = element_text(face = "bold")) + 
  theme_classic(base_size = 10) + 
plot_layout(nrow = 2, ncol = 2) +
plot_annotation(tag_levels = "a")
dev.off()
