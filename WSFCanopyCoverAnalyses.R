## Author: David Kottelenberg
## Institute: Wageningen University & Research
## Last edited: 13/11/2024
## Summary: Analyses of 2022-2024 field experiments, canopy cover

rm(list = ls())

library(tidyverse)
library(scales)
library(readxl)
library(ggpubr)
library(lubridate)
library(emdbook)
library(emmeans)
library(RColorBrewer)
library(xtable)
library(patchwork)

# Set directory
#setwd() # Set working directory to data location
dir.create("canopy_cover")

colourVector <- c(brewer.pal(n=8,"Set1"), "#FDC086")

modelColors <- c("Cereal" = colourVector[1],
                 "Legume" = colourVector[2],
                 "Faba" = colourVector[2],
                 "Triticale" = colourVector[1],
                 "Wheat" = colourVector[1],
                 "Wheat_Faba" = colourVector[3],
                 "Weed" = colourVector[4],
                 "T" = colourVector[1],
                 "F" = colourVector[2],
                 "T-375" = colourVector[1],
                 "F-375" = colourVector[2],
                 "TA" = colourVector[5],
                 "TM" = colourVector[9],
                 "FA" = colourVector[7],
                 "FM" =  colourVector[8],
                 "Intercrop" = colourVector[3],
                 "Sole cereal" = colourVector[1],
                 "Sole legume" = colourVector[2],
                 "1T:3F" = colourVector[3],
                 "3T:1F" = colourVector[3])

##################################################################
##################################################################
######################                      ######################
######################  Canopy cover  ######################
######################                      ######################
##################################################################
##################################################################

######################
##                  ##
## Experiment 2022  ##
##                  ##
######################

# Set sowing date
sowingDate2022 <- as.Date("2022-04-19")

# Load graph generation and model fitting code
source("WSFModelFitting.R")
source("WSFExploratoryGraphs.R")

# Load temperature data
tempData2022 <- read_delim("TempData2022.txt") %>% 
  mutate(Time = as.integer(as.Date(Date) - sowingDate2022))

cumulativeTemp2022 <- ((tempData2022$TMin + tempData2022$TMax) / 2)
cumulativeTemp2022[which(cumulativeTemp2022 < 0.0)] <- 0
cumulativeTemp2022 <- cumsum(cumulativeTemp2022)

tempData2022 <- tempData2022 %>% 
  mutate(CumulTemp = cumulativeTemp2022)

# Read data
data2022 <- read_xlsx("FE1_data.xlsx", sheet = "CanopyCover", range = "A1:F533", col_names = TRUE) %>% 
  rename("CanopyCover" = "PCover") %>% 
  mutate(DAS = as.integer(as.Date(Date) - sowingDate2022)) %>% 
  mutate(Time = tempData2022$CumulTemp[match(DAS, tempData2022$Time)]) %>% 
  dplyr::select(Plot, Treatment, TreatmentN, Block, Time, CanopyCover) %>% 
  mutate(CanopyCover = as.numeric(CanopyCover)) %>% 
  filter(!grepl("Lupine", Treatment))

plotFit <- function(data, par, fun){
  plot <- ggplot() +
    geom_point(data = data, aes(x = Time, y = CanopyCover)) +
    geom_function(fun = function(x)allDeterministicFunctions[[fun]](par, x),
                  size = 1)
  return(plot)
}

## Manual fits

dataTreatment2022 <- split(data2022, data2022$Treatment)
dataTreatment2022 <- lapply(dataTreatment2022, function(dat){
  datNew <- dat %>% 
    dplyr::select(Time, CanopyCover)
  return(datNew)
})




##################################
###                            ###
### Fit logistic normal curves ###
###                            ###
##################################     

fitsLogisticNormal2022 <- rep(list(NA), 14)

# Fit of Barley
fitsLogisticNormal2022[[1]] <- fitDataToModel(data = dataTreatment2022[[1]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[1]], fitsLogisticNormal2022[[1]]$par, "logistic")

# Fit of Barley_Faba
fitsLogisticNormal2022[[2]] <- fitDataToModel(data = dataTreatment2022[[2]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[2]], fitsLogisticNormal2022[[2]]$par, "logistic")

# Fit of Barley_Pea
fitsLogisticNormal2022[[3]] <- fitDataToModel(data = dataTreatment2022[[3]], model = logisticNormal, startParameters = c("r" = 0.005, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[3]], fitsLogisticNormal2022[[3]]$par, "logistic")

# Fit of Faba
fitsLogisticNormal2022[[4]] <- fitDataToModel(data = dataTreatment2022[[4]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[4]], fitsLogisticNormal2022[[4]]$par, "logistic")

# Fit of Pea
fitsLogisticNormal2022[[5]] <- fitDataToModel(data = dataTreatment2022[[5]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[5]], fitsLogisticNormal2022[[5]]$par, "logistic")

# Fit of Rye
fitsLogisticNormal2022[[6]] <- fitDataToModel(data = dataTreatment2022[[6]], model = logisticNormal, startParameters = c("r" = 0.012, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[6]], fitsLogisticNormal2022[[6]]$par, "logistic")

# Fit Rye_Faba
fitsLogisticNormal2022[[7]] <- fitDataToModel(data = dataTreatment2022[[7]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[7]], fitsLogisticNormal2022[[7]]$par, "logistic")

# Fit Rye_Pea
fitsLogisticNormal2022[[8]] <- fitDataToModel(data = dataTreatment2022[[8]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[8]], fitsLogisticNormal2022[[8]]$par, "logistic")

# Fit Triticale
fitsLogisticNormal2022[[9]] <- fitDataToModel(data = dataTreatment2022[[9]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[9]], fitsLogisticNormal2022[[9]]$par, "logistic")

# Fit Triticale_Faba
fitsLogisticNormal2022[[10]] <- fitDataToModel(data = dataTreatment2022[[10]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[10]], fitsLogisticNormal2022[[10]]$par, "logistic")

# Fit Triticale_Pea
fitsLogisticNormal2022[[11]] <- fitDataToModel(data = dataTreatment2022[[11]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[11]], fitsLogisticNormal2022[[11]]$par, "logistic")

# Fit Wheat
fitsLogisticNormal2022[[12]] <- fitDataToModel(data = dataTreatment2022[[12]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[12]], fitsLogisticNormal2022[[12]]$par, "logistic")

# Fit Wheat_Faba
fitsLogisticNormal2022[[13]] <- fitDataToModel(data = dataTreatment2022[[13]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[13]], fitsLogisticNormal2022[[13]]$par, "logistic")

# Fit Wheat_Pea
fitsLogisticNormal2022[[14]] <- fitDataToModel(data = dataTreatment2022[[14]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2022[[14]], fitsLogisticNormal2022[[14]]$par, "logistic")

AICLogisticNormal2022 <- sum(unlist(lapply(fitsLogisticNormal2022, function(fits)2 * 3 + 2 * fits$value)))


##############################################
###                                        ###
### Fit transformed logistic normal curves ###
###                                        ###
##############################################      

fitsTLogisticNormal2022 <- rep(list(NA), 14)

# Fit of Barley
fitsTLogisticNormal2022[[1]] <- fitDataToModel(data = dataTreatment2022[[1]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 300, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[1]], fitsTLogisticNormal2022[[1]]$par, "logistic")

# Fit of Barley_Faba
fitsTLogisticNormal2022[[2]] <- fitDataToModel(data = dataTreatment2022[[2]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[2]], fitsTLogisticNormal2022[[2]]$par, "logistic")

# Fit of Barley_Pea
fitsTLogisticNormal2022[[3]] <- fitDataToModel(data = dataTreatment2022[[3]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 350, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[3]], fitsTLogisticNormal2022[[3]]$par, "logistic")

# Fit of Faba
fitsTLogisticNormal2022[[4]] <- fitDataToModel(data = dataTreatment2022[[4]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 300, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[4]], fitsTLogisticNormal2022[[4]]$par, "logistic")

# Fit of Pea
fitsTLogisticNormal2022[[5]] <- fitDataToModel(data = dataTreatment2022[[5]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 300, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[5]], fitsTLogisticNormal2022[[5]]$par, "logistic")

# Fit of Rye
fitsTLogisticNormal2022[[6]] <- fitDataToModel(data = dataTreatment2022[[6]], model = tLogisticNormal, startParameters = c("r" = 0.03, "h" = 280, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[6]], fitsTLogisticNormal2022[[6]]$par, "logistic")

# Fit Rye_Faba
fitsTLogisticNormal2022[[7]] <- fitDataToModel(data = dataTreatment2022[[7]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 300, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[7]], fitsTLogisticNormal2022[[7]]$par, "logistic")

# Fit Rye_Pea
fitsTLogisticNormal2022[[8]] <- fitDataToModel(data = dataTreatment2022[[8]], model = tLogisticNormal, startParameters = c("r" = 0.02, "h" = 280, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[8]], fitsTLogisticNormal2022[[8]]$par, "logistic")

# Fit Triticale
fitsTLogisticNormal2022[[9]] <- fitDataToModel(data = dataTreatment2022[[9]], model = tLogisticNormal, startParameters = c("r" = 0.02, "h" = 280, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[9]], fitsTLogisticNormal2022[[9]]$par, "logistic")

# Fit Triticale_Faba
fitsTLogisticNormal2022[[10]] <- fitDataToModel(data = dataTreatment2022[[10]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 300, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[10]], fitsTLogisticNormal2022[[10]]$par, "logistic")

# Fit Triticale_Pea
fitsTLogisticNormal2022[[11]] <- fitDataToModel(data = dataTreatment2022[[11]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 300, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[11]], fitsTLogisticNormal2022[[11]]$par, "logistic")

# Fit Wheat
fitsTLogisticNormal2022[[12]] <- fitDataToModel(data = dataTreatment2022[[12]], model = tLogisticNormal, startParameters = c("r" = 0.02, "h" = 280, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[12]], fitsTLogisticNormal2022[[12]]$par, "logistic")

# Fit Wheat_Faba
fitsTLogisticNormal2022[[13]] <- fitDataToModel(data = dataTreatment2022[[13]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 300, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[13]], fitsTLogisticNormal2022[[13]]$par, "logistic")

# Fit Wheat_Pea
fitsTLogisticNormal2022[[14]] <- fitDataToModel(data = dataTreatment2022[[14]], model = tLogisticNormal, startParameters = c("r" = 0.02, "h" = 280, "sd" = 0.1, "d" = 0.8), dMax = 1.0)
plotFit(dataTreatment2022[[14]], fitsTLogisticNormal2022[[14]]$par, "logistic")

AICTLogisticNormal2022 <- sum(unlist(lapply(fitsTLogisticNormal2022, function(fits)2 * 3 + 2 * fits$value)))

#################################
###                           ###
### Fit logistic gamma curves ###
###                           ###
#################################     

fitsLogisticGamma2022 <- rep(list(NA), 14)

# Fit of Barley
fitsLogisticGamma2022[[1]] <- fitDataToModel(data = dataTreatment2022[[1]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0))
plotFit(dataTreatment2022[[1]], fitsLogisticGamma2022[[1]]$par, "logistic")

# Fit of Barley_Faba
fitsLogisticGamma2022[[2]] <- fitDataToModel(data = dataTreatment2022[[2]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 350, "shape" = 1.0))
plotFit(dataTreatment2022[[2]], fitsLogisticGamma2022[[2]]$par, "logistic")

# Fit of Barley_Pea
fitsLogisticGamma2022[[3]] <- fitDataToModel(data = dataTreatment2022[[3]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 350, "shape" = 1.0))
plotFit(dataTreatment2022[[3]], fitsLogisticGamma2022[[3]]$par, "logistic")

# Fit of Faba
fitsLogisticGamma2022[[4]] <- fitDataToModel(data = dataTreatment2022[[4]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0))
plotFit(dataTreatment2022[[4]], fitsLogisticGamma2022[[4]]$par, "logistic")

# Fit of Pea
fitsLogisticGamma2022[[5]] <- fitDataToModel(data = dataTreatment2022[[5]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0))
plotFit(dataTreatment2022[[5]], fitsLogisticGamma2022[[5]]$par, "logistic")

# Fit of Rye
fitsLogisticGamma2022[[6]] <- fitDataToModel(data = dataTreatment2022[[6]], model = logisticGamma, startParameters = c("r" = 0.03, "h" = 280, "shape" = 1.0))
plotFit(dataTreatment2022[[6]], fitsLogisticGamma2022[[6]]$par, "logistic")

# Fit Rye_Faba
fitsLogisticGamma2022[[7]] <- fitDataToModel(data = dataTreatment2022[[7]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0))
plotFit(dataTreatment2022[[7]], fitsLogisticGamma2022[[7]]$par, "logistic")

# Fit Rye_Pea
fitsLogisticGamma2022[[8]] <- fitDataToModel(data = dataTreatment2022[[8]], model = logisticGamma, startParameters = c("r" = 0.02, "h" = 280, "shape" = 1.0))
plotFit(dataTreatment2022[[8]], fitsLogisticGamma2022[[8]]$par, "logistic")

# Fit Triticale
fitsLogisticGamma2022[[9]] <- fitDataToModel(data = dataTreatment2022[[9]], model = logisticGamma, startParameters = c("r" = 0.02, "h" = 280, "shape" = 1.0))
plotFit(dataTreatment2022[[9]], fitsLogisticGamma2022[[9]]$par, "logistic")

# Fit Triticale_Faba
fitsLogisticGamma2022[[10]] <- fitDataToModel(data = dataTreatment2022[[10]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0))
plotFit(dataTreatment2022[[10]], fitsLogisticGamma2022[[10]]$par, "logistic")

# Fit Triticale_Pea
fitsLogisticGamma2022[[11]] <- fitDataToModel(data = dataTreatment2022[[11]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0))
plotFit(dataTreatment2022[[11]], fitsLogisticGamma2022[[11]]$par, "logistic")

# Fit Wheat
fitsLogisticGamma2022[[12]] <- fitDataToModel(data = dataTreatment2022[[12]], model = logisticGamma, startParameters = c("r" = 0.02, "h" = 280, "shape" = 1.0))
plotFit(dataTreatment2022[[12]], fitsLogisticGamma2022[[12]]$par, "logistic")

# Fit Wheat_Faba
fitsLogisticGamma2022[[13]] <- fitDataToModel(data = dataTreatment2022[[13]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0))
plotFit(dataTreatment2022[[13]], fitsLogisticGamma2022[[13]]$par, "logistic")

# Fit Wheat_Pea
fitsLogisticGamma2022[[14]] <- fitDataToModel(data = dataTreatment2022[[14]], model = logisticGamma, startParameters = c("r" = 0.02, "h" = 280, "shape" = 1.0))
plotFit(dataTreatment2022[[14]], fitsLogisticGamma2022[[14]]$par, "logistic")

AICLogisticGamma2022 <- sum(unlist(lapply(fitsLogisticGamma2022, function(fits)2 * 3 + 2 * fits$value)))


#############################################
###                                       ###
### Fit transformed logistic gamma curves ###
###                                       ###
#############################################      

fitsTLogisticGamma2022 <- rep(list(NA), 14)

# Fit of Barley
fitsTLogisticGamma2022[[1]] <- fitDataToModel(data = dataTreatment2022[[1]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[1]], fitsTLogisticGamma2022[[1]]$par, "logistic")

# Fit of Barley_Faba
fitsTLogisticGamma2022[[2]] <- fitDataToModel(data = dataTreatment2022[[2]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 350, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[2]], fitsTLogisticGamma2022[[2]]$par, "logistic")

# Fit of Barley_Pea
fitsTLogisticGamma2022[[3]] <- fitDataToModel(data = dataTreatment2022[[3]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 350, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[3]], fitsTLogisticGamma2022[[3]]$par, "logistic")

# Fit of Faba
fitsTLogisticGamma2022[[4]] <- fitDataToModel(data = dataTreatment2022[[4]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[4]], fitsTLogisticGamma2022[[4]]$par, "logistic")

# Fit of Pea
fitsTLogisticGamma2022[[5]] <- fitDataToModel(data = dataTreatment2022[[5]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[5]], fitsTLogisticGamma2022[[5]]$par, "logistic")

# Fit of Rye
fitsTLogisticGamma2022[[6]] <- fitDataToModel(data = dataTreatment2022[[6]], model = tLogisticGamma, startParameters = c("r" = 0.03, "h" = 280, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[6]], fitsTLogisticGamma2022[[6]]$par, "logistic")

# Fit Rye_Faba
fitsTLogisticGamma2022[[7]] <- fitDataToModel(data = dataTreatment2022[[7]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[7]], fitsTLogisticGamma2022[[7]]$par, "logistic")

# Fit Rye_Pea
fitsTLogisticGamma2022[[8]] <- fitDataToModel(data = dataTreatment2022[[8]], model = tLogisticGamma, startParameters = c("r" = 0.02, "h" = 280, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[8]], fitsTLogisticGamma2022[[8]]$par, "logistic")

# Fit Triticale
fitsTLogisticGamma2022[[9]] <- fitDataToModel(data = dataTreatment2022[[9]], model = tLogisticGamma, startParameters = c("r" = 0.02, "h" = 250, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[9]], fitsTLogisticGamma2022[[9]]$par, "logistic")

# Fit Triticale_Faba
fitsTLogisticGamma2022[[10]] <- fitDataToModel(data = dataTreatment2022[[10]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[10]], fitsTLogisticGamma2022[[10]]$par, "logistic")

# Fit Triticale_Pea
fitsTLogisticGamma2022[[11]] <- fitDataToModel(data = dataTreatment2022[[11]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[11]], fitsTLogisticGamma2022[[11]]$par, "logistic")

# Fit Wheat
fitsTLogisticGamma2022[[12]] <- fitDataToModel(data = dataTreatment2022[[12]], model = tLogisticGamma, startParameters = c("r" = 0.02, "h" = 280, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[12]], fitsTLogisticGamma2022[[12]]$par, "logistic")

# Fit Wheat_Faba
fitsTLogisticGamma2022[[13]] <- fitDataToModel(data = dataTreatment2022[[13]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 300, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[13]], fitsTLogisticGamma2022[[13]]$par, "logistic")

# Fit Wheat_Pea
fitsTLogisticGamma2022[[14]] <- fitDataToModel(data = dataTreatment2022[[14]], model = tLogisticGamma, startParameters = c("r" = 0.02, "h" = 280, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2022[[14]], fitsTLogisticGamma2022[[14]]$par, "logistic")

AICTLogisticGamma2022 <- sum(unlist(lapply(fitsTLogisticGamma2022, function(fits)2 * 3 + 2 * fits$value)))

AICLogisticNormal2022
AICTLogisticNormal2022
AICLogisticGamma2022
AICTLogisticGamma2022

# Transformed logistic normal is the best fitting model

# Calculate 95% CIs
parameters2022 <- names(fitsTLogisticNormal2022[[1]]$par)
parameterCIs2022 <- rep(NA, times = 3 * length(parameters2022))
names(parameterCIs2022) <- unlist(lapply(parameters2022, function(par)c(par, paste(par, "Lower", sep = ""), paste(par, "Upper", sep = ""))))
parameterCIs2022 <- rep(list(parameterCIs2022), length(fitsTLogisticNormal2022))

cI2022 <- function(par, x, z, opt, dMax){
  pars <- par
  pars[names(opt)] <- opt
  nll <- tLogisticNormal(pars, x, z, dMax)
  return(nll)
}

for(i in 1:length(fitsTLogisticNormal2022)){
  parameterVecs <- lapply(fitsTLogisticNormal2022, function(fit){
    parStart <- fit$par / 100
    parEnd <- fit$par * 20
    parVec <- list(seq(parStart[1], parEnd[1], length.out = 400), seq(parStart[2], parEnd[2], length.out = 400), seq(parStart[3], parEnd[3], length.out = 400), seq(parStart[4], parEnd[4], length.out = 400))
    names(parVec) <- names(fit$par)
    return(parVec)
  })
  optPar <- fitsTLogisticNormal2022[[i]]$par
  parameterCIs2022[[i]][names(optPar)] <- optPar
  parameterCIs2022[[i]][unlist((lapply(names(optPar), function(n)c(paste(n, "Lower", sep = ""), paste(n, "Upper", sep = "")))))] <-
    unlist(lapply(1:length(parameterVecs[[i]]), function(j){
      t(tibble(calculate95CI(parStartVec = optPar[-which(names(optPar) == names(fitsTLogisticNormal2022[[i]]$par)[[j]])],
                             optVec = parameterVecs[[i]][j],
                             model = cI2022,
                             optNLL = fitsTLogisticNormal2022[[i]]$value,
                             x = dataTreatment2022[[i]][,1][[1]],
                             z = dataTreatment2022[[i]][,2][[1]],
                             dMax = 1.0)))
    }))
  
}
names(parameterCIs2022) <- names(dataTreatment2022)

# Write to table
CIOut2022 <- tibble(Treatment = names(parameterCIs2022)) %>% 
  bind_cols(bind_rows(parameterCIs2022))
write.table(CIOut2022, "canopy_cover/canopy_cover_tLogisticNormal_parameters_2022.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

CIOut2022Table <- CIOut2022 %>% 
  mutate(r = paste0(round(r, 4), " (", round(rLower, 4), "-", round(rUpper, 4), ")"),
         h = paste0(round(h, 2), " (", round(hLower, 2), "-", round(hUpper, 2), ")"),
         sd = paste0(round(sd, 4), " (", round(sdLower, 4), "-", round(sdUpper, 4), ")"),
         d = paste0(round(d, 3), " (", round(dLower, 3), "-", round(dUpper, 3), ")")) %>% 
  dplyr::select(Treatment, r, h, sd, d)

print(xtable(CIOut2022Table, type = "latex"), include.rownames = FALSE)

dataTreatment2022 <- lapply(1:length(dataTreatment2022), function(i){
  dat <- dataTreatment2022[[i]]
  treatment <- names(dataTreatment2022)[i]
  dat <- dat %>% 
    mutate(Treatment = treatment)
})

dataTreatment2022 <- lapply(dataTreatment2022, function(Data){
  Average <- aggregate(Data$CanopyCover, list(Data$Time), function(x)mean(x, na.rm = TRUE))
  Data$Average <- Average$x[match(Data$Time, Average$Group.1)]
  return(Data)
})

# Plot triple combination graphs
geomPoints <- lapply(dataTreatment2022, function(d){
  return(geom_point(data = d, aes(x = Time, y = Average, col = as.factor(Treatment[1])), size = 2, shape = 1))
})

geomFunctions <- lapply(1:length(fitsTLogisticNormal2022), function(i){
  o <- fitsTLogisticNormal2022[[i]]
  return(geom_function(fun = function(x)o$par["d"] * (exp(o$par["r"]*(x - o$par["h"])))/(1 + exp(o$par["r"]*(x - o$par["h"]))),
                       aes(col = as.factor(dataTreatment2022[[i]]$Treatment[1])),
                       size = 1.5))
})

### Get all triple combinations
treatmentsToCompare <- unique(data2022$Treatment)

triples <- unlist(lapply(c("Barley", "Rye","Triticale", "Wheat"), function(c){
  lapply(c("Pea", "Faba"), function(l){
    return(c(c, paste(c, l, sep = "_"), l))
  })
}), recursive = FALSE)
names(dataTreatment2022) <- lapply(dataTreatment2022, function(l)l$Treatment[1])
names(geomPoints) <- names(dataTreatment2022)
names(geomFunctions) <- names(dataTreatment2022)

triplePlots <- lapply(triples, function(triple){
  plotMultipleGraphs(geomPoints[triple], geomFunctions[triple], triple, colourVector[c(1, 3, 2)], XLab = bquote(T[c]), YLab = bquote(p[cc]), Title = "")
})

# Create a combined plot for intercrops and normal density sole crops
combinedPlot <- plotMultipleGraphs(geomPoints[treatmentsToCompare], geomFunctions[treatmentsToCompare], treatmentsToCompare, c(colourVector[1], colourVector[2], colourVector[5], colourVector[4], colourVector[3], colourVector[6]), XLab = bquote(T[c]), YLab = bquote(p[cc]), Title = "")

### Combine data for AIC comparisons
combinations2022 <- lapply(triples, function(t)combn(t, 2))


# Fit of Barley-Pea
startParametersBarleyPea <- list(c("r" = 0.008, "h" = 400, "sd" = 0.1, "d" = 0.9),
                                 c("r" = 0.01, "h" = 380, "sd" = 0.1, "d" = 0.8),
                                 c("r" = 0.008, "h" = 380, "sd" = 0.1, "d" = 0.9),
                                 c("r" = 0.009, "h" = 420, "sd" = 0.1, "d" = 0.9),
                                 c("r" = 0.009, "h" = 410, "sd" = 0.1, "d" = 0.9))
combinationAICBarleyPea <- round(fitCombinations(data = dataTreatment2022, triple = triples[[1]], model = tLogisticNormal, startParameters = startParametersBarleyPea), 2)

# Fit of Barley-Faba
startParametersBarleyFaba <- list(c("r" = 0.008, "h" = 400, "sd" = 0.1, "d" = 0.9),
                                  c("r" = 0.01, "h" = 380, "sd" = 0.1, "d" = 0.8),
                                  c("r" = 0.008, "h" = 380, "sd" = 0.1, "d" = 0.9),
                                  c("r" = 0.009, "h" = 420, "sd" = 0.1, "d" = 0.9),
                                  c("r" = 0.009, "h" = 410, "sd" = 0.1, "d" = 0.9))
combinationAICBarleyFaba <- round(fitCombinations(data = dataTreatment2022, triple = triples[[2]], model = tLogisticNormal, startParameters = startParametersBarleyFaba), 2)

# Fit of Rye-Pea
startParametersRyePea <- list(c("r" = 0.01, "h" = 280, "sd" = 0.1, "d" = 0.95),
                              c("r" = 0.008, "h" = 320, "sd" = 0.1, "d" = 0.90),
                              c("r" = 0.008, "h" = 310, "sd" = 0.1, "d" = 0.9),
                              c("r" = 0.015, "h" = 320, "sd" = 0.1, "d" = 0.9),
                              c("r" = 0.01, "h" = 320, "sd" = 0.1, "d" = 0.9))
combinationAICRyePea <- round(fitCombinations(data = dataTreatment2022, triple = triples[[3]], model = tLogisticNormal, startParameters = startParametersRyePea), 2)

# Fit of Rye-Faba
startParametersRyeFaba <- list(c("r" = 0.01, "h" = 280, "sd" = 0.1, "d" = 0.95),
                               c("r" = 0.008, "h" = 320, "sd" = 0.1, "d" = 0.90),
                               c("r" = 0.008, "h" = 380, "sd" = 0.1, "d" = 0.9),
                               c("r" = 0.01, "h" = 320, "sd" = 0.1, "d" = 0.9),
                               c("r" = 0.01, "h" = 320, "sd" = 0.1, "d" = 0.9))
combinationAICRyeFaba <- round(fitCombinations(data = dataTreatment2022, triple = triples[[4]], model = tLogisticNormal, startParameters = startParametersRyeFaba), 2)

# Fit of Triticale-Pea
startParametersTriticalePea <- list(c("r" = 0.01, "h" = 280, "sd" = 0.1, "d" = 0.95),
                                    c("r" = 0.008, "h" = 320, "sd" = 0.1, "d" = 0.90),
                                    c("r" = 0.008, "h" = 380, "sd" = 0.1, "d" = 0.9),
                                    c("r" = 0.01, "h" = 320, "sd" = 0.1, "d" = 0.9),
                                    c("r" = 0.01, "h" = 320, "sd" = 0.1, "d" = 0.9))
combinationAICTriticalePea <- round(fitCombinations(data = dataTreatment2022, triple = triples[[5]], model = tLogisticNormal, startParameters = startParametersTriticalePea), 2)

# Fit of Triticale-Faba
startParametersTriticaleFaba <- list(c("r" = 0.01, "h" = 280, "sd" = 0.1, "d" = 0.95),
                                     c("r" = 0.008, "h" = 320, "sd" = 0.1, "d" = 0.90),
                                     c("r" = 0.008, "h" = 380, "sd" = 0.1, "d" = 0.9),
                                     c("r" = 0.01, "h" = 320, "sd" = 0.1, "d" = 0.9),
                                     c("r" = 0.01, "h" = 320, "sd" = 0.1, "d" = 0.9))
combinationAICTriticaleFaba <- round(fitCombinations(data = dataTreatment2022, triple = triples[[6]], model = tLogisticNormal, startParameters = startParametersTriticaleFaba), 2)

# Fit of Wheat-Pea
startParametersWheatPea <- list(c("r" = 0.01, "h" = 350, "sd" = 0.1, "d" = 0.95),
                                c("r" = 0.008, "h" = 320, "sd" = 0.1, "d" = 0.90),
                                c("r" = 0.008, "h" = 380, "sd" = 0.1, "d" = 0.9),
                                c("r" = 0.01, "h" = 320, "sd" = 0.1, "d" = 0.9),
                                c("r" = 0.01, "h" = 320, "sd" = 0.1, "d" = 0.9))
combinationAICWheatPea <- round(fitCombinations(data = dataTreatment2022, triple = triples[[7]], model = tLogisticNormal, startParameters = startParametersWheatPea), 2)

# Fit of Wheat-Faba
startParametersWheatFaba <- list(c("r" = 0.01, "h" = 420, "sd" = 0.1, "d" = 0.95),
                                 c("r" = 0.008, "h" = 320, "sd" = 0.1, "d" = 0.90),
                                 c("r" = 0.008, "h" = 380, "sd" = 0.1, "d" = 0.9),
                                 c("r" = 0.01, "h" = 320, "sd" = 0.1, "d" = 0.9),
                                 c("r" = 0.01, "h" = 320, "sd" = 0.1, "d" = 0.9))
combinationAICWheatFaba <- round(fitCombinations(data = dataTreatment2022, triple = triples[[8]], model = tLogisticNormal, startParameters = startParametersWheatFaba), 2)

combinationAIC2022 <- tibble(Treatment = rep(unique(data2022$Treatment)[1:8], each = 3),
                             Combination = rep(c("IC & C & L", "IC + C & L", "IC + L & C"), times = 8),
                             AIC = c(combinationAICBarleyPea,
                                     combinationAICBarleyFaba,
                                     combinationAICRyePea,
                                     combinationAICRyeFaba,
                                     combinationAICTriticalePea,
                                     combinationAICTriticaleFaba,
                                     combinationAICWheatPea,
                                     combinationAICWheatFaba),
                             Best = c("",
                                      "*",
                                      "",
                                      "",
                                      "*",
                                      "",
                                      "*",
                                      "",
                                      "",
                                      "*",
                                      "",
                                      "",
                                      "*",
                                      "*",
                                      "",
                                      "",
                                      "*",
                                      "",
                                      "*",
                                      "",
                                      "",
                                      "*",
                                      "",
                                      ""))

print(xtable(combinationAIC2022, type = "latex"), include.rownames = FALSE)

# Print all plots to PDF
pdf("canopy_cover/canopy_cover_tLogisticNormal_graphs_2022.pdf")
print(triplePlots)
print(combinedPlot)
dev.off()

# Print all plots to jpeg
dir.create("figures/canopy_cover")

# Combined plot

jpeg("figures/canopy_cover/plotCanopyCover2022TripleCombined.jpg", units = "px", width = 900, height = 1600, quality = 600)
triplePlots[[1]] + theme(legend.position = "none", plot.subtitle = element_text(hjust = -0.4)) + labs(y = "", subtitle = bquote({P[cc]}~~~~~~~~~~~"Barley_Pea")) + 
  ylim(c(0, 1.0)) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  geom_segment(aes(x = 525, y = 0.7, xend = 360, yend = 0.70),
               arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(aes(x = 525, y = 0.55, xend = 370, yend = 0.55),
               arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(aes(x = 525, y = 0.4, xend = 380, yend = 0.40),
               arrow = arrow(length = unit(0.5, "cm"))) +
  annotate("text", x = 525, y = c(0.7, 0.55, 0.4), label = c("Cereal", "Intercrop", "Legume"), size = 10, hjust = 0.0, col = modelColors[c("Cereal", "Intercrop", "Legume")]) +
triplePlots[[2]] + theme(legend.position = "none", plot.subtitle = element_text(hjust = -0.4)) + labs(y = "", subtitle = bquote(~~~~~~~~~~~~~~~~~~"Barley_Faba")) + 
  ylim(c(0, 1.0)) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
triplePlots[[3]] + theme(legend.position = "none", plot.subtitle = element_text(hjust = -0.4)) + labs(y = "", subtitle = bquote({P[cc]}~~~~~~~~~~~"Rye_Pea")) + 
  ylim(c(0, 1.0)) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
triplePlots[[4]] + theme(legend.position = "none", plot.subtitle = element_text(hjust = -0.4)) + labs(y = "", subtitle = bquote(~~~~~~~~~~~~~~~~~~"Rye_Faba")) + 
  ylim(c(0, 1.0)) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
triplePlots[[5]] + theme(legend.position = "none", plot.subtitle = element_text(hjust = -0.4)) + labs(y = "", subtitle = bquote({P[cc]}~~~~~~~~~~~"Triticale_Pea")) + 
  ylim(c(0, 1.0)) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
triplePlots[[6]] + theme(legend.position = "none", plot.subtitle = element_text(hjust = -0.4)) + labs(y = "", subtitle = bquote(~~~~~~~~~~~~~~~~~~"Triticale_Faba")) + 
  ylim(c(0, 1.0)) + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
triplePlots[[7]] + theme(legend.position = "none", plot.subtitle = element_text(hjust = -0.4)) + labs(y = "", subtitle = bquote({P[cc]}~~~~~~~~~~~"Wheat_Pea")) + 
  ylim(c(0, 1.0)) + 
triplePlots[[8]] + theme(plot.subtitle = element_text(hjust = -0.4)) +
  scale_color_manual(values = modelColors,
                     labels = c("Cereal_Legume", "Cereal", "Legume")) + labs(y = "", subtitle = bquote(~~~~~~~~~~~~~~~~~~"Wheat_Faba")) + 
  ylim(c(0, 1.0)) + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none") + 
plot_layout(design = "
  12
  34
  56
  78
", axis_titles = "collect") +
plot_annotation(tag_levels = "a")
dev.off()


##
## Segmented canopy cover
## 
Treatments2022 <- read_delim("PlotTreatments2022.txt", delim = "\t", col_names = TRUE)

dataSegCover2022 <- read_delim("cover_data_2022.txt", delim = "\t") %>% 
  mutate(relativeFaba = faba / all,
         relativeCereal = cereal / all,
         relativeWeed = weed / all,
         ground = 1 - all) %>% 
  mutate(Time = tempData2022$CumulTemp[match(time, tempData2022$Time)]) %>% 
  dplyr::select(plot, Time, faba, weed, cereal) %>% 
  rename(Plot = plot,
         Faba = faba,
         Weed = weed,
         Cereal = cereal) %>% 
  mutate(Treatment = str_replace(Treatments2022$Treatment[match(Plot, Treatments2022$Plot)], "Bean", "Faba")) %>% 
  dplyr::select(!Plot)

dataSegCover2022Split <- split(dataSegCover2022, dataSegCover2022$Treatment)

dataSeg2022BarleyFaba <- dataSegCover2022Split[[1]]
dataSeg2022BarleyFabaLong <- dataSeg2022BarleyFaba %>% 
  pivot_longer(cols = 2:4, names_to = "Type", values_to = "PCover")

(plotSeg2022BarleyFaba <- ggplot(data = dataSeg2022BarleyFabaLong, aes(x = Time, y = PCover, col = Type)) +
    stat_smooth(method = "loess", aes(fill = Type), alpha = 0.2, method.args = c(degree = 2)) +
    geom_point(shape = 1, aes(fill = Type), size = 2) +  # Adjusted to include fill mapping
    scale_color_manual(values = modelColors) +  # Set line colors
    scale_fill_manual(values = modelColors) +   # Set fill colors
    theme_classic(base_size = 40) +
    labs(title = "Segmented canopy cover",
         subtitle = "2022, Barley-Faba",
         x = bquote(T[c]),
         y = bquote(p[cc]),
         col = "Type"))

dataSeg2022RyeFaba <- dataSegCover2022Split[[2]]
dataSeg2022RyeFabaLong <- dataSeg2022RyeFaba %>% 
  pivot_longer(cols = 2:4, names_to = "Type", values_to = "PCover")

(plotSeg2022RyeFaba <- ggplot(data = dataSeg2022RyeFabaLong, aes(x = Time, y = PCover, col = Type)) +
    stat_smooth(method = "loess", aes(fill = Type), alpha = 0.2, method.args = c(degree = 2)) +
    geom_point(shape = 1, aes(fill = Type), size = 2) +  # Adjusted to include fill mapping
    scale_color_manual(values = modelColors) +  # Set line colors
    scale_fill_manual(values = modelColors) +   # Set fill colors
    theme_classic(base_size = 40) +
    labs(title = "Segmented canopy cover",
         subtitle = "2022, Rye-Faba",
         x = bquote(T[c]),
         y = bquote(p[cc]),
         col = "Type"))

dataSeg2022TriticaleFaba <- dataSegCover2022Split[[3]]
dataSeg2022TriticaleFabaLong <- dataSeg2022TriticaleFaba %>% 
  pivot_longer(cols = 2:4, names_to = "Type", values_to = "PCover")

(plotSeg2022TriticaleFaba <- ggplot(data = dataSeg2022TriticaleFabaLong, aes(x = Time, y = PCover, col = Type)) +
    stat_smooth(method = "loess", aes(fill = Type), alpha = 0.2, method.args = c(degree = 2)) +
    geom_point(shape = 1, aes(fill = Type), size = 2) +  # Adjusted to include fill mapping
    scale_color_manual(values = modelColors) +  # Set line colors
    scale_fill_manual(values = modelColors) +   # Set fill colors
    theme_classic(base_size = 40) +
    labs(title = "Segmented canopy cover",
         subtitle = "2022, Triticale-Faba",
         x = bquote(T[c]),
         y = bquote(p[cc]),
         col = "Type"))

dataSeg2022WheatFaba <- dataSegCover2022Split[[4]]
dataSeg2022WheatFabaLong <- dataSeg2022WheatFaba %>% 
  pivot_longer(cols = 2:4, names_to = "Type", values_to = "PCover")

(plotSeg2022WheatFaba <- ggplot(data = dataSeg2022WheatFabaLong, aes(x = Time, y = PCover, col = Type)) +
    stat_smooth(method = "loess", aes(fill = Type), alpha = 0.2, method.args = c(degree = 2)) +
    geom_point(shape = 1, aes(fill = Type), size = 2) +  # Adjusted to include fill mapping
    scale_color_manual(values = modelColors) +  # Set line colors
    scale_fill_manual(values = modelColors) +   # Set fill colors
    theme_classic(base_size = 40) +
    labs(title = "Segmented canopy cover",
         subtitle = "2022, Wheat-Faba",
         x = bquote(T[c]),
         y = bquote(p[cc]),
         col = "Type"))


########################
##                    ##
## Experiment 2023 A  ##
##                    ##
########################

sowingDate2023 <- as.Date("2023-03-02")

# Read data
data2023A <- read_xlsx("WSF_2023_data.xlsx", sheet = "CanopyCoverA", range = "A1:G721", col_names = TRUE) %>% 
  mutate(DAS = as.integer(as.Date(Date) - sowingDate2023)) %>% 
  rename(Time = DAS, Val = Pcover) %>% 
  filter(Time != 90) %>% 
  dplyr::select(Treatment, Weeds, Time, Val)

tempData2023 <- read_delim("TempData2023.txt", delim = "\t") %>% 
  mutate(Time = as.integer(as.Date(Date) - sowingDate2023))

cumulativeTemp2023 <- ((tempData2023$TMin + tempData2023$TMax) / 2)
cumulativeTemp2023[which(cumulativeTemp2023 < 0.0)] <- 0
cumulativeTemp2023 <- cumsum(cumulativeTemp2023)

tempData2023 <- tempData2023 %>% 
  mutate(cumulativeTemp = cumulativeTemp2023)

dataWeeds2023A <- data2023A %>% 
  filter(Weeds == "Y") %>% 
  dplyr::select(-Weeds) %>% 
  filter(Time != 91) %>% 
  mutate(Time = tempData2023$cumulativeTemp[match(Time, tempData2023$Time)])
dataNoWeeds2023A <- data2023A %>% 
  filter(Weeds == "N") %>% 
  dplyr::select(-Weeds) %>% 
  filter(Time != 91) %>% 
  mutate(Time = tempData2023$cumulativeTemp[match(Time, tempData2023$Time)])

## Manual fits
plotFit <- function(data, par, fun){
  plot <- ggplot() +
    geom_point(data = data, aes(x = Time, y = Val), shape = 1, size = 2) +
    geom_function(fun = function(x)allDeterministicFunctions[[fun]](par, x),
                  size = 1.5)
  return(plot)
}


##################
##################
####          ####
#### No weeds ####
####          ####
##################
##################

dataTreatment2023ANW <- split(dataNoWeeds2023A, dataNoWeeds2023A$Treatment)
dataTreatment2023ANW <- lapply(dataTreatment2023ANW, function(dat){
  datNew <- dat %>% 
    dplyr::select(Time, Val)
  return(datNew)
})




##################################
###                            ###
### Fit logistic normal curves ###
###                            ###
##################################     

fitsNWLogisticNormal2023A <- rep(list(NA), 8)

# Fit of 1T:1F
fitsNWLogisticNormal2023A[[1]] <- fitDataToModel(data = dataTreatment2023ANW[[1]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 520, "sd" = 0.1))
plotFit(dataTreatment2023ANW[[1]], fitsNWLogisticNormal2023A[[1]]$par, "logistic")

# Fit of 1T:3F
fitsNWLogisticNormal2023A[[2]] <- fitDataToModel(data = dataTreatment2023ANW[[2]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 520, "sd" = 0.1))
plotFit(dataTreatment2023ANW[[2]], fitsNWLogisticNormal2023A[[2]]$par, "logistic")

# Fit of 3T:1F
fitsNWLogisticNormal2023A[[3]] <- fitDataToModel(data = dataTreatment2023ANW[[3]], model = logisticNormal, startParameters = c("r" = 0.012, "h" = 475, "sd" = 0.1))
plotFit(dataTreatment2023ANW[[3]], fitsNWLogisticNormal2023A[[3]]$par, "logistic")

# Fit of F
fitsNWLogisticNormal2023A[[4]] <- fitDataToModel(data = dataTreatment2023ANW[[4]], model = logisticNormal, startParameters = c("r" = 0.005, "h" = 650, "sd" = 0.1))
plotFit(dataTreatment2023ANW[[4]], fitsNWLogisticNormal2023A[[4]]$par, "logistic")

# Fit of F+
fitsNWLogisticNormal2023A[[5]] <- fitDataToModel(data = dataTreatment2023ANW[[5]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 575, "sd" = 0.1))
plotFit(dataTreatment2023ANW[[5]], fitsNWLogisticNormal2023A[[5]]$par, "logistic")

# Fit of T
fitsNWLogisticNormal2023A[[6]] <- fitDataToModel(data = dataTreatment2023ANW[[6]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 400, "sd" = 0.1))
plotFit(dataTreatment2023ANW[[6]], fitsNWLogisticNormal2023A[[6]]$par, "logistic")

# Fit of T+
fitsNWLogisticNormal2023A[[7]] <- fitDataToModel(data = dataTreatment2023ANW[[7]], model = logisticNormal, startParameters = c("r" = 0.012, "h" = 350, "sd" = 0.1))
plotFit(dataTreatment2023ANW[[7]], fitsNWLogisticNormal2023A[[7]]$par, "logistic")

# Fit TF-M
fitsNWLogisticNormal2023A[[8]] <- fitDataToModel(data = dataTreatment2023ANW[[8]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 520, "sd" = 0.1))
plotFit(dataTreatment2023ANW[[8]], fitsNWLogisticNormal2023A[[8]]$par, "logistic")

AICLogisticNormal2023ANW <- sum(unlist(lapply(fitsNWLogisticNormal2023A, function(fits)2 * 3 + 2 * fits$value)))


##############################################
###                                        ###
### Fit transformed logistic normal curves ###
###                                        ###
##############################################      

fitsNWTLogisticNormal2023A <- rep(list(NA), 8)

# Fit of 1T:1F
fitsNWTLogisticNormal2023A[[1]] <- fitDataToModel(data = dataTreatment2023ANW[[1]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 480, "sd" = 0.1, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2023ANW[[1]], fitsNWTLogisticNormal2023A[[1]]$par, "tLogistic")

# Fit of 1T:3F
fitsNWTLogisticNormal2023A[[2]] <- fitDataToModel(data = dataTreatment2023ANW[[2]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 550, "sd" = 0.1, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2023ANW[[2]], fitsNWTLogisticNormal2023A[[2]]$par, "tLogistic")

# Fit of 3T:1F
fitsNWTLogisticNormal2023A[[3]] <- fitDataToModel(data = dataTreatment2023ANW[[3]], model = tLogisticNormal, startParameters = c("r" = 0.012, "h" = 475, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023ANW[[3]], fitsNWTLogisticNormal2023A[[3]]$par, "tLogistic")

# Fit of F
fitsNWTLogisticNormal2023A[[4]] <- fitDataToModel(data = dataTreatment2023ANW[[4]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 640, "sd" = 0.1, "d" = 0.95), dMax = 1.0)
plotFit(dataTreatment2023ANW[[4]], fitsNWTLogisticNormal2023A[[4]]$par, "tLogistic")

# Fit of F+
fitsNWTLogisticNormal2023A[[5]] <- fitDataToModel(data = dataTreatment2023ANW[[5]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 575, "sd" = 0.1, "d" = 0.95), dMax = 1.0)
plotFit(dataTreatment2023ANW[[5]], fitsNWTLogisticNormal2023A[[5]]$par, "tLogistic")

# Fit of T
fitsNWTLogisticNormal2023A[[6]] <- fitDataToModel(data = dataTreatment2023ANW[[6]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 400, "sd" = 0.1, "d" = 0.90), dMax = 1.0)
plotFit(dataTreatment2023ANW[[6]], fitsNWTLogisticNormal2023A[[6]]$par, "tLogistic")

# Fit of T+
fitsNWTLogisticNormal2023A[[7]] <- fitDataToModel(data = dataTreatment2023ANW[[7]], model = tLogisticNormal, startParameters = c("r" = 0.012, "h" = 350, "sd" = 0.1, "d" = 0.95), dMax = 1.0)
plotFit(dataTreatment2023ANW[[7]], fitsNWTLogisticNormal2023A[[7]]$par, "tLogistic")

# Fit TF-M
fitsNWTLogisticNormal2023A[[8]] <- fitDataToModel(data = dataTreatment2023ANW[[8]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 520, "sd" = 0.1, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2023ANW[[8]], fitsNWTLogisticNormal2023A[[8]]$par, "tLogistic")

AICTLogisticNormal2023ANW <- sum(unlist(lapply(fitsNWTLogisticNormal2023A, function(fits)2 * 4 + 2 * fits$value)))


#################################
###                           ###
### Fit logistic gamma curves ###
###                           ###
#################################     

fitsNWLogisticGamma2023A <- rep(list(NA), 8)

# Fit of 1T:1F
fitsNWLogisticGamma2023A[[1]] <- fitDataToModel(data = dataTreatment2023ANW[[1]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0))
plotFit(dataTreatment2023ANW[[1]], fitsNWLogisticGamma2023A[[1]]$par, "logistic")

# Fit of 1T:3F
fitsNWLogisticGamma2023A[[2]] <- fitDataToModel(data = dataTreatment2023ANW[[2]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0))
plotFit(dataTreatment2023ANW[[2]], fitsNWLogisticGamma2023A[[2]]$par, "logistic")

# Fit of 3T:1F
fitsNWLogisticGamma2023A[[3]] <- fitDataToModel(data = dataTreatment2023ANW[[3]], model = logisticGamma, startParameters = c("r" = 0.012, "h" = 475, "shape" = 1.0))
plotFit(dataTreatment2023ANW[[3]], fitsNWLogisticGamma2023A[[3]]$par, "logistic")

# Fit of F
fitsNWLogisticGamma2023A[[4]] <- fitDataToModel(data = dataTreatment2023ANW[[4]], model = logisticGamma, startParameters = c("r" = 0.005, "h" = 650, "shape" = 1.0))
plotFit(dataTreatment2023ANW[[4]], fitsNWLogisticGamma2023A[[4]]$par, "logistic")

# Fit of F+
fitsNWLogisticGamma2023A[[5]] <- fitDataToModel(data = dataTreatment2023ANW[[5]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 575, "shape" = 1.0))
plotFit(dataTreatment2023ANW[[5]], fitsNWLogisticGamma2023A[[5]]$par, "logistic")

# Fit of T
fitsNWLogisticGamma2023A[[6]] <- fitDataToModel(data = dataTreatment2023ANW[[6]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 400, "shape" = 1.0))
plotFit(dataTreatment2023ANW[[6]], fitsNWLogisticGamma2023A[[6]]$par, "logistic")

# Fit of T+
fitsNWLogisticGamma2023A[[7]] <- fitDataToModel(data = dataTreatment2023ANW[[7]], model = logisticGamma, startParameters = c("r" = 0.012, "h" = 350, "shape" = 1.0))
plotFit(dataTreatment2023ANW[[7]], fitsNWLogisticGamma2023A[[7]]$par, "logistic")

# Fit TF-M
fitsNWLogisticGamma2023A[[8]] <- fitDataToModel(data = dataTreatment2023ANW[[8]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0))
plotFit(dataTreatment2023ANW[[8]], fitsNWLogisticGamma2023A[[8]]$par, "logistic")

AICLogisticGamma2023ANW <- sum(unlist(lapply(fitsNWLogisticGamma2023A, function(fits)2 * 3 + 2 * fits$value)))


#############################################
###                                       ###
### Fit transformed logistic gamma curves ###
###                                       ###
#############################################      

fitsNWTLogisticGamma2023A <- rep(list(NA), 8)

# Fit of 1T:1F
fitsNWTLogisticGamma2023A[[1]] <- fitDataToModel(data = dataTreatment2023ANW[[1]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023ANW[[1]], fitsNWTLogisticGamma2023A[[1]]$par, "tLogistic")

# Fit of 1T:3F
fitsNWTLogisticGamma2023A[[2]] <- fitDataToModel(data = dataTreatment2023ANW[[2]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023ANW[[2]], fitsNWTLogisticGamma2023A[[2]]$par, "tLogistic")

# Fit of 3T:1F
fitsNWTLogisticGamma2023A[[3]] <- fitDataToModel(data = dataTreatment2023ANW[[3]], model = tLogisticGamma, startParameters = c("r" = 0.012, "h" = 475, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023ANW[[3]], fitsNWTLogisticGamma2023A[[3]]$par, "tLogistic")

# Fit of F
fitsNWTLogisticGamma2023A[[4]] <- fitDataToModel(data = dataTreatment2023ANW[[4]], model = tLogisticGamma, startParameters = c("r" = 0.012, "h" = 650, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023ANW[[4]], fitsNWTLogisticGamma2023A[[4]]$par, "tLogistic")

# Fit of F+
fitsNWTLogisticGamma2023A[[5]] <- fitDataToModel(data = dataTreatment2023ANW[[5]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 575, "shape" = 1.0, "d" = 1.0), dMax = 1.0)
plotFit(dataTreatment2023ANW[[5]], fitsNWTLogisticGamma2023A[[5]]$par, "tLogistic")

# Fit of T
fitsNWTLogisticGamma2023A[[6]] <- fitDataToModel(data = dataTreatment2023ANW[[6]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 400, "shape" = 1.0, "d" = 0.90), dMax = 1.0)
plotFit(dataTreatment2023ANW[[6]], fitsNWTLogisticGamma2023A[[6]]$par, "tLogistic")

# Fit of T+
fitsNWTLogisticGamma2023A[[7]] <- fitDataToModel(data = dataTreatment2023ANW[[7]], model = tLogisticGamma, startParameters = c("r" = 0.012, "h" = 350, "shape" = 1.0, "d" = 0.95), dMax = 1.0)
plotFit(dataTreatment2023ANW[[7]], fitsNWTLogisticGamma2023A[[7]]$par, "tLogistic")

# Fit TF-M
fitsNWTLogisticGamma2023A[[8]] <- fitDataToModel(data = dataTreatment2023ANW[[8]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2023ANW[[8]], fitsNWTLogisticGamma2023A[[8]]$par, "tLogistic")

AICTLogisticGamma2023ANW <- sum(unlist(lapply(fitsNWTLogisticGamma2023A, function(fits)2 * 4 + 2 * fits$value)))

AICLogisticNormal2023ANW
AICTLogisticNormal2023ANW
AICLogisticGamma2023ANW
AICTLogisticGamma2023ANW

# Transformed logistic normal is the best fitting model

# Calculate 95% CIs
parameters2023A <- names(fitsNWTLogisticNormal2023A[[1]]$par)
parameterCIs2023A <- rep(NA, times = 3 * length(parameters2023A))
names(parameterCIs2023A) <- unlist(lapply(parameters2023A, function(par)c(par, paste(par, "Lower", sep = ""), paste(par, "Upper", sep = ""))))
parameterCIs2023A <- rep(list(parameterCIs2023A), length(fitsNWTLogisticNormal2023A))

cI2023A <- function(par, x, z, opt, dMax){
  pars <- par
  pars[names(opt)] <- opt
  nll <- tLogisticNormal(pars, x, z, dMax)
  return(nll)
}

for(i in 1:length(fitsNWTLogisticNormal2023A)){
  parameterVecs <- lapply(fitsNWTLogisticNormal2023A, function(fit){
    parStart <- fit$par / 100
    parEnd <- fit$par * 20
    parVec <- list(seq(parStart[1], parEnd[1], length.out = 400), seq(parStart[2], parEnd[2], length.out = 400), seq(parStart[3], parEnd[3], length.out = 400), seq(parStart[4], parEnd[4], length.out = 400))
    names(parVec) <- names(fit$par)
    return(parVec)
  })
  optPar <- fitsNWTLogisticNormal2023A[[i]]$par
  parameterCIs2023A[[i]][names(optPar)] <- optPar
  parameterCIs2023A[[i]][unlist((lapply(names(optPar), function(n)c(paste(n, "Lower", sep = ""), paste(n, "Upper", sep = "")))))] <-
    unlist(lapply(1:length(parameterVecs[[i]]), function(j){
      t(tibble(calculate95CI(parStartVec = optPar[-which(names(optPar) == names(fitsNWTLogisticNormal2023A[[i]]$par)[[j]])],
                             optVec = parameterVecs[[i]][j],
                             model = cI2023A,
                             optNLL = fitsNWTLogisticNormal2023A[[i]]$value,
                             x = dataTreatment2023ANW[[i]][,1][[1]],
                             z = dataTreatment2023ANW[[i]][,2][[1]],
                             dMax = 1.0)))
    }))
  
}
names(parameterCIs2023A) <- names(dataTreatment2023ANW)

# Write to table
CIOut2023A <- tibble(Treatment = names(parameterCIs2023A)) %>% 
  bind_cols(bind_rows(parameterCIs2023A))
write.table(CIOut2023A, "canopy_cover/canopy_cover_NW_tLogisticNormal_parameters_2023A.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

CIOut2023ATable <- CIOut2023A %>% 
  mutate(r = paste0(round(r, 4), " (", round(rLower, 4), "-", round(rUpper, 4), ")"),
         h = paste0(round(h, 2), " (", round(hLower, 2), "-", round(hUpper, 2), ")"),
         sd = paste0(round(sd, 4), " (", round(sdLower, 4), "-", round(sdUpper, 4), ")"),
         d = paste0(round(d, 3), " (", round(dLower, 3), "-", round(dUpper, 3), ")")) %>% 
  dplyr::select(Treatment, r, h, sd, d)

print(xtable(CIOut2023ATable, type = "latex"), include.rownames = FALSE)

dataTreatment2023ANWM <- lapply(1:length(dataTreatment2023ANW), function(i){
  dat <- dataTreatment2023ANW[[i]]
  treatment <- names(dataTreatment2023ANW)[i]
  dat <- dat %>% 
    mutate(Treatment = treatment)
})

dataTreatment2023ANWM <- lapply(dataTreatment2023ANWM, function(Data){
  Average <- aggregate(Data$Val, list(Data$Time), function(x)mean(x, na.rm = TRUE))
  Data$Average <- Average$x[match(Data$Time, Average$Group.1)]
  return(Data)
})

# Plot triple combination graphs
geomPoints <- lapply(dataTreatment2023ANWM, function(d){
  return(geom_point(data = d, aes(x = Time, y = Average, col = as.factor(Treatment[1])), size = 2, shape = 1))
})

geomFunctions <- lapply(1:length(fitsNWTLogisticNormal2023A), function(i){
  o <- fitsNWTLogisticNormal2023A[[i]]
  return(geom_function(fun = function(x)o$par["d"] * (exp(o$par["r"]*(x - o$par["h"])))/(1 + exp(o$par["r"]*(x - o$par["h"]))),
                       aes(col = as.factor(dataTreatment2023ANWM[[i]]$Treatment[1])),
                       size = 1.5))
})

### Get all triple combinations
treatmentsToCompare <- c("1T:1F", "TF-M", "1T:3F", "3T:1F", "T", "F")
triples <- list(c("T", "1T:1F", "F"),
                c("T", "TF-M", "F"),
                c("T", "1T:3F", "F"),
                c("T", "3T:1F", "F"))
names(dataTreatment2023ANWM) <- lapply(dataTreatment2023ANWM, function(l)l$Treatment[1])
names(geomPoints) <- names(dataTreatment2023ANWM)
names(geomFunctions) <- names(dataTreatment2023ANWM)

triplePlotsNW <- lapply(triples, function(triple){
  plotMultipleGraphs(geomPoints[triple], geomFunctions[triple], triple, colourVector[c(1, 3, 2)], XLab = bquote(T[c]), YLab = bquote(p[cc]), Title = "")
})

# Create double comparison plots for the sole crops
doubles <- list(c("T", "T+"), c("F", "F+"))

doublePlotsNW <- lapply(doubles, function(double){
  plotMultipleGraphs(geomPoints[double], geomFunctions[double], double, c(colourVector[2], colourVector[1]), XLab = bquote(T[c]), YLab = bquote(p[cc]), Title = "")
})

# Create a combined plot for intercrops and normal density sole crops
combinedPlotNW <- plotMultipleGraphs(geomPoints[treatmentsToCompare], geomFunctions[treatmentsToCompare], treatmentsToCompare, c(colourVector[1], colourVector[2], colourVector[5], colourVector[4], colourVector[3], colourVector[6]), XLab = bquote(T[c]), YLab = bquote(p[cc]), Title = "")

### Combine data for AIC comparisons
# Fit of 1T:1F
startParameters1T1F <- list(c("r" = 0.008, "h" = 480, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 500, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 550, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75))
combinationAIC1T1F <- round(fitCombinations(data = dataTreatment2023ANW, triple = triples[[1]], model = tLogisticNormal, startParameters = startParameters1T1F), 2)

# Fit of 1T:3F
startParameters1T3F <- list(c("r" = 0.008, "h" = 480, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 500, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 550, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75))
combinationAIC1T3F <- round(fitCombinations(data = dataTreatment2023ANW, triple = triples[[2]], model = tLogisticNormal, startParameters = startParameters1T3F), 2)

# Fit of 3T:1F
startParameters3T1F <- list(c("r" = 0.008, "h" = 480, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 500, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 550, "sd" = 0.1, "d" = 0.75),
                            c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75))
combinationAIC3T1F <- round(fitCombinations(data = dataTreatment2023ANW, triple = triples[[3]], model = tLogisticNormal, startParameters = startParameters3T1F), 2)

# Fit of TF-M
startParametersTFM <- list(c("r" = 0.008, "h" = 480, "sd" = 0.1, "d" = 0.75),
                           c("r" = 0.008, "h" = 500, "sd" = 0.1, "d" = 0.75),
                           c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75),
                           c("r" = 0.008, "h" = 550, "sd" = 0.1, "d" = 0.75),
                           c("r" = 0.008, "h" = 600, "sd" = 0.1, "d" = 0.75))
combinationAICTFM <- round(fitCombinations(data = dataTreatment2023ANW, triple = triples[[4]], model = tLogisticNormal, startParameters = startParametersTFM), 2)

combinationAIC2023ANW <- tibble(Treatment = rep(unique(data2023A$Treatment)[c(1:3, 8)], each = 3),
                                Combination = rep(c("IC & C & L", "IC + C & L", "IC + L & C"), times = 4),
                                AIC = c(combinationAIC1T1F,
                                        combinationAIC1T3F,
                                        combinationAIC3T1F,
                                        combinationAICTFM),
                                Best = c("*",
                                         "",
                                         "",
                                         "*",
                                         "",
                                         "",
                                         "*",
                                         "",
                                         "",
                                         "*",
                                         "",
                                         ""))

print(xtable(combinationAIC2023ANW, type = "latex"), include.rownames = FALSE)

# Print all plots to PDF
pdf("canopy_cover/canopy_cover_NW_tLogisticNormal_graphs_2023A_NoWeeds.pdf")
print(triplePlotsNW)
print(combinedPlotNW)
print(doublePlotsNW)
dev.off()

# Combined plots
jpeg("figures/canopy_cover/plotCanopyCover2023ATripleCombined.jpg", units = "px", width = 900, height = 1000, quality = 600)
triplePlotsNW[[1]] + theme(legend.position = "none", plot.subtitle = element_text(hjust = -0.4)) + labs(y = "", subtitle = bquote({P[cc]}~~~~~~~~~~~"1T:1F")) + 
  ylim(c(0, 1.0)) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_segment(aes(x = 650, y = 0.45, xend = 425, yend = 0.55),
               arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(aes(x = 600, y = 0.30, xend = 450, yend = 0.39),
               arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(aes(x = 550, y = 0.15, xend = 500, yend = 0.20),
               arrow = arrow(length = unit(0.5, "cm"))) +
  annotate("text", x = c(650, 600, 550), y = c(0.45, 0.30, 0.15), label = c("Cereal", "Intercrop", "Legume"), size = 8, hjust = 0.0, col = modelColors[c("Cereal", "Intercrop", "Legume")]) + 
triplePlotsNW[[2]] + theme(legend.position = "none", plot.subtitle = element_text(hjust = -0.4)) + labs(y = "", subtitle = bquote(~~~~~~~~~~~~~~~~~~"TF-M")) + 
  ylim(c(0, 1.0)) + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
triplePlotsNW[[3]] + theme(legend.position = "none", plot.subtitle = element_text(hjust = -0.4)) + labs(y = "", subtitle = bquote({P[cc]}~~~~~~~~~~~"1T:3F")) + 
  ylim(c(0, 1.0)) + 
triplePlotsNW[[4]] + labs(y = "", subtitle = bquote(~~~~~~~~~~~~~~~~~~"3T:1F")) + 
  scale_color_manual(values = modelColors,
  labels = c("Cereal_Legume", "Cereal", "Legume")) + 
  ylim(c(0, 1.0)) + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none", plot.subtitle = element_text(hjust = -0.4)) + 
plot_layout(ncol = 2, nrow = 2, axis_titles = "collect") +
plot_annotation(tag_levels = "a")
dev.off()

##
## Segmented cover
##

Treatments2023A <- read_delim("PlotTreatments2023A.txt", delim = "\t", col_names = TRUE)

dataSegCover2023A <- read_delim("cover_data_2023.txt", delim = "\t") %>% 
  mutate(relativeFaba = faba / all,
         relativeCereal = cereal / all,
         relativeWeed = weed / all,
         ground = 1 - all) %>% 
  mutate(Time = tempData2023$cumulativeTemp[match(time, tempData2023$Time)]) %>% 
  dplyr::select(plot, Time, faba, weed, cereal) %>% 
  rename(Plot = plot,
         Faba = faba,
         Weed = weed,
         Cereal = cereal) %>% 
  mutate(Treatment = str_replace(Treatments2023A$Treatment[match(Plot, Treatments2023A$Plot)], "Bean", "Faba")) %>% 
  dplyr::select(!Plot)

dataSegCover2023ASplit <- split(dataSegCover2023A, dataSegCover2023A$Treatment)

dataSeg2023A1T1F <- dataSegCover2023ASplit[[1]]
dataSeg2023A1T1FLong <- dataSeg2023A1T1F %>% 
  pivot_longer(cols = 2:4, names_to = "Type", values_to = "PCover")

(plotSeg2023A1T1F <- ggplot(data = dataSeg2023A1T1FLong, aes(x = Time, y = PCover, col = Type)) +
    stat_smooth(method = "loess", aes(fill = Type), alpha = 0.2, method.args = c(degree = 2)) +
    geom_point(shape = 1, aes(fill = Type), size = 2) +  # Adjusted to include fill mapping
    scale_color_manual(values = modelColors) +  # Set line colors
    scale_fill_manual(values = modelColors) +   # Set fill colors
    theme_classic(base_size = 40) +
    labs(title = "Segmented canopy cover",
         subtitle = "2023A, 1T:1F",
         x = bquote(T[c]),
         y = bquote(p[cc]),
         col = "Type"))

dataSeg2023A1T3F <- dataSegCover2023ASplit[[2]]
dataSeg2023A1T3FLong <- dataSeg2023A1T3F %>% 
  pivot_longer(cols = 2:4, names_to = "Type", values_to = "PCover")

(plotSeg2023A1T3F <- ggplot(data = dataSeg2023A1T3FLong, aes(x = Time, y = PCover, col = Type)) +
    stat_smooth(method = "loess", aes(fill = Type), alpha = 0.2, method.args = c(degree = 2)) +
    geom_point(shape = 1, aes(fill = Type), size = 2) +  # Adjusted to include fill mapping
    scale_color_manual(values = modelColors) +  # Set line colors
    scale_fill_manual(values = modelColors) +   # Set fill colors
    theme_classic(base_size = 40) +
    labs(title = "Segmented canopy cover",
         subtitle = "2023A, 1T:3F",
         x = bquote(T[c]),
         y = bquote(p[cc]),
         col = "Type"))

dataSeg2023A3T1F <- dataSegCover2023ASplit[[3]]
dataSeg2023A3T1FLong <- dataSeg2023A3T1F %>% 
  pivot_longer(cols = 2:4, names_to = "Type", values_to = "PCover")

(plotSeg2023A3T1F <- ggplot(data = dataSeg2023A3T1FLong, aes(x = Time, y = PCover, col = Type)) +
    stat_smooth(method = "loess", aes(fill = Type), alpha = 0.2, method.args = c(degree = 2)) +
    geom_point(shape = 1, aes(fill = Type), size = 2) +  # Adjusted to include fill mapping
    scale_color_manual(values = modelColors) +  # Set line colors
    scale_fill_manual(values = modelColors) +   # Set fill colors
    theme_classic(base_size = 40) +
    labs(title = "Segmented canopy cover",
         subtitle = "2023A, 3T:1F",
         x = bquote(T[c]),
         y = bquote(p[cc]),
         col = "Type"))

dataSeg2023AF <- dataSegCover2023ASplit[[4]]
dataSeg2023AFLong <- dataSeg2023AF %>% 
  pivot_longer(cols = 2:4, names_to = "Type", values_to = "PCover")

(plotSeg2023AF <- ggplot(data = dataSeg2023AFLong, aes(x = Time, y = PCover, col = Type)) +
    stat_smooth(method = "loess", aes(fill = Type), alpha = 0.2, method.args = c(degree = 2)) +
    geom_point(shape = 1, aes(fill = Type), size = 2) +  # Adjusted to include fill mapping
    scale_color_manual(values = modelColors) +  # Set line colors
    scale_fill_manual(values = modelColors) +   # Set fill colors
    theme_classic(base_size = 40) +
    labs(title = "Segmented canopy cover",
         subtitle = "2023A, F",
         x = bquote(T[c]),
         y = bquote(p[cc]),
         col = "Type"))

dataSeg2023AFF <- dataSegCover2023ASplit[[5]]
dataSeg2023AFFLong <- dataSeg2023AFF %>% 
  pivot_longer(cols = 2:4, names_to = "Type", values_to = "PCover")

(plotSeg2023AFF <- ggplot(data = dataSeg2023AFFLong, aes(x = Time, y = PCover, col = Type)) +
    stat_smooth(method = "loess", aes(fill = Type), alpha = 0.2, method.args = c(degree = 2)) +
    geom_point(shape = 1, aes(fill = Type), size = 2) +  # Adjusted to include fill mapping
    scale_color_manual(values = modelColors) +  # Set line colors
    scale_fill_manual(values = modelColors) +   # Set fill colors
    theme_classic(base_size = 40) +
    labs(title = "Segmented canopy cover",
         subtitle = "2023A, F+",
         x = bquote(T[c]),
         y = bquote(p[cc]),
         col = "Type"))

dataSeg2023AT <- dataSegCover2023ASplit[[6]]
dataSeg2023ATLong <- dataSeg2023AT %>% 
  pivot_longer(cols = 2:4, names_to = "Type", values_to = "PCover")

(plotSeg2023AT <- ggplot(data = dataSeg2023ATLong, aes(x = Time, y = PCover, col = Type)) +
    stat_smooth(method = "loess", aes(fill = Type), alpha = 0.2, method.args = c(degree = 2)) +
    geom_point(shape = 1, aes(fill = Type), size = 2) +  # Adjusted to include fill mapping
    scale_color_manual(values = modelColors) +  # Set line colors
    scale_fill_manual(values = modelColors) +   # Set fill colors
    theme_classic(base_size = 40) +
    labs(title = "Segmented canopy cover",
         subtitle = "2023A, T",
         x = bquote(T[c]),
         y = bquote(p[cc]),
         col = "Type"))

dataSeg2023ATT <- dataSegCover2023ASplit[[7]]
dataSeg2023ATTLong <- dataSeg2023ATT %>% 
  pivot_longer(cols = 2:4, names_to = "Type", values_to = "PCover")

(plotSeg2023ATT <- ggplot(data = dataSeg2023ATTLong, aes(x = Time, y = PCover, col = Type)) +
    stat_smooth(method = "loess", aes(fill = Type), alpha = 0.2, method.args = c(degree = 2)) +
    geom_point(shape = 1, aes(fill = Type), size = 2) +  # Adjusted to include fill mapping
    scale_color_manual(values = modelColors) +  # Set line colors
    scale_fill_manual(values = modelColors) +   # Set fill colors
    theme_classic(base_size = 40) +
    labs(title = "Segmented canopy cover",
         subtitle = "2023A, T+",
         x = bquote(T[c]),
         y = bquote(p[cc]),
         col = "Type"))

dataSeg2023ATFM <- dataSegCover2023ASplit[[8]]
dataSeg2023ATFMLong <- dataSeg2023ATFM %>% 
  pivot_longer(cols = 2:4, names_to = "Type", values_to = "PCover")

(plotSeg2023ATFM <- ggplot(data = dataSeg2023ATFMLong, aes(x = Time, y = PCover, col = Type)) +
    stat_smooth(method = "loess", aes(fill = Type), alpha = 0.2, method.args = c(degree = 2)) +
    geom_point(shape = 1, aes(fill = Type), size = 2) +  # Adjusted to include fill mapping
    scale_color_manual(values = modelColors) +  # Set line colors
    scale_fill_manual(values = modelColors) +   # Set fill colors
    theme_classic(base_size = 40) +
    labs(title = "Segmented canopy cover",
         subtitle = "2023A, TF-M",
         x = bquote(T[c]),
         y = bquote(p[cc]),
         col = "Type"))


########################
##                    ##
## Experiment 2023 B  ##
##                    ##
########################

sowingDate2023 <- as.Date("2023-03-02")

# Read data
data2023B <- read_xlsx("WSF_2023_data.xlsx", sheet = "CanopyCoverB", range = "A1:G793", col_names = TRUE) %>% 
  dplyr::select(Treatment, Weeds, Date, PCover) %>% 
  mutate(Date = as.character(Date)) %>% 
  mutate(Time = as.integer(as.Date(as.character(Date)) - sowingDate2023)) %>% 
  rename(Val = PCover) %>% 
  dplyr::select(-Date) %>% 
  filter(!(Treatment %in% c("T-25", "F-25", "1T:1F-25", "TF-M-25")))

tempData2023 <- read_delim("TempData2023.txt", delim = "\t") %>% 
  mutate(Time = as.integer(as.Date(Date) - sowingDate2023))

cumulativeTemp2023 <- ((tempData2023$TMin + tempData2023$TMax) / 2)
cumulativeTemp2023[which(cumulativeTemp2023 < 0.0)] <- 0
cumulativeTemp2023 <- cumsum(cumulativeTemp2023)

tempData2023 <- tempData2023 %>% 
  mutate(cumulativeTemp = cumulativeTemp2023)

data2023B <- data2023B %>% 
  mutate(Treatment = factor(Treatment, levels = c("F", "F-375", "TF-M", "TF-M-375", "T", "T-375", "1T:1F", "1T:1F-375"))) %>% 
  arrange(Treatment)

dataWeeds2023B <- data2023B %>% 
  filter(Weeds == "Y") %>% 
  dplyr::select(-Weeds) %>% 
  filter(Time != 91) %>%
  mutate(Time = tempData2023$cumulativeTemp[match(Time, tempData2023$Time)]) %>% 
  dplyr::select(Treatment, Time, Val)
dataNoWeeds2023B <- data2023B %>% 
  filter(Weeds == "N") %>% 
  dplyr::select(-Weeds) %>% 
  filter(Time != 91) %>%
  mutate(Time = tempData2023$cumulativeTemp[match(Time, tempData2023$Time)]) %>% 
  dplyr::select(Treatment, Time, Val)

## Manual fits

##################
##################
####          ####
#### No Weeds ####
####          ####
##################
##################

dataTreatment2023BNW <- split(dataNoWeeds2023B, dataNoWeeds2023B$Treatment)
dataTreatment2023BNW <- lapply(dataTreatment2023BNW, function(dat){
  datNew <- dat %>% 
    dplyr::select(Time, Val)
  return(datNew)
})


##################################
###                            ###
### Fit logistic normal curves ###
###                            ###
##################################     

fitsNWLogisticNormal2023B <- rep(list(NA), 8)

# Fit of F
fitsNWLogisticNormal2023B[[1]] <- fitDataToModel(data = dataTreatment2023BNW[[1]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 700, "sd" = 0.1))
plotFit(dataTreatment2023BNW[[1]], fitsNWLogisticNormal2023B[[1]]$par, "logistic")

# Fit of F-375
fitsNWLogisticNormal2023B[[2]] <- fitDataToModel(data = dataTreatment2023BNW[[2]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 700, "sd" = 0.1))
plotFit(dataTreatment2023BNW[[2]], fitsNWLogisticNormal2023B[[2]]$par, "logistic")

# Fit of TF-M
fitsNWLogisticNormal2023B[[3]] <- fitDataToModel(data = dataTreatment2023BNW[[3]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 550, "sd" = 0.1))
plotFit(dataTreatment2023BNW[[3]], fitsNWLogisticNormal2023B[[3]]$par, "logistic")

# Fit of TF-M-375
fitsNWLogisticNormal2023B[[4]] <- fitDataToModel(data = dataTreatment2023BNW[[4]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 575, "sd" = 0.1))
plotFit(dataTreatment2023BNW[[4]], fitsNWLogisticNormal2023B[[4]]$par, "logistic")

# Fit of T
fitsNWLogisticNormal2023B[[5]] <- fitDataToModel(data = dataTreatment2023BNW[[5]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 450, "sd" = 0.1))
plotFit(dataTreatment2023BNW[[5]], fitsNWLogisticNormal2023B[[5]]$par, "logistic")

# Fit T-375
fitsNWLogisticNormal2023B[[6]] <- fitDataToModel(data = dataTreatment2023BNW[[6]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 550, "sd" = 0.1))
plotFit(dataTreatment2023BNW[[6]], fitsNWLogisticNormal2023B[[6]]$par, "logistic")

# Fit 1T:1F
fitsNWLogisticNormal2023B[[7]] <- fitDataToModel(data = dataTreatment2023BNW[[7]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 550, "sd" = 0.1))
plotFit(dataTreatment2023BNW[[7]], fitsNWLogisticNormal2023B[[7]]$par, "logistic")

# Fit 1T:1F-375
fitsNWLogisticNormal2023B[[8]] <- fitDataToModel(data = dataTreatment2023BNW[[8]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 600, "sd" = 0.1))
plotFit(dataTreatment2023BNW[[8]], fitsNWLogisticNormal2023B[[8]]$par, "logistic")


AICLogisticNormal2023BNW <- sum(unlist(lapply(fitsNWLogisticNormal2023B, function(fits)2 * 3 + 2 * fits$value)))


##############################################
###                                        ###
### Fit transformed logistic normal curves ###
###                                        ###
##############################################      

fitsNWTLogisticNormal2023B <- rep(list(NA), 8)

# Fit of F
fitsNWTLogisticNormal2023B[[1]] <- fitDataToModel(data = dataTreatment2023BNW[[1]], model = tLogisticNormal, startParameters = c("r" = 0.005, "h" = 700, "sd" = 0.1, "d" = 0.95), dMax = 1.0)
plotFit(dataTreatment2023BNW[[1]], fitsNWTLogisticNormal2023B[[1]]$par, "tLogistic")

# Fit of F-375
fitsNWTLogisticNormal2023B[[2]] <- fitDataToModel(data = dataTreatment2023BNW[[2]], model = tLogisticNormal, startParameters = c("r" = 0.03, "h" = 700, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[2]], fitsNWTLogisticNormal2023B[[2]]$par, "tLogistic")

# Fit of TF-M
fitsNWTLogisticNormal2023B[[3]] <- fitDataToModel(data = dataTreatment2023BNW[[3]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 550, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[3]], fitsNWTLogisticNormal2023B[[3]]$par, "tLogistic")

# Fit of TF-M-375
fitsNWTLogisticNormal2023B[[4]] <- fitDataToModel(data = dataTreatment2023BNW[[4]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 575, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[4]], fitsNWTLogisticNormal2023B[[4]]$par, "tLogistic")

# Fit of T
fitsNWTLogisticNormal2023B[[5]] <- fitDataToModel(data = dataTreatment2023BNW[[5]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 450, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[5]], fitsNWTLogisticNormal2023B[[5]]$par, "tLogistic")

# Fit T-375
fitsNWTLogisticNormal2023B[[6]] <- fitDataToModel(data = dataTreatment2023BNW[[6]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 550, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[6]], fitsNWTLogisticNormal2023B[[6]]$par, "tLogistic")

# Fit 1T:1F
fitsNWTLogisticNormal2023B[[7]] <- fitDataToModel(data = dataTreatment2023BNW[[7]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 550, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[7]], fitsNWTLogisticNormal2023B[[7]]$par, "tLogistic")

# Fit 1T:1F-375
fitsNWTLogisticNormal2023B[[8]] <- fitDataToModel(data = dataTreatment2023BNW[[8]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 600, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[8]], fitsNWTLogisticNormal2023B[[8]]$par, "tLogistic")

AICTLogisticNormal2023BNW <- sum(unlist(lapply(fitsNWTLogisticNormal2023B, function(fits)2 * 4 + 2 * fits$value)))


#################################
###                           ###
### Fit logistic gamma curves ###
###                           ###
#################################     

fitsNWLogisticGamma2023B <- rep(list(NA), 8)

# Fit of F
fitsNWLogisticGamma2023B[[1]] <- fitDataToModel(data = dataTreatment2023BNW[[1]], model = logisticGamma, startParameters = c("r" = 0.03, "h" = 700, "shape" = 1.0))
plotFit(dataTreatment2023BNW[[1]], fitsNWLogisticGamma2023B[[1]]$par, "logistic")

# Fit of F-375
fitsNWLogisticGamma2023B[[2]] <- fitDataToModel(data = dataTreatment2023BNW[[2]], model = logisticGamma, startParameters = c("r" = 0.02, "h" = 700, "shape" = 1.0))
plotFit(dataTreatment2023BNW[[2]], fitsNWLogisticGamma2023B[[2]]$par, "logistic")

# Fit of TF-M
fitsNWLogisticGamma2023B[[3]] <- fitDataToModel(data = dataTreatment2023BNW[[3]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 550, "shape" = 1.0))
plotFit(dataTreatment2023BNW[[3]], fitsNWLogisticGamma2023B[[3]]$par, "logistic")

# Fit of TF-M-375
fitsNWLogisticGamma2023B[[4]] <- fitDataToModel(data = dataTreatment2023BNW[[4]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 575, "shape" = 1.0))
plotFit(dataTreatment2023BNW[[4]], fitsNWLogisticGamma2023B[[4]]$par, "logistic")

# Fit of T
fitsNWLogisticGamma2023B[[5]] <- fitDataToModel(data = dataTreatment2023BNW[[5]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 450, "shape" = 1.0))
plotFit(dataTreatment2023BNW[[5]], fitsNWLogisticGamma2023B[[5]]$par, "logistic")

# Fit T-375
fitsNWLogisticGamma2023B[[6]] <- fitDataToModel(data = dataTreatment2023BNW[[6]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 550, "shape" = 1.0))
plotFit(dataTreatment2023BNW[[6]], fitsNWLogisticGamma2023B[[6]]$par, "logistic")

# Fit 1T:1F
fitsNWLogisticGamma2023B[[7]] <- fitDataToModel(data = dataTreatment2023BNW[[7]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 550, "shape" = 1.0))
plotFit(dataTreatment2023BNW[[7]], fitsNWLogisticGamma2023B[[7]]$par, "logistic")

# Fit 1T:1F-375
fitsNWLogisticGamma2023B[[8]] <- fitDataToModel(data = dataTreatment2023BNW[[8]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 600, "shape" = 1.0))
plotFit(dataTreatment2023BNW[[8]], fitsNWLogisticGamma2023B[[8]]$par, "logistic")

AICLogisticGamma2023BNW <- sum(unlist(lapply(fitsNWLogisticGamma2023B, function(fits)2 * 3 + 2 * fits$value)))


#############################################
###                                       ###
### Fit transformed logistic gamma curves ###
###                                       ###
#############################################      

fitsNWTLogisticGamma2023B <- rep(list(NA), 8)

# Fit of F
fitsNWTLogisticGamma2023B[[1]] <- fitDataToModel(data = dataTreatment2023BNW[[1]], model = tLogisticGamma, startParameters = c("r" = 0.02, "h" = 700, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[1]], fitsNWTLogisticGamma2023B[[1]]$par, "tLogistic")

# Fit of F-375
fitsNWTLogisticGamma2023B[[2]] <- fitDataToModel(data = dataTreatment2023BNW[[2]], model = tLogisticGamma, startParameters = c("r" = 0.04, "h" = 700, "shape" = 1.0, "d" = 0.90), dMax = 1.0)
plotFit(dataTreatment2023BNW[[2]], fitsNWTLogisticGamma2023B[[2]]$par, "tLogistic")

# Fit of TF-M
fitsNWTLogisticGamma2023B[[3]] <- fitDataToModel(data = dataTreatment2023BNW[[3]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 550, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[3]], fitsNWTLogisticGamma2023B[[3]]$par, "tLogistic")

# Fit of TF-M-375
fitsNWTLogisticGamma2023B[[4]] <- fitDataToModel(data = dataTreatment2023BNW[[4]], model = tLogisticGamma, startParameters = c("r" = 0.02, "h" = 575, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[4]], fitsNWTLogisticGamma2023B[[4]]$par, "tLogistic")

# Fit of T
fitsNWTLogisticGamma2023B[[5]] <- fitDataToModel(data = dataTreatment2023BNW[[5]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 450, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[5]], fitsNWTLogisticGamma2023B[[5]]$par, "tLogistic")

# Fit T-375
fitsNWTLogisticGamma2023B[[6]] <- fitDataToModel(data = dataTreatment2023BNW[[6]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 550, "shape" = 1.0, "d" = 0.75), dMax = 1.0)
plotFit(dataTreatment2023BNW[[6]], fitsNWTLogisticGamma2023B[[6]]$par, "tLogistic")

# Fit 1T:1F
fitsNWTLogisticGamma2023B[[7]] <- fitDataToModel(data = dataTreatment2023BNW[[7]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 550, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[7]], fitsNWTLogisticGamma2023B[[7]]$par, "tLogistic")

# Fit 1T:1F-375
fitsNWTLogisticGamma2023B[[8]] <- fitDataToModel(data = dataTreatment2023BNW[[8]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 600, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2023BNW[[8]], fitsNWTLogisticGamma2023B[[8]]$par, "tLogistic")


AICTLogisticGamma2023BNW <- sum(unlist(lapply(fitsNWTLogisticGamma2023B, function(fits)2 * 4 + 2 * fits$value)))

AICLogisticNormal2023BNW
AICTLogisticNormal2023BNW
AICLogisticGamma2023BNW
AICTLogisticGamma2023BNW

# Logistic gamma is the best fitting model

# Calculate 95% CIs
parameters2023BNW <- names(fitsNWLogisticGamma2023B[[1]]$par)
parameterCIs2023BNW <- rep(NA, times = 3 * length(parameters2023BNW))
names(parameterCIs2023BNW) <- unlist(lapply(parameters2023BNW, function(par)c(par, paste(par, "Lower", sep = ""), paste(par, "Upper", sep = ""))))
parameterCIs2023BNW <- rep(list(parameterCIs2023BNW), length(fitsNWLogisticGamma2023B))

cI2023BNW <- function(par, x, z, opt){
  pars <- par
  pars[names(opt)] <- opt
  nll <- logisticGamma(pars, x, z)
  return(nll)
}

for(i in 1:length(fitsNWLogisticGamma2023B)){
  parameterVecs <- lapply(fitsNWLogisticGamma2023B, function(fit){
    parStart <- fit$par / 100
    parEnd <- fit$par * 20
    parVec <- list(seq(parStart[1], parEnd[1], length.out = 400), seq(parStart[2], parEnd[2], length.out = 400), seq(parStart[3], parEnd[3], length.out = 400))
    names(parVec) <- names(fit$par)
    return(parVec)
  })
  optPar <- fitsNWLogisticGamma2023B[[i]]$par
  parameterCIs2023BNW[[i]][names(optPar)] <- optPar
  parameterCIs2023BNW[[i]][unlist((lapply(names(optPar), function(n)c(paste(n, "Lower", sep = ""), paste(n, "Upper", sep = "")))))] <-
    unlist(lapply(1:length(parameterVecs[[i]]), function(j){
      t(tibble(calculate95CI(parStartVec = optPar[-which(names(optPar) == names(fitsNWLogisticGamma2023B[[i]]$par)[[j]])],
                             optVec = parameterVecs[[i]][j],
                             model = cI2023BNW,
                             optNLL = fitsNWLogisticGamma2023B[[i]]$value,
                             x = dataTreatment2023BNW[[i]][,1][[1]],
                             z = dataTreatment2023BNW[[i]][,2][[1]])))
    }))
  
}
names(parameterCIs2023BNW) <- names(dataTreatment2023BNW)

# Write to table
CIOut2023BNW <- tibble(Treatment = names(parameterCIs2023BNW)) %>% 
  bind_cols(bind_rows(parameterCIs2023BNW))
write.table(CIOut2023BNW, "canopy_cover/canopy_cover_logisticGamma_parameters_2023BNW.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

CIOut2023BNWTable <- CIOut2023BNW %>% 
  mutate(r = paste0(round(r, 4), " (", round(rLower, 4), "-", round(rUpper, 4), ")"),
         h = paste0(round(h, 2), " (", round(hLower, 2), "-", round(hUpper, 2), ")"),
         shape = paste0(round(shape, 2), " (", round(shapeLower, 2), "-", round(shapeUpper, 2), ")")) %>% 
  dplyr::select(Treatment, r, h, shape)

print(xtable(CIOut2023BNWTable, type = "latex"), include.rownames = FALSE)


dataTreatment2023BNWM <- lapply(1:length(dataTreatment2023BNW), function(i){
  dat <- dataTreatment2023BNW[[i]]
  treatment <- names(dataTreatment2023BNW)[i]
  dat <- dat %>% 
    mutate(Treatment = treatment)
})

dataTreatment2023BNWM <- lapply(dataTreatment2023BNWM, function(Data){
  Average <- aggregate(Data$Val, list(Data$Time), function(x)mean(x, na.rm = TRUE))
  Data$Average <- Average$x[match(Data$Time, Average$Group.1)]
  return(Data)
})

# Plot triple combination graphs
geomPoints2023BNW <- lapply(dataTreatment2023BNWM, function(d){
  return(geom_point(data = d, aes(x = Time, y = Average, col = as.factor(Treatment[1])), size = 2, shape = 1))
})

geomFunctions2023BNW <- lapply(1:length(fitsNWLogisticGamma2023B), function(i){
  o <- fitsNWLogisticGamma2023B[[i]]
  return(geom_function(fun = function(x)(exp(o$par["r"]*(x - o$par["h"])))/(1 + exp(o$par["r"]*(x - o$par["h"]))),
                       aes(col = as.factor(dataTreatment2023BNWM[[i]]$Treatment[1])),
                       size = 1.5))
})

values2023BNW <- colourVector[1:length(fitsNWLogisticGamma2023B)]
names(values2023BNW) <- unname(sapply(dataTreatment2023BNWM, function(d)d$Treatment[1]))
labels2023BNW <- unname(sapply(dataTreatment2023BNWM, function(d)d$Treatment[1]))

plots2023BNW <- lapply(1:length(geomPoints2023BNW), function(i){
  plot <- ggplot() +
    geomPoints2023BNW[[i]] +
    geomFunctions2023BNW[[i]] +
    theme_classic(base_size = 40) +
    scale_color_manual(labels = labels2023BNW[i],
                       values = "#7FC97F") +
    labs(title = "",
         y = bquote(p[cc]),
         x = bquote(T[c]),
         color = "Treatment") 
})

# Create a combined plot for intercrops and normal density sole crops
combinedPlot2023BNW <- plotMultipleGraphs(geomPoints2023BNW, geomFunctions2023BNW, names(dataTreatment2023BNW), c(colourVector[c(2, 9, 8, 4, 1, 5, 3, 7)]), XLab = bquote(T[c]), YLab = bquote(p[cc]), Title = "")

combinedFaba2023BNW <- plotMultipleGraphs(geomPoints2023BNW[c(1, 2)], geomFunctions2023BNW[c(1, 2)], c("F", "F-375"), c(colourVector[c(2, 9)]), XLab = bquote(T[c]), YLab = bquote(p[cc]), Title = "")
combinedTriticale2023BNW <- plotMultipleGraphs(geomPoints2023BNW[c(5, 6)], geomFunctions2023BNW[c(5, 6)], c("T", "T-375"), c(colourVector[c(1, 5)]), XLab = bquote(T[c]), YLab = bquote(p[cc]), Title = "")
combinedIC2023BNW <- plotMultipleGraphs(geomPoints2023BNW[c(3, 4, 7, 8)], geomFunctions2023BNW[c(3, 4, 7, 8)], c("TF-M", "TF-M-375", "1T:1F", "1T:1F-375"), c(colourVector[c(8, 4, 3, 7)]), XLab = bquote(T[c]), YLab = bquote(p[cc]), Title = "")
combined1252023BNW <- plotMultipleGraphs(geomPoints2023BNW[c(1, 3, 5, 7)], geomFunctions2023BNW[c(1, 3, 5, 7)], c("F", "TF-M", "T", "1T:1F"), c(colourVector[c(2, 8, 1, 3)]), XLab = bquote(T[c]), YLab = bquote(p[cc]), Title = "")
combined3752023BNW <- plotMultipleGraphs(geomPoints2023BNW[c(2, 4, 6, 8)], geomFunctions2023BNW[c(2, 4, 6, 8)], c("F-375", "TF-M-375", "T-375", "1T:1F-375"), c(colourVector[c(9, 4, 5, 7)]), XLab = bquote(T[c]), YLab = bquote(p[cc]), Title = "")

### Combine data for AIC comparisons
# Fit of 1T:1F
startParameters1T1F <- list(c("r" = 0.008, "h" = 520, "shape" = 1.0),
                             c("r" = 0.008, "h" = 500, "shape" = 1.0),
                             c("r" = 0.008, "h" = 600, "shape" = 1.0),
                             c("r" = 0.008, "h" = 550, "shape" = 1.0),
                             c("r" = 0.008, "h" = 600, "shape" = 1.0))
combinationAIC1T1F <- round(fitCombinations(data = dataTreatment2023BNW, triple = c("1T:1F", "T", "F"), model = logisticGamma, startParameters = startParameters1T1F, dMax = NA), 2)

# Fit of 1T:1F-375
startParameters1T1F375 <- list(c("r" = 0.008, "h" = 480, "shape" = 1.0),
                             c("r" = 0.008, "h" = 500, "shape" = 1.0),
                             c("r" = 0.008, "h" = 600, "shape" = 1.0),
                             c("r" = 0.008, "h" = 550, "shape" = 1.0),
                             c("r" = 0.008, "h" = 600, "shape" = 1.0))
combinationAIC1T1F375 <- round(fitCombinations(data = dataTreatment2023BNW, triple = c("1T:1F-375", "T-375", "F-375"), model = logisticGamma, startParameters = startParameters1T1F375, dMax = NA), 2)

# Fit of TF-M
startParametersTFM <- list(c("r" = 0.008, "h" = 520, "shape" = 1.0),
                            c("r" = 0.008, "h" = 500, "shape" = 1.0),
                            c("r" = 0.008, "h" = 600, "shape" = 1.0),
                            c("r" = 0.008, "h" = 550, "shape" = 1.0),
                            c("r" = 0.008, "h" = 600, "shape" = 1.0))
combinationAICTFM <- round(fitCombinations(data = dataTreatment2023BNW, triple = c("TF-M", "T", "F"), model = logisticGamma, startParameters = startParametersTFM, dMax = NA), 2)

# Fit of TF-M-375
startParametersTFM375 <- list(c("r" = 0.008, "h" = 480, "shape" = 1.0),
                            c("r" = 0.008, "h" = 500, "shape" = 1.0),
                            c("r" = 0.008, "h" = 600, "shape" = 1.0),
                            c("r" = 0.008, "h" = 550, "shape" = 1.0),
                            c("r" = 0.008, "h" = 600, "shape" = 1.0))
combinationAICTFM375 <- round(fitCombinations(data = dataTreatment2023BNW, triple = c("TF-M-375", "T-375", "F-375"), model = logisticGamma, startParameters = startParametersTFM375, dMax = NA), 2)

combinationAIC2023BNW <- tibble(Treatment = rep(unique(data2023B$Treatment)[c(1:3, 8)], each = 3),
                                Combination = rep(c("IC & C & L", "IC + C & L", "IC + L & C"), times = 4),
                                AIC = c(combinationAIC1T1F,
                                        combinationAIC1T1F375,
                                        combinationAICTFM,
                                        combinationAICTFM375),
                                Best = c("*",
                                         "*",
                                         "",
                                         "*",
                                         "",
                                         "",
                                         "*",
                                         "",
                                         "",
                                         "*",
                                         "",
                                         ""))

print(xtable(combinationAIC2023BNW, type = "latex"), include.rownames = FALSE)



# Print all plots to PDF
pdf("canopy_cover/canopy_cover_logisticGamma_graphs_2023BNW.pdf")
print(plots2023BNW)
combinedPlot2023BNW
combinedFaba2023BNW
combinedTriticale2023BNW
combinedIC2023BNW
combined1252023BNW
combined3752023BNW
dev.off()

jpeg("figures/canopy_cover/plotCanopyCover2023BCombinedNW.jpg", units = "px", width = 1250, height = 2000, quality = 600)
combinedFaba2023BNW + labs(y = "", subtitle = bquote({P[cc]})) + 
  ylim(c(0, 1.0)) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none", plot.subtitle = element_text(hjust = -0.05)) +
  geom_segment(aes(x = 650, y = 0.75, xend = 730, yend = 0.58),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1.0) +
  geom_segment(aes(x = 650, y = 0.2, xend = 620, yend = 0.31),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1.0) +
  annotate("text", x = c(625, 650), y = c(0.75, 0.20), label = c("F", "F-375"), size = 10, hjust = 0.0, col = modelColors[c("F", "TM")]) + 
combinedTriticale2023BNW + labs(y = "", subtitle = "", x = "") + 
  ylim(c(0, 1.0)) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none") + 
  geom_segment(aes(x = 450, y = 0.9, xend = 540, yend = 0.78),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1.0) +
  geom_segment(aes(x = 650, y = 0.5, xend = 600, yend = 0.6),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1.0) +
  annotate("text", x = c(425, 650), y = c(0.9, 0.5), label = c("T", "T-375"), size = 10, hjust = 0.0, col = modelColors[c("T", "TA")]) + 
combinedIC2023BNW + labs(y = "", subtitle = bquote({P[cc]})) + 
  ylim(c(0, 1.0)) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none", plot.subtitle = element_text(hjust = -0.05)) + 
  geom_segment(aes(x = 500, y = 0.85, xend = 625, yend = 0.75),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1.0) +
  geom_segment(aes(x = 500, y = 0.75, xend = 600, yend = 0.65),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1.0) +
  geom_segment(aes(x = 600, y = 0.30, xend = 550, yend = 0.40),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1.0) +
  geom_segment(aes(x = 550, y = 0.15, xend = 475, yend = 0.30),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1.0) +
  annotate("text", x = c(400, 400, 600, 550), y = c(0.85, 0.75, 0.30, 0.15), label = c("TF-M", "1T:1F", "TF-M-375", "1T:1F-375"), size = 10, hjust = 0.0, col = modelColors[c("FM", "Intercrop", "Weed", "FA")]) + 
combined1252023BNW + labs(y = "", subtitle = "") + 
  ylim(c(0, 1.0)) + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none") + 
  geom_segment(aes(x = 475, y = 0.95, xend = 580, yend = 0.85),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1.0) +
  geom_segment(aes(x = 400, y = 0.75, xend = 575, yend = 0.65),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1.0) +
  geom_segment(aes(x = 675, y = 0.30, xend = 525, yend = 0.50),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1.0) +
  geom_segment(aes(x = 600, y = 0.15, xend = 550, yend = 0.24),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1.0) +
  annotate("text", x = c(450, 320, 675, 600), y = c(0.95, 0.75, 0.30, 0.15), label = c("T", "TF-M", "1T:1F", "F"), size = 10, hjust = 0.0, col = modelColors[c("T", "FM", "Intercrop", "F")]) + 
combined3752023BNW + labs(y = "", subtitle = bquote({P[cc]})) + ylim(c(0, 1.0)) +
  theme(legend.position = "none",
        plot.subtitle = element_text(hjust = -0.05)) +
  geom_segment(aes(x = 500, y = 0.85, xend = 650, yend = 0.70),
             arrow = arrow(length = unit(0.5, "cm")),
             size = 1.0) +
  geom_segment(aes(x = 500, y = 0.70, xend = 623, yend = 0.56),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1.0) +
  geom_segment(aes(x = 675, y = 0.30, xend = 590, yend = 0.47),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1.0) +
  geom_segment(aes(x = 550, y = 0.10, xend = 525, yend = 0.18),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1.0) +
  annotate("text", x = c(400, 325, 675, 550), y = c(0.85, 0.70, 0.30, 0.10), label = c("T-375", "1T:1F-375", "TF-M-375", "F-375"), size = 10, hjust = 0.0, col = modelColors[c("TA", "FA", "Weed", "TM")]) + 
  plot_layout(design = "
  12
  34
  5#
", axis_titles = "collect") +
plot_annotation(tag_levels = "a")
dev.off()


jpeg("figures/canopy_cover/plotSegCombined.jpg", units = "px", width = 2000, height = 3200, quality = 600)
plotSeg2022BarleyFaba + labs(title = "", subtitle = bquote({P[cc]}~~~~~~~~~~~"Barley-Faba"), y = "", x = bquote(T[c])) + 
  ylim(c(0, 1.0)) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), legend.position = "none",
        plot.subtitle = element_text(hjust = -0.1)) + 
  geom_segment(aes(x = 525, y = 0.92, xend = 560, yend = 0.82),
               arrow = arrow(length = unit(1.0, "cm")),
               size = 2,
               col = "black") +
  geom_segment(aes(x = 500, y = 0.28, xend = 525, yend = 0.17),
               arrow = arrow(length = unit(1.0, "cm")),
               size = 2,
               col = "black") +
  geom_segment(aes(x = 650, y = 0.08, xend = 600, yend = 0.02),
               arrow = arrow(length = unit(1.0, "cm")),
               size = 2,
               col = "black") +
  annotate("text", x = c(450, 450, 650), y = c(0.97, 0.34, 0.10), label = c("Cereal", "Faba", "Weed"), size = 16, hjust = 0.0, col = modelColors[c("Cereal", "Faba", "Weed")]) + 
plotSeg2022RyeFaba + labs(title = "", subtitle = bquote(~~~~~~~~~~~~~~~~~~"Rye-Faba"), y = "", x = bquote(T[c])) + 
  ylim(c(0, 1.0)) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), legend.position = "none",
        plot.subtitle = element_text(hjust = -0.1)) + 
plotSeg2022TriticaleFaba + labs(title = "", subtitle = bquote({P[cc]}~~~~~~~~~~~"Triticale-Faba"), y = "", x = bquote(T[cc])) + 
  ylim(c(0, 1.0)) + 
  theme(legend.position = "none",
        plot.subtitle = element_text(hjust = -0.1)) + 
plotSeg2022WheatFaba + labs(title = "", subtitle = bquote(~~~~~~~~~~~~~~~~~~"Wheat-Faba"), y = "", x = bquote(T[c])) + 
  ylim(c(0, 1.0)) + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), legend.position = "none",
        plot.subtitle = element_text(hjust = -0.1)) + 
plotSeg2023A1T1F + labs(title = "", subtitle = bquote({P[cc]}~~~~~~~~~~~"1T:1F"), y = "", x = bquote(T[c])) + 
  ylim(c(0, 1.0)) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), legend.position = "none",
        plot.subtitle = element_text(hjust = -0.1)) + 
plotSeg2023ATFM + labs(title = "", subtitle = bquote(~~~~~~~~~~~~~~~~~~"TF-M"), y = "", x = bquote(T[c])) + 
  ylim(c(0, 1.0)) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), legend.position = "none",
        plot.subtitle = element_text(hjust = -0.1)) + 
plotSeg2023A1T3F + labs(title = "", subtitle = bquote({P[cc]}~~~~~~~~~~~"1T:3F"), y = "", x = bquote(T[c])) + 
  ylim(c(0, 1.0)) + theme(legend.position = "none",
                          plot.subtitle = element_text(hjust = -0.1)) +
plotSeg2023A3T1F + labs(title = "", subtitle = bquote(~~~~~~~~~~~~~~~~~~"3T:1F"), y = "", x = bquote(T[c])) + 
  ylim(c(0, 1.0)) + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.subtitle = element_text(hjust = -0.1),
        legend.position = "none") + 
plot_layout(nrow = 4, ncol = 2, axis_titles = "collect") + 
plot_annotation(tag_levels = "a")
dev.off()


######################
##                  ##
## Experiment 2024  ##
##                  ##
######################

sowingDate2024 <- as.Date("2023-04-12")

# Read data
data2024 <- read_xlsx("WSF_2024_data.xlsx", sheet = "CanopyCover", range = "A1:E181", col_names = TRUE) %>% 
  mutate(DAS = as.numeric(as.Date(Date) - sowingDate2024)) %>% 
  rename(Time = DAS, Val = PCover) %>% 
  dplyr::select(Plot, Treatment, Time, Weeds, Val)

data2024$Val[which(data2024$Val < 0.001)] <- 0.001
data2024$Val[which(data2024$Val > 0.999)] <- 0.999

tempData2024 <- read_delim("TempData2024.txt", delim = "\t") %>% 
  mutate(Time = as.integer(as.Date(Date) - sowingDate2024))

cumulativeTemp2024 <- ((tempData2024$TMin + tempData2024$TMax) / 2)
cumulativeTemp2024[which(cumulativeTemp2024 < 0.0)] <- 0
cumulativeTemp2024 <- cumsum(cumulativeTemp2024)

tempData2024 <- tempData2024 %>% 
  mutate(cumulativeTemp = cumulativeTemp2024)

dataWeeds2024 <- data2024 %>% 
  filter(Weeds == "Y") %>% 
  dplyr::select(-Weeds) %>% 
  mutate(Time = tempData2024$cumulativeTemp[match(Time, tempData2024$Time)])
dataNoWeeds2024 <- data2024 %>% 
  filter(Weeds == "N") %>% 
  dplyr::select(-Weeds) %>% 
  mutate(Time = tempData2024$cumulativeTemp[match(Time, tempData2024$Time)])

## Manual fits
plotFit <- function(data, par, fun){
  plot <- ggplot() +
    geom_point(data = data, aes(x = Time, y = Val), size = 2, shape = 1) +
    geom_function(fun = function(x)allDeterministicFunctions[[fun]](par, x),
                  size = 1.5)
  return(plot)
}


##################
##################
####          ####
#### No weeds ####
####          ####
##################
##################

dataTreatment2024NW <- split(dataNoWeeds2024, dataNoWeeds2024$Treatment)
dataTreatment2024NW <- lapply(dataTreatment2024NW, function(dat){
  datNew <- dat %>% 
    dplyr::select(Time, Val)
  return(datNew)
})




##################################
###                            ###
### Fit logistic normal curves ###
###                            ###
##################################     

fitsNWLogisticNormal2024 <- rep(list(NA), 3)

# Fit of 1T:1F
fitsNWLogisticNormal2024[[1]] <- fitDataToModel(data = dataTreatment2024NW[[1]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 450, "sd" = 0.1))
plotFit(dataTreatment2024NW[[1]], fitsNWLogisticNormal2024[[1]]$par, "logistic")

# Fit of T
fitsNWLogisticNormal2024[[2]] <- fitDataToModel(data = dataTreatment2024NW[[2]], model = logisticNormal, startParameters = c("r" = 0.01, "h" = 520, "sd" = 0.1))
plotFit(dataTreatment2024NW[[2]], fitsNWLogisticNormal2024[[2]]$par, "logistic")

# Fit of F
fitsNWLogisticNormal2024[[3]] <- fitDataToModel(data = dataTreatment2024NW[[3]], model = logisticNormal, startParameters = c("r" = 0.012, "h" = 320, "sd" = 0.1))
plotFit(dataTreatment2024NW[[3]], fitsNWLogisticNormal2024[[3]]$par, "logistic")

AICLogisticNormal2024NW <- sum(unlist(lapply(fitsNWLogisticNormal2024, function(fits)2 * 3 + 2 * fits$value)))


##############################################
###                                        ###
### Fit transformed logistic normal curves ###
###                                        ###
##############################################      

fitsNWTLogisticNormal2024 <- rep(list(NA), 3)

# Fit of 1T:1F
fitsNWTLogisticNormal2024[[1]] <- fitDataToModel(data = dataTreatment2024NW[[1]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 375, "sd" = 0.1, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2024NW[[1]], fitsNWTLogisticNormal2024[[1]]$par, "tLogistic")

# Fit of T
fitsNWTLogisticNormal2024[[2]] <- fitDataToModel(data = dataTreatment2024NW[[2]], model = tLogisticNormal, startParameters = c("r" = 0.01, "h" = 375, "sd" = 0.1, "d" = 0.9), dMax = 1.0)
plotFit(dataTreatment2024NW[[2]], fitsNWTLogisticNormal2024[[2]]$par, "tLogistic")

# Fit of F
fitsNWTLogisticNormal2024[[3]] <- fitDataToModel(data = dataTreatment2024NW[[3]], model = tLogisticNormal, startParameters = c("r" = 0.02, "h" = 325, "sd" = 0.1, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2024NW[[3]], fitsNWTLogisticNormal2024[[3]]$par, "tLogistic")

AICTLogisticNormal2024NW <- sum(unlist(lapply(fitsNWTLogisticNormal2024, function(fits)2 * 4 + 2 * fits$value)))


#################################
###                           ###
### Fit logistic gamma curves ###
###                           ###
#################################     

fitsNWLogisticGamma2024 <- rep(list(NA), 3)

# Fit of 1T:1F
fitsNWLogisticGamma2024[[1]] <- fitDataToModel(data = dataTreatment2024NW[[1]], model = logisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0))
plotFit(dataTreatment2024NW[[1]], fitsNWLogisticGamma2024[[1]]$par, "logistic")

# Fit of T
fitsNWLogisticGamma2024[[2]] <- fitDataToModel(data = dataTreatment2024NW[[2]], model = logisticGamma, startParameters = c("r" = 0.02, "h" = 500, "shape" = 1.0))
plotFit(dataTreatment2024NW[[2]], fitsNWLogisticGamma2024[[2]]$par, "logistic")

# Fit of F
fitsNWLogisticGamma2024[[3]] <- fitDataToModel(data = dataTreatment2024NW[[3]], model = logisticGamma, startParameters = c("r" = 0.012, "h" = 475, "shape" = 1.0))
plotFit(dataTreatment2024NW[[3]], fitsNWLogisticGamma2024[[3]]$par, "logistic")

AICLogisticGamma2024NW <- sum(unlist(lapply(fitsNWLogisticGamma2024, function(fits)2 * 3 + 2 * fits$value)))


#############################################
###                                       ###
### Fit transformed logistic gamma curves ###
###                                       ###
#############################################      

fitsNWTLogisticGamma2024 <- rep(list(NA), 3)

# Fit of 1T:1F
fitsNWTLogisticGamma2024[[1]] <- fitDataToModel(data = dataTreatment2024NW[[1]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 420, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2024NW[[1]], fitsNWTLogisticGamma2024[[1]]$par, "tLogistic")

# Fit of T
fitsNWTLogisticGamma2024[[2]] <- fitDataToModel(data = dataTreatment2024NW[[2]], model = tLogisticGamma, startParameters = c("r" = 0.01, "h" = 520, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2024NW[[2]], fitsNWTLogisticGamma2024[[2]]$par, "tLogistic")

# Fit of F
fitsNWTLogisticGamma2024[[3]] <- fitDataToModel(data = dataTreatment2024NW[[3]], model = tLogisticGamma, startParameters = c("r" = 0.012, "h" = 350, "shape" = 1.0, "d" = 0.85), dMax = 1.0)
plotFit(dataTreatment2024NW[[3]], fitsNWTLogisticGamma2024[[3]]$par, "tLogistic")

AICTLogisticGamma2024NW <- sum(unlist(lapply(fitsNWTLogisticGamma2024, function(fits)2 * 4 + 2 * fits$value)))

AICLogisticNormal2024NW
AICTLogisticNormal2024NW
AICLogisticGamma2024NW
AICTLogisticGamma2024NW

# Logistic normal is the best fitting model

# Calculate 95% CIs
parameters2024 <- names(fitsNWLogisticNormal2024[[1]]$par)
parameterCIs2024 <- rep(NA, times = 3 * length(parameters2024))
names(parameterCIs2024) <- unlist(lapply(parameters2024, function(par)c(par, paste(par, "Lower", sep = ""), paste(par, "Upper", sep = ""))))
parameterCIs2024 <- rep(list(parameterCIs2024), length(fitsNWLogisticNormal2024))

cI2024 <- function(par, x, z, opt){
  pars <- par
  pars[names(opt)] <- opt
  nll <- logisticNormal(pars, x, z)
  return(nll)
}

for(i in 1:length(fitsNWLogisticNormal2024)){
  parameterVecs <- lapply(fitsNWLogisticNormal2024, function(fit){
    parStart <- fit$par / 100
    parEnd <- fit$par * 20
    parVec <- list(seq(parStart[1], parEnd[1], length.out = 400), seq(parStart[2], parEnd[2], length.out = 400), seq(parStart[3], parEnd[3], length.out = 400))
    names(parVec) <- names(fit$par)
    return(parVec)
  })
  optPar <- fitsNWLogisticNormal2024[[i]]$par
  parameterCIs2024[[i]][names(optPar)] <- optPar
  parameterCIs2024[[i]][unlist((lapply(names(optPar), function(n)c(paste(n, "Lower", sep = ""), paste(n, "Upper", sep = "")))))] <-
    unlist(lapply(1:length(parameterVecs[[i]]), function(j){
      t(tibble(calculate95CI(parStartVec = optPar[-which(names(optPar) == names(fitsNWLogisticNormal2024[[i]]$par)[[j]])],
                             optVec = parameterVecs[[i]][j],
                             model = cI2024,
                             optNLL = fitsNWLogisticNormal2024[[i]]$value,
                             x = dataTreatment2024NW[[i]][,1][[1]],
                             z = dataTreatment2024NW[[i]][,2][[1]])))
    }))
  
}
names(parameterCIs2024) <- names(dataTreatment2024NW)

# Write to table
CIOut2024 <- tibble(Treatment = names(parameterCIs2024)) %>% 
  bind_cols(bind_rows(parameterCIs2024))
write.table(CIOut2024, "light_interception/light_interception_NW_logisticNormal_parameters_2024.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

CIOut2024Table <- CIOut2024 %>% 
  mutate(r = paste0(round(r, 4), " (", round(rLower, 4), "-", round(rUpper, 4), ")"),
         h = paste0(round(h, 2), " (", round(hLower, 2), "-", round(hUpper, 2), ")"),
         sd = paste0(round(sd, 4), " (", round(sdLower, 4), "-", round(sdUpper, 4), ")")) %>% 
  dplyr::select(Treatment, r, h, sd)

print(xtable(CIOut2024Table, type = "latex"), include.rownames = FALSE)

dataTreatment2024NWM <- lapply(1:length(dataTreatment2024NW), function(i){
  dat <- dataTreatment2024NW[[i]]
  treatment <- names(dataTreatment2024NW)[i]
  dat <- dat %>% 
    mutate(Treatment = treatment)
})

dataTreatment2024NWM <- lapply(dataTreatment2024NWM, function(Data){
  Average <- aggregate(Data$Val, list(Data$Time), function(x)mean(x, na.rm = TRUE))
  Data$Average <- Average$x[match(Data$Time, Average$Group.1)]
  return(Data)
})

# Plot triple combination graphs
geomPoints <- lapply(dataTreatment2024NWM, function(d){
  return(geom_point(data = d, aes(x = Time, y = Average, col = as.factor(Treatment[1])), size = 2, shape = 1))
})

geomFunctions <- lapply(1:length(fitsNWTLogisticNormal2024), function(i){
  o <- fitsNWTLogisticNormal2024[[i]]
  return(geom_function(fun = function(x)(exp(o$par["r"]*(x - o$par["h"])))/(1 + exp(o$par["r"]*(x - o$par["h"]))),
                       aes(col = as.factor(dataTreatment2024NWM[[i]]$Treatment[1])),
                       size = 1.5))
})

### Get all triple combinations
treatmentsToCompare <- unique(data2024$Treatment)
triples <- list(c("T", "1T:1F", "F"))
names(dataTreatment2024NWM) <- lapply(dataTreatment2024NWM, function(l)l$Treatment[1])
names(geomPoints) <- names(dataTreatment2024NWM)
names(geomFunctions) <- names(dataTreatment2024NWM)

triplePlotsNW <- lapply(triples, function(triple){
  plotMultipleGraphs(geomPoints[triple], geomFunctions[triple], triple, colourVector[c(1, 3, 2)], XLab = bquote(T[c]), YLab = bquote(p[cc]), Title = "")
})

### Combine data for AIC comparisons
# Fit of 1T:1F
startParameters1T1F <- list(c("r" = 0.01, "h" = 480, "sd" = 0.1, "d" = 0.9),
                            c("r" = 0.01, "h" = 350, "sd" = 0.1, "d" = 0.9),
                            c("r" = 0.012, "h" = 550, "sd" = 0.1, "d" = 0.85),
                            c("r" = 0.01, "h" = 420, "sd" = 0.1, "d" = 0.9),
                            c("r" = 0.01, "h" = 480, "sd" = 0.1, "d" = 0.9))
combinationAIC1T1F <- round(fitCombinations(data = dataTreatment2024NW, triple = triples[[1]], model = tLogisticNormal, startParameters = startParameters1T1F), 2)

combinationAIC2024NW <- tibble(Treatment = rep(unique(data2024$Treatment)[c(3)], each = 3),
                               Combination = rep(c("IC & C & L", "IC + C & L", "IC + L & C"), times = 1),
                               AIC = c(combinationAIC1T1F),
                               Best = c("*",
                                        "",
                                        ""))

print(xtable(combinationAIC2024NW, type = "latex"), include.rownames = FALSE)

# Print all plots to PDF
pdf("light_interception/light_interception_NW_tLogisticNormal_graphs_2024_NoWeeds.pdf")
print(triplePlotsNW)
dev.off()

jpeg("figures/canopy_cover/plotCanopyCover2024Combined.jpg", units = "px", width = 850, height = 800, quality = 600)
triplePlotsNW[[1]] + labs(y = "", subtitle = bquote({P[cc]})) +
  theme(legend.position = "none",
        plot.subtitle = element_text(hjust = -0.05)) +
  geom_segment(aes(x = 375, y = 0.9, xend = 420, yend = 0.85),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1) +
  geom_segment(aes(x = 470, y = 0.70, xend = 440, yend = 0.75),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1) +
  geom_segment(aes(x = 500, y = 0.3, xend = 460, yend = 0.34),
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1) +
  annotate("text", x = c(300, 470, 500), y = c(0.95, 0.70, 0.3), label = c("Cereal", "Intercrop", "Legume"), size = 8, hjust = 0.0, col = modelColors[c("Cereal", "Intercrop", "Legume")])
  
dev.off()

AICsTable <- tibble(Year = c(rep("2022", 4), rep("2023A", 4), rep("2023B", 4), rep("2024", 4)),
                    Model = rep(c("Logistic normal",
                                  "Transformed logistic normal",
                                  "Logistic gamma",
                                  "Transformed logistic gamma"), times = 4),
                    AIC = c(AICLogisticNormal2022,
                              AICTLogisticNormal2022,
                              AICLogisticGamma2022,
                              AICTLogisticGamma2022,
                              AICLogisticNormal2023ANW,
                              AICTLogisticNormal2023ANW,
                              AICLogisticGamma2023ANW,
                              AICTLogisticGamma2023ANW,
                              AICLogisticNormal2023BNW,
                              AICTLogisticNormal2023BNW,
                              AICLogisticGamma2023BNW,
                              AICTLogisticGamma2023BNW,
                              AICLogisticNormal2024NW,
                              AICTLogisticNormal2024NW,
                              AICLogisticGamma2024NW,
                              AICTLogisticGamma2024NW
                              )) %>% 
  mutate(AIC = round(AIC, 2))
print(xtable(AICsTable, type = "latex"), include.rownames = FALSE)


combinationAICAll <- bind_rows(combinationAIC2022,
                               combinationAIC2023ANW,
                               combinationAIC2023BNW,
                               combinationAIC2024NW
                               ) %>% 
  mutate(Experiment = c(rep("2022", nrow(combinationAIC2022)),
                        rep("2023A", nrow(combinationAIC2023ANW)),
                        rep("2023B", nrow(combinationAIC2023BNW)),
                        rep("2024", nrow(combinationAIC2024NW)))) %>% 
  dplyr::select(Experiment, Treatment, Combination, AIC)

print(xtable(combinationAICAll, type = "latex"), include.rownames = FALSE)
