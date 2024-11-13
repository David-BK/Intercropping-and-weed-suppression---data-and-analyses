## Author: David Kottelenberg
## Institute: Wageningen University & Research
## Last edited: 13/11/2024
## Summary: Functions to aid data analysis of weed biomass

library(tidyverse)
library(emmeans)
library(multcomp)
library(stringr)

source("WSFTreatments.R")


calcArithmeticHarmonic <- function(data, dataSC){
  soleC2022 <- c("Barley", "Rye", "Triticale", "Wheat")
  soleL2022 <- c("Pea", "Lupine", "Faba")
  soleC2023B <- c("T", "T-25", "T-375")
  soleL2023B <- c("F", "F-25", "F-375")
  inter2023B <- c("1T:1F", "1T:1F-25", "1T:1F-375", "TF-M", "TF-M-25", "TF-M-375")
  dataComb <- lapply(data, function(dat){
    datSCCSel <- dataSC[which(dataSC$Treatment %in% c(soleC2022, soleC2023B, "T")),]
    datSCCSelW <- datSCCSel[which(datSCCSel$Treatment %in% c(strsplit(dat$Treatment[1], "_")[[1]][1], "T", paste("T-", strsplit(dat$Treatment[1], "-")[[1]][length(strsplit(dat$Treatment[1], "-")[[1]])], sep = ""))),]
    if(nrow(datSCCSelW) > nrow(dat)){
      datSCCSelWW <- filter(datSCCSelW, Treatment == "T-375")$W
    }else{
      datSCCSelWW <- filter(datSCCSelW, Treatment %in% c("T", datSCCSelW$Treatment[1]))$W
    }
    dat$SCC <- datSCCSelWW
    dat$SCC[which(is.na(dat$SCC))] <- mean(dat$SCC, na.rm = TRUE)
    datSCLSel <- dataSC[which(dataSC$Treatment %in% c(soleL2022, soleL2023B, "F")),]
    datSCLSelW <- datSCLSel[which(datSCLSel$Treatment %in% c(strsplit(dat$Treatment[1], "_")[[1]][2], "F", paste("F-", strsplit(dat$Treatment[1], "-")[[1]][length(strsplit(dat$Treatment[1], "-")[[1]])], sep = ""))),]
    if(nrow(datSCLSelW) > nrow(dat)){
      datSCLSelWW <- filter(datSCLSelW, Treatment == "F-375")$W
    }else{
      datSCLSelWW <- filter(datSCLSelW, Treatment %in% c("F", datSCLSelW$Treatment[1]))$W
    }
    dat$SCL <- datSCLSelWW
    dat$SCL[which(is.na(dat$SCL))] <- mean(dat$SCL, na.rm = TRUE)
    if(dat$Treatment[1] %in% c("1T:1F", "TF-M", inter2023B) | grepl("_", dat$Treatment[1])){
      dat <- mutate(dat, Arithmetic = (dat$SCC + dat$SCL) / 2)
      dat$Harmonic <- 1 / (0.5 * (1 / dat$SCC) + 0.5 * (1 / dat$SCL))
    }else if(dat$Treatment[1] == "1T:3F"){
      dat$Arithmetic <- (0.5 * dat$SCC + 1.5 * dat$SCL) / 2
      dat$Harmonic <- 1 / (0.25 * (1 / dat$SCC) + 0.75 * (1 / dat$SCL))
    }else{
      dat$Arithmetic <- (0.5 * dat$SCL + 1.5 * dat$SCC) / 2
      dat$Harmonic <- 1 / (0.25 * (1 / dat$SCL) + 0.75 * (1 / dat$SCC))
    }
    return(dat)
  })
  return(dataComb)
}


calcWeeds <- function(data){
  # modWeedsA2 <- glmer(W ~ Treatment + (1|Block), data = dataWeedsA2, family = "Gamma")
  modWeeds <- lmer(W ~ Treatment + (1|Block), data = data)
  
  # variable <- ifelse(relative, "RelativeBiomass", "Biomass")
  # Check for normality and heteroscedasticity
  ## T
  # modbiomass <- lmer(W ~ (1|Block) + Treatment, data = data)
  
  nonNormal <- shapiro.test(residuals(modWeeds))$p.value < 0.05
  unequalVariance <- bartlett.test(residuals(modWeeds) ~ data$Treatment)$p.value < 0.05
  
  newDat <- data
  if(nonNormal | unequalVariance){
    bc <- boxcox(data$W ~ data$Treatment, lambda = seq(-5, 5, 1/10))
    lambda <- bc$x[which.max(bc$y)]
    if(lambda > -0.3 & lambda < 0.3){
      newDat$W <- log(newDat$W)
    }else{
      newDat$W <- newDat$W^(lambda)  
    }
    
    modWeeds <- lmer(newDat$W ~ (1|Block) + Treatment, data = newDat)
    nonNormal <- shapiro.test(residuals(modWeeds))$p.value < 0.05
    unequalVariance <- bartlett.test(residuals(modWeeds) ~ newDat$Treatment)$p.value < 0.05
    if(nonNormal | unequalVariance){
      cat("Cannot guarantee normal distribution and/or equal variances. Please manually check for these assumptions. \n")
    }
  }
  
  cat(paste("p-value = ", summary(aov(newDat$W ~ (1|Block) + Treatment, data = newDat))[[1]]$`Pr(>F)`[1], "\n", sep = ""))
  
  PHWeeds <- emmeans(modWeeds, list(pairwise ~ Treatment), adjust = "tukey")
  CLDWeeds <- cld(PHWeeds$emmeans,
                    Letters = letters,
                    decreasing = TRUE)
  
  return(CLDWeeds)
  
  # dataWeedsA2 <- mutate(dataWeedsA2, group = sapply(dataWeedsA2$Treatment, function(t)ifelse(t == "T" | t == "T+", "Sole triticale", ifelse(t == "F" | t == "F+", "Sole faba", "Intercrop"))))
  
}

calcWeedsGrouped <- function(data){
  modWeeds <- lmer(W ~ Group + (1|Block), data = data)
  nonNormal <- shapiro.test(residuals(modWeeds))$p.value < 0.05
  unequalVariance <- bartlett.test(residuals(modWeeds) ~ data$Group)$p.value < 0.05
  
  newDat <- data
  if(nonNormal | unequalVariance){
    bc <- boxcox(data$W ~ data$Group, lambda = seq(-5, 5, 1/10))
    newDat$WeedBiomass <- newDat$W^(bc$x[which.max(bc$y)])
    modWeeds <- lmer(newDat$W ~ (1|Block) + Group, data = newDat)
    nonNormal <- shapiro.test(residuals(modWeeds))$p.value < 0.05
    unequalVariance <- bartlett.test(residuals(modWeeds) ~ newDat$Group)$p.value < 0.05
    if(nonNormal | unequalVariance){
      cat("Cannot guarantee normal distribution and/or equal variances. Please manually check for these assumptions. \n")
    }
  }
  
  cat(paste("p-value = ", summary(aov(newDat$W ~ (1|Block) + Group, data = newDat))[[1]]$`Pr(>F)`[1], "\n", sep = ""))
  
  PHWeeds <- emmeans(modWeeds, list(pairwise ~ Group), adjust = "tukey")
  CLDWeeds <- cld(PHWeeds$emmeans,
                  Letters = letters,
                  decreasing = TRUE)
  return(CLDWeeds)
}

predictWeeds <- function(dataWeeds){
  sole2022 <- c("Barley", "Rye", "Triticale", "Wheat", "Lupine", "Pea", "Faba")
  sole2023A <- c("T", "T+", "F", "F+")
  sole2023B <- c("T", "T-25", "T-375", "F", "F-25", "F-375")
  sole2024 <- c("T", "F")
  dataWeedsSC <- dataWeeds %>% 
    filter(Treatment %in% c(sole2022, sole2023A, sole2023B, sole2024))
  dataWeedsIC <- dataWeeds %>% 
    filter(!Treatment %in% c(sole2022, sole2023A, sole2023B, sole2024))

  dataWeedsICSplit <- split(dataWeedsIC, dataWeedsIC$Treatment)
  dataWeedsICSplitPred <- calcArithmeticHarmonic(dataWeedsICSplit, dataWeedsSC)
  dataWeedsICPred <- bind_rows(dataWeedsICSplitPred)

  lmArithmetic <- lm(W ~ Arithmetic, data = dataWeedsICPred)
  lmHarmonic <- lm(W ~ Harmonic, data = dataWeedsICPred)
  
  dataWeedsICPredLong <- dataWeedsICPred %>% 
  pivot_longer(cols = c(Arithmetic, Harmonic), names_to = "Model", values_to = "Predicted")
  dataWeedsICPredLong <- mutate(dataWeedsICPredLong, Treatment = factor(Treatment))
  plotWeedPred <- ggplot(data = dataWeedsICPredLong, aes(x = Predicted, y = W, group_by = Model, color = Model)) +
    geom_point(size = 6, shape = 1, stroke = 2) + 
    geom_abline(intercept = 0, slope = 1,
                size = 2, 
                linetype = "dashed") +
    geom_function(fun = function(x)return(summary(lmArithmetic)$coefficients[1,1] + summary(lmArithmetic)$coefficients[2,1] * x),
                  col = "#F8766D",
                  size = 2) +
    geom_function(fun = function(x)return(summary(lmHarmonic)$coefficients[1,1] + summary(lmHarmonic)$coefficients[2,1] * x),
                  col = "#00BFC4",
                  size = 2) +
    theme_bw(base_size = 30) +
    labs(y = bquote("Observed weed biomass (g "~m^-2~")"),
         x = bquote("Predicted weed biomass (g "~m^-2~")"),
         title = "Predicted and observed weed biomass") +
    scale_shape_manual(values = 1:nlevels(dataWeedsICPredLong$Treatment))
  
  predWeeds2022Long <- dataWeedsICPred %>% 
    pivot_longer(cols = c("W", "Arithmetic", "Harmonic"), names_to = "Type", values_to = "WeedBiomass")
  
  modWeedPred <- lmer(WeedBiomass ~ Type + (1|Block), data = predWeeds2022Long)
  
  nonNormal <- shapiro.test(residuals(modWeedPred))$p.value < 0.05
  unequalVariance <- bartlett.test(residuals(modWeedPred) ~ na.omit(predWeeds2022Long)$Type)$p.value < 0.05
  
  newDat <- predWeeds2022Long
  if(nonNormal | unequalVariance){
    bc <- boxcox(predWeeds2022Long$WeedBiomass ~ predWeeds2022Long$Type, lambda = seq(-5, 5, 1/10))
    modVal <- bc$x[which.max(bc$y)]
    if(modVal < 0.2 & modVal > -0.2){
      newDat$WeedBiomass <- log(newDat$WeedBiomass)
    }else{
      newDat$WeedBiomass <- newDat$WeedBiomass^(bc$x[which.max(bc$y)])
    }
    modWeedPred <- lmer(newDat$WeedBiomass ~ (1|Block) + Type, data = newDat)
    nonNormal <- shapiro.test(residuals(modWeedPred))$p.value < 0.05
    unequalVariance <- bartlett.test(residuals(modWeedPred) ~ na.omit(newDat)$Type)$p.value < 0.05
    if(nonNormal | unequalVariance){
      cat("Cannot guarantee normal distribution and/or equal variances. Please manually check for these assumptions. \n")
    }
  }

  
    return(list(dataWeedsICPred, plotWeedPred))
}
