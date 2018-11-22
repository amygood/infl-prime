library(ggplot2)
library(lme4)
library(simr)
library(LMERConvenienceFunctions)

# Read in and prepare data ----
pilot_data <- read.csv("input/pilot_data.csv")
nrow(pilot_data)
# Code variables
pilot_data$isMorph <- 0
pilot_data[pilot_data$Prime.Condition == "z_pl",]$isMorph <- 1
pilot_data$isPhon <- 1
pilot_data[pilot_data$Prime.Condition == "z_sing",]$isPhon <- 0

# Scale all numeric variables
pilot_data$zTrial <- scale(pilot_data$Trial)[,]
pilot_data$zFreq <- scale(pilot_data$Freq)[,]
pilot_data$zLSA <- scale(log(pilot_data$LSA+1))[,]
pilot_data$zPrime.ISI <- scale(pilot_data$Prime.ISI)[,]
pilot_data$zDur <- scale(pilot_data$Dur)[,]
pilot_data$zPrime.Dur <- scale(pilot_data$Prime.Dur)[,]
pilot_data$Prime.logRT <- as.numeric(pilot_data$Prime.logRT)
pilot_data$zPrime.logRT <- scale(pilot_data$Prime.logRT)[,]
pilot_data$zPrime.Freq <- scale(pilot_data$Prime.Freq)[,]
pilot_data$zTransProb <- scale(pilot_data$TransProb)[,]

# center by Prime.Condition for select numeric variables
pilot_data$c.zPrime.Freq <- ave(pilot_data$Prime.Freq,
                                pilot_data$Prime.Condition,
                                FUN = scale)
pilot_data$c.zPrime.logRT <- ave(pilot_data$Prime.logRT,
                                 pilot_data$Prime.Condition,
                                 FUN = scale)
pilot_data$c.zPrime.Dur <- ave(pilot_data$Prime.Dur,
                                        pilot_data$Prime.Condition,
                                        FUN = scale)
pilot_data$c.zLevenshteinDist <- ave(pilot_data$LevenshteinDist,
                                     pilot_data$Prime.Condition,
                                     FUN = scale)
pilot_data$c.zMFCC <- ave(pilot_data$mfcc,
                          pilot_data$Prime.Condition,
                          FUN = scale)

# Change the names and levels of the Subj column
pilot_data$Subj <- as.factor(pilot_data$Subj)
pilot_data$Subj <- droplevels(pilot_data$Subj)
levels(pilot_data$Subj) <- c(1:length(levels(pilot_data$Subj)))
pilot_data$Subj <- as.numeric(pilot_data$Subj)
unique(pilot_data$Subj)

# Code prime condition
pilot_data$Prime.Condition <- relevel(pilot_data$Prime.Condition, ref = "z_pl")
contrasts(pilot_data$Prime.Condition)

# Run the maximal model
Mod <- lmer(logRT ~ Prime.Condition +
              zTrial + zLSA + zPrime.ISI + zDur + zFreq +
              c.zMFCC + c.zLevenshteinDist + c.zPrime.Freq + c.zPrime.logRT +
              (1|Subj) +
              (1|Word),
            data=pilot_data)
pilot_data_trim <- romr.fnc(Mod, pilot_data, trim = 2.5)$data
FinalMod <- lmer(logRT ~ Prime.Condition +
                   zTrial + zLSA + zPrime.ISI + zDur + zFreq +
                   c.zMFCC + c.zLevenshteinDist + c.zPrime.Freq + c.zPrime.logRT +
                   (1|Subj) +
                   (1|Word),
                 data=pilot_data_trim)
nrow(pilot_data_trim) == nrow(FinalMod@frame)
summary(FinalMod)
cor(fitted(FinalMod),pilot_data_trim$logRT)^2
str(pilot_data)

# Power calculations ----
NumSims <- 1000
InflPrimePowMod<- extend(FinalMod, along = "Subj", n = 200)
InflPrimePowModData <- getData(InflPrimePowMod)

# 15ms analysis ----
InflPrimePowMod15 <- InflPrimePowMod
betaAdjustment <- 15
newInflPrimeBeta <- log(exp(fixef(FinalMod)[["(Intercept)"]]) - betaAdjustment) -  fixef(FinalMod)[["(Intercept)"]]

fixef(InflPrimePowMod15)["Prime.Conditionz_sing"] <- newInflPrimeBeta

PC15 <- powerCurve(InflPrimePowMod15, along="Subj", fixed("Prime.Conditionz_sing","t"), nsim=NumSims, seed=18)

exp(fixef(InflPrimePowMod15)[["(Intercept)"]]) - exp(fixef(InflPrimePowMod15)[["(Intercept)"]] - fixef(InflPrimePowMod15)[["Prime.Conditionz_sing"]])

PC15$errors
print(PC15)
plot(PC15)
save(PC15, file = "output/PC15.RData")