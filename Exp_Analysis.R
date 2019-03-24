# Packages ----
library(dplyr)
library(ggplot2)
library(ggridges)
library(lme4)
library(LMERConvenienceFunctions)
library(lmerTest)
library(plyr)
library(RePsychLing)
library(sets)
library(MuMIn)

# Customised ggplot theme
theme_AGD <- function () { 
  theme_bw(base_size = 16) %+replace% 
    theme(
      axis.title.x = element_text(face="bold", margin = margin(5, 0, 5, 0)),
      axis.title.y  = element_text(face="bold", angle=90, margin = margin(0, 5, 0, 5)),
      axis.text.y  = element_text(),
      plot.title = element_text(face="bold", margin = margin(5, 0, 5, 0)),
      legend.text = element_text(margin = margin(5, 0, 5, 0)),
      legend.title = element_text(face="bold", margin = margin(5, 0, 5, 0)),
      strip.text.x = element_text(margin = margin(10, 0, 10, 0)),
      strip.text.y = element_text(angle=90, margin = margin(0, 10, 0, 10))
    )
}

# Colour palette
IPPalette <- c("#56B4E9", # z_pl
               "#009E73", # z_o
               "#CC79A7") # z_sing

# Function to subset and drop levels
dsubset <- function (...) { return(droplevels(subset(...))) }

# Function to check for NAs in any column
nacols <- function(df) {
  colnames(df)[unlist(lapply(df, function(x) anyNA(x)))]
}

# Turn off scientific notation (can turn on again with options(scipen=0))
options(scipen=999)

# Load data ----
infl_prime_data <- read.csv("input/infl_prime_data.csv") # trimmed data
nrow(infl_prime_data) # 51072
nacols(infl_prime_data)
str(infl_prime_data)

# Stimuli properties ----
# Table of mean/SD frequency and LSA
infl_prime_design <- ddply(infl_prime_data,
                           .(Prime.Word, Word),
                           summarise,
                           Prime.Condition = unique(Prime.Condition),
                           Condition = unique(Condition),
                           LSA = unique(LSA),
                           Prime.Freq = unique(Prime.Freq))
ddply(infl_prime_design, .(Prime.Condition, Condition), summarise,
      meanLSA = round(mean(LSA),3),
      sdLSA = round(sd(LSA),3),
      meanFreq = round(mean(Prime.Freq),2),
      sdFreq = round(sd(Prime.Freq),2),
      count = length(Word))

# Remove subjs with <70% accuracy ----

infl_prime_data_pre_low_acc_subj <- infl_prime_data
subject_accuracy_pre <-  ddply(infl_prime_data, .(Subj), summarise,
      overall_acc = mean((as.numeric(as.character(Acc)) + as.numeric(as.character(Prime.Acc)))/2))
mean_acc <- mean(subject_accuracy_pre$overall_acc)
sd_acc <- sd(subject_accuracy_pre$overall_acc)
lowerbound <- mean(subject_accuracy_pre$overall_acc) - 2.5*sd(subject_accuracy_pre$overall_acc)
upperbound <- mean(subject_accuracy_pre$overall_acc) + 2.5*sd(subject_accuracy_pre$overall_acc)
InflPrime3AccuracyDistPre <- ggplot(subject_accuracy_pre, aes(x = overall_acc)) +
  geom_density(fill = "gray") +
  geom_vline(xintercept = lowerbound, colour = "blue") +
  geom_vline(xintercept = upperbound, colour = "blue") +
  geom_vline(xintercept = .7, colour = "red")
InflPrime3AccuracyDistPre
ggsave(plot=InflPrime3AccuracyDistPre, filename="output/InflPrime3AccuracyDistPre.png", width = 7, height = 7, dpi=300, device = "png")

ggplot(subject_accuracy_pre, aes(x = Subj, y = overall_acc)) +
         geom_point() +
  geom_hline(yintercept = 0.7, colour = "red")
infl_prime_data <- merge(infl_prime_data, subject_accuracy_pre, all.x = TRUE)
infl_prime_data <- subset(infl_prime_data, infl_prime_data$overall_acc > 0.7)

subject_accuracy_post <-  ddply(infl_prime_data, .(Subj), summarise,
                           overall_acc = mean((as.numeric(as.character(Acc)) + as.numeric(as.character(Prime.Acc)))/2))
InflPrime3AccuracyDistPost <- ggplot(subject_accuracy_post, aes(x = overall_acc)) +
  geom_density(fill = "gray") +
  geom_vline(xintercept = lowerbound, colour = "blue") +
  geom_vline(xintercept = upperbound, colour = "blue") +
  geom_vline(xintercept = .7, colour = "red")
InflPrime3AccuracyDistPost
ggsave(plot=InflPrime3AccuracyDistPost, filename="output/InflPrime3AccuracyDistPost.png", width = 7, height = 7, dpi=300, device = "png")

length(unique(infl_prime_data_pre_low_acc_subj$Subj)) - length(unique(infl_prime_data$Subj))

# Remove filler ----
infl_prime_data <- subset(infl_prime_data, infl_prime_data$ExpFiller == "exp1")

# Accuracy info ----
ddply(infl_prime_data, .(Acc), summarise, count = length(Word))
infl_prime_data$isAccNumeric <- as.numeric(as.character(infl_prime_data$Acc))
ddply(infl_prime_data, .(Prime.Condition), summarise,
      count = length(Word),
      meanAcc= round(mean(isAccNumeric),2))

# Subset to accurate data ----
infl_prime_data_pre_acc <- infl_prime_data
infl_prime_data <- subset(infl_prime_data, Acc == 1 & Prime.Acc == 1)

# Remove absurd RT outliers ----
infl_prime_data_pre_absurd <- infl_prime_data
ggplot(infl_prime_data, aes(x = Trial, y = RT)) +
  geom_point() +
  ylim(0, 6000) +
  geom_hline(yintercept = 300, colour = "red") +
  geom_hline(yintercept = 3000, colour = "red")
infl_prime_data <- subset(infl_prime_data, !(RT < 300 | RT > 3000 |
                                               Prime.RT < 300 | Prime.RT > 3000))

# Remove ISI outliers ----
infl_prime_data_pre_isi <- infl_prime_data
ggplot(infl_prime_data, aes(x = Trial, y = Prime.ISI)) +
  geom_point() +
  ylim(0, 2000) +
  geom_hline(yintercept = 900, colour = "red")
infl_prime_data <- subset(infl_prime_data, infl_prime_data$Prime.ISI < 900)

# Trimming by Subj ----
infl_prime_data_pre_subj_trim <- infl_prime_data

# Load bounds and remove outliers
subj_bounds <- read.csv("input/subj_rt_bounds.csv")
infl_prime_data <- merge(infl_prime_data, subj_bounds, by='Subj', all.x=TRUE)
infl_prime_data <- subset(infl_prime_data, RT < S_highBound & RT > S_lowBound)

# Trimming by target ----
infl_prime_data_pre_target_trim <- infl_prime_data
ggplot(infl_prime_data, aes(sample = RT)) +
geom_qq() +
  facet_wrap(~ Word)

# Load bounds and remove outliers
target_bounds <- read.csv("input/target_rt_bounds.csv")
infl_prime_data <- merge(infl_prime_data, target_bounds, by='Word', all.x=TRUE)
infl_prime_data <- subset(infl_prime_data, RT < T_highBound & RT > T_lowBound)

# Transform RT ----

ggplot(infl_prime_data, aes(sample = RT)) + 
  geom_qq() +
  geom_qq_line()
infl_prime_data$logRT <- log2(infl_prime_data$RT)
ggplot(infl_prime_data, aes(sample = logRT)) +
  geom_qq() +
  geom_qq_line()

ggplot(infl_prime_data, aes(sample = Prime.RT)) + 
  geom_qq() +
  geom_qq_line()
infl_prime_data$Prime.logRT <- log2(infl_prime_data$Prime.RT)
ggplot(infl_prime_data, aes(sample = Prime.logRT)) +
  geom_qq() +
  geom_qq_line()

# Summary of response times ----

ddply(infl_prime_data, .(Prime.Condition), summarise,
      meanRT = round(mean(RT),1),
      sdRT = round(sd(RT),1))

# Coding variables ----

# drop levels
infl_prime_data <- droplevels(infl_prime_data)

# Code prime condition
infl_prime_data$Prime.Condition <- relevel(infl_prime_data$Prime.Condition, ref = "z_pl")
contrasts(infl_prime_data$Prime.Condition)

# Code isMorph
infl_prime_data$isMorph <- 0
infl_prime_data[infl_prime_data$Prime.Condition=="z_pl",]$isMorph <- 1
infl_prime_data$isMorph <- as.factor(infl_prime_data$isMorph)
contrasts(infl_prime_data$isMorph)

# Code isPhonCntrlCntrl
infl_prime_data$isPhonCntrl <- 0
infl_prime_data$isPhonCntrl[infl_prime_data$Prime.Condition=="z_o"] <- 1
infl_prime_data$isPhonCntrl <- as.factor(infl_prime_data$isPhonCntrl)
contrasts(infl_prime_data$isPhonCntrl)

# Code isSingCntrl
infl_prime_data$isSingCntrl <- 0
infl_prime_data$isSingCntrl[infl_prime_data$Prime.Condition=="z_sing"] <- 1
infl_prime_data$isSingCntrl <- as.factor(infl_prime_data$isSingCntrl)
contrasts(infl_prime_data$isSingCntrl)

# Centre trial
infl_prime_data$cTrial <- scale(infl_prime_data$Trial, center = TRUE, scale = FALSE)[,]

# Scale all other numeric variables
infl_prime_data$zFreq <- scale(infl_prime_data$Freq)[,]
infl_prime_data$zLSA <- scale(infl_prime_data$LSA)[,]
infl_prime_data$zPrime.ISI <- scale(infl_prime_data$Prime.ISI)[,]
infl_prime_data$zDur <- scale(infl_prime_data$Dur)[,]
infl_prime_data$zPrime.Dur <- scale(infl_prime_data$Prime.Dur)[,]
infl_prime_data$Prime.logRT <- as.numeric(infl_prime_data$Prime.logRT)
infl_prime_data$zPrime.logRT <- scale(infl_prime_data$Prime.logRT)[,]
infl_prime_data$zPrime.Freq <- scale(infl_prime_data$Prime.Freq)[,]
infl_prime_data$zTransProb <- scale(infl_prime_data$TransProb)[,]

# Center by Prime.Condition for select numeric variables
infl_prime_data$c.zPrime.Freq <- ave(infl_prime_data$Prime.Freq, infl_prime_data$Prime.Condition, FUN = scale)
infl_prime_data$c.zPrime.logRT <- ave(infl_prime_data$Prime.logRT, infl_prime_data$Prime.Condition, FUN = scale)
infl_prime_data$c.zPrime.Dur <- ave(infl_prime_data$Prime.Dur, infl_prime_data$Prime.Condition, FUN = scale)
infl_prime_data$c.zLevenshteinDist <- ave(infl_prime_data$LevenshteinDist, infl_prime_data$Prime.Condition, FUN = scale)
infl_prime_data$c.zMFCC <- ave(infl_prime_data$mfcc, infl_prime_data$Prime.Condition, FUN = scale)

# mod0 <- lmer(logRT ~ Prime.Condition +
#                               cTrial + zLSA + zPrime.ISI + zDur + zFreq +
#                               c.zMFCC + c.zLevenshteinDist + c.zPrime.Freq + c.zPrime.logRT +
#                               (isPhonCntrl + isSingCntrl|Subj) +
#                               (isPhonCntrl + isSingCntrl|Word) +
#                               (1|Prime.Word),
#                             data = infl_prime_data)
# failed to converge

mod1 <- lmer(logRT ~ Prime.Condition +
                     cTrial + zLSA + zPrime.ISI + zDur + zFreq +
                     c.zMFCC + c.zLevenshteinDist + c.zPrime.Freq + c.zPrime.logRT +
                     (isPhonCntrl|Subj) +
                     (isPhonCntrl + isSingCntrl|Word) +
                     (1|Prime.Word),
                   data = infl_prime_data)

untrimmed_mod <- lmer(logRT ~ Prime.Condition +
              cTrial + zLSA + zPrime.ISI + zDur + zFreq +
              c.zMFCC + c.zLevenshteinDist + c.zPrime.Freq + c.zPrime.logRT +
              (1|Subj) +
              (1|Word) +
              (1|Prime.Word),
            data = infl_prime_data)

anova(mod1, untrimmed_mod) # n.s.

infl_prime_data_trim <- romr.fnc(untrimmed_mod, infl_prime_data, trim = 2.5)$data

# Re-centre trial
infl_prime_data_trim$cTrial <- scale(infl_prime_data_trim$Trial, scale = FALSE)[,]

# Re-scale all numeric variables
infl_prime_data_trim$zFreq <- scale(infl_prime_data_trim$Freq)[,]
infl_prime_data_trim$zLSA <- scale(infl_prime_data_trim$LSA)[,]
infl_prime_data_trim$zPrime.ISI <- scale(infl_prime_data_trim$Prime.ISI)[,]
infl_prime_data_trim$zDur <- scale(infl_prime_data_trim$Dur)[,]
infl_prime_data_trim$zPrime.Dur <- scale(infl_prime_data_trim$Prime.Dur)[,]
infl_prime_data_trim$Prime.logRT <- as.numeric(infl_prime_data_trim$Prime.logRT)
infl_prime_data_trim$zPrime.logRT <- scale(infl_prime_data_trim$Prime.logRT)[,]
infl_prime_data_trim$zTransProb <- scale(infl_prime_data_trim$TransProb)[,]
infl_prime_data_trim$zMFCC <- scale(infl_prime_data_trim$mfcc)[,]
infl_prime_data_trim$zLevenshteinDist <- scale(infl_prime_data_trim$LevenshteinDist)[,]

# Center by Prime.Condition for select numeric variables
infl_prime_data_trim$c.zPrime.Freq <- ave(infl_prime_data_trim$Prime.Freq, infl_prime_data_trim$Prime.Condition, FUN = scale)
infl_prime_data_trim$c.zPrime.logRT <- ave(infl_prime_data_trim$Prime.logRT, infl_prime_data_trim$Prime.Condition, FUN = scale)
infl_prime_data_trim$c.zPrime.Dur <- ave(infl_prime_data_trim$Prime.Dur, infl_prime_data_trim$Prime.Condition, FUN = scale)
infl_prime_data_trim$c.zLevenshteinDist <- ave(infl_prime_data_trim$LevenshteinDist, infl_prime_data_trim$Prime.Condition, FUN = scale)
infl_prime_data_trim$c.zMFCC <- ave(infl_prime_data_trim$mfcc, infl_prime_data_trim$Prime.Condition, FUN = scale)

final_mod <- lmer(logRT~ Prime.Condition +
                   cTrial + zLSA + zPrime.ISI + zDur + zFreq +
                   c.zMFCC + c.zLevenshteinDist + c.zPrime.Freq + c.zPrime.logRT +
                   (1|Subj) +
                   (1|Word) +
                   (1|Prime.Word),
                 data=infl_prime_data_trim)
summary(final_mod)

r.squaredGLMM(final_mod) # revised statistics based on Nakagawa et al. (2013) paper

# Interpretting coefficients ----
pred <- expand.grid(Prime.Condition=c("z_o","z_pl","z_sing"),
                    cTrial=0, zLSA=0, zPrime.ISI=0, zDur=0, zFreq=0,
                    c.zMFCC=0, c.zLevenshteinDist=0, c.zPrime.Freq=0, c.zPrime.logRT=0)
pred$logRT <- predict(final_mod, newdata=pred, re.form=NA)
pred$RT <- 2^pred$logRT
predTable <- ddply(pred, .(Prime.Condition), summarise, meanlogRT=mean(logRT),
                   meanRT=mean(RT))

round(fixef(final_mod)["Prime.Conditionz_o"],3)
100-round((2^fixef(final_mod)["Prime.Conditionz_o"])*100,2)

predTable[predTable$Prime.Condition=="z_o",]$meanlogRT - predTable[predTable$Prime.Condition=="z_pl",]$meanlogRT

round((predTable[predTable$Prime.Condition=="z_o",]$meanRT - 
    predTable[predTable$Prime.Condition=="z_pl",]$meanRT)*100/
  predTable[predTable$Prime.Condition=="z_pl",]$meanRT,2)
round(predTable[predTable$Prime.Condition=="z_o",]$meanRT - 
  predTable[predTable$Prime.Condition=="z_pl",]$meanRT,2)

round(fixef(final_mod)["Prime.Conditionz_sing"],3)
round((2^fixef(final_mod)["Prime.Conditionz_sing"])*100,2)

predTable[predTable$Prime.Condition=="z_sing",]$meanlogRT - predTable[predTable$Prime.Condition=="z_pl",]$meanlogRT

(predTable[predTable$Prime.Condition=="z_sing",]$meanRT - 
    predTable[predTable$Prime.Condition=="z_pl",]$meanRT)*100/
  predTable[predTable$Prime.Condition=="z_pl",]$meanRT
round(predTable[predTable$Prime.Condition=="z_sing",]$meanRT - 
  predTable[predTable$Prime.Condition=="z_pl",]$meanRT,2)

# Visualise effect sizes and confidence intervals ----

confint_final_mod <- as.data.frame(confint(final_mod))
confint_final_mod$name <- rownames(confint_final_mod)
confint_final_mod$estimate <- NA
confint_final_mod[confint_final_mod$name == "(Intercept)",]$estimate <- 
  fixef(final_mod)[["(Intercept)"]]
confint_final_mod[confint_final_mod$name == "Prime.Conditionz_o",]$estimate <- 
 fixef(final_mod)[["Prime.Conditionz_o"]]
confint_final_mod[confint_final_mod$name == "Prime.Conditionz_sing",]$estimate <- 
  fixef(final_mod)[["Prime.Conditionz_sing"]]
confint_final_mod[confint_final_mod$name == "cTrial",]$estimate <- 
  fixef(final_mod)[["cTrial"]]
confint_final_mod[confint_final_mod$name == "zLSA",]$estimate <- 
  fixef(final_mod)[["zLSA"]]
confint_final_mod[confint_final_mod$name == "zPrime.ISI",]$estimate <- 
  fixef(final_mod)[["zPrime.ISI"]]
confint_final_mod[confint_final_mod$name == "zDur",]$estimate <- 
  fixef(final_mod)[["zDur"]]
confint_final_mod[confint_final_mod$name == "zFreq",]$estimate <- 
  fixef(final_mod)[["zFreq"]]
confint_final_mod[confint_final_mod$name == "c.zMFCC",]$estimate <- 
  fixef(final_mod)[["c.zMFCC"]]
confint_final_mod[confint_final_mod$name == "c.zLevenshteinDist",]$estimate <- 
  fixef(final_mod)[["c.zLevenshteinDist"]]
confint_final_mod[confint_final_mod$name == "c.zPrime.Freq",]$estimate <- 
  fixef(final_mod)[["c.zPrime.Freq"]]
confint_final_mod[confint_final_mod$name == "c.zPrime.logRT",]$estimate <- 
  fixef(final_mod)[["c.zPrime.logRT"]]

confint_final_mod$estimate <- round(confint_final_mod$estimate,3)
confint_final_mod$`2.5 %` <- round(confint_final_mod$`2.5 %`,3)
confint_final_mod$`97.5 %` <- round(confint_final_mod$`97.5 %`,3)

confint_final_mod <- subset(confint_final_mod, !(confint_final_mod$name %in% c("(Intercept)", ".sigma", ".sig03", ".sig02", ".sig01")))


fe_labels <- c("Prime.Conditionz_o" = "Phon. cntrl",
              "Prime.Conditionz_sing" = "Sing. cntrl",
              "cTrial" = "Trial",                     
              "zLSA" = "LSA",                       
              "zPrime.ISI" = "ISI",                   
              "zDur" = "Duration",               
              "zFreq" = "Target freq.",                     
              "c.zMFCC" = "MFCC",                    
              "c.zLevenshteinDist" = "Levenshtein Dist.",          
              "c.zPrime.Freq" = "Prime freq.",             
              "c.zPrime.logRT" = "Prime RT")
InflPrimeFEplot <- ggplot(confint_final_mod, aes(x=name, y=estimate)) +
  geom_point(size=2) + 
  geom_errorbar(aes(ymin=`2.5 %`, ymax=`97.5 %`), width=0.2, size=1) +
  geom_hline(yintercept=0,colour="red", size=1, linetype=2) +
  ylab("Estimate") +
  xlab("Fixed effect") +
  scale_x_discrete(labels=fe_labels) +
  theme_AGD() +
  coord_flip() +
  theme(axis.text.y = element_text(hjust = 1))
InflPrimeFEplot
ggsave(plot=InflPrimeFEplot, filename="output/fixedeffectplot.pdf", width = 9, height = 7, dpi=300, device = "pdf")

# Plot random effects ----

# Summary of data removal
nrow(infl_prime_data_pre_acc)
nrow(infl_prime_data_pre_acc) - nrow(infl_prime_data_pre_absurd)
round(100 * (nrow(infl_prime_data_pre_acc) - nrow(infl_prime_data_pre_absurd)) / nrow(infl_prime_data_pre_acc),2)
nrow(infl_prime_data_pre_absurd) - nrow(infl_prime_data_pre_isi)
round(100 * (nrow(infl_prime_data_pre_absurd) - nrow(infl_prime_data_pre_isi)) / nrow(infl_prime_data_pre_acc),2)
nrow(infl_prime_data_pre_isi) - nrow(infl_prime_data_pre_subj_trim)
round(100 * (nrow(infl_prime_data_pre_isi) - nrow(infl_prime_data_pre_subj_trim)) / nrow(infl_prime_data_pre_acc),2)
nrow(infl_prime_data_pre_subj_trim) - nrow(infl_prime_data_pre_target_trim)
round(100 * (nrow(infl_prime_data_pre_subj_trim) - nrow(infl_prime_data_pre_target_trim)) / nrow(infl_prime_data_pre_acc),2)
nrow(infl_prime_data) - nrow(infl_prime_data_trim)
round(100 * (nrow(infl_prime_data) - nrow(infl_prime_data_trim)) / nrow(infl_prime_data_pre_acc),2)
nrow(infl_prime_data_trim)

# Plots ----
infl_prime_data_trim$Prime.Condition <- factor(infl_prime_data_trim$Prime.Condition, levels = c("z_o", "z_sing", "z_pl"))
InflPrime3Boxplot <- ggplot() +
  geom_boxplot(data=infl_prime_data_trim, aes(x=Prime.Condition, y=logRT, fill=Prime.Condition))+
  geom_hline(yintercept = median(infl_prime_data_trim[infl_prime_data_trim$Prime.Condition=="z_pl",]$logRT), colour="black", linetype=2, size = 1) +
  xlab("Prime condition") +
  ylab("Log2-transformed RT") +
  scale_fill_manual(values=alpha(IPPalette,0.9)) +
  guides(fill=FALSE) +
  theme_AGD() +
  theme_bw(base_size = 20) +
  scale_x_discrete(name="Prime Condition",
                   breaks=c("z_o", "z_pl", "z_sing"),
                   labels=c("Phon.\nControl", "Plural", "Sing.\nControl")) +
  scale_y_continuous(limits=c(9,11))
InflPrime3Boxplot
ggsave(plot=InflPrime3Boxplot, filename="output/boxplot.pdf", width = 7, height = 7, dpi=300, device = "pdf")

InflPrime3Ridge <- ggplot(data=infl_prime_data_trim, aes(y=Prime.Condition, x=logRT, fill=Prime.Condition)) + geom_density_ridges(scale = 0.5) +
  ylab("Prime condition") +
  xlab("Log2-transformed RT") +
  guides(fill=FALSE) +
  theme_AGD() +
  theme_bw(base_size = 20) +
  scale_fill_manual(values=alpha(IPPalette,0.9)) +
  scale_y_discrete(name="Prime Condition",
                   breaks=c("z_o", "z_pl", "z_sing"),
                   labels=c("Phon.\nControl", "Plural", "Sing.\nControl")) +
  scale_x_continuous(limits=c(9,11)) +
  geom_vline(xintercept = median(infl_prime_data_trim[infl_prime_data_trim$Prime.Condition=="z_pl",]$logRT), colour="black", linetype=2, size=1) +
  coord_flip()
InflPrime3Ridge
ggsave(plot=InflPrime3Ridge, filename="output/ridgeplot.pdf", width = 7, height = 7, dpi=300, device = "pdf")