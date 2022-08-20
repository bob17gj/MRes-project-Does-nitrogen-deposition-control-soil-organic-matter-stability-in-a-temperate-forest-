library(tidyverse) #for all data wrangling
library(car) # used to test for outliers
library(cowplot) #for manuscript ready figures
library(car)# for bonferroni outlier test
library(multcomp) # for ancova tukey test


rm(list=ls())
setwd("~/General Waring/Project")

# save.image(file ="Figure 3 - per gram MBC.RData")
load("Figure 3 - per gram MBC.RData")

hydrolytic_enzymes <- read.csv("Hydrolytic enzyme master sheet for R (2).csv")
str(hydrolytic_enzymes)


# create data frame with new variable relating to enzyme activity per unit microbial biomass C
AP.mbc <- hydrolytic_enzymes
AP.mbc <- mutate(hydrolytic_enzymes, AP.mbc = enzyme.activity/microbial.biomass.C)
AP.mbc <-AP.mbc%>%filter(enzyme=="AP")

BG.mbc <- hydrolytic_enzymes
BG.mbc <- mutate(hydrolytic_enzymes, BG.mbc = enzyme.activity/microbial.biomass.C)
BG.mbc <-BG.mbc%>%filter(enzyme=="BG")

BX.mbc <- hydrolytic_enzymes
BX.mbc <- mutate(hydrolytic_enzymes, BX.mbc = enzyme.activity/microbial.biomass.C)
BX.mbc <-BX.mbc%>%filter(enzyme=="BX")

PPO.mbc <- hydrolytic_enzymes
PPO.mbc <- mutate(hydrolytic_enzymes, PPO.mbc = absorbance/microbial.biomass.C)
PPO.mbc <-PPO.mbc%>%filter(enzyme=="PPO")

NAG.mbc <- hydrolytic_enzymes
NAG.mbc <- mutate(hydrolytic_enzymes, NAG.mbc = enzyme.activity/microbial.biomass.C)
NAG.mbc <-NAG.mbc%>%filter(enzyme=="NAG")

LAP.mbc <- hydrolytic_enzymes
LAP.mbc <- mutate(hydrolytic_enzymes, LAP.mbc = enzyme.activity/microbial.biomass.C)
LAP.mbc <-LAP.mbc%>%filter(enzyme=="LAP")

# convert all "NaN values in DT to NA
AP.mbc <- AP.mbc %>% mutate_all(~ifelse(is.nan(.), NA, .))

BG.mbc <- BG.mbc %>% mutate_all(~ifelse(is.nan(.), NA, .))

BX.mbc <- BX.mbc %>% mutate_all(~ifelse(is.nan(.), NA, .))

PPO.mbc <- PPO.mbc %>% mutate_all(~ifelse(is.nan(.), NA, .))

NAG.mbc <- NAG.mbc %>% mutate_all(~ifelse(is.nan(.), NA, .))

LAP.mbc <- LAP.mbc %>% mutate_all(~ifelse(is.nan(.), NA, .))

# change "inf" values to NA
AP.mbc <- AP.mbc %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf)))

BG.mbc <-BG.mbc %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf)))

BX.mbc <-BX.mbc %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf)))

PPO.mbc <-PPO.mbc %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf)))

NAG.mbc <-NAG.mbc %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf)))

LAP.mbc <-LAP.mbc %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf)))

#### Per gram MBC

######### AP

## visualising data
ggplot(AP.mbc, aes(x = AP.mbc, y = TOC, colour = horizon, group = forest.type))+
  geom_point(aes(shape=forest.type))



# simple lm
AP.model <- lm(TOC ~ AP.mbc + horizon + forest.type + NOy.NHx.combined.28.year + soil.type
               , data = AP.mbc)

print(AP.model)
summary(AP.model)


# remove outliers
outlierTest(AP.model,cutoff=Inf, TOC~1)

ggplot(AP.model, aes(x= TOC , y=AP.mbc, group = horizon)) + 
  geom_boxplot()


AP.mbc_1=AP.mbc

AP.mbc_1 <- AP.mbc[-c(25,129,155,103,51,77,98,150, 141, 35, 104, 61, 147, 159, 9, 207, 17, 145, 47, 142, 42,125,21, 64, 12, 221),] # also removing values that later on cause non-normality

AP.mbc_1["AP.mbc"][AP.mbc_1["AP.mbc"]>0.02]<-NA # remove additional X-axis outlier - lone value no where near any other
AP.mbc_1["AP.mbc"][AP.mbc_1["AP.mbc"]>4e-05 & AP.mbc_1["horizon"] == "A"]<-NA # remove additional X-axis outlier - lone value no where near any other
AP.mbc_1["TOC"][AP.mbc_1["TOC"]<5 & AP.mbc_1["horizon"] == "A"]<-NA 


AP.model_1 <- lm(TOC ~ AP.mbc + horizon + forest.type + NOy.NHx.combined.28.year + soil.type
                 , data = AP.mbc_1)

# check that outliers have been removed
ggplot(AP.model_1, aes(x= TOC , y=AP.mbc, group = horizon)) + 
  geom_boxplot()

ggplot(AP.model_1, aes(x = AP.mbc, y = TOC, colour = horizon))+
  geom_point(aes(shape=forest.type))

summary(AP.model_1)

# simple ANOVA
anova(AP.model_1)


## Assumptions
# The linearity assumption states that the general relationship between the 
# response and predictor variable should look like a straight line. 
# We can evaluate this assumption by constructing a residuals vs. fitted values plot.

## linearity - don't want a pattern. E.g., residuals get larger/smaller as X gets larger
plot(AP.model_1, which = 1, add.smooth = T)

## normality 
plot(AP.model_1, which = 2)

## variance
plot(AP.model_1, add.smooth = T, which = 3)


##### diagnostics - Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption


## plot ANCOVA

AP.ancova.model <- lm(TOC ~ AP.mbc*horizon*forest.type*NOy.NHx.combined.28.year*soil.type,
                      data = AP.mbc_1)

par(mfrow=c(2,2))
plot(AP.ancova.model)
par(mfrow=c(1,1))

## accounting for leverage - e.g., one extreme value shift the trend line

AP.high <- as.data.frame(cooks.distance(AP.ancova.model))
colnames(AP.high) <- 'x'
AP.high <- subset(AP.high, x>4/nrow(AP.ancova.model))
AP.mbc_1_1 <- subset(AP.mbc_1,!rownames(AP.mbc_1)%in% rownames(AP.high))

AP.ancova.model_1 <- lm(TOC ~ AP.mbc*horizon*forest.type*NOy.NHx.combined.28.year*soil.type,
                        data = AP.mbc_1)

par(mfrow=c(2,2))
plot(AP.ancova.model_1)
par(mfrow=c(1,1))



AP.mbc.TOC.sqrt_1 <- mutate(AP.mbc_1_1, sqrtTOC = sqrt(TOC))
AP.mbc.TOC.ancova.mod.sqrt <- lm(sqrtTOC ~ AP.mbc*horizon*forest.type*NOy.NHx.combined.28.year*soil.type,
                                 data = AP.mbc_1)

par(mfrow=c(2,2))
plot(AP.mbc.TOC.ancova.mod.sqrt)
par(mfrow=c(1,1))

### DO NOT use logged data - equal variance assumption not improved - 
# Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption

## ANOVA
AP.mbc=as.character(AP.ancova.model_1[["effects"]][["AP.mbc"]])
horizon=factor(AP.ancova.model_1[["effects"]][["horizonO"]])

summary(AP.ancova.model_1)
anova(AP.ancova.model_1)

# Tukey HSD

AP.ancova.model_1_aov <- aov(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*NOy.NHx.combined.28.year*soil.type,
                             data = AP.mbc_1)
AP.ancova.model_1_aov_tt <- TukeyHSD(AP.ancova.model_1_aov, which = 'horizon:forest.type') # p-values too small for tukey to recognize
AP.ancova.model_1_aov_tt                                                                                                 # values < 10-3 can't be easily distinguished wtih aov function

str(AP.ancova.model_1_aov_tt)

print(AP.ancova.model_1_aov_tt,digits=17)


### make a split plot

AP.mbc.plt <- ggplot(AP.ancova.model_1, aes(x = AP.mbc, y = TOC, colour = horizon))+
  geom_point(aes(shape=forest.type)) +stat_smooth(method="lm", se=T, colour = "black")+
  labs(shape = "Forest Type") + labs(colour = "Horizon") +
  labs(title="Effect of AP activity on SOC")+
  ylab("SOC (%)")+
  xlab(expression(paste('Enzyme activity (',mu,'mol g'^-1," microbial biomass C)")))+
  theme_classic()+
  facet_grid(. ~ horizon, scales='free')

AP.mbc.plt


## presenting results with logged axes for clear representation

AP.mbc.2=AP.mbc_1
AP.mbc.2["AP.mbc"][AP.mbc.2["AP.mbc"]<7.681545e-19]<-NA # removing Inf values

AP.ancova.model.2 <- lm(TOC ~ AP.mbc*horizon*forest.type*NOy.NHx.combined.28.year*soil.type,
                        data = AP.mbc.2)

AP.mbc.plt.log <- ggplot(AP.ancova.model.2, aes(x = AP.mbc, y = TOC, colour = horizon))+
  geom_point(aes(shape=forest.type)) +stat_smooth(method="lm", se=T, colour = "black")+
  labs(shape = "Forest Type") + labs(colour = "Horizon") +
  labs(title="Effect of AP activity on SOM",
       x ="ln(Enzyme activity)", y = "ln(SOM)")+
  scale_x_continuous(trans = scales::log_trans(),
                     breaks = scales::log_breaks())+
  scale_y_continuous(trans = scales::log_trans(),
                     breaks = scales::log_breaks())+
  theme_classic()

AP.mbc.plt.log



######### BG

## visualising data
ggplot(BG.mbc, aes(x = BG.mbc, y = TOC, colour = horizon, group = forest.type))+
  geom_point(aes(shape=forest.type))


# simple lm
BG.model <- lm(TOC ~ BG.mbc + horizon + forest.type + NOy.NHx.combined.28.year+ soil.type
               , data = BG.mbc)

print(BG.model)
summary(BG.model)


# remove outliers
outlierTest(BG.model,cutoff=Inf, TOC~1)

ggplot(BG.model, aes(x= TOC , y=BG.mbc, group = horizon)) + 
  geom_boxplot()

BG.mbc_1=BG.mbc

BG.mbc_1<-BG.mbc_1[-c(25,51,77,103,129,155,98,150,104,147,141,17,61,9,156,30,35, 159, 185, 221, 38, 64, 142, 7, 69, 95, 12, 151, 21, 47,73, 99,125),] 

BG.mbc_1["BG.mbc"][BG.mbc_1["BG.mbc"]>0.002& BG.mbc_1["horizon"] == "O"]<-NA # remove additional X-axis outlier - lone value no where near any other
BG.mbc_1["BG.mbc"][BG.mbc_1["BG.mbc"]>0.00012 & BG.mbc_1["horizon"] == "A"]<-NA
BG.mbc_1["BG.mbc"][BG.mbc_1["BG.mbc"]>5.0e-05 & BG.mbc_1["TOC"]>40 & BG.mbc_1["horizon"] == "A"]<-NA # random, unexplainable outlier


BG.model_1 <- lm(TOC ~ BG.mbc + horizon + forest.type + NOy.NHx.combined.28.year + soil.type
                 , data = BG.mbc_1)

# check that outliers have been removed
ggplot(BG.model_1, aes(x= TOC , y=BG.mbc, group = horizon)) + 
  geom_boxplot()

ggplot(BG.mbc_1, aes(x = BG.mbc, y = TOC, colour = horizon))+
  geom_point(aes(shape=forest.type))

# simple ANOVA
anova(BG.model_1)



## Assumptions
# The linearity assumption states that the general relationship between the 
# response and predictor variable should look like a straight line. 
# We can evaluate this assumption by constructing a residuals vs. fitted values plot.

## linearity - don't want a pattern. E.g., residuals get larger/smaller as X gets larger
par(mfrow=c(2,2))
plot(BG.model_1)
par(mfrow=c(1,1))


##### diagnostics - Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption


## plot ANCOVA

BG.ancova.model <- lm(TOC ~ BG.mbc *horizon *forest.type *NOy.NHx.combined.28.year *soil.type,
                      data = BG.mbc_1)

par(mfrow=c(2,2))
plot(BG.ancova.model)
par(mfrow=c(1,1))

## accounting for leverage - e.g., one extreme value shift the trend line

BG.high <- as.data.frame(cooks.distance(BG.ancova.model))
colnames(BG.high) <- 'x'
BG.high <- subset(BG.high, x>4/nrow(BG.ancova.model))
BG.mbc_1_1 <- subset(BG.mbc_1,!rownames(BG.mbc_1)%in% rownames(BG.high))

BG.ancova.model_1 <- lm(TOC ~ BG.mbc *horizon *forest.type *NOy.NHx.combined.28.year *soil.type,
                        data = BG.mbc_1_1)

par(mfrow=c(2,2))
plot(BG.ancova.model_1)
par(mfrow=c(1,1))


# transformation
BG.mbc.TOC.ln_1 <- mutate(BG.mbc_1_1, lnTOC = log(TOC))
BG.mbc.TOC.ancova.mod.log <- lm(lnTOC ~  BG.mbc *horizon *forest.type *NOy.NHx.combined.28.year *soil.type,
                                data = BG.mbc.TOC.ln_1)

par(mfrow=c(2,2))
plot(BG.mbc.TOC.ancova.mod.log)
par(mfrow=c(1,1))

BG.mbc.TOC.sqrt_1 <- mutate(BG.mbc_1_1, sqrtTOC = sqrt(TOC))
BG.mbc.TOC.ancova.mod.sqrt <- lm(sqrtTOC ~ BG.mbc *horizon *forest.type *NOy.NHx.combined.28.year *soil.type,
                                 data = BG.mbc.TOC.sqrt_1)

par(mfrow=c(2,2))
plot(BG.mbc.TOC.ancova.mod.sqrt)
par(mfrow=c(1,1))

### DO NOT use logged data - equal variance assumption not improved - 
# Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption

## ANOVA
BG.mbc=as.character(BG.ancova.model_1[["effects"]][["BG.mbc"]])
horizon=factor(BG.ancova.model_1[["effects"]][["horizonO"]])

summary(BG.ancova.model_1)
anova(BG.ancova.model_1)

# Tukey HSD

BG.ancova.model_1_aov <- aov(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*NOy.NHx.combined.28.year*soil.type,
                             data = BG.mbc_1)
BG.ancova.model_1_aov_tt <- TukeyHSD(BG.ancova.model_1_aov, which = 'horizon:forest.type') # p-values too small for tukey to recognize
BG.ancova.model_1_aov_tt                                                                                                 # values < 10-3 can't be easily distinguished wtih aov function

str(BG.ancova.model_1_aov_tt)

print(BG.ancova.model_1_aov_tt,digits=17)


### make a split plot

BG.mbc.plt <- ggplot(BG.ancova.model_1, aes(x = BG.mbc, y = TOC, colour = horizon))+
  geom_point(aes(shape=forest.type)) +stat_smooth(method="lm", se=T, colour = "black")+
  labs(shape = "Forest Type") + labs(colour = "Horizon") +
  labs(title="Effect of BG activity on SOC")+
  ylab("SOC (%)")+
  xlab(expression(paste('Enzyme activity (',mu,'mol g'^-1," microbial biomass C)")))+
  theme_classic()+
  facet_grid(. ~ horizon, scales='free')

BG.mbc.plt



######### BX

## visualising data
ggplot(BX.mbc, aes(x = BX.mbc, y = TOC, colour = horizon, group = forest.type))+
  geom_point(aes(shape=forest.type))



# simple lm
BX.model <- lm(TOC ~ BX.mbc + horizon + forest.type + NOy.NHx.combined.28.year + soil.type 
               , data = BX.mbc)

print(BX.model)
summary(BX.model)


# remove outliers
outlierTest(BX.model,cutoff=Inf, TOC~1)

ggplot(BX.model, aes(x= TOC , y=BX.mbc, group = horizon)) + 
  geom_boxplot()

BX.mbc_1=BX.mbc


BX.mbc_1<-BX.mbc_1[-c(25,51,77,103,129,155,150,98, 141, 104, 109, 59, 147,17, 183, 61, 221, 142,38, 
                      151, 12, 73, 63, 64, 116, 185, 35, 209, 21, 99, 125, 47),] # removed point 71 as it is skewing the distribution + 42 + 155 as it was causing high leverage

BX.mbc_1["BX.mbc"][BX.mbc_1["BX.mbc"]>0.004& BX.mbc_1["horizon"] == "O"]<-NA # remove additional X-axis outlier - lone value no where near any other
BX.mbc_1["BX.mbc"][BX.mbc_1["BX.mbc"]>1e-05 & BX.mbc_1["horizon"] == "A"]<-NA

BX.model_1 <- lm(TOC ~ BX.mbc + horizon + forest.type + NOy.NHx.combined.28.year + soil.type
                 , data = BX.mbc_1)


# check that outliers have been removed
ggplot(BX.model_1, aes(x= TOC , y=BX.mbc, group = horizon)) + 
  geom_boxplot()

# check that outliers have been removed
ggplot(BX.mbc_1, aes(x = BX.mbc, y = TOC, colour = horizon))+
  geom_point(aes(shape=forest.type))

# simple ANOVA
anova(BX.model_1)


## Assumptions
# The linearity assumption states that the general relationship between the 
# response and predictor variable should look like a straight line. 
# We can evaluate this assumption by constructing a residuals vs. fitted values plot.

## linearity - don't want a pattern. E.g., residuals get larger/smaller as X gets larger
par(mfrow=c(2,2))
plot(BX.model_1)
par(mfrow=c(1,1))


##### diagnostics - Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption


## plot ANCOVA

BX.ancova.model <- lm(TOC ~ BX.mbc * horizon * forest.type * NOy.NHx.combined.28.year * soil.type,
                      data = BX.mbc_1)

par(mfrow=c(2,2))
plot(BX.ancova.model)
par(mfrow=c(1,1))

## accounting for leverage - e.g., one extreme value shift the trend line

BX.high <- as.data.frame(cooks.distance(BX.ancova.model))
colnames(BX.high) <- 'x'
BX.high <- subset(BX.high, x>4/nrow(BX.ancova.model))
BX.mbc_1_1 <- subset(BX.mbc_1,!rownames(BX.mbc_1)%in% rownames(BX.high))

BX.ancova.model_1 <- lm(TOC ~ BX.mbc * horizon * forest.type * NOy.NHx.combined.28.year * soil.type,
                        data = BX.mbc_1_1)

par(mfrow=c(2,2))
plot(BX.ancova.model_1)
par(mfrow=c(1,1))


# transformation
BX.mbc.TOC.ln_1 <- mutate(BX.mbc_1_1, lnTOC = log(TOC))
BX.mbc.TOC.ancova.mod.log <- lm(lnTOC ~ BX.mbc * horizon * forest.type * NOy.NHx.combined.28.year * soil.type,
                                data = BX.mbc.TOC.ln_1)

par(mfrow=c(2,2))
plot(BX.mbc.TOC.ancova.mod.log)
par(mfrow=c(1,1))

BX.mbc.TOC.sqrt_1 <- mutate(BX.mbc_1_1, sqrtTOC = sqrt(TOC))
BX.mbc.TOC.ancova.mod.sqrt <- lm(sqrtTOC ~ BX.mbc * horizon * forest.type * NOy.NHx.combined.28.year * soil.type,
                                 data = BX.mbc.TOC.sqrt_1)

par(mfrow=c(2,2))
plot(BX.mbc.TOC.ancova.mod.log)
par(mfrow=c(1,1))

### DO NOT use logged data - equal variance assumption not improved - 
# Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption

## ANOVA
BX.mbc=as.character(BX.ancova.model_1[["effects"]][["BX.mbc"]])
horizon=factor(BX.ancova.model_1[["effects"]][["horizonO"]])

summary(BX.ancova.model_1)
anova(BX.ancova.model_1)

# Tukey HSD

BX.ancova.model_1_aov <- aov(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*NOy.NHx.combined.28.year*soil.type,
                             data = BX.mbc_1)
BX.ancova.model_1_aov_tt <- TukeyHSD(BX.ancova.model_1_aov, which = 'horizon:forest.type') # p-values too small for tukey to recognize
BX.ancova.model_1_aov_tt                                                                                                 # values < 10-3 can't be easily distinguished wtih aov function

str(BX.ancova.model_1_aov_tt)

print(BX.ancova.model_1_aov_tt,digits=17)


### make a split plot

BX.mbc.plt <- ggplot(BX.ancova.model_1, aes(x = BX.mbc, y = TOC, colour = horizon))+
  geom_point(aes(shape=forest.type)) +stat_smooth(method="lm", se=T, colour = "black")+
  labs(shape = "Forest Type") + labs(colour = "Horizon") +
  labs(title="Effect of BX activity on SOC")+
  ylab("SOC (%)")+
  xlab(expression(paste('Enzyme activity (',mu,'mol g'^-1," microbial biomass C)")))+
  theme_classic()+
  facet_grid(. ~ horizon, scales='free')

BX.mbc.plt



######### PPO ####

## visualising data
ggplot(PPO.mbc, aes(x = PPO.mbc, y = TOC, colour = horizon, group = forest.type))+
  geom_point(aes(shape=forest.type))



# simple lm
PPO.model <- lm(TOC ~ PPO.mbc + horizon + forest.type + NOy.NHx.combined.28.year + soil.type
               , data = PPO.mbc)

print(PPO.model)
summary(PPO.model)


# remove outliers
outlierTest(PPO.model,cutoff=Inf, TOC~1)

ggplot(PPO.model, aes(x= TOC , y=PPO.mbc, group = horizon)) + 
  geom_boxplot()


PPO.mbc_1=PPO.mbc

PPO.mbc_1 <- PPO.mbc[-c(208, 205, 209, 210, 206, 207,
                        160,162),] # also removing values that later on cause non-normality

PPO.mbc_1["PPO.mbc"][PPO.mbc_1["PPO.mbc"]>0.178 & PPO.mbc_1["horizon"] == "O"]<-NA # remove additional X-axis outlier - lone value no where near any other


PPO.model_1 <- lm(TOC ~ PPO.mbc + horizon + forest.type + NOy.NHx.combined.28.year + soil.type
                 , data = PPO.mbc_1)

# check that outliers have been removed
ggplot(PPO.model_1, aes(x= TOC , y=PPO.mbc, group = horizon)) + 
  geom_boxplot()

ggplot(PPO.model_1, aes(x = PPO.mbc, y = TOC, colour = horizon))+
  geom_point(aes(shape=forest.type))

summary(PPO.model_1)

# simple ANOVA
anova(PPO.model_1)


## Assumptions
# The linearity assumption states that the general relationship between the 
# response and predictor variable should look like a straight line. 
# We can evaluate this assumption by constructing a residuals vs. fitted values plot.

## linearity - don't want a pattern. E.g., residuals get larger/smaller as X gets larger
plot(PPO.model_1, which = 1, add.smooth = T)

## normality 
plot(PPO.model_1, which = 2)

## variance
plot(PPO.model_1, add.smooth = T, which = 3)


##### diagnostics - Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption


## plot ANCOVA

PPO.ancova.model <- lm(TOC ~ PPO.mbc*horizon*forest.type*NOy.NHx.combined.28.year*soil.type,
                      data = PPO.mbc_1)

par(mfrow=c(2,2))
plot(PPO.ancova.model)
par(mfrow=c(1,1))

## accounting for leverage - e.g., one extreme value shift the trend line

PPO.high <- as.data.frame(cooks.distance(PPO.ancova.model))
colnames(PPO.high) <- 'x'
PPO.high <- subset(PPO.high, x>4/nrow(PPO.ancova.model))
PPO.mbc_1_1 <- subset(PPO.mbc_1,!rownames(PPO.mbc_1)%in% rownames(PPO.high))

PPO.ancova.model_1 <- lm(TOC ~ PPO.mbc*horizon*forest.type*NOy.NHx.combined.28.year*soil.type*pH.H2O,
                        data = PPO.mbc_1_1)

par(mfrow=c(2,2))
plot(PPO.ancova.model_1)
par(mfrow=c(1,1))



PPO.mbc.TOC.ln_1 <- mutate(PPO.mbc_1_1, lnTOC = log(TOC))
PPO.mbc.TOC.ancova.mod.ln <- lm(lnTOC ~ PPO.mbc*horizon*forest.type*NOy.NHx.combined.28.year*soil.type*pH.H2O,
                                 data = PPO.mbc.TOC.ln_1)

par(mfrow=c(2,2))
plot(PPO.mbc.TOC.ancova.mod.ln)
par(mfrow=c(1,1))

### DO NOT use logged data - equal variance assumption not improved - 
# Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption

## ANOVA
PPO.mbc=as.character(PPO.ancova.model_1[["effects"]][["PPO.mbc"]])
horizon=factor(PPO.ancova.model_1[["effects"]][["horizonO"]])

summary(PPO.ancova.model_1)
anova(PPO.ancova.model_1)

# Tukey HSD

PPO.ancova.model_1_aov <- aov(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*NOy.NHx.combined.28.year*soil.type*pH.H2O,
                             data = PPO.mbc_1)
PPO.ancova.model_1_aov_tt <- TukeyHSD(PPO.ancova.model_1_aov, which = 'horizon:forest.type') # p-values too small for tukey to recognize
PPO.ancova.model_1_aov_tt                                                                                                 # values < 10-3 can't be easily distinguished wtih aov function

str(PPO.ancova.model_1_aov_tt)

print(PPO.ancova.model_1_aov_tt,digits=17)


### make a split plot

PPO.mbc.plt <- ggplot(PPO.ancova.model_1, aes(x = PPO.mbc, y = TOC, colour = horizon))+
  geom_point(aes(shape=forest.type)) +stat_smooth(method="lm", se=T, colour = "black")+
  labs(shape = "Forest Type") + labs(colour = "Horizon") +
  labs(title="Effect of PPO activity on SOC")+
  ylab("SOC (%)")+
  xlab(expression(paste('Enzyme activity (',mu,'mol g'^-1," microbial biomass C)")))+
  theme_classic()+
  facet_grid(. ~ horizon, scales='free')

PPO.mbc.plt


## presenting results with logged axes for clear representation

PPO.mbc.2=PPO.mbc_1
PPO.mbc.2["PPO.mbc"][PPO.mbc.2["PPO.mbc"]<7.681545e-19]<-NA # removing Inf values

PPO.ancova.model.2 <- lm(TOC ~ PPO.mbc*horizon*forest.type*NOy.NHx.combined.28.year*soil.type,
                        data = PPO.mbc.2)

PPO.mbc.plt.log <- ggplot(PPO.ancova.model.2, aes(x = PPO.mbc, y = TOC, colour = horizon))+
  geom_point(aes(shape=forest.type)) +stat_smooth(method="lm", se=T, colour = "black")+
  labs(shape = "Forest Type") + labs(colour = "Horizon") +
  labs(title="Effect of PPO activity on SOM",
       x ="ln(Enzyme activity)", y = "ln(SOM)")+
  scale_x_continuous(trans = scales::log_trans(),
                     breaks = scales::log_breaks())+
  scale_y_continuous(trans = scales::log_trans(),
                     breaks = scales::log_breaks())+
  theme_classic()

PPO.mbc.plt.log



##### NAG

## visualising data
ggplot(NAG.mbc, aes(x = NAG.mbc, y = TOC, colour = horizon, group = forest.type))+
  geom_point(aes(shape=forest.type))



# simple lm
NAG.model <- lm(TOC ~ NAG.mbc + horizon + forest.type + NOy.NHx.combined.28.year + soil.type 
                , data = NAG.mbc)

print(NAG.model)
summary(NAG.model)


# remove outliers
outlierTest(NAG.model,cutoff=Inf, TOC~1)

ggplot(NAG.model, aes(x= TOC , y=NAG.mbc, group = horizon)) + 
  geom_boxplot()

NAG.mbc_1=NAG.mbc

NAG.mbc_1<-NAG.mbc_1[-c(77,103,129,155,25,51,150,98, 147, 17, 104, 41, 61, 9, 67, 151, 35, 38, 142, 64, 141, 156, 122, 73, 120),] 

NAG.mbc_1["NAG.mbc"][NAG.mbc_1["NAG.mbc"]>0.004& NAG.mbc_1["horizon"] == "O"]<-NA # remove additional X-axis outlier - lone value no where near any other
NAG.mbc_1["NAG.mbc"][NAG.mbc_1["NAG.mbc"]>0.00002 & NAG.mbc_1["horizon"] == "A"]<-NA

NAG.model_1 <- lm(TOC ~ NAG.mbc + horizon + forest.type + NOy.NHx.combined.28.year 
                  , data = NAG.mbc_1)


# check that outliers have been removed
ggplot(NAG.model_1, aes(x= TOC , y=NAG.mbc, group = horizon)) + 
  geom_boxplot()

NAG.model_1 <- lm(TOC ~ NAG.mbc + horizon + forest.type + NOy.NHx.combined.28.year + soil.type
                  , data = NAG.mbc_1)

# check that outliers have been removed
ggplot(NAG.mbc_1, aes(x = NAG.mbc, y = TOC, colour = horizon))+
  geom_point(aes(shape=forest.type))

# simple ANOVA
anova(NAG.model_1)


## Assumptions
# The linearity assumption states that the general relationship between the 
# response and predictor variable should look like a straight line. 
# We can evaluate this assumption by constructing a residuals vs. fitted values plot.

## linearity - don't want a pattern. E.g., residuals get larger/smaller as X gets larger
par(mfrow=c(2,2))
plot(NAG.model_1)
par(mfrow=c(1,1))


##### diagnostics - Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption


## plot ANCOVA

NAG.ancova.model <- lm(TOC ~ NAG.mbc * horizon * forest.type * NOy.NHx.combined.28.year * soil.type,
                       data = NAG.mbc_1)

par(mfrow=c(2,2))
plot(NAG.ancova.model)
par(mfrow=c(1,1))

## accounting for leverage - e.g., one extreme value shift the trend line

NAG.high <- as.data.frame(cooks.distance(NAG.ancova.model))
colnames(NAG.high) <- 'x'
NAG.high <- subset(NAG.high, x>4/nrow(NAG.ancova.model))
NAG.mbc_1_1 <- subset(NAG.mbc_1,!rownames(NAG.mbc_1)%in% rownames(NAG.high))

NAG.ancova.model_1 <- lm(TOC ~ NAG.mbc * horizon * forest.type * NOy.NHx.combined.28.year * soil.type,
                         data = NAG.mbc_1_1)

par(mfrow=c(2,2))
plot(NAG.ancova.model_1)
par(mfrow=c(1,1))


# transformation
NAG.mbc.TOC.ln_1 <- mutate(NAG.mbc_1_1, lnTOC = log(TOC))
NAG.mbc.TOC.ancova.mod.log <- lm(lnTOC ~ NAG.mbc * horizon * forest.type * NOy.NHx.combined.28.year * soil.type,
                                 data = NAG.mbc.TOC.ln_1)

par(mfrow=c(2,2))
plot(NAG.mbc.TOC.ancova.mod.log)
par(mfrow=c(1,1))

NAG.mbc.TOC.sqrt_1 <- mutate(NAG.mbc_1_1, sqrtTOC = sqrt(TOC))
NAG.mbc.TOC.ancova.mod.sqrt <- lm(sqrtTOC ~ NAG.mbc * horizon * forest.type * NOy.NHx.combined.28.year * soil.type,
                                  data = NAG.mbc.TOC.sqrt_1)

par(mfrow=c(2,2))
plot(NAG.mbc.TOC.ancova.mod.log)
par(mfrow=c(1,1))

### DO NOT use logged data - equal variance assumption not improved - 
# Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption

## ANOVA
NAG.mbc=as.character(NAG.ancova.model_1[["effects"]][["NAG.mbc"]])
horizon=factor(NAG.ancova.model_1[["effects"]][["horizonO"]])

summary(NAG.ancova.model_1)
anova(NAG.ancova.model_1)

# Tukey HSD

NAG.ancova.model_1_aov <- aov(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*NOy.NHx.combined.28.year*soil.type,
                             data = NAG.mbc_1)
NAG.ancova.model_1_aov_tt <- TukeyHSD(NAG.ancova.model_1_aov, which = 'horizon:forest.type') # p-values too small for tukey to recognize
NAG.ancova.model_1_aov_tt                                                                                                 # values < 10-3 can't be easily distinguished wtih aov function

str(NAG.ancova.model_1_aov_tt)

print(NAG.ancova.model_1_aov_tt,digits=17)


### make a split plot

NAG.mbc.plt <- ggplot(NAG.ancova.model_1, aes(x = NAG.mbc, y = TOC, colour = horizon))+
  geom_point(aes(shape=forest.type)) +stat_smooth(method="lm", se=T, colour = "black")+
  labs(shape = "Forest Type") + labs(colour = "Horizon") +
  labs(title="Effect of NAG activity on SOC")+
  ylab("SOC (%)")+
  xlab(expression(paste('Enzyme activity (',mu,'mol g'^-1," microbial biomass C)")))+
  theme_classic()+
  facet_grid(. ~ horizon, scales='free')

NAG.mbc.plt


##### LAP

## visualising data
ggplot(LAP.mbc, aes(x = LAP.mbc, y = TOC, colour = horizon, group = forest.type))+
  geom_point(aes(shape=forest.type))



# simple lm
LAP.model <- lm(TOC ~ LAP.mbc + horizon + forest.type + NOy.NHx.combined.28.year + soil.type
                ,data = LAP.mbc)

print(LAP.model)
summary(LAP.model)


# remove outliers
outlierTest(LAP.model,cutoff=Inf, TOC~1)

ggplot(NAG.model, aes(x= TOC , y=NAG.mbc, group = horizon)) + 
  geom_boxplot()


LAP.mbc_1=LAP.mbc

LAP.mbc_1<-LAP.mbc_1[-c(145,146,147,148,149,150,118,120),] 

LAP.mbc_1["LAP.mbc"][LAP.mbc_1["LAP.mbc"] > 0.015] <- NA # remove additional X-axis outliers - lone values no where near any other
LAP.mbc_1["LAP.mbc"][LAP.mbc_1["LAP.mbc"] > 1e-04 & LAP.mbc_1["horizon"] == "A"]<-NA
LAP.mbc_1["LAP.mbc"][LAP.mbc_1["LAP.mbc"] > 0.0025 & LAP.mbc_1["horizon"] == "O"]<-NA
LAP.mbc_1["TOC"][LAP.mbc_1["TOC"] > 40 & LAP.mbc_1["horizon"] == "A"]<-NA


LAP.model_1 <- lm(TOC ~ LAP.mbc + horizon + forest.type + NOy.NHx.combined.28.year + soil.type, 
                  data = LAP.mbc_1)

# check that outliers have been removed
ggplot(LAP.model_1, aes(x = LAP.mbc, y = TOC, colour = horizon))+
  geom_point(aes(shape=forest.type))

# simple ANOVA
anova(LAP.model_1)


## Assumptions
# The linearity assumption states that the general relationship between the 
# response and predictor variable should look like a straight line. 
# We can evaluate this assumption by constructing a residuals vs. fitted values plot.

## linearity - don't want a pattern. E.g., residuals get larger/smaller as X gets larger
par(mfrow=c(2,2))
plot(LAP.model_1)
par(mfrow=c(1,1))


##### diagnostics - Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption


## plot ANCOVA

LAP.ancova.model <- lm(TOC ~ LAP.mbc*horizon*forest.type*NOy.NHx.combined.28.year*soil.type,
                       data = LAP.mbc_1)

par(mfrow=c(2,2))
plot(LAP.ancova.model)
par(mfrow=c(1,1))

## accounting for leverage - e.g., one extreme value shift the trend line

LAP.high <- as.data.frame(cooks.distance(LAP.ancova.model))
colnames(LAP.high) <- 'x'
LAP.high <- subset(LAP.high, x>4/nrow(LAP.ancova.model))
LAP.mbc_1_1 <- subset(LAP.mbc_1,!rownames(LAP.mbc_1)%in% rownames(LAP.high))

LAP.ancova.model_1 <- lm(TOC ~ LAP.mbc*horizon*forest.type*NOy.NHx.combined.28.year*soil.type,
                         data = LAP.mbc_1_1)

par(mfrow=c(2,2))
plot(LAP.ancova.model_1)
par(mfrow=c(1,1))


# transformation
LAP.mbc.TOC.ln_1 <- mutate(LAP.mbc_1_1, lnTOC = log(TOC))
LAP.mbc.TOC.ancova.mod.log <- lm(lnTOC ~ LAP.mbc*horizon*forest.type*NOy.NHx.combined.28.year*soil.type,
                                 data = LAP.mbc.TOC.ln_1)

par(mfrow=c(2,2))
plot(LAP.mbc.TOC.ancova.mod.log)
par(mfrow=c(1,1))

LAP.mbc.TOC.sqrt_1 <- mutate(LAP.mbc_1_1, sqrtTOC = sqrt(TOC))
LAP.mbc.TOC.ancova.mod.sqrt <- lm(sqrtTOC ~ LAP.mbc*horizon*forest.type*NOy.NHx.combined.28.year*soil.type,
                                  data = LAP.mbc.TOC.sqrt_1)

par(mfrow=c(2,2))
plot(LAP.mbc.TOC.ancova.mod.log)
par(mfrow=c(1,1))

### DO NOT use logged data - equal variance assumption not improved - 
# Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption

## ANOVA
LAP.mbc=as.character(LAP.ancova.model_1[["effects"]][["LAP.mbc"]])
horizon=factor(LAP.ancova.model_1[["effects"]][["horizonO"]])

summary(LAP.ancova.model_1)
anova(LAP.ancova.model_1)

# Tukey HSD

LAP.ancova.model_1_aov <- aov(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*NOy.NHx.combined.28.year*soil.type,
                              data = LAP.mbc_1)
LAP.ancova.model_1_aov_tt <- TukeyHSD(LAP.ancova.model_1_aov, which = 'horizon:forest.type') # p-values too small for tukey to recognize
LAP.ancova.model_1_aov_tt                                                                                                 # values < 10-3 can't be easily distinguished wtih aov function

str(LAP.ancova.model_1_aov_tt)

print(LAP.ancova.model_1_aov_tt,digits=17)


### make a split plot

LAP.mbc.plt <- ggplot(LAP.ancova.model_1, aes(x = LAP.mbc, y = TOC, colour = horizon))+
  geom_point(aes(shape=forest.type)) +stat_smooth(method="lm", se=T, colour = "black")+
  labs(shape = "Forest Type") + labs(colour = "Horizon") +
  labs(title="Effect of LAP activity on SOC")+
  ylab("SOC (%)")+
  xlab(expression(paste('Enzyme activity (',mu,'mol g'^-1," microbial biomass C)")))+
  theme_classic()+
  facet_grid(. ~ horizon, scales='free')

LAP.mbc.plt


#### make grid of plots
plot_grid(AP.mbc.plt,BG.mbc.plt, BX.mbc.plt, NAG.mbc.plt, LAP.mbc.plt, nrow = 3, labels = c("auto"), label_size = 10)



