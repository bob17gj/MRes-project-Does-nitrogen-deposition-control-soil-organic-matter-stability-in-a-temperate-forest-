library(tidyverse) #for all data wrangling
library(multcomp) # for ANCOVA post-hoc tests 
library(cowplot) #for manuscript ready figures
library(car)# for bonferroni outlier test
library(ggpmisc)
library(agricolae)


rm(list=ls())
setwd("~/General Waring/MRes Project")

hydrolytic_enzymes <- read.csv("Hydrolytic enzyme master sheet for R (2).csv")
str(hydrolytic_enzyames)

AP.enzyme <- filter(hydrolytic_enzymes, enzyme %in% c("AP"))
AP.enzyme

BG.enzyme <- filter(hydrolytic_enzymes, enzyme %in% c("BG"))
BG.enzyme

BX.enzyme <- filter(hydrolytic_enzymes, enzyme %in% c("BX"))
BX.enzyme

PPO.enzyme <- filter(hydrolytic_enzymes, enzyme %in% c("PPO"))
PPO.enzyme

LAP.enzyme <- filter(hydrolytic_enzymes, enzyme %in% c("LAP"))
LAP.enzyme

NAG.enzyme <- filter(hydrolytic_enzymes, enzyme %in% c("NAG"))
NAG.enzyme

# convert all "NaN values in DT to NA
AP.enzyme <- AP.enzyme %>% mutate_all(~ifelse(is.nan(.), NA, .))

BG.enzyme <- BG.enzyme %>% mutate_all(~ifelse(is.nan(.), NA, .))

BX.enzyme <- BX.enzyme %>% mutate_all(~ifelse(is.nan(.), NA, .))

PPO.enzyme <- PPO.enzyme %>% mutate_all(~ifelse(is.nan(.), NA, .))

NAG.enzyme <- NAG.enzyme %>% mutate_all(~ifelse(is.nan(.), NA, .))

LAP.enzyme <- LAP.enzyme %>% mutate_all(~ifelse(is.nan(.), NA, .))

# change "inf" values to NA
AP.enzyme <- AP.enzyme %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf)))

BG.enzyme <-BG.enzyme %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf)))

BX.enzyme <-BX.enzyme %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf)))

PPO.enzyme <-PPO.enzyme %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf)))

NAG.enzyme <-NAG.enzyme %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf)))

LAP.enzyme <-LAP.enzyme %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf)))

#### Per gram soil per hour

######### AP

# simple lm
AP.model <- lm(standardised.enzyme.activity.rate ~ NOy.NHx.combined.28.year*horizon*forest.type*soil.type*TOC,
               data = AP.enzyme)

print(AP.model)
summary(AP.model)


# remove outliers
outlierTest(AP.model,cutoff=Inf, standardised.enzyme.activity.rate~1)

ggplot(AP.model, aes(x= NOy.NHx.combined.28.year , y=standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

AP.enzyme_1=AP.enzyme


AP.enzyme_1 <- AP.enzyme_1[-c(24,61,49,207,
                              30,85,
                              45),] 

AP.enzyme_1["standardised.enzyme.activity.rate"][AP.enzyme_1["standardised.enzyme.activity.rate"]>1.6& AP.enzyme_1["horizon"] == "O"]<-NA # remove additional X-axis outlier - lone value no where near any other
AP.enzyme_1["standardised.enzyme.activity.rate"][AP.enzyme_1["standardised.enzyme.activity.rate"]>1.3 & AP.enzyme_1["horizon"] == "A"]<-NA
AP.enzyme_1["standardised.enzyme.activity.rate"][AP.enzyme_1["standardised.enzyme.activity.rate"]<0.1e-17]<-NA
AP.enzyme_1["standardised.enzyme.activity.rate"][AP.enzyme_1["standardised.enzyme.activity.rate"]>1.5& AP.enzyme_1["forest.type"] == "broadleaf"]<-NA # remove additional X-axis outlier - lone value no where near any other


AP.model_1 <- lm(standardised.enzyme.activity.rate ~ NOy.NHx.combined.28.year + horizon + forest.type + NOy.NHx.combined.28.year + soil.type + TOC
                 , data = AP.enzyme_1)


# check that outliers have been removed
ggplot(AP.model_1, aes(x= NOy.NHx.combined.28.year , y=standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

# check that outliers have been removed
ggplot(AP.enzyme_1, aes(x = NOy.NHx.combined.28.year , y = standardised.enzyme.activity.rate, colour = horizon))+
  geom_point(aes(shape=forest.type))

## ANOVA
print(AP.model_1)

anova(AP.model_1)

## Assumptions
# The linearity assumption states that the general relationship between the 
# response and predictor variable should look like a straight line. 
# We can evaluate this assumption by constructing a residuals vs. fitted values plot.

## linearity - don't want a pattern. E.g., residuals get larger/smaller as X gets larger
plot(AP.model, which = 1, add.smooth = T)

## normality 
plot(AP.model, which = 2)

## variance
plot(AP.model, add.smooth = T, which = 3)


##### diagnostics - use log transformed data - Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption


## plot ANCOVA

AP.ancova.model <- lm(standardised.enzyme.activity.rate ~ NOy.NHx.combined.28.year*horizon*forest.type*soil.type*TOC,
                      data = AP.enzyme_1)

par(mfrow=c(2,2))
plot(AP.ancova.model)
par(mfrow=c(1,1))

## accounting for leverage - e.g., one extreme value shift the trend line

AP.high <- as.data.frame(cooks.distance(AP.ancova.model))
colnames(AP.high) <- 'x'
AP.high <- subset(AP.high, x>4/nrow(AP.ancova.model))
AP.enzyme_1_1 <- subset(AP.enzyme_1,!rownames(AP.enzyme_1)%in% rownames(AP.high))

AP.ancova.model_1 <- lm(standardised.enzyme.activity.rate ~NOy.NHx.combined.28.year*horizon*forest.type*soil.type*TOC,
                        data = AP.enzyme_1_1)

par(mfrow=c(2,2))
plot(AP.ancova.model_1)
par(mfrow=c(1,1))


AP.deg.Ndepenzyme.TC.ln_1 <- mutate(AP.enzyme_1_1, lnstandardised.enzyme.activity.rate = log(standardised.enzyme.activity.rate))

AP.deg.Ndepenzyme.TC.ln_1 <- AP.deg.Ndepenzyme.TC.ln_1 %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf))) #convert Inf values to NA
AP.deg.Ndepenzyme.TC.ln_1 <- AP.deg.Ndepenzyme.TC.ln_1 %>% mutate_all(~ifelse(is.nan(.), NA, .)) # convert NaN values to NA



AP.deg.Ndepenzyme.TC.ancova.mod.log <- lm(lnstandardised.enzyme.activity.rate ~ horizon * forest.type * TOC * soil.type + NOy.NHx.combined.28.year,
                                          data = AP.deg.Ndepenzyme.TC.ln_1)

par(mfrow=c(2,2))
plot(AP.deg.Ndepenzyme.TC.ancova.mod.log)
par(mfrow=c(1,1))

### do not use logged data - equal variance assumption not improved though- 
# Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption

## ANOVA
horizon=factor(AP.ancova.model_1[["effects"]][["horizonO"]])

summary(AP.ancova.model_1)
anova(AP.ancova.model_1)

# tukey HSD - not needed - no interaction


## presenting results
AP.enzyme.plt <- ggplot(AP.ancova.model_1, aes(x = NOy.NHx.combined.28.year, y = standardised.enzyme.activity.rate))+
  stat_smooth(aes(group=1),method="lm", se=T, colour = "black")+
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")), size = 4) +
  geom_point(aes(shape=forest.type, colour = horizon), size = 2) +
  labs(shape = "Forest Type") + labs(colour = "Horizon") +
  labs(title="AP")+
  xlab(expression(paste('N depositon (kg ha'^-1," yr"^-1,")")))+ 
  ylab(expression(atop("Enzyme activity",paste(" (",mu,'mol g'^-1," soil h"^-1,")"))))+
  scale_x_continuous(breaks=c(0,5,10,15,20,25,30,35))+
  scale_y_continuous(breaks=c(0,0.3,0.6,0.9,1.2,1.5,1.8))+
  expand_limits(y=c(0,1.8), x =c(5,35))+
  scale_colour_manual(name="Horizon", # Legend label, use darker colors
                      breaks=c("A", "O"),
                      labels=c("A", "O"),
                      values=c("#CC9900","#990000"))+
  theme_classic() +   
  theme(axis.text=element_text(size = 12), axis.title=element_text(size = 12),legend.position="none",legend.title = element_text(size=12), legend.text=element_text(size=12))

AP.enzyme.plt



######### BG

# simple lm
BG.model <- lm(standardised.enzyme.activity.rate ~ NOy.NHx.combined.28.year*horizon*forest.type*soil.type*TOC,
               data = BG.enzyme)

print(BG.model)
summary(BG.model)


# remove outliers
outlierTest(BG.model,cutoff=Inf, standardised.enzyme.activity.rate~1)

ggplot(BG.model, aes(x= NOy.NHx.combined.28.year , y=standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

BG.enzyme_1=BG.enzyme


BG.enzyme_1 <- BG.enzyme[-c(61,101,48,187,86,126),] 

BG.enzyme_1["standardised.enzyme.activity.rate"][BG.enzyme_1["standardised.enzyme.activity.rate"]>0.5& BG.enzyme_1["horizon"] == "O"]<-NA # remove additional X-axis outlier - lone value no where near any other
BG.enzyme_1["standardised.enzyme.activity.rate"][BG.enzyme_1["standardised.enzyme.activity.rate"]>0.5 & BG.enzyme_1["horizon"] == "A"]<-NA
BG.enzyme_1["standardised.enzyme.activity.rate"][BG.enzyme_1["standardised.enzyme.activity.rate"]<0.1e-17]<-NA

BG.model_1 <- lm(standardised.enzyme.activity.rate ~ NOy.NHx.combined.28.year + horizon + forest.type + NOy.NHx.combined.28.year + soil.type + TOC
                 , data = BG.enzyme_1)


# check that outliers have been removed
ggplot(BG.model_1, aes(x= NOy.NHx.combined.28.year , y=standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

# check that outliers have been removed
ggplot(BG.enzyme_1, aes(x = NOy.NHx.combined.28.year , y = standardised.enzyme.activity.rate, colour = horizon))+
  geom_point(aes(shape=forest.type))

## ANOVA
print(BG.model_1)

anova(BG.model_1)

## Assumptions
# The linearity assumption states that the general relationship between the 
# response and predictor variable should look like a straight line. 
# We can evaluate this assumption by constructing a residuals vs. fitted values plot.

## linearity - don't want a pattern. E.g., residuals get larger/smaller as X gets larger
plot(BG.model, which = 1, add.smooth = T)

## normality 
plot(BG.model, which = 2)

## variance
plot(BG.model, add.smooth = T, which = 3)


##### diagnostics - use log transformed data - Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption


## plot ANCOVA

BG.ancova.model <- lm(standardised.enzyme.activity.rate ~ NOy.NHx.combined.28.year*horizon*forest.type*soil.type*TOC,
                      data = BG.enzyme_1)

par(mfrow=c(2,2))
plot(BG.ancova.model)
par(mfrow=c(1,1))

## accounting for leverage - e.g., one extreme value shift the trend line

BG.high <- as.data.frame(cooks.distance(BG.ancova.model))
colnames(BG.high) <- 'x'
BG.high <- subset(BG.high, x>4/nrow(BG.ancova.model))
BG.enzyme_1_1 <- subset(BG.enzyme_1,!rownames(BG.enzyme_1)%in% rownames(BG.high))

BG.ancova.model_1 <- lm(standardised.enzyme.activity.rate ~NOy.NHx.combined.28.year*horizon*forest.type*soil.type*TOC,
                        data = BG.enzyme_1_1)

par(mfrow=c(2,2))
plot(BG.ancova.model_1)
par(mfrow=c(1,1))

# transformation

BG.deg.Ndepenzyme.TC.ln_1 <- mutate(BG.enzyme_1_1, lnstandardised.enzyme.activity.rate = log(standardised.enzyme.activity.rate))

BG.deg.Ndepenzyme.TC.ln_1 <- BG.deg.Ndepenzyme.TC.ln_1 %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf))) #convert Inf values to NA
BG.deg.Ndepenzyme.TC.ln_1 <- BG.deg.Ndepenzyme.TC.ln_1 %>% mutate_all(~ifelse(is.nan(.), NA, .)) # convert NaN values to NA



BG.deg.Ndepenzyme.TC.ancova.mod.log <- lm(lnstandardised.enzyme.activity.rate ~ horizon * forest.type * TOC * soil.type + NOy.NHx.combined.28.year,
                                          data = BG.deg.Ndepenzyme.TC.ln_1)

par(mfrow=c(2,2))
plot(BG.deg.Ndepenzyme.TC.ancova.mod.log)
par(mfrow=c(1,1))

### DO NOT use logged data - equal variance assumption not improved - 
# Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption

## ANOVA
horizon=factor(BG.ancova.model_1[["effects"]][["horizonO"]])

summary(BG.ancova.model_1)
anova(BG.ancova.model_1)

# tukey HSD - not needed - no interaction


## presenting results
BG.enzyme.plt <- ggplot(BG.ancova.model_1, aes(x = NOy.NHx.combined.28.year, y = standardised.enzyme.activity.rate, colour = horizon))+
  geom_point(aes(shape=forest.type, colour=horizon), size = 2)+
  labs(title="BG")+
  xlab(expression(paste('N depositon (kg ha'^-1," yr"^-1,")")))+ 
  ylab(expression(atop("Enzyme activity",paste(" (",mu,'mol g'^-1," soil h"^-1,")"))))+
  scale_x_continuous(breaks=c(0,5,10,15,20,25,30,35))+
  expand_limits(x=c(5,35))+
  scale_y_continuous()+ 
  scale_colour_manual(name="Horizon", # Legend label, use darker colors
                      breaks=c("A", "O"),
                      labels=c("A", "O"),
                      values=c("#CC9900","#990000"))+
  scale_shape_discrete(labels=c("broadleaf","coniferous"),
                       name="Forest type")+
  theme_classic() +   
  theme(axis.text=element_text(size = 12), axis.title=element_text(size = 12),legend.position="none")

BG.enzyme.plt





#### BX # 

# simple lm
BX.model <- lm(standardised.enzyme.activity.rate ~ NOy.NHx.combined.28.year*horizon*forest.type*soil.type*TOC,
               data = BX.enzyme)

print(BX.model)
summary(BX.model)


# remove outliers
outlierTest(BX.model,cutoff=Inf, standardised.enzyme.activity.rate~1)

ggplot(BX.model, aes(x= NOy.NHx.combined.28.year , y=standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

BX.enzyme_1=BX.enzyme


BX.enzyme_1 <- BX.enzyme[-c(61,101,48,187,86,126,171,141),] 

BX.enzyme_1["standardised.enzyme.activity.rate"][BX.enzyme_1["standardised.enzyme.activity.rate"]>0.6& BX.enzyme_1["horizon"] == "O"]<-NA # remove additional X-axis outlier - lone value no where near any other
BX.enzyme_1["standardised.enzyme.activity.rate"][BX.enzyme_1["standardised.enzyme.activity.rate"]>0.75 & BX.enzyme_1["horizon"] == "A"]<-NA
BX.enzyme_1["standardised.enzyme.activity.rate"][BX.enzyme_1["standardised.enzyme.activity.rate"]<0.1e-17]<-NA

BX.model_1 <- lm(standardised.enzyme.activity.rate ~ NOy.NHx.combined.28.year + horizon + forest.type + NOy.NHx.combined.28.year + soil.type + TOC
                 , data = BX.enzyme_1)


# check that outliers have been removed
ggplot(BX.model_1, aes(x= NOy.NHx.combined.28.year , y=standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

# check that outliers have been removed
ggplot(BX.enzyme_1, aes(x = NOy.NHx.combined.28.year , y = standardised.enzyme.activity.rate, colour = horizon))+
  geom_point(aes(shape=forest.type))

## ANOVA
print(BX.model_1)

anova(BX.model_1)

## Assumptions
# The linearity assumption states that the general relationship between the 
# response and predictor variable should look like a straight line. 
# We can evaluate this assumption by constructing a residuals vs. fitted values plot.

## linearity - don't want a pattern. E.g., residuals get larger/smaller as X gets larger
plot(BX.model, which = 1, add.smooth = T)

## normality 
plot(BX.model, which = 2)

## variance
plot(BX.model, add.smooth = T, which = 3)


##### diagnostics - use log transformed data - Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption


## plot ANCOVA

BX.ancova.model <- lm(standardised.enzyme.activity.rate ~ NOy.NHx.combined.28.year*horizon*forest.type*soil.type*TOC,
                      data = BX.enzyme_1)

par(mfrow=c(2,2))
plot(BX.ancova.model)
par(mfrow=c(1,1))

## accounting for leverage - e.g., one extreme value shift the trend line

BX.high <- as.data.frame(cooks.distance(BX.ancova.model))
colnames(BX.high) <- 'x'
BX.high <- subset(BX.high, x>4/nrow(BX.ancova.model))
BX.enzyme_1_1 <- subset(BX.enzyme_1,!rownames(BX.enzyme_1)%in% rownames(BX.high))

BX.ancova.model_1 <- lm(standardised.enzyme.activity.rate ~NOy.NHx.combined.28.year*horizon*forest.type*soil.type*TOC,
                        data = BX.enzyme_1_1)

par(mfrow=c(2,2))
plot(BX.ancova.model_1)
par(mfrow=c(1,1))

# transformation

BX.deg.Ndepenzyme.TC.ln_1 <- mutate(BX.enzyme_1_1, lnstandardised.enzyme.activity.rate = log(standardised.enzyme.activity.rate))

BX.deg.Ndepenzyme.TC.ln_1 <- BX.deg.Ndepenzyme.TC.ln_1 %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf))) #convert Inf values to NA
BX.deg.Ndepenzyme.TC.ln_1 <- BX.deg.Ndepenzyme.TC.ln_1 %>% mutate_all(~ifelse(is.nan(.), NA, .)) # convert NaN values to NA



BX.deg.Ndepenzyme.TC.ancova.mod.log <- lm(lnstandardised.enzyme.activity.rate ~ horizon * forest.type * TOC * soil.type + NOy.NHx.combined.28.year,
                                          data = BX.deg.Ndepenzyme.TC.ln_1)

par(mfrow=c(2,2))
plot(BX.deg.Ndepenzyme.TC.ancova.mod.log)
par(mfrow=c(1,1))

### DO NOT use logged data - equal variance assumption not improved - 
# Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption

## ANOVA
horizon=factor(BX.ancova.model_1[["effects"]][["horizonO"]])

summary(BX.ancova.model_1)
anova(BX.ancova.model_1)

# tukey HSD - not needed - no interaction


## presenting results
BX.enzyme.plt <- ggplot(BX.ancova.model_1, aes(x = NOy.NHx.combined.28.year, y = standardised.enzyme.activity.rate, colour = horizon))+
  geom_point(aes(shape=forest.type), size = 2) +
  labs(shape = "Forest Type") + labs(colour = "Horizon") +
  labs(title="BX")+
  xlab(expression(paste('N depositon (kg ha'^-1," yr"^-1,")")))+
  ylab(expression(atop("Enzyme activity",paste(" (",mu,'mol g'^-1," soil h"^-1,")"))))+
  scale_x_continuous(breaks=c(0,5,10,15,20,25,30,35))+
  expand_limits(x=c(5,35))+
  scale_y_continuous()+
  expand_limits(y=c(0,0.8))+
  scale_colour_manual(name="Horizon", # Legend label, use darker colors
                      breaks=c("A", "O"),
                      labels=c("A", "O"),
                      values=c("#CC9900","#990000"))+
  scale_shape_discrete(labels=c("broadleaf","coniferous"),
                       name="Forest type")+
  theme_classic() +   
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12),
        legend.title = element_text(size=12), legend.text=element_text(size=12),
        legend.position = c(0.93, 0.72),
        legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="white"),
        legend.key.size = unit(0.1, "mm"))+
  guides(colour = guide_legend(override.aes = list(size=4))) # legend point size

BX.enzyme.plt



# PPO####

# simple lm
PPO.model <- lm(standardised.enzyme.activity.rate ~ NOy.NHx.combined.28.year*horizon*forest.type*soil.type*TOC,
                data = PPO.enzyme)

print(PPO.model)
summary(PPO.model)


# remove outliers
outlierTest(PPO.model,cutoff=Inf, standardised.enzyme.activity.rate~1)

ggplot(PPO.model, aes(x= NOy.NHx.combined.28.year , y=standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

PPO.enzyme_1=PPO.enzyme


PPO.enzyme_1 <- PPO.enzyme_1[-c(29,
                                175,
                                130,
                                220,158,166),] 

PPO.enzyme_1["standardised.enzyme.activity.rate"][PPO.enzyme_1["standardised.enzyme.activity.rate"]>0.0125& PPO.enzyme_1["horizon"] == "O"]<-NA # remove additional X-axis outlier - lone value no where near any other
#PPO.enzyme_1["standardised.enzyme.activity.rate"][PPO.enzyme_1["standardised.enzyme.activity.rate"]>1.3 & PPO.enzyme_1["horizon"] == "A"]<-NA
PPO.enzyme_1["standardised.enzyme.activity.rate"][PPO.enzyme_1["standardised.enzyme.activity.rate"]>0.028]<-NA


PPO.model_1 <- lm(standardised.enzyme.activity.rate ~ NOy.NHx.combined.28.year + horizon + forest.type + NOy.NHx.combined.28.year + soil.type + TOC
                  , data = PPO.enzyme_1)


# check that outliers have been removed
ggplot(PPO.model_1, aes(x= NOy.NHx.combined.28.year , y=standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

# check that outliers have been removed
ggplot(PPO.enzyme_1, aes(x = NOy.NHx.combined.28.year , y = standardised.enzyme.activity.rate, colour = horizon))+
  geom_point(aes(shape=forest.type))

## ANOVA
print(PPO.model_1)

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


##### diagnostics - use log transformed data - Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption


## plot ANCOVA

PPO.ancova.model <- lm(standardised.enzyme.activity.rate ~ NOy.NHx.combined.28.year*horizon*forest.type*soil.type*TOC,
                       data = PPO.enzyme_1)

par(mfrow=c(2,2))
plot(PPO.ancova.model)
par(mfrow=c(1,1))

## accounting for leverage - e.g., one extreme value shift the trend line

PPO.high <- as.data.frame(cooks.distance(PPO.ancova.model))
colnames(PPO.high) <- 'x'
PPO.high <- subset(PPO.high, x>4/nrow(PPO.ancova.model))
PPO.enzyme_1_1 <- subset(PPO.enzyme_1,!rownames(PPO.enzyme_1)%in% rownames(PPO.high))

PPO.ancova.model_1 <- lm(standardised.enzyme.activity.rate ~NOy.NHx.combined.28.year*horizon*forest.type*soil.type*TOC,
                         data = PPO.enzyme_1_1)

par(mfrow=c(2,2))
plot(PPO.ancova.model_1)
par(mfrow=c(1,1))


PPO.deg.Ndepenzyme.TC.ln_1 <- mutate(PPO.enzyme_1_1, lnstandardised.enzyme.activity.rate = log(standardised.enzyme.activity.rate))

PPO.deg.Ndepenzyme.TC.ln_1 <- PPO.deg.Ndepenzyme.TC.ln_1 %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf))) #convert Inf values to NA
PPO.deg.Ndepenzyme.TC.ln_1 <- PPO.deg.Ndepenzyme.TC.ln_1 %>% mutate_all(~ifelse(is.nan(.), NA, .)) # convert NaN values to NA



PPO.deg.Ndepenzyme.TC.ancova.mod.log <- lm(lnstandardised.enzyme.activity.rate ~ horizon * forest.type * TOC * soil.type + NOy.NHx.combined.28.year,
                                           data = PPO.deg.Ndepenzyme.TC.ln_1)

par(mfrow=c(2,2))
plot(PPO.deg.Ndepenzyme.TC.ancova.mod.log)
par(mfrow=c(1,1))

### do not use logged data - equal variance assumption not improved though- 
# Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption

## ANOVA
horizon=factor(PPO.ancova.model_1[["effects"]][["horizonO"]])

summary(PPO.ancova.model_1)
anova(PPO.ancova.model_1)

aov <- aov(PPO.ancova.model_1)

TukeyHSD(aov, which = 'forest.type')
HSD.test(aov, trt = c("forest.type"), console = TRUE)

PPO_sum.forest.type <- 
  PPO.enzyme_1_1 %>% 
  group_by(enzyme, forest.type) %>% 
  summarise(mean_activity = mean(standardised.enzyme.activity.rate, na.rm = T), n=n(), sd=sd(standardised.enzyme.activity.rate, na.rm = T), se=sd/sqrt(n))
PPO_sum.forest.type


## presenting PPO results ############################
PPO.enzyme.plt <- ggplot(PPO.ancova.model_1, aes(x = NOy.NHx.combined.28.year, y = standardised.enzyme.activity.rate))+
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")), size =4) +
  geom_point(aes(shape=forest.type, colour=horizon), size = 2)+
  stat_smooth(aes(group=1),method="lm", se=T, colour = "black")+
  labs(shape = "Forest Type") + labs(colour = "Horizon") +
  labs(title="PPO")+
  xlab(expression(paste('N depositon (kg ha'^-1," yr"^-1,")")))+ 
  ylab(expression(atop("Enzyme activity",paste(" (",mu,'mol g'^-1," soil h"^-1,")"))))+
  scale_x_continuous(breaks=c(0,5,10,15,20,25,30,35))+
  scale_y_continuous()+ 
  expand_limits(y=c(0,0.03), x =c(5,35))+
  scale_colour_manual(name="Horizon", # Legend label, use darker colors
                      breaks=c("A", "O"),
                      labels=c("A", "O"),
                      values=c("#CC9900","#990000"))+
  theme_classic() +   
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12),legend.position="none",legend.title = element_text(size=12), legend.text=element_text(size=12))

PPO.enzyme.plt





#### NAG #

# simple lm
NAG.model <- lm(standardised.enzyme.activity.rate ~ NOy.NHx.combined.28.year*horizon*forest.type*soil.type*TOC,
                data = NAG.enzyme)

print(NAG.model)
summary(NAG.model)


# remove outliers
outlierTest(NAG.model,cutoff=Inf, standardised.enzyme.activity.rate~1)

ggplot(NAG.model, aes(x= NOy.NHx.combined.28.year , y=standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

NAG.enzyme_1=NAG.enzyme


NAG.enzyme_1 <- NAG.enzyme[-c(61,168),] 

NAG.enzyme_1["standardised.enzyme.activity.rate"][NAG.enzyme_1["standardised.enzyme.activity.rate"]>0.8& NAG.enzyme_1["horizon"] == "O"]<-NA # remove additional X-axis outlier - lone value no where near any other
NAG.enzyme_1["standardised.enzyme.activity.rate"][NAG.enzyme_1["standardised.enzyme.activity.rate"]>1.3 & NAG.enzyme_1["horizon"] == "A"]<-NA
NAG.enzyme_1["standardised.enzyme.activity.rate"][NAG.enzyme_1["standardised.enzyme.activity.rate"]<0.1e-17]<-NA

NAG.model_1 <- lm(standardised.enzyme.activity.rate ~ NOy.NHx.combined.28.year + horizon + forest.type + NOy.NHx.combined.28.year + soil.type + TOC
                  , data = NAG.enzyme_1)


# check that outliers have been removed
ggplot(NAG.model_1, aes(x= NOy.NHx.combined.28.year , y=standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

# check that outliers have been removed
ggplot(NAG.enzyme_1, aes(x = NOy.NHx.combined.28.year , y = standardised.enzyme.activity.rate, colour = horizon))+
  geom_point(aes(shape=forest.type))

## ANOVA
print(NAG.model_1)

anova(NAG.model_1)

## Assumptions
# The linearity assumption states that the general relationship between the 
# response and predictor variable should look like a straight line. 
# We can evaluate this assumption by constructing a residuals vs. fitted values plot.

## linearity - don't want a pattern. E.g., residuals get larger/smaller as X gets larger
plot(NAG.model, which = 1, add.smooth = T)

## normality 
plot(NAG.model, which = 2)

## variance
plot(NAG.model, add.smooth = T, which = 3)


##### diagnostics - use log transformed data - Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption


## plot ANCOVA

NAG.ancova.model <- lm(standardised.enzyme.activity.rate ~ NOy.NHx.combined.28.year*horizon*forest.type*soil.type*TOC,
                       data = NAG.enzyme_1)

par(mfrow=c(2,2))
plot(NAG.ancova.model)
par(mfrow=c(1,1))

## accounting for leverage - e.g., one extreme value shift the trend line

NAG.high <- as.data.frame(cooks.distance(NAG.ancova.model))
colnames(NAG.high) <- 'x'
NAG.high <- subset(NAG.high, x>4/nrow(NAG.ancova.model))
NAG.enzyme_1_1 <- subset(NAG.enzyme_1,!rownames(NAG.enzyme_1)%in% rownames(NAG.high))

NAG.ancova.model_1 <- lm(standardised.enzyme.activity.rate ~NOy.NHx.combined.28.year*horizon*forest.type*soil.type*TOC,
                         data = NAG.enzyme_1_1)

par(mfrow=c(2,2))
plot(NAG.ancova.model_1)
par(mfrow=c(1,1))

# transformation

NAG.deg.Ndepenzyme.TC.ln_1 <- mutate(NAG.enzyme_1_1, lnstandardised.enzyme.activity.rate = log(standardised.enzyme.activity.rate))

NAG.deg.Ndepenzyme.TC.ln_1 <- NAG.deg.Ndepenzyme.TC.ln_1 %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf))) #convert Inf values to NA
NAG.deg.Ndepenzyme.TC.ln_1 <- NAG.deg.Ndepenzyme.TC.ln_1 %>% mutate_all(~ifelse(is.nan(.), NA, .)) # convert NaN values to NA



NAG.deg.Ndepenzyme.TC.ancova.mod.log <- lm(lnstandardised.enzyme.activity.rate ~ horizon * forest.type * TOC * soil.type + NOy.NHx.combined.28.year,
                                           data = NAG.deg.Ndepenzyme.TC.ln_1)

par(mfrow=c(2,2))
plot(NAG.deg.Ndepenzyme.TC.ancova.mod.log)
par(mfrow=c(1,1))

### DO NOT use logged data - equal variance assumption not improved - 
# Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption

## ANOVA
horizon=factor(NAG.ancova.model_1[["effects"]][["horizonO"]])

summary(NAG.ancova.model_1)
anova(NAG.ancova.model_1)

# tukey HSD - not needed - no interaction


## presenting results
NAG.enzyme.plt <- ggplot(NAG.ancova.model_1, aes(x = NOy.NHx.combined.28.year, y = standardised.enzyme.activity.rate, colour = horizon))+
  geom_point(aes(shape=forest.type), size = 2) +
  labs(shape = "Forest Type") + labs(colour = "Horizon") +
  labs(title="NAG")+
  xlab(expression(paste('N depositon (kg ha'^-1," yr"^-1,")")))+
  ylab(expression(atop("Enzyme activity",paste(" (",mu,'mol g'^-1," soil h"^-1,")"))))+
  scale_x_continuous(breaks=c(0,5,10,15,20,25,30,35))+
  expand_limits(x=c(5,35))+
  scale_y_continuous()+ 
  expand_limits(y=c(0,1.5))+
  scale_colour_manual(name="Horizon", # Legend label, use darker colors
                      breaks=c("A", "O"),
                      labels=c("A", "O"),
                      values=c("#CC9900","#990000"))+
  theme_classic() +   theme(axis.text=element_text(size=12), axis.title=element_text(size=12),legend.position="none",legend.title = element_text(size=12), legend.text=element_text(size=12))

NAG.enzyme.plt



#### LAP #

# simple lm
LAP.model <- lm(standardised.enzyme.activity.rate ~ NOy.NHx.combined.28.year*horizon*forest.type*soil.type*TOC,
                data = LAP.enzyme)

print(LAP.model)
summary(LAP.model)


# remove outliers
outlierTest(LAP.model,cutoff=Inf, standardised.enzyme.activity.rate~1)

ggplot(LAP.model, aes(x= NOy.NHx.combined.28.year , y=standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

LAP.enzyme_1=LAP.enzyme


LAP.enzyme_1 <- LAP.enzyme[-c(24,168),] 

LAP.enzyme_1["standardised.enzyme.activity.rate"][LAP.enzyme_1["standardised.enzyme.activity.rate"]>0.75 & LAP.enzyme_1["horizon"] == "O"]<-NA # remove additional X-axis outlier - lone value no where near any other
LAP.enzyme_1["standardised.enzyme.activity.rate"][LAP.enzyme_1["standardised.enzyme.activity.rate"]>0.75 & LAP.enzyme_1["horizon"] == "A"]<-NA
LAP.enzyme_1["standardised.enzyme.activity.rate"][LAP.enzyme_1["standardised.enzyme.activity.rate"]<0.1e-17]<-NA
LAP.enzyme_1["NOy.NHx.combined.28.year"][LAP.enzyme_1["NOy.NHx.combined.28.year"]<10]<-NA


LAP.model_1 <- lm(standardised.enzyme.activity.rate ~ NOy.NHx.combined.28.year + horizon + forest.type + NOy.NHx.combined.28.year + soil.type + TOC
                  , data = LAP.enzyme_1)


# check that outliers have been removed
ggplot(LAP.model_1, aes(x= NOy.NHx.combined.28.year , y=standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

# check that outliers have been removed
ggplot(LAP.enzyme_1, aes(x = NOy.NHx.combined.28.year , y = standardised.enzyme.activity.rate, colour = horizon))+
  geom_point(aes(shape=forest.type))

## ANOVA
print(LAP.model_1)

anova(LAP.model_1)

## Assumptions
# The linearity assumption states that the general relationship between the 
# response and predictor variable should look like a straight line. 
# We can evaluate this assumption by constructing a residuals vs. fitted values plot.

## linearity - don't want a pattern. E.g., residuals get larger/smaller as X gets larger
plot(LAP.model, which = 1, add.smooth = T)

## normality 
plot(LAP.model, which = 2)

## variance
plot(LAP.model, add.smooth = T, which = 3)


##### diagnostics - use log transformed data - Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption


## plot ANCOVA

LAP.ancova.model <- lm(standardised.enzyme.activity.rate ~ NOy.NHx.combined.28.year*horizon*forest.type*soil.type*TOC,
                       data = LAP.enzyme_1)

par(mfrow=c(2,2))
plot(LAP.ancova.model)
par(mfrow=c(1,1))

## accounting for leverage - e.g., one extreme value shift the trend line

LAP.high <- as.data.frame(cooks.distance(LAP.ancova.model))
colnames(LAP.high) <- 'x'
LAP.high <- subset(LAP.high, x>4/nrow(LAP.ancova.model))
LAP.enzyme_1_1 <- subset(LAP.enzyme_1,!rownames(LAP.enzyme_1)%in% rownames(LAP.high))

LAP.ancova.model_1 <- lm(standardised.enzyme.activity.rate ~NOy.NHx.combined.28.year*horizon*forest.type*soil.type*TOC,
                         data = LAP.enzyme_1_1)

par(mfrow=c(2,2))
plot(LAP.ancova.model_1)
par(mfrow=c(1,1))

# transformation

LAP.deg.Ndepenzyme.TC.ln_1 <- mutate(LAP.enzyme_1_1, lnstandardised.enzyme.activity.rate = log(standardised.enzyme.activity.rate))

LAP.deg.Ndepenzyme.TC.ln_1 <- LAP.deg.Ndepenzyme.TC.ln_1 %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf))) #convert Inf values to NA
LAP.deg.Ndepenzyme.TC.ln_1 <- LAP.deg.Ndepenzyme.TC.ln_1 %>% mutate_all(~ifelse(is.nan(.), NA, .)) # convert NaN values to NA



LAP.deg.Ndepenzyme.TC.ancova.mod.log <- lm(lnstandardised.enzyme.activity.rate ~ horizon * forest.type * TOC * soil.type + NOy.NHx.combined.28.year,
                                           data = LAP.deg.Ndepenzyme.TC.ln_1)

par(mfrow=c(2,2))
plot(LAP.deg.Ndepenzyme.TC.ancova.mod.log)
par(mfrow=c(1,1))

### DO NOT use logged data - equal variance assumption not improved - 
# Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption

## ANOVA
horizon=factor(LAP.ancova.model_1[["effects"]][["horizonO"]])

summary(LAP.ancova.model_1)
anova(LAP.ancova.model_1)

# tukey HSD - not needed - no interaction


## presenting results
LAP.enzyme.plt <- ggplot(LAP.ancova.model_1, aes(x = NOy.NHx.combined.28.year, y = standardised.enzyme.activity.rate))+
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")), size = 4) +
  geom_point(aes(shape=forest.type, colour=horizon), size = 2)+
  stat_smooth(aes(group=1),method="lm", se=T, colour = "black")+
  labs(shape = "Forest Type") + labs(colour = "Horizon") +
  labs(title="LAP")+
  xlab(expression(paste('N depositon (kg ha'^-1," yr"^-1,")")))+ 
  ylab(expression(atop("Enzyme activity",paste(" (",mu,'mol g'^-1," soil h"^-1,")"))))+
  scale_x_continuous(breaks=c(10,15,20,25,30,35))+
  expand_limits(x=c(10,35))+
  scale_y_continuous()+ 
  expand_limits(y=c(0,0.8))+
  scale_colour_manual(name="Horizon", 
                      breaks=c("A", "O"),
                      labels=c("A", "O"),
                      values=c("#CC9900","#990000"))+
  theme_classic() +   
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12),legend.position="none",legend.title = element_text(size=12), legend.text=element_text(size=12))

LAP.enzyme.plt


### make grid of plots

plot_grid(BG.enzyme.plt, BX.enzyme.plt, NAG.enzyme.plt, LAP.enzyme.plt, AP.enzyme.plt, PPO.enzyme.plt, 
          nrow = 3, labels = c("a   ", "b    ","c    ", "d    ", "e    ", "f    "), label_size = 14)


ggsave("EEA~pH.pdf", height = 12, width = 13)
