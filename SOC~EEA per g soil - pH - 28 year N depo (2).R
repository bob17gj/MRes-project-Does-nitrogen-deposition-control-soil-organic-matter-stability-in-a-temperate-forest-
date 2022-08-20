library(tidyverse) #for all data wrangling
library(multcomp) # for ANCOVA post-hoc tests 
library(cowplot) #for manuscript ready figures
library(car)# for bonferroni outlier test
library(ggpmisc) # for regression lable


rm(list=ls())
setwd("~/General Waring/Project")

hydrolytic_enzymes <- read.csv("Hydrolytic enzyme master sheet for R (2).csv")
str(hydrolytic_enzymes)

AP.enzyme <- filter(hydrolytic_enzymes, enzyme %in% c("AP"))
AP.enzyme

BG.enzyme <- filter(hydrolytic_enzymes, enzyme %in% c("BG"))
BG.enzymecx c 

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

## visualising data
ggplot(AP.enzyme, aes(x = standardised.enzyme.activity.rate, y = TOC, colour = horizon, group = forest.type))+
  geom_point(aes(shape=forest.type))

# simple lm
AP.model <- lm(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*NOy.NHx.combined.28.year*soil.type*pH.H2O,
               data = AP.enzyme)

print(AP.model)
summary(AP.model)


# remove outliers
outlierTest(AP.model,cutoff=Inf, TOC~1)

ggplot(AP.model, aes(x= TOC , y=standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

AP.enzyme_1=AP.enzyme


AP.enzyme_1 <- AP.enzyme[-c(129,155,25,134, 82, 73, 151,99, 21, 108, 30, 125, 47, 4, 57),] 

AP.enzyme_1["standardised.enzyme.activity.rate"][AP.enzyme_1["standardised.enzyme.activity.rate"]>1.6& AP.enzyme_1["horizon"] == "O"]<-NA # remove additional X-axis outlier - lone value no where near any other
AP.enzyme_1["standardised.enzyme.activity.rate"][AP.enzyme_1["standardised.enzyme.activity.rate"]>1.7 & AP.enzyme_1["horizon"] == "A"]<-NA

AP.model_1 <- lm(TOC ~ standardised.enzyme.activity.rate + horizon + forest.type + NOy.NHx.combined.28.year + soil.type
                 , data = AP.enzyme_1)


# check that outliers have been removed
ggplot(AP.model_1, aes(x= TOC , y=standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

# check that outliers have been removed
ggplot(AP.enzyme_1, aes(x = standardised.enzyme.activity.rate, y = TOC, colour = horizon))+
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

AP.ancova.model <- lm(TOC ~standardised.enzyme.activity.rate*horizon*forest.type*soil.type+NOy.NHx.combined.28.year+pH.H2O,
                      data = AP.enzyme_1)

par(mfrow=c(2,2))
plot(AP.ancova.model)
par(mfrow=c(1,1))

## accounting for leverage - e.g., one extreme value shift the trend line

AP.high <- as.data.frame(cooks.distance(AP.ancova.model))
colnames(AP.high) <- 'x'
AP.high <- subset(AP.high, x>4/nrow(AP.ancova.model))
AP.enzyme_1_1 <- subset(AP.enzyme_1,!rownames(AP.enzyme_1)%in% rownames(AP.high))

AP.ancova.model_1 <- lm(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*soil.type+NOy.NHx.combined.28.year+pH.H2O,
                        data = AP.enzyme_1_1)

par(mfrow=c(2,2))
plot(AP.ancova.model_1)
par(mfrow=c(1,1))

# transformation
AP.enzyme.TOC.ln_1 <- mutate(AP.enzyme_1, lnTOC = log(TOC))
AP.enzyme.TOC.ancova.mod.log <- lm(lnTOC ~ standardised.enzyme.activity.rate*horizon*forest.type*soil.type+NOy.NHx.combined.28.year+pH.H2O,
                                   data = AP.enzyme.TOC.ln_1)

par(mfrow=c(2,2))
plot(AP.enzyme.TOC.ancova.mod.log)
par(mfrow=c(1,1))

### DO NOT use logged data - equal variance assumption not improved - 
# Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption

## ANOVA
standardised.enzyme.activity.rate=as.character(AP.ancova.model_1[["effects"]][["standardised.enzyme.activity.rate"]])
horizon=factor(AP.ancova.model_1[["effects"]][["horizonO"]])

summary(AP.ancova.model_1)
anova(AP.ancova.model_1)

# Tukey HSD

AP.ancova.model_1_aov <- aov(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*soil.type+NOy.NHx.combined.28.year+pH.H2O,
                             data = AP.enzyme_1_1)
AP.ancova.model_1_aov_tt <- TukeyHSD(AP.ancova.model_1_aov, which = 'horizon:forest.type') # p-values too small for tukey to recognize
AP.ancova.model_1_aov_tt                                                                                                 # values < 10-3 can't be easily distinguished wtih aov function

str(AP.ancova.model_1_aov_tt)

print(AP.ancova.model_1_aov_tt,digits=17)

## presenting results
AP.enzyme.plt <- ggplot(AP.ancova.model_1, aes(x = standardised.enzyme.activity.rate, y = TOC))+
 # stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                # after_stat(rr.label), sep = "*\", \"*")), size = 3) +
  geom_point(aes(shape=forest.type, colour = horizon), size = 2) +
#  stat_smooth(method="lm", se=T, colour = "black")+
  labs(shape = "Forest Type") + labs(colour = "Horizon") +
  labs(title="AP ")+
  ylab("SOC (%)")+
  xlab(expression(paste('Enzyme activity (',mu,'mol g'^-1," soil h"^-1,")")))+
  scale_y_continuous(limits=c(0, 60))+
  scale_x_continuous(limits=c(0, 2))+
  scale_colour_manual(name="Horizon", # Legend label, use darker colors
                      breaks=c("A", "O"),
                      labels=c("A", "O"),
                      values=c("#CC9900","#990000"))+
  theme_classic() +   theme(legend.position="none",axis.text=element_text(size=15), axis.title=element_text(size=20))

AP.enzyme.plt



######### BG

## visualising data
ggplot(BG.enzyme, aes(x = standardised.enzyme.activity.rate, y = TOC, colour = horizon, group = forest.type))+
  geom_point(aes(shape=forest.type))

# simple lm
BG.model <- lm(TOC ~ standardised.enzyme.activity.rate + horizon + forest.type + NOy.NHx.combined.28.year + soil.type
               , data = BG.enzyme)

print(BG.model)
summary(BG.model)


# remove outliers
outlierTest(BG.model,cutoff=Inf, TOC~1)

ggplot(BG.model, aes(x= TOC , y=standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

BG.enzyme_1=BG.enzyme


BG.enzyme_1 <- BG.enzyme[-c(151,134,46, 137, 189, 25,12,162,207, 18, 56, 150, 47, 73, 99, 21, 125, 98),] # removing points to improve normality


BG.enzyme_1["standardised.enzyme.activity.rate"][BG.enzyme_1["standardised.enzyme.activity.rate"]>1 & BG.enzyme_1["horizon"] == "O"]<-NA # remove additional X-axis outlier - lone value no where near any other
BG.enzyme_1["standardised.enzyme.activity.rate"][BG.enzyme_1["standardised.enzyme.activity.rate"]>0.55 & BG.enzyme_1["horizon"] == "A"]<-NA

BG.model_1 <- lm(TOC ~ standardised.enzyme.activity.rate + horizon + forest.type + NOy.NHx.combined.28.year +soil.type 
                 , data = BG.enzyme_1)


# check that outliers have been removed
ggplot(BG.model_1, aes(x= TOC , y=standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

# check that outliers have been removed
ggplot(BG.enzyme_1, aes(x = standardised.enzyme.activity.rate, y = TOC, colour = horizon))+
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


##### diagnostics - do not use log transformed data - Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption


## plot ANCOVA

BG.ancova.model <- lm(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*soil.type+NOy.NHx.combined.28.year+pH.H2O,
                      data = BG.enzyme_1)

par(mfrow=c(2,2))
plot(BG.ancova.model)
par(mfrow=c(1,1))

## accounting for leverage - e.g., one extreme value shift the trend line

BG.high <- as.data.frame(cooks.distance(BG.ancova.model))
colnames(BG.high) <- 'x'
BG.high <- subset(BG.high, x>4/nrow(BG.ancova.model))
BG.enzyme_1_1 <- subset(BG.enzyme_1,!rownames(BG.enzyme_1)%in% rownames(BG.high))

BG.ancova.model_1 <- lm(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*soil.type+NOy.NHx.combined.28.year+pH.H2O,
                        data = BG.enzyme_1_1)

par(mfrow=c(2,2))
plot(BG.ancova.model_1)
par(mfrow=c(1,1))

# transformation
BG.enzyme.TOC.ln_1 <- mutate(BG.enzyme_1, lnTOC = log(TOC))
BG.enzyme.TOC.ancova.mod.log <- lm(lnTOC ~ standardised.enzyme.activity.rate*horizon*forest.type*soil.type+NOy.NHx.combined.28.year+pH.H2O,
                                   data = BG.enzyme.TOC.ln_1)

par(mfrow=c(2,2))
plot(BG.enzyme.TOC.ancova.mod.log)
par(mfrow=c(1,1))


BG.enzyme.TOC.sqrt_1 <- mutate(BG.enzyme_1, sqrtTOC = sqrt(TOC))
BG.enzyme.TOC.ancova.mod.sqrt <- lm(sqrtTOC ~ standardised.enzyme.activity.rate*horizon*forest.type*soil.type+NOy.NHx.combined.28.year+pH.H2O,
                                    data = BG.enzyme.TOC.sqrt_1)

par(mfrow=c(2,2))
plot(BG.enzyme.TOC.ancova.mod.sqrt)
par(mfrow=c(1,1))

### do not use logged data - equal variance assumption not improved - 
# Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption

## ANOVA
standardised.enzyme.activity.rate=as.character(BG.ancova.model_1[["effects"]][["standardised.enzyme.activity.rate"]])
horizon=factor(BG.ancova.model_1[["effects"]][["horizonO"]])

summary(BG.ancova.model_1)
anova(BG.ancova.model_1)


# Tukey HSD

BG.ancova.model_1_aov <- aov(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*soil.type+NOy.NHx.combined.28.year+pH.H2O,
                             data = BG.enzyme_1_1)
                                                                                                 # values < 10-3 can't be easily distinguished wtih aov function

TukeyHSD(BG.ancova.model_1_aov, which = 'standardised.enzyme.activity.rate:horizon')
HSD.test(BG.ancova.model_1_aov, trt = c("standardised.enzyme.activity.rate:horizon"), console = TRUE)

str(BG.ancova.model_1_aov_tt)

print(BG.ancova.model_1_aov_tt,digits=17)

## presenting results
BG.enzyme.plt <- ggplot(BG.ancova.model_1, aes(x = standardised.enzyme.activity.rate, y = TOC, colour = horizon))+
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")), size = 5) +
  geom_point(aes(shape=forest.type, colour = horizon), size = 2) +
  stat_smooth(method="lm", se=T)+
  labs(shape = "Forest Type") + labs(colour = "Horizon") +
  labs(title="BG ") +
  ylab("SOC (%)")+
  xlab(expression(paste("")))+
  scale_y_continuous(limits=c(0, 80))+
  scale_x_continuous(limits=c(0, 0.8))+
  scale_colour_manual(name="Horizon", # Legend label, use darker colors
                      breaks=c("A", "O"),
                      labels=c("A", "O"),
                      values=c("#CC9900","#990000"))+
  theme_classic()+   
  theme(legend.position="none",axis.text=element_text(size=15), axis.title=element_text(size=20))

BG.enzyme.plt


######### BX

## visualising data
ggplot(BX.enzyme, aes(x = standardised.enzyme.activity.rate, y = TOC, colour = horizon, group = forest.type))+
  geom_point(aes(shape=forest.type))

# simple lm
BX.model <- lm(TOC ~ standardised.enzyme.activity.rate + horizon + forest.type + NOy.NHx.combined.28.year +soil.type 
               , data = BX.enzyme)

print(BX.model)
summary(BX.model)


# remove outliers
outlierTest(BX.model,cutoff=Inf, TOC~1)

ggplot(BX.model, aes(x= TOC , y=standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

BX.enzyme_1=BX.enzyme


BX.enzyme_1 <- BX.enzyme[-c(19,108,116,118,151, 148, 4 ,137, 130, 139, 61),] 


BX.enzyme_1["standardised.enzyme.activity.rate"][BX.enzyme_1["standardised.enzyme.activity.rate"]>0.15 & BX.enzyme_1["horizon"] == "O"]<-NA # remove additional X-axis outlier - lone value no where near any other
BX.enzyme_1["standardised.enzyme.activity.rate"][BX.enzyme_1["standardised.enzyme.activity.rate"]>0.29 & BX.enzyme_1["horizon"] == "A"]<-NA
BX.enzyme_1["TOC"][BX.enzyme_1["TOC"]>40 & BX.enzyme_1["horizon"] == "A"]<-NA


BX.model_1 <- lm(TOC ~ standardised.enzyme.activity.rate + horizon + forest.type + NOy.NHx.combined.28.year + soil.type
                 , data = BX.enzyme_1)


# check that outliers have been removed
ggplot(BX.model_1, aes(x= TOC , y=standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

# check that outliers have been removed
ggplot(BX.enzyme_1, aes(x = standardised.enzyme.activity.rate, y = TOC, colour = horizon))+
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

BX.ancova.model <- lm(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*soil.type+NOy.NHx.combined.28.year+pH.H2O,
                      data = BX.enzyme_1)

par(mfrow=c(2,2))
plot(BX.ancova.model)
par(mfrow=c(1,1))

## accounting for leverage - e.g., one extreme value shift the trend line

BX.high <- as.data.frame(cooks.distance(BX.ancova.model))
colnames(BX.high) <- 'x'
BX.high <- subset(BX.high, x>4/nrow(BX.ancova.model))
BX.enzyme_1_1 <- subset(BX.enzyme_1,!rownames(BX.enzyme_1)%in% rownames(BX.high))


BX.ancova.model_1 <- lm(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*soil.type+NOy.NHx.combined.28.year+pH.H2O,
                        data = BX.enzyme_1_1)

par(mfrow=c(2,2))
plot(BX.ancova.model_1)
par(mfrow=c(1,1))

# transformation
BX.enzyme.TOC.ln_1 <- mutate(BX.enzyme_1, lnTOC = log(TOC))
BX.enzyme.TOC.ancova.mod.log <- lm(lnTOC ~ standardised.enzyme.activity.rate + 
                                     horizon + forest.type + horizon:standardised.enzyme.activity.rate + NOy.NHx.combined.28.year +
                                     horizon:forest.type + forest.type:standardised.enzyme.activity.rate+ 
                                     NOy.NHx.combined.28.year:standardised.enzyme.activity.rate+
                                     NOy.NHx.combined.28.year:horizon+ 
                                     NOy.NHx.combined.28.year:forest.type + 
                                     NOy.NHx.combined.28.year:standardised.enzyme.activity.rate:horizon+
                                     NOy.NHx.combined.28.year:standardised.enzyme.activity.rate:forest.type+
                                     NOy.NHx.combined.28.year:forest.type:horizon+
                                     NOy.NHx.combined.28.year:standardised.enzyme.activity.rate:horizon:forest.type,
                                   data = BX.enzyme.TOC.ln_1)

par(mfrow=c(2,2))
plot(BX.enzyme.TOC.ancova.mod.log)
par(mfrow=c(1,1))

### DO NOT use transformed data - equal variance assumption not improved - 
# Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption

summary(BX.ancova.model_1)
anova(BX.ancova.model_1)

# tukey HSD

BX.ancova.model_1_aov <- aov(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*soil.type+NOy.NHx.combined.28.year+pH.H2O,
                             data = BX.enzyme_1_1)
BX.ancova.model_1_aov_tt <- TukeyHSD(BX.ancova.model_1_aov, which = 'horizon:forest.type') # p-values too small for tukey to recognize
BX.ancova.model_1_aov_tt                                                                                                 # values < 10-3 can't be easily distinguished wtih aov function

str(BX.ancova.model_1_aov_tt)

print(BX.ancova.model_1_aov_tt,digits=17)

## presenting results
BX.enzyme.plt <- ggplot(BX.ancova.model_1, aes(x = standardised.enzyme.activity.rate, y = TOC))+
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")), size = 5) +
  geom_point(aes(shape=forest.type, colour = horizon), size = 2) +
  stat_smooth(method="lm", se=T, colour = "black")+
  labs(shape = "Forest Type") + labs(colour = "Horizon") +
  labs(title="BX ")+
  ylab("")+
  xlab(expression(paste("")))+
  scale_y_continuous(limits=c(0, 70))+
  scale_x_continuous(limits=c(0, 0.3))+
  scale_colour_manual(name="Horizon", # Legend label, use darker colors
                      breaks=c("A", "O"),
                      labels=c("A", "O"),
                      values=c("#CC9900","#990000"))+
  scale_shape_discrete(labels=c("broadleaf","coniferous"),
                       name="Forest type")+
  theme_classic()+
  theme(legend.position="none",axis.text=element_text(size=15), axis.title=element_text(size=20))
      
         # legend.position = c(0.85, 0.75), 
       # legend.background = element_rect(size=0.5, linetype="solid", 
                                        # colour ="white"),
        #legend.text=element_text(size=15),
      #  legend.key.size = unit(0.3, 'cm'))+ #change legend key size
  # guides(colour = guide_legend(override.aes = list(size=4))) # legend point size
                                                                                                                                        

BX.enzyme.plt


######### PPO ####

## visualising data
ggplot(PPO.enzyme, aes(x = standardised.enzyme.activity.rate, y = TOC, colour = horizon, group = forest.type))+
  geom_point(aes(shape=forest.type))

# simple lm
PPO.model <- lm(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*NOy.NHx.combined.28.year*soil.type*pH.H2O,
               data = PPO.enzyme)

print(PPO.model)
summary(PPO.model)


# remove outliers
outlierTest(PPO.model,cutoff=Inf, TOC~1)

ggplot(PPO.model, aes(x= TOC , y=standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

PPO.enzyme_1=PPO.enzyme


PPO.enzyme_1 <- PPO.enzyme[-c(222,221,210,159,209,207,29),] 

PPO.enzyme_1["standardised.enzyme.activity.rate"][PPO.enzyme_1["standardised.enzyme.activity.rate"]>0.025]<-NA # remove additional X-axis outlier - lone value no where near any other


PPO.model_1 <- lm(TOC ~ standardised.enzyme.activity.rate + horizon + forest.type + NOy.NHx.combined.28.year + soil.type +pH.H2O
                 , data = PPO.enzyme_1)


# check that outliers have been removed
ggplot(PPO.model_1, aes(x= TOC , y=standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

# check that outliers have been removed
ggplot(PPO.enzyme_1, aes(x = standardised.enzyme.activity.rate, y = TOC, colour = horizon))+
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

PPO.ancova.model <- lm(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*soil.type+NOy.NHx.combined.28.year+pH.H2O,
                      data = PPO.enzyme_1)

par(mfrow=c(2,2))
plot(PPO.ancova.model)
par(mfrow=c(1,1))

## accounting for leverage - e.g., one extreme value shift the trend line

PPO.high <- as.data.frame(cooks.distance(PPO.ancova.model))
colnames(PPO.high) <- 'x'
PPO.high <- subset(PPO.high, x>4/nrow(PPO.ancova.model))
PPO.enzyme_1_1 <- subset(PPO.enzyme_1,!rownames(PPO.enzyme_1)%in% rownames(PPO.high))

PPO.ancova.model_1 <- lm(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*soil.type+NOy.NHx.combined.28.year+pH.H2O,
                        data = PPO.enzyme_1_1)

par(mfrow=c(2,2))
plot(PPO.ancova.model_1)
par(mfrow=c(1,1))

# transformation
PPO.enzyme.TOC.ln_1 <- mutate(PPO.enzyme_1, lnTOC = log(TOC))
PPO.enzyme.TOC.ancova.mod.log <- lm(lnTOC ~standardised.enzyme.activity.rate*horizon*forest.type*soil.type+NOy.NHx.combined.28.year+pH.H2O,
                                   data = PPO.enzyme.TOC.ln_1)

par(mfrow=c(2,2))
plot(PPO.enzyme.TOC.ancova.mod.log)
par(mfrow=c(1,1))

### DO NOT use logged data - equal variance assumption not improved - 
# Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption

## ANOVA
standardised.enzyme.activity.rate=as.character(PPO.ancova.model_1[["effects"]][["standardised.enzyme.activity.rate"]])
horizon=factor(PPO.ancova.model_1[["effects"]][["horizonO"]])

summary(PPO.ancova.model_1)
anova(PPO.ancova.model_1)

# Tukey HSD

PPO.ancova.model_1_aov <- aov(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*soil.type+NOy.NHx.combined.28.year+pH.H2O,
                             data = PPO.enzyme_1_1)
PPO.ancova.model_1_aov_tt <- TukeyHSD(PPO.ancova.model_1_aov, which = 'horizon') # p-values too small for tukey to recognize
PPO.ancova.model_1_aov_tt                                                                                                 # values < 10-3 can't be easily distinguished wtih aov function

str(PPO.ancova.model_1_aov_tt)

print(PPO.ancova.model_1_aov_tt,digits=17)

## presenting results
PPO.enzyme.plt <- ggplot(PPO.ancova.model_1, aes(x = standardised.enzyme.activity.rate, y = TOC, colour = horizon))+
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")), size = 5) +
  geom_point(aes(shape=forest.type, colour = horizon), size = 2) +
  stat_smooth(fullrange=TRUE, method="lm",se=T)+
  labs(shape = "Forest Type") + labs(colour = "Horizon") +
  labs(title="PPO ")+
  ylab("")+
  xlab(expression(paste('Enzyme activity (',mu,'mol g'^-1," soil h"^-1,")")))+
  scale_x_continuous(expand=c(0,0), limits=c(0,0.024873635)) +
  scale_y_continuous(expand=c(0,0), limits=c(-10,80)) +
  coord_cartesian(xlim=c(0,0.024873635), ylim=c(-10,80)) +
  scale_x_continuous(breaks=seq(0,0.025,0.005)) +
  scale_colour_manual(name="Horizon", # Legend label, use darker colors
                      breaks=c("A", "O"),
                      labels=c("A", "O"),
                      values=c("#CC9900","#990000"))+
  theme_classic() +   
  theme(legend.position="none",axis.text=element_text(size=15), axis.title=element_text(size=20))

PPO.enzyme.plt



####### NAG

## visualising data
ggplot(NAG.enzyme, aes(x = standardised.enzyme.activity.rate, y = TOC, colour = horizon, group = forest.type))+
  geom_point(aes(shape=forest.type))

# simple lm
NAG.model <- lm(TOC ~ standardised.enzyme.activity.rate + horizon + forest.type + NOy.NHx.combined.28.year + soil.type 
                , data = NAG.enzyme)

print(NAG.model)
summary(NAG.model)


# remove outliers
outlierTest(NAG.model,cutoff=Inf, TOC~1)

ggplot(NAG.model, aes(x= TOC , y=standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

NAG.enzyme_1=NAG.enzyme


NAG.enzyme_1 <- NAG.enzyme[-c(76, 77, 103, 46, 129, 155, 498, 150, 148, 139, 18, 96, 137, 183, 161,
                              144,142,141,
                              60),] # removing values to improve normality

NAG.enzyme_1["standardised.enzyme.activity.rate"][NAG.enzyme_1["standardised.enzyme.activity.rate"]> 0.1 & NAG.enzyme_1["horizon"] == "O"]<-NA # remove additional X-axis outlier - lone value no where near any other
NAG.enzyme_1["standardised.enzyme.activity.rate"][NAG.enzyme_1["standardised.enzyme.activity.rate"]> 0.8 & NAG.enzyme_1["horizon"] == "A"]<-NA
NAG.enzyme_1["TOC"][NAG.enzyme_1["TOC"]> 40 & NAG.enzyme_1["horizon"] == "A"]<-NA

NAG.model_1 <- lm(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*soil.type+NOy.NHx.combined.28.year+pH.H2O,
                  data = NAG.enzyme_1)


# check that outliers have been removed
ggplot(NAG.model_1, aes(x= TOC , y=standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

# check that outliers have been removed
ggplot(NAG.enzyme_1, aes(x = standardised.enzyme.activity.rate, y = TOC, colour = horizon))+
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

NAG.ancova.model <- lm(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*soil.type+NOy.NHx.combined.28.year+pH.H2O, 
                       data = NAG.enzyme_1)

par(mfrow=c(2,2))
plot(NAG.ancova.model)
par(mfrow=c(1,1))

## accounting for leverage - e.g., one extreme value shift the trend line

NAG.high <- as.data.frame(cooks.distance(NAG.ancova.model))
colnames(NAG.high) <- 'x'
NAG.high <- subset(NAG.high, x>4/nrow(NAG.ancova.model))
NAG.enzyme_1_1 <- subset(NAG.enzyme_1,!rownames(NAG.enzyme_1)%in% rownames(NAG.high))

NAG.ancova.model_1 <- lm(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*soil.type+NOy.NHx.combined.28.year+pH.H2O,
                         data = NAG.enzyme_1_1)

par(mfrow=c(2,2))
plot(NAG.ancova.model)
par(mfrow=c(1,1))


# transformation
NAG.enzyme.TOC.ln_1 <- mutate(NAG.enzyme_1, lnTOC = log(TOC))
NAG.enzyme.TOC.ancova.mod.log <- lm(lnTOC ~ standardised.enzyme.activity.rate*horizon*forest.type*soil.type+NOy.NHx.combined.28.year+pH.H2O,
                                    data = NAG.enzyme.TOC.ln_1)

par(mfrow=c(2,2))
plot(NAG.enzyme.TOC.ancova.mod.log)
par(mfrow=c(1,1))

### DO NOT use logged data - equal variance assumption not improved - 
# Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption

## ANOVA
standardised.enzyme.activity.rate=as.character(NAG.ancova.model_1[["effects"]][["standardised.enzyme.activity.rate"]])
horizon=factor(NAG.ancova.model_1[["effects"]][["horizonO"]])

summary(NAG.ancova.model_1)
anova(NAG.ancova.model_1)

# tukey HSD

NAG.ancova.model_1_aov <- aov(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*soil.type+NOy.NHx.combined.28.year+pH.H2O,
                              data = NAG.enzyme_1_1)
NAG.ancova.model_1_aov_tt <- TukeyHSD(NAG.ancova.model_1_aov, which = 'horizon:forest.type') # p-values too small for tukey to recognize
NAG.ancova.model_1_aov_tt                                                                                                 # values < 10-3 can't be easily distinguished wtih aov function

str(NAG.ancova.model_1_aov_tt)

print(NAG.ancova.model_1_aov_tt,digits=17)

## presenting results
NAG.enzyme.plt <- ggplot(NAG.ancova.model_1, aes(x = standardised.enzyme.activity.rate, y = TOC))+
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*")), size = 5) +
  geom_point(aes(shape=forest.type, colour = horizon), size = 2) +
  stat_smooth(method="lm", se=T, colour = "black")+
  labs(shape = "Forest Type") + labs(colour = "Horizon") +
  labs(title="NAG ") +
  ylab("SOC (%)")+
  xlab(expression(paste("")))+
  # scale_x_continuous(expand=c(0,0), limits=c(0,0.8)) +
  scale_x_continuous(expand=c(0,0), limits=c(0,0.8)) +
  scale_y_continuous(expand=c(0,0), limits=c(-10,70)) +
  coord_cartesian(xlim=c(0,0.8), ylim=c(-10,70)) +
  expand_limits(x=c(0,0.8))+
  scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8)) +
  
  scale_colour_manual(name="Horizon", # Legend label, use darker colors
                      breaks=c("A", "O"),
                      labels=c("A", "O"),
                      values=c("#CC9900","#990000"))+
  theme_classic()+   
  theme(legend.position="none",axis.text=element_text(size=15), axis.title=element_text(size=20))

NAG.enzyme.plt



####### LAP

## visualising data
ggplot(LAP.enzyme, aes(x = standardised.enzyme.activity.rate, y = TOC, colour = horizon, group = forest.type))+
  geom_point(aes(shape=forest.type))

# simple lm
LAP.model <- lm(TOC ~ standardised.enzyme.activity.rate + horizon + forest.type + NOy.NHx.combined.28.year + soil.type 
                , data = LAP.enzyme)

print(LAP.model)
summary(LAP.model)


# remove outliers
outlierTest(LAP.model,cutoff=Inf, TOC~1)

ggplot(LAP.model, aes(x= TOC , y=standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

LAP.enzyme_1=LAP.enzyme


LAP.enzyme_1 <- LAP.enzyme[-c(),] # not needed, but kept for completeness

LAP.enzyme_1["standardised.enzyme.activity.rate"][LAP.enzyme_1["standardised.enzyme.activity.rate"]> 0.6 & LAP.enzyme_1["horizon"] == "O"]<-NA # remove additional X-axis outlier - lone value no where near any other
LAP.enzyme_1["standardised.enzyme.activity.rate"][LAP.enzyme_1["standardised.enzyme.activity.rate"]> 4 & LAP.enzyme_1["horizon"] == "A"]<-NA
LAP.enzyme_1["standardised.enzyme.activity.rate"][LAP.enzyme_1["standardised.enzyme.activity.rate"]> 1 & LAP.enzyme_1["forest.type"] == "deciduous" & LAP.enzyme_1["horizon"] == "A"]<-NA


LAP.model_1 <- lm(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*NOy.NHx.combined.28.year*soil.type*pH.H2O, 
                  data = LAP.enzyme_1)


# check that outliers have been removed
ggplot(LAP.model_1, aes(x= TOC , y=standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

# check that outliers have been removed
ggplot(LAP.enzyme_1, aes(x = standardised.enzyme.activity.rate, y = TOC, colour = horizon))+
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

LAP.ancova.model <- lm(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*soil.type+NOy.NHx.combined.28.year+pH.H2O, 
                       data = LAP.enzyme_1)

par(mfrow=c(2,2))
plot(LAP.ancova.model)
par(mfrow=c(1,1))

## accounting for leverage - e.g., one extreme value shift the trend line

LAP.high <- as.data.frame(cooks.distance(LAP.ancova.model))
colnames(LAP.high) <- 'x'
LAP.high <- subset(LAP.high, x>4/nrow(LAP.ancova.model))
LAP.enzyme_1_1 <- subset(LAP.enzyme_1,!rownames(LAP.enzyme_1)%in% rownames(LAP.high))

LAP.ancova.model_1 <- lm(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*soil.type+NOy.NHx.combined.28.year+pH.H2O, 
                         data = LAP.enzyme_1_1)

par(mfrow=c(2,2))
plot(LAP.ancova.model)
par(mfrow=c(1,1))

# transformation
LAP.enzyme.TOC.ln_1 <- mutate(LAP.enzyme_1, lnTOC = log(TOC))
LAP.enzyme.TOC.ancova.mod.log <- lm(lnTOC ~ standardised.enzyme.activity.rate*horizon*forest.type*NOy.NHx.combined.28.year*soil.type*pH.H2O, 
                                    data = LAP.enzyme.TOC.ln_1)

plot(LAP.enzyme.TOC.ancova.mod.log, add.smooth = T, which = 1)
plot(LAP.enzyme.TOC.ancova.mod.log, which = 2)
plot(LAP.enzyme.TOC.ancova.mod.log, add.smooth = T, which = 3)

### DO NOT use transformed data - equal variance assumption not improved - 
# Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption

## ANOVA
standardised.enzyme.activity.rate=as.character(LAP.ancova.model_1[["effects"]][["standardised.enzyme.activity.rate"]])
horizon=factor(LAP.ancova.model_1[["effects"]][["horizonO"]])

summary(LAP.ancova.model_1)
anova(LAP.ancova.model_1)

# tukey HSD

LAP.ancova.model_1_aov <- aov(TOC ~ standardised.enzyme.activity.rate*horizon*forest.type*soil.type+NOy.NHx.combined.28.year+pH.H2O,
                              data = LAP.enzyme_1_1)
LAP.ancova.model_1_aov_tt <- TukeyHSD(LAP.ancova.model_1_aov, which = 'horizon:forest.type') # p-values too small for tukey to recognize
LAP.ancova.model_1_aov_tt                                                                                                 # values < 10-3 can't be easily distinguished wtih aov function

str(LAP.ancova.model_1_aov_tt)

print(LAP.ancova.model_1_aov_tt,digits=17)

## presenting results
LAP.enzyme.plt <- ggplot(LAP.ancova.model_1, aes(x = standardised.enzyme.activity.rate, y = TOC, colour = horizon))+
  geom_point(aes(shape=forest.type), size = 2) +
  labs(shape = "Forest Type") + labs(colour = "Horizon") +
  labs(title="LAP ") +
  ylab("")+
  xlab(expression(paste("")))+
  scale_y_continuous(limits=c(0, 60))+
  scale_x_continuous(limits=c(0, 4))+
  scale_colour_manual(name="Horizon", # Legend label, use darker colors
                      breaks=c("A", "O"),
                      labels=c("A", "O"),
                      values=c("#CC9900","#990000"))+
  theme_classic()+
  theme(legend.position="none",axis.text=element_text(size=15), axis.title=element_text(size=20))


LAP.enzyme.plt


#### make grid of plots

SOC.plot <- 
  plot_grid(BG.enzyme.plt, BX.enzyme.plt, NAG.enzyme.plt, LAP.enzyme.plt, AP.enzyme.plt, PPO.enzyme.plt, nrow = 3, labels = c("auto"), label_size = 17)

SOC.plot

ggsave("SOC~EEA____.pdf", height = 12, width = 13)

