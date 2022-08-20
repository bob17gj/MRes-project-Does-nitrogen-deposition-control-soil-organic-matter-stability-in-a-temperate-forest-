library(tidyverse) #for all data wrangling
library(multcomp) # for ANCOVA post-hoc tests 
library(multcompView)
library(cowplot) #for manuscript ready figures
library(car)# for bonferroni outlier test
library(agricolae)

rm(list=ls())
setwd("~/General Waring/MRes Project")

hydrolytic_enzymes <- read.csv("Hydrolytic enzyme master sheet for R (2).csv")
str(hydrolytic_enzymes)


# convert all "NaN values in DT to NA
EEA <- hydrolytic_enzymes %>% mutate_all(~ifelse(is.nan(.), NA, .))



# change "inf" values to NA
EEA <- EEA %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf)))


# per g soil #####

## visualising data
EEA.per.g.soil.plt <- ggplot(EEA)+
  geom_col(aes(x=enzyme, y = standardised.enzyme.activity.rate, fill = horizon))
EEA.per.g.soil.plt

## carrying on without fill changing for enzyme type


# simple lm
EEA.per.g.soil.model <- lm(standardised.enzyme.activity.rate ~ enzyme + horizon + forest.type + NOy.NHx.combined.28.year +soil.type
                           , data = EEA)

print(EEA.per.g.soil.model)
summary(EEA.per.g.soil.model)


# remove outliers
outlierTest(EEA.per.g.soil.model,cutoff=Inf, standardised.enzyme.activity.rate~1)

ggplot(EEA.per.g.soil.model, aes(y= standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

EEA_1=EEA

EEA_1 <- EEA_1[-c(1008,492,1080, 
                  638,644,640,
                  1105,
                  835,
                  881,892,
                  725,
                  310,459,
                  942,
                  307,308,
                  373,
                  1101,290,746,974),] # removing points to improve normality


EEA_1["standardised.enzyme.activity.rate"][EEA_1["standardised.enzyme.activity.rate"]>10]<-NA # remove additional X-axis outlier - lone value no where near any other



EEA.per.g.soil.model_1 <- lm(standardised.enzyme.activity.rate ~ enzyme + horizon + forest.type + NOy.NHx.combined.28.year +soil.type
                             , data = EEA_1)

# check that outliers have been removed
ggplot(EEA.per.g.soil.model_1, aes(y= standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

# check that outliers have been removed
ggplot(EEA.per.g.soil.model_1)+
  geom_col(aes(x=enzyme, y = standardised.enzyme.activity.rate, fill = horizon))

## ANOVA
print(EEA.per.g.soil.model_1)

anova(EEA.per.g.soil.model_1)

## Assumptions
# The linearity assumption states that the general relationship between the 
# response and predictor variable should look like a straight line. 
# We can evaluate this assumption by constructing a residuals vs. fitted values plot.

## linearity - don't want a pattern. E.g., residuals get larger/smaller as X gets larger
plot(EEA.per.g.soil.model_1, which = 1, add.smooth = T)

## normality 
plot(EEA.per.g.soil.model_1, which = 2)

## variance
plot(EEA.per.g.soil.model_1, add.smooth = T, which = 3)


##### diagnostics - use transformed data - Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption


## plot ANCOVA
EEA.per.g.soil.ANCOVA.model<- lm(standardised.enzyme.activity.rate ~ enzyme * horizon * forest.type * soil.type + NOy.NHx.combined.28.year
                                 , data = EEA_1)

par(mfrow=c(2,2))
plot(EEA.per.g.soil.ANCOVA.model)
par(mfrow=c(1,1))


# transformation
EEA.ln_1 <- mutate(EEA_1, lnstandardised.enzyme.activity.rate = log(standardised.enzyme.activity.rate))

EEA.ln_1 <- EEA.ln_1 %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf))) #convert Inf values to NA
EEA.ln_1 <- EEA.ln_1 %>% mutate_all(~ifelse(is.nan(.), NA, .)) # convert NaN values to NA



EEA.ln_1.mod.log <- lm(lnstandardised.enzyme.activity.rate ~ enzyme * horizon * forest.type * soil.type+NOy.NHx.combined.28.year,
                       EEA.ln_1)

par(mfrow=c(2,2))
plot(EEA.ln_1.mod.log)
par(mfrow=c(1,1))

### use transformed data - equal variance assumption not improved - 
# Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption

## ANOVA
standardised.enzyme.activity.rate=as.character(EEA.ln_1.mod.log[["effects"]][["standardised.enzyme.activity.rate"]])
horizon=factor(EEA.ln_1.mod.log[["effects"]][["horizonO"]])

EEA.ln_1$enzyme = as.factor(EEA.ln_1$enzyme)
EEA.ln_1$horizon = as.character(EEA.ln_1$horizon)
EEA.ln_1$forest.type = as.factor(EEA.ln_1$forest.type)
EEA.ln_1$soil.type = as.factor(EEA.ln_1$soil.type)
EEA.ln_1$NOy.NHx.combined.28.year=as.factor(EEA.ln_1$NOy.NHx.combined.28.year)
str(EEA.ln_1)

anova(EEA.ln_1.mod.log)

aov <- aov(EEA.ln_1.mod.log)

aov.2 <- aov(EEA.per.g.soil.ANCOVA.model)

TukeyHSD(aov.2, which = 'enzyme')
HSD.test(aov.2, trt = c("enzyme"), console = TRUE)

TukeyHSD(aov.2, which = 'horizon:enzyme')
HSD.test(aov.2, trt = c("horizon", "enzyme"), console = TRUE)


TukeyHSD(aov.2, which = 'forest.type:enzyme')
HSD.test(aov.2, trt = c("forest.type", "enzyme"), console = TRUE)

TukeyHSD(aov.2, which = 'forest.type')
HSD.test(aov.2, trt = c("forest.type"), console = TRUE)


# PPO per g soil analysis #### 

PPO <-hydrolytic_enzymes%>%filter(enzyme %in% c("PPO"))

# convert all "NaN values in DT to NA
PPO <- hydrolytic_enzymes %>% mutate_all(~ifelse(is.nan(.), NA, .))


# change "inf" values to NA
PPO <- PPO %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf)))


## visualising data
PPO.per.g.soil.plt <- ggplot(PPO)+
  geom_col(aes(x=enzyme, y = standardised.enzyme.activity.rate, fill = horizon))
PPO.per.g.soil.plt

## carrying on without fill changing for enzyme type


# simple lm
PPO.per.g.soil.model <- lm(standardised.enzyme.activity.rate ~ horizon + forest.type + NOy.NHx.combined.28.year +soil.type
                           , data = PPO)

print(PPO.per.g.soil.model)
summary(PPO.per.g.soil.model)


# remove outliers
outlierTest(PPO.per.g.soil.model,cutoff=Inf, ab~1)

ggplot(PPO.per.g.soil.model, aes(y= standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()

PPO_1=PPO

PPO_1 <- PPO_1[-c(29),] # removing points to improve normality


PPO_1["standardised.enzyme.activity.rate"][PPO_1["standardised.enzyme.activity.rate"]>0.12]<-NA # remove additional X-axis outlier - lone value no where near any other



PPO.per.g.soil.model_1 <- lm(standardised.enzyme.activity.rate ~ horizon + forest.type + NOy.NHx.combined.28.year +soil.type
                             , data = PPO_1)

# check that outliers have been removed
ggplot(PPO.per.g.soil.model_1, aes(y= standardised.enzyme.activity.rate, group = horizon)) + 
  geom_boxplot()


## ANOVA
print(PPO.per.g.soil.model_1)

anova(PPO.per.g.soil.model_1)

## Assumptions
# The linearity assumption states that the general relationship between the 
# response and predictor variable should look like a straight line. 
# We can evaluate this assumption by constructing a residuals vs. fitted values plot.

## linearity - don't want a pattern. E.g., residuals get larger/smaller as X gets larger
plot(PPO.per.g.soil.model_1, which = 1, add.smooth = T)

## normality 
plot(PPO.per.g.soil.model_1, which = 2)

## variance
plot(PPO.per.g.soil.model_1, add.smooth = T, which = 3)


##### diagnostics - use transformed data - Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption


## plot ANCOVA
PPO.per.g.soil.ANCOVA.model<- lm(standardised.enzyme.activity.rate ~ horizon * forest.type * NOy.NHx.combined.28.year *soil.type
                                 , data = PPO_1)

par(mfrow=c(2,2))
plot(PPO.per.g.soil.ANCOVA.model)
par(mfrow=c(1,1))


# transformation
PPO.ln_1 <- mutate(PPO_1, lnabsorbance = log(standardised.enzyme.activity.rate))

PPO.ln_1 <- PPO.ln_1 %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf))) #convert Inf values to NA
PPO.ln_1 <- PPO.ln_1 %>% mutate_all(~ifelse(is.nan(.), NA, .)) # convert NaN values to NA



PPO.ln_1.mod.log <- lm(lnabsorbance ~  horizon * forest.type * soil.type+ NOy.NHx.combined.28.year,
                       PPO.ln_1)

par(mfrow=c(2,2))
plot(PPO.ln_1.mod.log)
par(mfrow=c(1,1))

### use transformed data - equal variance assumption not improved - 
# Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption

## ANOVA
standardised.enzyme.activity.rate=as.character(PPO.ln_1.mod.log[["effects"]][["standardised.enzyme.activity.rate"]])
horizon=factor(PPO.ln_1.mod.log[["effects"]][["horizonO"]])

PPO.ln_1$enzyme = as.factor(PPO.ln_1$enzyme)
PPO.ln_1$horizon = as.character(PPO.ln_1$horizon)
PPO.ln_1$forest.type = as.factor(PPO.ln_1$forest.type)
PPO.ln_1$soil.type = as.factor(PPO.ln_1$soil.type)
PPO.ln_1$NOy.NHx.combined.28.year=as.factor(PPO.ln_1$NOy.NHx.combined.28.year)
str(PPO.ln_1)

anova(PPO.ln_1.mod.log)

aov <- aov(PPO.ln_1.mod.log)

TukeyHSD(aov, which = 'horizon')
HSD.test(aov, trt = c("horizon"), console = TRUE)

TukeyHSD(aov, which = 'forest.type')
HSD.test(aov, trt = c("forest.type"), console = TRUE)



## hydrolytic graphs

EEA_2 <- EEA_1

EEA_2 <-
  EEA_2 %>%
  filter(enzyme !='PPO')

enzyme_sum.horizon <- 
  EEA_2 %>% 
  group_by(enzyme, horizon) %>% 
  summarise(mean_activity = mean(standardised.enzyme.activity.rate, na.rm = T), n=n(), sd=sd(standardised.enzyme.activity.rate, na.rm = T), se=sd/sqrt(n))
enzyme_sum.horizon

# Compact letter display
EEA.per.g.soil.ANCOVA.model<- lm(standardised.enzyme.activity.rate ~ enzyme * horizon * forest.type * soil.type + NOy.NHx.combined.28.year
                                 , data = EEA_1)

aov.2 <- aov(EEA.per.g.soil.ANCOVA.model)

TukeyHSD(aov.2, which = 'enzyme')
HSD.test(aov.2, trt = c("enzyme"), console = TRUE)

TukeyHSD(aov.2, which = 'horizon:enzyme')
HSD.test(aov.2, trt = c("horizon", "enzyme"), console = TRUE)

enzyme_sum.horizon$cld <- c("a","ab","bc","cd","cd","cd","cd","d","bc","cd")

# plot

EEA.per.g.soil.horizon.plt <- ggplot(enzyme_sum.horizon, 
                                     aes(x=enzyme, y=mean_activity, fill = horizon)) + 
  geom_bar(position=position_dodge(), stat="identity", colour = "black") +
  geom_bar(position=position_dodge(), stat="identity") + 
  geom_errorbar(data = enzyme_sum.horizon,
                aes(ymin = mean_activity-se, ymax= mean_activity+se), 
                width = .2, position=position_dodge(width=0.9))+
  geom_text(data = enzyme_sum.horizon, aes(x = enzyme, label = cld, y= mean_activity + se + 0.02),
            position = position_dodge(0.9))+
#  coord_flip() +
  labs(x ="Enzyme")+
  ylab(expression(paste('Enzyme activity (',mu,'mol g'^-1,"soil h"^-1,")")))+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75),
                     #sec.axis = sec_axis(trans=~./1000,name=expression(paste('PPO activity (',mu,'mol g'^-1,"soil h"^-1,")")))
                     )+
  expand_limits(y=c(0,0.75))+
  scale_fill_manual(name="Horizon", # Legend label, use darker colors
                 breaks=c("A", "O"),
                 labels=c("A", "O"),
                 values=c("#CC9900","#990000")) +
    theme_classic()+
    theme(legend.position="none")
    


EEA.per.g.soil.horizon.plt


# inc;ude PPO in same graph
#enzyme_sum.horizon.2 <- 
#  EEA_1 %>% 
#  group_by(enzyme, horizon) %>% 
#  summarise(mean_activity = mean(standardised.enzyme.activity.rate, na.rm = T), n=n(), sd=sd(standardised.enzyme.activity.rate, na.rm = T), se=sd/sqrt(n))
# enzyme_sum.horizon.2


#per.g.soil.horizon.plt <- ggplot(enzyme_sum.horizon.2, 
#                                     aes(x=enzyme, y=mean_activity, fill = horizon)) + 
#  geom_bar(position=position_dodge(), stat="identity", colour = "black") +
#  geom_bar(position=position_dodge(), stat="identity") + 
#  geom_errorbar(data = enzyme_sum.horizon.2,
#                aes(ymin = mean_activity-se, ymax= mean_activity+se), 
#                width = .2, position=position_dodge(width=0.9))+
#  geom_text(data = enzyme_sum.horizon, aes(x = enzyme, label = cld, y= mean_activity + se + 0.02),
#            position = position_dodge(0.9))+
  #  coord_flip() +
#  labs(x ="Enzyme")+
#  ylab(expression(paste('Enzyme activity (',mu,'mol g'^-1,"soil h"^-1,")")))+
#  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8),
#                     sec.axis = sec_axis(trans=~.*10,name=expression(paste('PPO activity (',mu,'mol g'^-1,"soil h"^-1,")")))
 # )+
#  expand_limits(y=c(0,0.75))+
 # scale_fill_manual(name="Horizon", # Legend label, use darker colors
 #                   breaks=c("A", "O"),
#                    labels=c("A", "O"),
 #                   values=c("#CC9900","#990000")) +
 # theme_classic()+
  #theme(legend.position="none")

#per.g.soil.horizon.plt


# forest type plot

TukeyHSD(aov.2, which = 'forest.type:enzyme')
HSD.test(aov.2, trt = c("forest.type", "enzyme"), console = TRUE)

enzyme_sum.forest.type <- 
  EEA_2 %>% 
  group_by(enzyme, forest.type) %>% 
  summarise(mean_activity = mean(standardised.enzyme.activity.rate, na.rm = T), n=n(), sd=sd(standardised.enzyme.activity.rate, na.rm = T), se=sd/sqrt(n))
enzyme_sum.forest.type


enzyme_sum.forest.type$cld <- c("bc","a","d","b","d","b","cd","bc","cd","b")

EEA.per.g.soil.forest.type.plt <- 
  ggplot(enzyme_sum.forest.type, aes(x=enzyme, y=mean_activity, fill=forest.type)) + 
  geom_bar(position=position_dodge(), stat="identity", colour = "black") +
  geom_bar(position=position_dodge(), stat="identity") + 
  geom_errorbar(data = enzyme_sum.forest.type,
                aes(ymin = mean_activity-se, ymax= mean_activity+se), 
                width = .2, position=position_dodge(width=0.9))+
  geom_text(data = enzyme_sum.forest.type, aes(x = enzyme, label = cld, y= mean_activity + se + 0.03),
            position = position_dodge(0.9))+
  xlab("Enzyme") +
  ylab(expression(paste('Enzyme activity (',mu,'mol g'^-1,"soil h"^-1,")")))+
#  coord_flip() +
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),)+
  expand_limits(y=c(0,1))+
  theme(axis.title.x=element_blank())+
  scale_fill_manual(name="Forest type", # Legend label, use darker colors
                 breaks=c("broadleaf", "evergreen"),
                 labels=c("broadleaf", "coniferous"),values=c("#56B4E9","#0072B2")) +
  theme_classic() +
  theme(legend.position="none")

EEA.per.g.soil.forest.type.plt


#### soil type ####

enzyme_sum.soil.type <- 
  EEA_2 %>% 
  group_by(enzyme, soil.type) %>% 
  summarize(mean_activity = mean(lnstandardised.enzyme.activity.rate, na.rm = T), n=n(), sd=sd(lnstandardised.enzyme.activity.rate, na.rm = T), se=sd/sqrt(n))
enzyme_sum.soil.type

EEA.per.g.soil.soil.type.plt <- 
  ggplot(enzyme_sum.soil.type, aes(x=enzyme, y=mean_activity, fill=soil.type)) + 
  geom_bar(position=position_dodge(), stat="identity", colour = "black") +
  geom_bar(position=position_dodge(), stat="identity") +  
  xlab("enzyme") +
  ylab("ln(Enzyme activity)") +
#  coord_flip() +
  labs(x ="Enzyme", y = "ln(Enzyme activity)")+
  theme(axis.title.x=element_blank())+
  scale_fill_manual(name="Soil type", # Legend label, use darker colors
                    breaks=c("peaty gley", "brown earth"),
                    labels=c("peaty gley", "brown earth"),values=c("#003333","#99CCCC")) +
  scale_y_reverse(lim=c(0,-2))+
  theme_classic() +   theme(axis.text=element_text(size=10), axis.title=element_text(size=9))


EEA.per.g.soil.soil.type.plt


## PPO grpahs

EEA_3 <- EEA_1

EEA_3 <-
  EEA_3 %>%
  filter(enzyme =="PPO")
EEA_3

PPO_sum.horizon <- 
  EEA_3 %>% 
  group_by(enzyme, horizon) %>% 
  summarise(mean_activity = mean(standardised.enzyme.activity.rate, na.rm = T), n=n(), sd=sd(standardised.enzyme.activity.rate, na.rm = T), se=sd/sqrt(n))
PPO_sum.horizon

# Compact letter display

# plot


PPO.per.g.soil.horizon.plt <- ggplot(PPO_sum.horizon, 
                                     aes(x=enzyme, y=mean_activity, fill=horizon)) + 
  geom_bar(position=position_dodge(), stat="identity", colour = "black", width = 0.5) +
  geom_bar(position=position_dodge(), stat="identity", width = 0.5) + 
  geom_errorbar(data = PPO_sum.horizon,
                aes(ymin = mean_activity-se, ymax= mean_activity+se), 
                width = .2, position=position_dodge(width=0.5))+
#  coord_flip() +
  ylab(expression(paste('Enzyme activity (',mu,'mol g'^-1,"soil h"^-1,")")))+
  labs(x ="")+
  scale_y_continuous(breaks=c(0,0.004,0.008,0.012,0.016))+
  expand_limits(y=c(0,0.016))+
  scale_fill_manual(name="Horizon", # Legend label, use darker colors
                    breaks=c("A", "O"),
                    labels=c("A", "O"),
                    values=c("#CC9900","#990000")) +
  theme_classic() +
  theme(legend.background = element_rect(size=0.3, linetype="solid", 
                                          colour ="black"))

PPO.per.g.soil.horizon.plt
# forest type plot

TukeyHSD(aov, which = 'forest.type:enzyme')
HSD.test(aov, trt = c("forest.type", "enzyme"), console = TRUE)

PPO_sum.forest.type <- 
  EEA_3 %>% 
  group_by(enzyme, forest.type) %>% 
  summarise(mean_activity = mean(standardised.enzyme.activity.rate, na.rm = T), n=n(), sd=sd(standardised.enzyme.activity.rate, na.rm = T), se=sd/sqrt(n))
PPO_sum.forest.type



PPO.per.g.soil.forest.type.plt <- 
  ggplot(PPO_sum.forest.type, aes(x=enzyme, y=mean_activity, fill=forest.type)) + 
  geom_bar(position=position_dodge(), stat="identity", colour = "black", width = 0.5) +
  geom_bar(position=position_dodge(), stat="identity", width = 0.5) + 
  geom_errorbar(data = PPO_sum.forest.type,
                aes(ymin = mean_activity-se, ymax= mean_activity+se), 
                width = .2, position=position_dodge(width=0.5))+
  #  coord_flip() +
  scale_y_continuous(breaks=c(0,0.004,0.008,0.012))+
  expand_limits(y=c(0,0.012))+
  ylab(expression(paste('Enzyme activity (',mu,'mol g'^-1,"soil h"^-1,")")))+
  labs(x ="")+
  theme(axis.title.x=element_blank())+
  scale_fill_manual(name="Forest type", # Legend label, use darker colors
                    breaks=c("broadleaf", "evergreen"),
                    labels=c("broadleaf", "coniferous"),values=c("#56B4E9","#0072B2")) +
  theme_classic() +
  theme(legend.background = element_rect(size=0.3, linetype="solid", 
                                          colour ="black"))


PPO.per.g.soil.forest.type.plt


#### make grid of plots

plot_grid(PPO.per.g.soil.horizon.plt, PPO.per.g.soil.forest.type.plt,
          EEA.per.g.soil.horizon.plt, EEA.per.g.soil.forest.type.plt,
          nrow = 2, labels = c("a","b","",""), label_size = 14)




# per g mbc ####

rm(list=ls())
setwd("~/General Waring/Project")

hydrolytic_enzymes <- read.csv("Hydrolytic enzyme master sheet for R (2).csv")
str(hydrolytic_enzymes)

EEA.mbc <- hydrolytic_enzymes
EEA.mbc <- mutate(hydrolytic_enzymes, EEA.mbc = enzyme.activity/microbial.biomass.C)
PPO.mbc <- mutate(hydrolytic_enzymes, PPO.microbioc = standardised.enzyme.activity.rate/microbial.biomass.C)


EEA.mbc <-EEA.mbc%>%filter(enzyme %in% c("BG", "BX","LAP","NAG","AP"))
PPO.mbc <-PPO.mbc%>%filter(enzyme %in% c("PPO"))


# convert all "NaN values in DT to NA
EEA.mbc <- EEA.mbc %>% mutate_all(~ifelse(is.nan(.), NA, .))
PPO.mbc <- PPO.mbc %>% mutate_all(~ifelse(is.nan(.), NA, .))

# change "inf" values to NA
EEA.mbc <- EEA.mbc %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf)))
PPO.mbc <- PPO.mbc %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf)))


# replace NA in PPO.microbioc column with 0
PPO.mbc <- PPO.mbc %>% mutate_at(vars("PPO.microbioc"), ~replace_na(.,0))




## visualising data
EEA.per.g.soil.plt <- ggplot(EEA.mbc)+
  geom_col(aes(x=enzyme, y = EEA.mbc, fill = horizon))
EEA.per.g.soil.plt

## carrying on without fill changing for enzyme type


# simple lm
EEA.per.g.soil.model <- lm(EEA.mbc ~ enzyme + horizon + forest.type + NOy.NHx.combined.28.year +soil.type
                           , data = EEA.mbc)

print(EEA.per.g.soil.model)
summary(EEA.per.g.soil.model)


# remove outliers
outlierTest(EEA.per.g.soil.model,cutoff=Inf, EEA.mbc~1)

ggplot(EEA.per.g.soil.model, aes(y= EEA.mbc, group = horizon)) + 
  geom_boxplot()

EEA_1.mbc=EEA.mbc

EEA_1.mbc <- EEA_1.mbc[-c(1008,1020,1078,934,813,815,814,
                          1105,
                          835,
                          725,
                          459,155,61,
                          942,
                          1057,
                          413,
                          145),] # removing points to improve normality

EEA_1.mbc["EEA.mbc"][EEA_1.mbc["EEA.mbc"]>0.04]<-NA # remove additional X-axis outlier - lone value no where near any other
EEA_1.mbc["EEA.mbc"][EEA_1.mbc["EEA.mbc"]>0.003 & EEA_1.mbc["horizon"] == "O"]<-NA

EEA.per.g.mbc.model_1 <- lm(EEA.mbc ~ enzyme + horizon + forest.type + NOy.NHx.combined.28.year +soil.type
                             , data = EEA_1.mbc)

# check that outliers have been removed
ggplot(EEA.per.g.mbc.model_1, aes(y= EEA.mbc, group = horizon)) + 
  geom_boxplot()

# check that outliers have been removed
ggplot(EEA.per.g.mbc.model_1)+
  geom_col(aes(x=enzyme, y = EEA.mbc, fill = horizon))

## ANOVA
print(EEA.per.g.mbc.model_1)

anova(EEA.per.g.mbc.model_1)

## Assumptions
# The linearity assumption states that the general relationship between the 
# response and predictor variable should look like a straight line. 
# We can evaluate this assumption by constructing a residuals vs. fitted values plot.

## linearity - don't want a pattern. E.g., residuals get larger/smaller as X gets larger
plot(EEA.per.g.mbc.model_1, which = 1, add.smooth = T)

## normality 
plot(EEA.per.g.mbc.model_1, which = 2)

## variance
plot(EEA.per.g.mbc.model_1, add.smooth = T, which = 3)


##### diagnostics - use transformed data - Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption


## plot ANCOVA
EEA.per.g.soil.ANCOVA.model<- lm(EEA.mbc ~  enzyme * horizon * forest.type * soil.type + NOy.NHx.combined.28.year
                                , data = EEA_1.mbc)

par(mfrow=c(2,2))
plot(EEA.per.g.soil.ANCOVA.model)
par(mfrow=c(1,1))


# transformation
EEA.mbc.ln_1 <- mutate(EEA_1.mbc, lnEEA.mbc = log(EEA.mbc))

EEA.mbc.ln_1 <- EEA.mbc.ln_1 %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf))) #convert Inf values to NA
EEA.mbc.ln_1 <- EEA.mbc.ln_1 %>% mutate_all(~ifelse(is.nan(.), NA, .)) # convert NaN values to NA



EEA.mbc.ln_1.mod.log <- lm(lnEEA.mbc ~ enzyme * horizon * forest.type * soil.type+NOy.NHx.combined.28.year ,
                       EEA.mbc.ln_1)

par(mfrow=c(2,2))
plot(EEA.mbc.ln_1.mod.log)
par(mfrow=c(1,1))

### use transformed data - equal variance assumption not improved - 
# Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption

## ANOVA
lnEEA.mbc=as.character(EEA.mbc.ln_1.mod.log[["effects"]][["lnEEA.mbc"]])
horizon=factor(EEA.mbc.ln_1.mod.log[["effects"]][["horizonO"]])

summary(EEA.mbc.ln_1.mod.log)
anova(EEA.mbc.ln_1.mod.log)

aov.mbc <- aov(EEA.mbc.ln_1.mod.log)

TukeyHSD(aov.mbc, which = 'enzyme')
HSD.test(aov.mbc, trt = c("enzyme"), console = TRUE)

TukeyHSD(aov.mbc, which = 'horizon')
HSD.test(aov.mbc, trt = c("horizon"), console = TRUE)

TukeyHSD(aov.mbc, which = 'horizon:enzyme')
HSD.test(aov.mbc, trt = c("horizon", "enzyme"), console = TRUE)

TukeyHSD(aov.mbc, which = 'forest.type')
HSD.test(aov.mbc, trt = c("forest.type"), console = TRUE)

TukeyHSD(aov.mbc, which = 'forest.type:enzyme')
HSD.test(aov.mbc, trt = c("forest.type", "enzyme"), console = TRUE)

# PPO data ####

## visualising data
PPO.per.g.mbc.plt <- ggplot(PPO.mbc)+
  geom_col(aes(x=enzyme, y = PPO.microbioc, fill = forest.type))
PPO.per.g.mbc.plt

## carrying on without fill changing for enzyme type


# simple lm
PPO.per.g.mbc.model <- lm(PPO.microbioc ~ horizon + forest.type + NOy.NHx.combined.28.year +soil.type
                           , data = PPO.mbc)

print(PPO.per.g.mbc.model)
summary(PPO.per.g.mbc.model)


# remove outliers
outlierTest(PPO.per.g.mbc.model,cutoff=Inf, PPO.microbioc~1)

ggplot(PPO.per.g.mbc.model, aes(y= PPO.microbioc, group = forest.type)) + 
  geom_boxplot()

PPO_1.mbc=PPO.mbc

PPO_1.mbc <- PPO_1.mbc[-c(54,51,53,52,49),] # removing points to improve normality

PPO_1.mbc["PPO.microbioc"][PPO_1.mbc["PPO.microbioc"]>0.08]<-NA # remove additional X-axis outlier - lone value no where near any other
PPO_1.mbc["PPO.microbioc"][PPO_1.mbc["PPO.microbioc"]>0.03 & PPO_1.mbc["forest.type"]== "evergreen"]<-NA # remove additional X-axis outlier - lone value no where near any other


PPO.per.g.mbc.model_1 <- lm(PPO.microbioc ~ horizon + forest.type + NOy.NHx.combined.28.year +soil.type
                            , data = PPO_1.mbc)

# check that outliers have been removed
ggplot(PPO.per.g.mbc.model_1, aes(y= PPO.microbioc, group = forest.type)) + 
  geom_boxplot()

## ANOVA
print(PPO.per.g.mbc.model_1)

anova(PPO.per.g.mbc.model_1)

## Assumptions
# The linearity assumption states that the general relationship between the 
# response and predictor variable should look like a straight line. 
# We can evaluate this assumption by constructing a residuals vs. fitted values plot.

## linearity - don't want a pattern. E.g., residuals get larger/smaller as X gets larger
plot(PPO.per.g.mbc.model_1, which = 1, add.smooth = T)

## normality 
plot(PPO.per.g.mbc.model_1, which = 2)

## variance
plot(PPO.per.g.mbc.model_1, add.smooth = T, which = 3)


##### diagnostics - use transformed data - Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption


## plot ANCOVA
PPO.per.g.mbc.ANCOVA.model<- lm(PPO.microbioc ~ horizon * forest.type *soil.type+ NOy.NHx.combined.28.year,
                                  data = PPO_1.mbc)

par(mfrow=c(2,2))
plot(PPO.per.g.mbc.ANCOVA.model)
par(mfrow=c(1,1))


# transformation
PPO.mbc.ln_1 <- mutate(PPO_1.mbc, lnPPO.microbioc = log(PPO.microbioc))

PPO.mbc.ln_1 <- PPO.mbc.ln_1 %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf))) #convert Inf values to NA
PPO.mbc.ln_1 <- PPO.mbc.ln_1 %>% mutate_all(~ifelse(is.nan(.), NA, .)) # convert NaN values to NA



PPO.mbc.ln_1.mod.log <- lm(lnPPO.microbioc ~ horizon * forest.type * soil.type+NOy.NHx.combined.28.year,
                           PPO.mbc.ln_1)

par(mfrow=c(2,2))
plot(PPO.mbc.ln_1.mod.log)
par(mfrow=c(1,1))

### use transformed data - equal variance assumption not improved - 
# Bonnie says having residuals the
# same distance from the fitted line (homoscedastic), 
# more important than meeting the normality assumption

## ANOVA
lnPPO.microbioc=as.character(PPO.mbc.ln_1.mod.log[["effects"]][["lnPPO.microbioc"]])
forest.type=factor(PPO.mbc.ln_1.mod.log[["effects"]][["horizonO"]])

summary(PPO.mbc.ln_1.mod.log)
anova(PPO.mbc.ln_1.mod.log)

aov.mbc <- aov(PPO.mbc.ln_1.mod.log)

TukeyHSD(aov.mbc, which = 'forest.type')
HSD.test(aov.mbc, trt = c("forest.type"), console = TRUE)

TukeyHSD(aov.mbc, which = 'horizon')
HSD.test(aov.mbc, trt = c("horizon"), console = TRUE)



## presenting results


# horizon plot

EEA_2.mbc <- EEA_1.mbc

EEA_2.mbc <-
  EEA_2.mbc %>%
  filter(enzyme !='PPO')

enzyme_sum.mbc<- 
  EEA_2.mbc %>% 
  group_by(enzyme) %>% 
  summarise(mean_activity = mean(EEA.mbc, na.rm = T), n=n(), sd=sd(EEA.mbc, na.rm = T), se=sd/sqrt(n))
enzyme_sum.mbc

enzyme_sum.mbc.horizon <- 
  EEA_2.mbc %>% 
  group_by(enzyme, horizon) %>% 
  summarise(mean_activity = mean(EEA.mbc, na.rm = T), n=n(), sd=sd(EEA.mbc, na.rm = T), se=sd/sqrt(n))
enzyme_sum.mbc.horizon


enzyme_sum.mbc.forest.type <- 
  EEA_2.mbc %>% 
  group_by(enzyme, forest.type) %>% 
  summarise(mean_activity = mean(EEA.mbc, na.rm = T), n=n(), sd=sd(EEA.mbc, na.rm = T), se=sd/sqrt(n))
enzyme_sum.mbc.forest.type

aov.mbc <- aov(EEA.per.g.soil.ANCOVA.model)

TukeyHSD(aov.mbc, which = 'enzyme')
HSD.test(aov.mbc, trt = c("enzyme"), console = TRUE)

TukeyHSD(aov.mbc, which = 'horizon:enzyme')
HSD.test(aov.mbc, trt = c("horizon", "enzyme"), console = TRUE)


enzyme_sum.mbc.horizon$cld <- c("c","a","c","b","c","b","c","c","c","ab")

# make new function for powers on y-axis

scientific_10 <- 
  function(x) {   parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x))) }

enzyme.per.g.mbc.horizon.plt <- ggplot(enzyme_sum.mbc.horizon, aes(x=enzyme, y=mean_activity, fill=horizon, group = horizon)) + 
  geom_bar(position=position_dodge(), stat="identity", colour = "black") +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(data = enzyme_sum.mbc.horizon,
                aes(ymin = mean_activity-se, ymax= mean_activity+se), 
                width = .2, position=position_dodge(width=0.9))+
 geom_text(data = enzyme_sum.mbc.horizon, aes(x = enzyme, y = mean_activity+se+0.00004, label = cld, group = horizon),
                  position = position_dodge(0.9),
            inherit.aes = TRUE)+
  #coord_flip() +
  xlab(expression(paste("Enzyme")))+ 
  ylab(expression(paste('Enzyme activity (',mu,'mol g'^-1," MBC h"^-1,")")))+
  scale_fill_manual(name="Horizon", # Legend label, use darker colors
                    breaks=c("A", "O"),
                    labels=c("A", "O"),
                    values=c("#CC9900","#990000")) +
  scale_y_continuous(label=scientific_10)+
  expand_limits(y=c(0,0.001))+
  theme_classic() +   
  theme(legend.background = element_rect(size=0.3, linetype="solid", 
        colour ="black"))+
  theme(legend.position="none")


enzyme.per.g.mbc.horizon.plt

# forest type plot

aov.mbc <- aov(EEA.per.g.soil.ANCOVA.model)


TukeyHSD(aov.mbc, which = 'forest.type:enzyme')
HSD.test(aov.mbc, trt = c("forest.type", "enzyme"), console = TRUE)

enzyme_sum.mbc.forest.type <- 
  EEA_2.mbc %>% 
  group_by(enzyme, forest.type) %>% 
  summarise(mean_activity = mean(EEA.mbc, na.rm = T), n=n(), sd=sd(EEA.mbc, na.rm = T), se=sd/sqrt(n))
enzyme_sum.mbc.forest.type

enzyme_sum.mbc.forest.type$cld <- c("de","a","de","abc","e","abcd","e","cde","bcde","ab")

# make new function for powers on y-axis

scientific_10 <- 
  function(x) {   parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x))) }

enzyme.per.g.mbc.forest.type.plt <- 
  ggplot(enzyme_sum.mbc.forest.type, aes(x=enzyme, y=mean_activity, fill=forest.type)) + 
  geom_bar(position=position_dodge(), stat="identity", colour = "black") +
  geom_bar(position=position_dodge(), stat="identity") +  
  geom_errorbar(data = enzyme_sum.mbc.forest.type,
                aes(ymin = mean_activity-se, ymax= mean_activity+se), 
                width = .2, position=position_dodge(width=0.9))+
  geom_text(data = enzyme_sum.mbc.forest.type,  aes(x = enzyme, y = mean_activity+se+0.00003, label = cld, group = forest.type),
           position = position_dodge(0.9),
              inherit.aes = TRUE)+
  xlab(expression(paste("Enzyme")))+ 
  ylab(expression(paste('Enzyme activity (',mu,'mol g'^-1," MBC h"^-1,")")))+
  #coord_flip() +
  expand_limits(y=c(0,5e-4))+
  scale_y_continuous(label=scientific_10)+
  theme(axis.title.x=element_blank())+
  theme(legend.position="none")+
  scale_fill_manual(name="Forest type", # Legend label, use darker colors
                    breaks=c("broadleaf", "evergreen"),
                    labels=c("broadleaf", "coniferous"),values=c("#56B4E9","#0072B2")) +
  theme_classic() +   
  theme(legend.background = element_rect(size=0.3, linetype="solid", 
        colour ="black"))+
  theme(legend.position="none")


enzyme.per.g.mbc.forest.type.plt


## PPO graphs

PPO_2.mbc <- PPO_1.mbc

PPO_2.mbc <-
  PPO_2.mbc %>%
  filter(enzyme =='PPO')

PPO_sum.mbc.horizon <- 
  PPO_2.mbc %>% 
  group_by(enzyme, horizon) %>% 
  summarise(mean_activity = mean(PPO.microbioc, na.rm = T), n=n(), sd=sd(PPO.microbioc, na.rm = T), se=sd/sqrt(n))
PPO_sum.mbc.horizon

# enzyme_sum.mbc.horizon$cld <- c("b","a","bc","a","c","a","b","a","bc","a")

# make new function for powers on y-axis

scientific_10 <- 
  function(x) {   parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x))) }

PPO.per.g.mbc.horizon.plt <- ggplot(PPO_sum.mbc.horizon, aes(x=enzyme, y=mean_activity, fill=horizon, group = horizon)) + 
  geom_bar(position=position_dodge(), stat="identity", colour = "black", width = 0.5) +
  geom_bar(position=position_dodge(), stat="identity", width = 0.5) +
  geom_errorbar(data = PPO_sum.mbc.horizon,
                aes(ymin = mean_activity-se, ymax= mean_activity+se), 
                width = .2, position=position_dodge(width=0.5))+
  #coord_flip() +
  xlab(expression(paste("")))+ 
  ylab(expression(paste('Enzyme activity (',mu,'mol g'^-1," MBC h"^-1,")")))+  scale_fill_manual(name="Horizon", # Legend label, use darker colors
                    breaks=c("A", "O"),
                    labels=c("A", "O"),
                    values=c("#CC9900","#990000")) +
  scale_y_continuous(label=scientific_10,
                     breaks=seq(0, 7.5e-3, 2.5e-3))+
  expand_limits(y=c(0,7.5e-3))+
  theme_classic() +   
  theme(legend.position="none",legend.background = element_rect(size=0.3, linetype="solid", 
                                                              colour ="black"))


PPO.per.g.mbc.horizon.plt

# forest type plot

PPO_sum.mbc.forest.type <- 
  PPO_2.mbc %>% 
  group_by(enzyme, forest.type) %>% 
  summarise(mean_activity = mean(PPO.microbioc, na.rm = T), n=n(), sd=sd(PPO.microbioc, na.rm = T), se=sd/sqrt(n))
PPO_sum.mbc.forest.type

# PPO_sum.mbc.forest.type$cld <- c("de","ab","cde","bcd","e","bcd","de","a","cde","bc")

# make new function for powers on y-axis

scientific_10 <- 
  function(x) {   parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x))) }

PPO.per.g.mbc.forest.type.plt <- 
  ggplot(PPO_sum.mbc.forest.type, aes(x=enzyme, y=mean_activity, fill=forest.type)) + 
  geom_bar(position=position_dodge(), stat="identity", colour = "black", width = 0.5) +
  geom_bar(position=position_dodge(), stat="identity", width = 0.5) +  
  geom_errorbar(data = PPO_sum.mbc.forest.type,
                aes(ymin = mean_activity-se, ymax= mean_activity+se), 
                width = .2, position=position_dodge(width=0.5))+
  xlab(expression(paste("")))+ 
  ylab(expression(paste('Enzyme activity (',mu,'mol g'^-1," MBC h"^-1,")")))+
  expand_limits(y=c(0,4e-3))+
  scale_y_continuous(label=scientific_10)+
  theme(axis.title.x=element_blank())+
  scale_fill_manual(name="Forest type", # Legend label, use darker colors
                    breaks=c("broadleaf", "evergreen"),
                    labels=c("broadleaf", "coniferous"),values=c("#56B4E9","#0072B2")) +
  theme_classic() +   
  theme(legend.position="none",legend.background = element_rect(size=0.3, linetype="solid", 
                                                             colour ="black"))
                            
  

PPO.per.g.mbc.forest.type.plt


#### make grid of plots

plot_grid(PPO.per.g.mbc.horizon.plt, PPO.per.g.mbc.forest.type.plt,
          enzyme.per.g.mbc.horizon.plt, enzyme.per.g.mbc.forest.type.plt,
          nrow = 2, labels = c("c","d"), label_size = 14)


