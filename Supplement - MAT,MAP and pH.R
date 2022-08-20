library(tidyverse) #for all data wrangling
library(cowplot) #for manuscript ready figures
library(car)# for bonferroni outlier test
library(multcomp) # for acova tukey test



rm(list=ls())
setwd("~/General Waring/MRes Project")

hydrolytic_enzymes <- read.csv("Hydrolytic enzyme master sheet for R (2).csv")
str(hydrolytic_enzymes)




#### MAT ####

pH.H2O.MAT.model <- lm(pH.H2O~MAT, data = hydrolytic_enzymes)

# remove outliers

outlierTest(pH.H2O.MAT.model,cutoff=Inf, pH.H2O~1) # no anomalous response values


# diagnostics

par(mfrow=c(2,2))
plot(pH.H2O.MAT.model)
par(mfrow=c(1,1))

# normality assumption not met, but variance assumption, which is more important, is met

# ANCOVA

anova(pH.H2O.MAT.model)

# make plot

pH.H2O.MAT.plt <- ggplot(pH.H2O.MAT.model, aes(x = MAT, y = pH.H2O))+
  geom_point()+
  stat_smooth(method="lm", se=T, color = "red") +
  labs(title="",
       y ="pH", x = "MAT (\u00B0C)")+
  expand_limits(y=c(2,8))+
  expand_limits(x=c(6,11))+
  scale_y_continuous(breaks=c(2,3,4,5,6,7,8))+
  scale_x_continuous(breaks=c(6,7,8,9,10,11))+ 
  theme_classic()+
  theme(axis.title.y=element_blank())

pH.H2O.MAT.plt





#### MAP ####

pH.H2O.MAP.model <- lm(pH.H2O~MAP, data = hydrolytic_enzymes)

# remove outliers

ggplot(pH.H2O.MAP.model, aes(x= MAP , y=pH.H2O)) + 
  geom_boxplot()

outlierTest(pH.H2O.MAP.model,cutoff=Inf, MAP~1) 

hydrolytic_enzymes_1=hydrolytic_enzymes

hydrolytic_enzymes_1 <- hydrolytic_enzymes_1[-c(71,227,383,539,695),] 


pH.H2O.MAP.model_1 <- lm(pH.H2O~MAP, data = hydrolytic_enzymes_1)

# convert pH.H2O to a factor
hydrolytic_enzymes_1$pH.H2O <- factor(hydrolytic_enzymes_1$pH.H2O)

# diagnostics

par(mfrow=c(2,2))
plot(pH.H2O.MAP.model_1) # large amounst of colinearity
par(mfrow=c(1,1))

vif(pH.H2O.MAP.model_1)

# normality assumption not met, but variance assumption, which is more important, is met

# ANCOVA

anova(pH.H2O.MAP.model)

# make plot

pH.H2O.MAP.plt <- ggplot(pH.H2O.MAP.model_1, aes(x = MAP, y = pH.H2O))+
  geom_point()+
  geom_smooth(method=lm,colour = "red")+
  labs(y ="pH", x = "MAP (mm)")+ 
  expand_limits(y=c(2,8))+
  expand_limits(x=c(500,2500))+
  scale_x_continuous(breaks=seq(500,2500,500))+
  scale_y_continuous(breaks=c(2,3,4,5,6,7,8))+ 
  theme_classic()

pH.H2O.MAP.plt


plot_grid(pH.H2O.MAP.plt ,pH.H2O.MAT.plt, ncol = 2, labels = c("auto"), label_size = 10)

