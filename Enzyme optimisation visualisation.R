# 22/02/22

library(ggplot2)
library (dplyr)
library(readr)
library("tidyr")

enzyme.vis <- read_csv("General Waring/Project/enzyme optimisation visualisation.csv")
View(enzyme.vis)

### drop NAs
enzyme.vis <- 
  enzyme.vis %>% 
  drop_na()
View(enzyme.vis)

enzyme.vis$site_horizon=paste(enzyme.vis$site,enzyme.vis$horizon)

enzyme.sum <-
  enzyme.vis %>% 
  group_by(site,concentration,horizon,site_horizon) %>% 
  summarise(mean.absorbance = mean(absorbance), std = sd(absorbance)) # he warning message is a friendly one. We can specify the .groups with one of the values i.e. 'drop' removes the group attribute. By default, it removes the last group

dev.off() # remove all current plots

 site_list= unique(enzyme.sum$site)
 for (i in site_list ){
p=ggplot(data=subset(enzyme.sum,site==i),aes(x = concentration, y = mean.absorbance,fill = horizon))+
  geom_bar(position = "dodge",
               stat="identity")+
geom_errorbar(aes(ymin = mean.absorbance - std, ymax = mean.absorbance + std),
              position = position_dodge(width = 0.9, preserve = "single"), width = 0.3)+
  scale_x_discrete(limits=c(1,2,5,10))+
  labs(y="Mean Absorbance",x="Concentration")+
  ggtitle(i)+
  theme_classic()
print (p)
 }


########## 2-way ANOVA
enzyme.vis$concentration=factor(enzyme.vis$concentration) # convert concentration to a factor so R recognizes it as a discrete variable rather than continuous
soil_model <- lm(absorbance ~ horizon + concentration + horizon : concentration, data = enzyme.vis)

#### assumptions

  ## normality
plot(soil_model, which = 2, add.smooth = FALSE)

  ## variance
plot(soil_model, which = 3, add.smooth = FALSE)

### do ANOVA
anova(soil_model)

### multiple comparisons test
soil_aov <- aov(soil_model)

TukeyHSD(soil_aov, which = 'horizon:concentration')


# one-way ANOVA: concentration
enzyme.vis$site_horizon=paste(enzyme.vis$site,enzyme.vis$horizon)
conc_list <- unique(enzyme.vis$site_horizon)
for (j in conc_list){
  conc_model <- lm(absorbance ~ concentration, data=subset(enzyme.vis,site_horizon==j))
conc_anova <- 
  anova(conc_model)
print(conc_anova)
}
  

# one-way ANOVA: horizon

enzyme.vis$horizon <- as.character(enzyme.vis$horizon)

# subset for sitesthat ontain both o and a horizons????
a_sum=subset(enzyme.vis, horizon=="a")
o_sum=subset(enzyme.vis, horizon=="o")
common_site=intersect(a_sum$site,o_sum$site)
common_site
common=subset(enzyme.vis,site%in%common_site)

hor_list <- unique(common$site)
for (k in hor_list){
  hor_model <- lm(absorbance ~ horizon, data=subset(common,site==k))
  hor_anova <- anova(hor_model)
  print(hor_anova)
}

