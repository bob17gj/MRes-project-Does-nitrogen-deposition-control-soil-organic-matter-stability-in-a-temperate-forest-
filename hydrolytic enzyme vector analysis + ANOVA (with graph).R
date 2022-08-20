library(VecStatGraphs2D) # vector analysis package

rm(list=ls())
setwd("~/General Waring/Project")

hydrolytic_vector_enzymes <- read.csv("Vector analysis.csv")
str(hydrolytic_vector_enzymes)

#### fig. 4 vector plot

hydrolytic_vector_enzymes["y"][hydrolytic_vector_enzymes["y"]==1]<-NA # removed values = 1, as there are no relative enzyme activities to compare
hydrolytic_vector_enzymes["x"][hydrolytic_vector_enzymes["x"]==1]<-NA

vector.plot <- ggplot(hydrolytic_vector_enzymes, aes(x = x, y = y, colour = horizon)) + 
  geom_point(aes(shape=forest.type),alpha = 0.8) + 
  geom_abline()+
  xlab("[C/(C+P)]") + ylab("[C/(C+N)]")+
  theme_classic()

vector.plot


## Fig. 5 ANOVA

