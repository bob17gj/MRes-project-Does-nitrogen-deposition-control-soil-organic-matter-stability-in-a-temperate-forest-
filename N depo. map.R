library(tidyverse)
library(ggplot2)
library(maps)

rm(list=ls())
setwd("~/General Waring/MRes Project")

sites <- read_csv("Hydrolytic enzyme master sheet for R (2).csv")

worldmap = map_data('world')

ggplot() + 
  geom_polygon(data = worldmap, 
               aes(x = long, y = lat, group = group), 
               fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.3, xlim = c(-10,3), ylim = c(50, 59)) + 
  geom_point(data = sites, aes(x = longitude, y = latitude, color=soil.type, shape =forest.type), 
             size = 4)+
  labs(shape = "Forest type", size = 10) + labs(colour = "Soil type", size = 10) +
  scale_colour_manual(name="Soil type", # Legend label, use darker colors
                      breaks=c("peaty gley", "brown earth"),
                      labels=c("peaty gley", "brown earth"),
                      values=c("#CC9999","#660000"))+
  scale_shape_manual(name="Forest type", # Legend label, use darker colors
                      breaks=c("broadleaf", "evergreen"),
                      labels=c("broadleaf", "coniferous"),
                      values=c(16, 17))+
                    theme_void()+
  theme(legend.text=element_text(size=10),
        legend.position = c(0.9, 0.6),
        legend.background = element_rect(size=0.3, 
                                         linetype="solid", 
                                         colour ="black"),
        legend.margin = margin(6, 6, 6, 6))
               







