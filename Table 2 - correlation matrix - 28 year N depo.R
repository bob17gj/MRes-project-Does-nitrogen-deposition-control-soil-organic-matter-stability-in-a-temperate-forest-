library(tidyverse) #for all data wrangling
library(Hmisc) # for calculating significance levels
library(sjPlot) # for producing matrix tables for publication
library(webshot) # for saving matrices
library(stargazer)

rm(list=ls())
setwd("~/General Waring/Project")

hydrolytic_enzymes <- read.csv("Correlation - 28 year N depo - master sheet for R.csv")
str(hydrolytic_enzymes)

#### broadleaf correlation ####


correlation <- hydrolytic_enzymes
correlation.broadleaf <- correlation %>% filter(forest.type == "deciduous") 
names(correlation.broadleaf) <- gsub("\\.", " ", names(correlation.broadleaf)) # remove . from column names
correlation.broadleaf <- correlation.broadleaf[,c("BG","BX", "NAG", "LAP","AP", 
                                                  "vector length", "vector angle",
                                                  "N deposition","SOC","C","N","pH","MAT","MAP","GSM","site")]




hist(correlation.broadleaf) # looking at distribution of variables
# none are normally distributed - Spearman's correlation necessary

correlation.broadleaf$P <- as.numeric(as.character(correlation.broadleaf$P))

correlation.broadleaf$SOC <- as.numeric(as.character(correlation.broadleaf$SOC))

####convert NA values to 0####
correlation.broadleaf[is.na(correlation.broadleaf)] <- 0

glimpse(correlation.broadleaf)

map_broadleaf.spearman <-as.dist(round(cor(correlation.broadleaf[,1:16], method = "spearman"),2))
map_broadleaf.spearman

##### broadleaf significance levels ####

map_broadleaf_rcorr <- rcorr(as.matrix(correlation.broadleaf[,1:16])) # p-values = 0 when the p-value is very small
map_broadleaf_rcorr

#### broadleaf matrix table ####

cor.table.broadleaf <- 
  tab_corr(correlation.broadleaf,  
           na.deletion = c("listwise"),
           corr.method = c("spearman"),
           title = "Spearman's correlation coefficients 
        relating soil microbe and enzyme activities to edaphic variables 
        across temperate broadleaf sites",
           var.labels = NULL,
           wrap.labels = 40,
           show.p = TRUE,
           p.numeric = FALSE,
           fade.ns = TRUE,
           val.rm = FALSE,
           digits = 3,
           triangle = "lower",
           string.diag = NULL,
           CSS = NULL,
           encoding = NULL,
           file = NULL,
           use.viewer = TRUE,
           remove.spaces = TRUE
  )

cor.table.broadleaf



#### coniferous correlation ####

correlation <- hydrolytic_enzymes
correlation.coniferous <- correlation %>% filter(forest.type == "evergreen") 
names(correlation.coniferous) <- gsub("\\.", " ", names(correlation.coniferous)) # remove . from column names
correlation.coniferous <- correlation.coniferous[,c("BG","BX", "NAG", "LAP","AP", 
                                                    "vector length", "vector angle",
                                                    "N deposition","SOC","C","N","pH","MAT","MAP","GSM","site")]




hist(correlation.coniferous) # looking at distribution of variables
# none are normally distributed - Spearman's correlation necessary


####convert NA values to 0####
correlation.coniferous[is.na(correlation.coniferous)] <- 0

glimpse(correlation.coniferous)

map_coniferous.spearman <-as.dist(round(cor(correlation.coniferous[,1:16], method = "spearman"),2))
map_coniferous.spearman

##### coniferous significance levels ####

map_coniferous_rcorr <- rcorr(as.matrix(correlation.coniferous[,1:17])) # p-values = 0 when the p-value is very small
map_coniferous_rcorr

#### coniferous matrix table ####

cor.table.coniferous <- 
  tab_corr(correlation.coniferous,  
           na.deletion = c("listwise"),
           corr.method = c("spearman"),
           title = "Spearman's correlation coefficients 
        relating soil microbe and enzyme activities to edaphic variables 
        across coniferous sites",
           var.labels = NULL,
           wrap.labels = 40,
           show.p = TRUE,
           p.numeric = FALSE,
           fade.ns = TRUE,
           val.rm = FALSE,
           digits = 3,
           triangle = "lower",
           string.diag = NULL,
           CSS = NULL,
           encoding = NULL,
           file = NULL,
           use.viewer = TRUE,
           remove.spaces = TRUE
  )

cor.table.coniferous


# O horizon matrix table ####

correlation <- hydrolytic_enzymes
correlation.O <- correlation %>% filter(horizon == "O") 
names(correlation.O) <- gsub("\\.", " ", names(correlation.O)) # remove . from column names
correlation.O <- correlation.O[,c("BG","BX", "NAG", "LAP","AP", 
                                  "vector length", "vector angle",
                                  "N deposition","SOC","C","N","pH","MAT","MAP","GSM","site")]




hist(correlation.O) # looking at distribution of variables
# none are normally distributed - Spearman's correlation necessary


####convert NA values to 0####
correlation.O[is.na(correlation.O)] <- 0

glimpse(correlation.O)

map_O.spearman <-as.dist(round(cor(correlation.O[,1:17], method = "spearman"),2))
map_O.spearman

##### O hotizon significance levels ####

map_O_rcorr <- rcorr(as.matrix(correlation.O[,1:17])) # p-values = 0 when the p-value is very small
map_O_rcorr

#### O horizon matrix table ####

cor.table.O <- 
  tab_corr(correlation.O,  
           na.deletion = c("listwise"),
           corr.method = c("spearman"),
           title = "Spearman's correlation coefficients 
        relating soil microbe and enzyme activities to edaphic variables 
        in the O horizon across sites",
           var.labels = NULL,
           wrap.labels = 40,
           show.p = TRUE,
           p.numeric = FALSE,
           fade.ns = TRUE,
           val.rm = FALSE,
           digits = 3,
           triangle = "lower",
           string.diag = NULL,
           CSS = NULL,
           encoding = NULL,
           file = NULL,
           use.viewer = TRUE,
           remove.spaces = TRUE
  )

cor.table.O



# A horizon matrix table ####

correlation <- hydrolytic_enzymes
correlation.A <- correlation %>% filter(horizon == "A") 
names(correlation.A) <- gsub("\\.", " ", names(correlation.A)) # remove . from column names
correlation.A <- correlation.A[,c("BG","BX", "NAG", "LAP","AP", 
                                  "vector length", "vector angle",
                                  "N deposition","SOC","C","N","pH","MAT","MAP","GSM","site")]




hist(correlation.O) # looking at distribution of variables
# none are normally distributed - Spearman's correlation necessary


####convert NA values to 0####
correlation.A[is.na(correlation.A)] <- 0

glimpse(correlation.A)

map_A.spearman <-as.dist(round(cor(correlation.A[,1:17], method = "spearman"),2))
map_A.spearman

##### A hotizon significance levels ####

map_A_rcorr <- rcorr(as.matrix(correlation.A[,1:17])) # p-values = 0 when the p-value is very small
map_A_rcorr

#### A horizon matrix table ####

cor.table.A <- 
  tab_corr(correlation.A,  
           na.deletion = c("listwise"),
           corr.method = c("spearman"),
           title = "Spearman's correlation coefficients 
        relating soil microbe and enzyme activities to edaphic variables 
        in the A horizon across sites",
           var.labels = NULL,
           wrap.labels = 40,
           show.p = TRUE,
           p.numeric = FALSE,
           fade.ns = TRUE,
           val.rm = FALSE,
           digits = 3,
           triangle = "lower",
           string.diag = NULL,
           CSS = NULL,
           encoding = NULL,
           file = NULL,
           use.viewer = TRUE,
           remove.spaces = TRUE
  )

cor.table.A


# whole correlation table ####

correlation <- hydrolytic_enzymes
names(correlation) <- gsub("\\.", " ", names(correlation)) # remove . from column names
correlation <- correlation[,c("AP","BG","BX", "NAG", "LAP", "PPO", 
                              "vector length", "vector angle",
                                                      "MBC",
                                                  "N deposition","SOC","C","N","P",
                                                  "GSM","pH","MAP","MAT")]




hist(correlation) # looking at distribution of variables
# none are normally distributed - Spearman's correlation necessary



####convert NA values to 0####
#correlation[is.na(correlation)] <- 0

glimpse(correlation)

map.spearman <-as.dist(round(cor(correlation[,1:17], method = "spearman"),2))
map.spearman

##### significance levels ####

map_rcorr <- rcorr(as.matrix(correlation[,1:17])) # p-values = 0 when the p-value is very small
map_rcorr

#### matrix table ####

cor.table <- 
  tab_corr(correlation,  
           na.deletion = c("listwise"),
           corr.method = c("spearman"),
           title = NULL,
           var.labels = NULL,
           wrap.labels = 40,
           show.p = TRUE,
           p.numeric = FALSE,
           fade.ns = TRUE,
           val.rm = FALSE,
           digits = 3,
           triangle = "upper",
           string.diag = NULL,
           CSS = NULL,
           encoding = NULL,
           file = NULL,
           use.viewer = TRUE,
           remove.spaces = T
  )

cor.table

# Save it


# first save table to html file
tab_model(cor.table, file = "cor.table.html")

# then take this html file and make .png file
webshot("cor.table.html", "cor.table.png")


#save_plot("edaphic and climate variable Spearman's correlation matrix.png",
 # cor.table,
#  base_asp = 1.1)

# save_plot("edaphic and climate variable Spearman's correlation matrix.png", fig = ggplot2::last_plot(), width = 12, height = 9,
         # dpi = 300, theme = ggplot2::theme_get(), label.size = 2.4,
          # axis.textsize = 0.8, axis.titlesize = 0.75, legend.textsize = 0.6,
          # legend.titlesize = 0.65, legend.itemsize = 0.5)



# whole ratio correlation table ####

correlation <- hydrolytic_enzymes
names(correlation) <- gsub("\\.", " ", names(correlation)) # remove . from column names

C.N <- as.numeric("C.N")
C.P <- as.numeric("C.P")
N.P <- as.numeric("N.P")

ratio.correlation <- correlation[,c("AP","BG","BX", "NAG", "LAP", "PPO", 
                              "vector length", "vector angle",
                              "MBC",
                              "N deposition","SOC","C","N","P",
                              "C.N", "C.P", "N.P",
                              "pH","GSM","MAP","MAT")]




hist(ratio.correlation) # looking at distribution of variables
# none are normally distributed - Spearman's correlation necessary



####convert NA values to 0####
#correlation[is.na(correlation)] <- 0

glimpse(ratio.correlation)

map.ratio.spearman <-as.dist(round(cor(ratio.correlation[,1:21], method = "spearman"),2))
map.ratio.spearman

##### significance levels ####

map_rcorr <- rcorr(as.matrix(correlation[,1:21])) # p-values = 0 when the p-value is very small
map_rcorr

#### matrix table ####

cor.table <- 
  tab_corr(correlation,  
           na.deletion = c("listwise"),
           corr.method = c("spearman"),
           title = NULL,
           var.labels = NULL,
           wrap.labels = 40,
           show.p = TRUE,
           p.numeric = FALSE,
           fade.ns = TRUE,
           val.rm = FALSE,
           digits = 3,
           triangle = "upper",
           string.diag = NULL,
           CSS = NULL,
           encoding = NULL,
           file = NULL,
           use.viewer = TRUE,
           remove.spaces = T
  )

cor.table

# Save it
