library("ggplot2")
library("tidyverse")
library("dplyr")
library("GGally")

df <- read.csv("data.csv")
head(df,20)

df %>% count(Treatment)

cols <- c('WAKE', 'NREM', 'REM', 'qsec')
ggpairs(df,columns = cols ggplot2::aes(colour=as.factor(Treatment)))