library("brms")
library("ggplot2")
library("bayesplot")

df <- read.csv("data.csv")
df[,'Genotype']<-factor(df[,'Genotype'])
head(df,20)


fit <- brm(formula = WAKE ~ Treatment + Genotype + Treatment:Genotype + (1 | id),
            data = df,
            warmup = 1000, iter = 2000, chains = 4)