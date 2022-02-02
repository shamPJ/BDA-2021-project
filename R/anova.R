library("lme4")
library("ggplot2")

df <- read.csv("data.csv")
df[,'Genotype']<-factor(df[,'Genotype'])
head(df,20)

# Fixed effects regression (with interaction)
fixed <- lm(WAKE~Treatment + Genotype + Treatment:Genotype, data = df)
summary(fixed)

# Mixedeffects regression (with interaction)
mixed <- lmer(WAKE ~ Treatment + Genotype + Treatment:Genotype + (1 | id), df)
summary(mixed)

anova(mixed, fixed)
