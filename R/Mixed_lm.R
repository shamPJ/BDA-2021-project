library("lme4")
library("ggplot2")

df <- read.csv("data.csv")
df[,'Genotype']<-factor(df[,'Genotype'])
head(df,20)

# Mixedeffects regression (with interaction)
mixed <- lmer(WAKE ~ Treatment + Genotype + Treatment:Genotype + (1 | id), df)
summary(mixed)
# Diagnostic Plots for Linear Regression Analysis
par(mfrow = c(2,2))
plot(fixed, col = df$Treatment+1)

ggplot(df,aes(x=Treatment,y=WAKE, colour=Genotype, group=Genotype)) + 
  geom_point() + geom_smooth(method = "lm") 




