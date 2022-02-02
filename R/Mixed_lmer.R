library("lme4")
library("ggplot2")

df <- read.csv("../data.csv")
df[,'Genotype']<-factor(df[,'Genotype'])
head(df,20)

# Mixedeffects regression (with interaction)
mixed <- lmer(WAKE ~ Treatment + Genotype + Treatment:Genotype + (1 | id), df)
summary(mixed)
# Diagnostic Plots for Linear Regression Analysis
par(mfrow = c(2,2))
plot(mixed, col = df$Treatment+1)

# compute MSE (~44.66147)
y = df$WAKE
ypred = predict(mixed)
MSE = mean((ypred-y)^2)

# Scatter plot
plot(y, ypred,
     pch = 19,
     col = factor(df$Treatment))

# Legend
legend("topleft",
       legend = levels(factor(df$Treatment)),
       pch = 19,
       col = factor(levels(factor(df$Treatment))))





