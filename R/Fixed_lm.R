library("lme4")
library("ggplot2")

df <- read.csv("data.csv")
df[,'Genotype']<-factor(df[,'Genotype'])
head(df,20)

# Fixed effects regression (with interaction)
fixed <- lm(WAKE~Treatment + Genotype + Treatment:Genotype, data = df)
summary(fixed)
# Diagnostic Plots for Linear Regression Analysis
par(mfrow = c(2,2))
plot(fixed, col = df$Treatment+1)

ggplot(df,aes(x=Treatment,y=WAKE, colour=Genotype, group=Genotype)) + 
       geom_point() + geom_smooth(method = "lm") 

# compute MSE (~59.43669)
y = df$WAKE
ypred = predict(fixed)
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

