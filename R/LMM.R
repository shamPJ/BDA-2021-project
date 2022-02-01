library("lme4")
library("ggplot2")

df <- read.csv("data.csv")
df[,'Genotype']<-factor(df[,'Genotype'])
head(df,20)

# Fixed effects model
fixed <- lm(WAKE~Treatment + Genotype + Treatment:Genotype, data = df)
summary(fixed)
ggplot(df,aes(x=Treatment,y=WAKE, col=Genotype)) + geom_point() + geom_smooth(method = "lm")

# Plot data separately
ggplot(df,aes(x=Treatment,y=WAKE)) + geom_smooth(method = "lm",level = 0.95) + 
  geom_point() + facet_wrap(~id) + scale_x_discrete(name ="Dose (mg)", 
                                                      limits=c(0,3,10))

# Random effects
mixed <- lmer(WAKE ~ Treatment + Genotype + Treatment:Genotype + (1 | id), df)
summary(mixed)

sqrt(sum(residuals(fixed)^2)/(dim(df)[1]-2))
sqrt(sum(resid(mixed)^2)/(dim(df)[1]-2))

mixed <- lmer(WAKE ~ Treatment + Genotype + Treatment:Genotype + (1 | id), df, REML = FALSE)
anova(mixed, fixed)

# Estimate confidence intervals with bootstrapping
# fixef() returns intercept & slope

N_boot <- 1000
#DEFINE MATRICES OF POPULATION/INDIVIDUAL LEVEL DATA POINTS AND DO BOOTSTRAPPING
intercept_fixef <- vector(); slope_fixef <- vector()
newdata_pop_level <- data.frame(Treatment = c(0,3,10))
datapoints_fixef <- matrix(ncol = N_boot, nrow = dim(newdata_pop_level)[1])
newdata_individ_level<-data.frame(Treatment=df$Treatment, Genotype=df$Genotype, id=df$id)


datapoints_individ_level <- matrix(ncol=N_boot, nrow=dim(newdata_individ_level)[1])
for(i in 1:N_boot)
{
  # choose 25% of data to fit to LMM
  df_boot <- df[sample(1:dim(df)[1], dim(df)[1]*0.75, 
                                       replace = TRUE),]
  lmerfit_boot <- lmer(WAKE ~ Treatment + Genotype + Treatment:Genotype + (1 | id), df_boot)
  intercept_fixef <- append(intercept_fixef, as.numeric(fixef(lmerfit_boot)[1]))
  slope_fixef <- append(slope_fixef, as.numeric(fixef(lmerfit_boot)[2]))
  datapoints_individ_level[,i] <- predict(lmerfit_boot, newdata_individ_level, 
                                          allow.new.levels = TRUE)
}

fixef_df <- data.frame(Intercept = intercept_fixef, Treatment = slope_fixef)
for(i in 1:N_boot) 
{
  datapoints_fixef[,i] <- 
    # get Treatment 
    model.matrix(~ 1 + Treatment, data = newdata_pop_level) %*% as.matrix(fixef_df)[i,]
}
# %*% is matrix multiplication

bootstrap_conf_interval_pop_level <- data.frame(Treatment = c(0,3,10),
                                                lwr=apply(datapoints_fixef,1,quantile,prob=.05),
                                                fit=apply(datapoints_fixef,1,quantile,prob=.5),
                                                upr=apply(datapoints_fixef,1,quantile,prob=.95))
#PLOT POPULATION LEVEL CONFIDENCE INTERVALS
p1 <- ggplot(df, aes(x = Treatment, y = WAKE)) + geom_point(shape = 1) + 
  geom_abline(data=fixef_df,alpha=0.1,size=2,aes(intercept=Intercept,slope=Treatment)) + 
  geom_smooth(method = "lm", color = "blue", size = 0.5)
p2 <- ggplot(df, aes(x = Treatment, y = WAKE)) + geom_point(shape = 1) + 
  geom_smooth(method = "lm", color = "blue", size = 0.5) + 
  geom_line(data=bootstrap_conf_interval_pop_level,aes(y=fit),size=0.5,color="red") + 
  geom_line(data=bootstrap_conf_interval_pop_level,aes(y=lwr),lty=2, color = "red") + 
  geom_line(data=bootstrap_conf_interval_pop_level,aes(y=upr),lty=2, color = "red")
library("gridExtra")
grid.arrange(p1, p2, nrow = 1)

#PLOT HISTOGRAMS OF SLOPES AND INTERCEPTS OF FIXED EFFECTS
par(mfrow = c(1, 2))
hist(fixef_df$Intercept, breaks = 100, col = "darkgreen", xlab = "Intercept", main = "Bootstrapped Intercept Values")
hist(fixef_df$Treatment, breaks = 100, col = "darkred", xlab = "Slope / Treatment", main = "Bootstrapped Slope Values")

bootstrap_conf_interval_individ_level <- data.frame(Treatment = df$Treatment, 
                                                    id = df$id,
                                                    lwr = apply(datapoints_individ_level, 1, 
                                                                quantile, prob = 0.05),
                                                    fit = apply(datapoints_individ_level, 1, 
                                                                quantile, prob = 0.5),
                                                    upr = apply(datapoints_individ_level, 1, 
                                                                quantile, prob = 0.95))
#PLOT INDIVIDUAL LEVEL CONFIDENCE INTERVALS
ggplot(df, aes(x = Treatment, y = WAKE)) + geom_point(size = 1) + 
  geom_smooth(method = "lm", level = 0.95, size = 0.5) + 
  facet_wrap(~id) + 
  geom_line(data = bootstrap_conf_interval_individ_level,aes(y=fit), size = 0.5, 
            color = "red") + 
  geom_line(data = bootstrap_conf_interval_individ_level,aes(y=lwr),lty=2,size=0.5, 
            color = "red") + 
  geom_line(data = bootstrap_conf_interval_individ_level,aes(y=upr),lty=2,size=0.5, 
            color = "red")

