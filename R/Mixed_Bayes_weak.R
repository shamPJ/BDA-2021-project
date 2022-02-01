library("ggplot2")
library("rstan")
library(aaltobda)
library(tidyr)
library(gridExtra)
library(bayesplot)
require(GGally)

df <- read.csv("data.csv")
head(df,20)

data_list <- list(N=length(df$WAKE),
                  J = length(unique(df$id)),
                  y=df$WAKE,
                  id = df$id,
                  treat=df$Treatment,
                  gen=df$Genotype)
# flat priors
stanmodel <- "
data {
  int<lower=0> N;                         // n.o. samples
  int<lower=0> J;                         // n.o. subj
  real y[N];                              // measurements, response var
  int id[N];                              // subject id
  real treat[N];                          // predictor treatment
  int<lower=0,upper=1> gen[N];            // predictor genotype
}

parameters {
  real<lower=0> sigma;    // sheared standard deviation
  vector[4] beta;         // intercept and slopes
  vector[J] u;            // subj intercepts
}

model {
  real mu;
  sigma ~ normal(0,10);     // weakly informative prior
  beta ~ normal(0, 100);    // weakly informative prior
  u ~ normal(0,10);         // subj intercepts
  
  for (i in 1:N){
      mu = beta[1] + u[id[i]] + beta[2]*treat[i] + beta[3]*gen[i] + beta[4]*treat[i]*gen[i];
      y[i] ~ normal(mu, sigma);
      }
}

generated quantities {
  real ypred[N];
  for(i in 1:N){
      ypred[i] = beta[1] + u[id[i]] + beta[2]*treat[i] + beta[3]*gen[i] + beta[4]*treat[i]*gen[i];
  }
}
"

# fit data
fit <- stan(model_code = stanmodel, 
            data = data_list,
            iter = 2000,
            refresh=0,
            seed=42,
            verbose = FALSE)
print(fit, 'beta')

traceplot(fit, 'beta', inc_warmup = TRUE)
plot(fit)

posterior <- as.matrix(fit)
plot_title <- ggtitle("Posterior distributions",
                      "with medians and 90% intervals")
color_scheme_set(scheme = "blue")
mcmc_areas(posterior,
           pars = c("beta[1]", "beta[2]", "beta[3]", "beta[4]"),
           prob = 0.9) + plot_title
sumr = summary(fit, pars=c("ypred"))

plot(df$WAKE, sumr$summary[,1], pch = 19, col=factor(df$Treatment))

df[,'Genotype']<-factor(df[,'Genotype'])
ggplot(aes(x=df$WAKE, y=sumr$summary[,1], shape=df$Genotype, color=df$Treatment)) +
  geom_point()

