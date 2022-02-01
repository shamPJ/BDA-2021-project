library("lme4")
library("ggplot2")
library("rstan")
library(aaltobda)
library(tidyr)
library(gridExtra)
library("bayesplot")

df <- read.csv("data.csv")
head(df,20)

data_list <- list(N=length(df$WAKE),
                  y=df$WAKE,
                  treat=df$Treatment,
                  gen=df$Genotype)
# flat priors
stanmodel <- "
data {
  int<lower=0> N;                         // n.o. samples
  real y[N];                              // measurements, response var
  real treat[N];                          // predictor treatment
  int<lower=0,upper=1> gen[N];            // predictor genotype
}

parameters {
  real<lower=0> sigma;    // sheared standard deviation, (0,+inf)
  vector[4] beta;         // intercept and slopes, (-inf,+inf)
}

model {
  real mu;
  for (i in 1:N){
      mu = beta[1] + beta[2]*treat[i] + beta[3]*gen[i] + beta[4]*treat[i]*gen[i];
      y[i] ~ normal(mu, sigma);
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
