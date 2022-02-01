data {
  int<lower=0> N;                         // n.o. samples
  real y[N];                              // measurements, response var
  real treat[N];                          // predictor treatment
  int<lower=0,upper=1> gen[N];            // predictor genotype
}


parameters {
  real<lower=0> sigma;    // sheared standard deviation
  vector[4] beta;         // intercept and slopes
}


model {
  real mu;
  sigma ~ normal(0,10);     // weakly informative prior
  beta ~ normal(0, 100);    // weakly informative prior
  
  for (i in 1:N){
    mu = beta[1] + beta[2]*treat[i] + beta[3]*gen[i] + beta[4]*treat[i]*gen[i];
    y[i] ~ normal(mu, sigma);
  }
}


generated quantities {
  // the log-likelihood values of each observation for every posterior draw
  
  real log_lik[N];
  real mu;
  
  for (i in 1:N) {
    mu = beta[1] + beta[2]*treat[i] + beta[3]*gen[i] + beta[4]*treat[i]*gen[i];
    log_lik[i] = normal_lpdf(y[i] | mu, sigma);
  }
}
