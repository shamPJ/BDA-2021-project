library(rstan)
library(tidyr)
library(aaltobda)
library(loo)

df <- read.csv("../data.csv")
data_list <- list(N=length(df$WAKE),
                  J = length(unique(df$id)),
                  y=df$WAKE,
                  id = df$id,
                  treat=df$Treatment,
                  gen=df$Genotype)

fixed <- stan(
  file = "fixed.stan",    # Stan program
  data = data_list ,      # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 1,              # number of cores (could use one per chain)
  refresh = 0,             # no progress shown
  seed=42
)

mixed <- stan(
  file = "mixed.stan",    # Stan program
  data = data_list ,      # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 1,              # number of cores (could use one per chain)
  refresh = 0,             # no progress shown
  seed=42
)

# Extract pointwise log-likelihood
ll_fixed <- extract_log_lik(fixed, merge_chains = FALSE)
ll_mixed <- extract_log_lik(mixed, merge_chains = FALSE)

# relative effective sample sizes
r_eff_fixed <- relative_eff(exp(ll_fixed), cores = 2)
r_eff_mixed <- relative_eff(exp(ll_mixed), cores = 2)

# loo
loo_fixed <- loo(ll_fixed, r_eff = r_eff_fixed, cores = 2)
loo_mixed <- loo(ll_mixed, r_eff = r_eff_mixed, cores = 2)

print(loo_fixed)
print(loo_mixed)

comp <- loo_compare(loo_fixed, loo_mixed)
comp

