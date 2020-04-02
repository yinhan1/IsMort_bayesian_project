
library(tidyverse)
library(data.table)
library(magrittr)
library(ggsci)
library(kableExtra)

source("./scripts/functions.R")

#### ACO =========================== ####

#### step 1: filter out ACO and dummy categorical features ####

is_dummy <- 
  is_data %>% 
  filter(Population == "ACO") %>% 
  mutate(
    intercept = 1,
    infected = ifelse(Treatment == "Infected", 1, 0),
    female = ifelse(Sex == "Female", 1, 0),
    age_14 = ifelse(Age == 14, 1, 0),
    age_28 = ifelse(Age == 28, 1, 0),
    inf_female = infected*female, 
    inf_14 = infected*age_14,
    inf_28 = infected*age_28,
    female_14 = female*age_14,
    female_28 = female*age_28
  ) %>% 
  dplyr::select(c(Day, intercept, infected, female, age_14, age_28,
                  inf_female, inf_14, inf_28, female_14, female_28))

#### step 2: initial values and draws from MH algorithm ####

B <- 5
a1 <- a2 <- 0

outcomes <- matrix(0, nrow = B + 1, ncol = ncol(is_dummy)-1+2)
colnames(outcomes) <- c("alpha", "theta", colnames(is_dummy)[-1])
outcomes[1,] = c(1.8,1.7,1.8,rep(0, ncol(outcomes)-3))

#### step 3: for loop of MH algorithm

for(iter in 2:(B+1)){
  
  #### ---------------------- MH for shape ---------------------- ####
  # propose new shape
  log_shape_proposed = log_shape_sampler(log_alpha = log(outcomes[iter-1,1]),
                                         log_theta = log(outcomes[iter-1,2]),
                                         multiplier = 0.05
  )
  loglike = log_likelihood(data = is_dummy,
                           alpha = exp(log_shape_proposed[1]),
                           theta = exp(log_shape_proposed[2]),
                           coefs = outcomes[iter-1,3:ncol(outcomes)])
  if(loglike == 0){
    q1 = 0
  }else{
    # ratio of posterior 
    q1_post = 
      (
        loglike + 
          log(prior(alpha = exp(log_shape_proposed[1]),
                    theta = exp(log_shape_proposed[2]),
                    coefs = outcomes[iter-1,3:ncol(outcomes)])) -
          log_likelihood(data = is_dummy,
                         alpha = outcomes[iter-1,1],
                         theta = outcomes[iter-1,2],
                         coefs = outcomes[iter-1,3:ncol(outcomes)]) -
          log(prior(alpha = exp(log_shape_proposed[1]),
                    theta = exp(log_shape_proposed[2]),
                    coefs = outcomes[iter-1,3:ncol(outcomes)]))
      ) %>% 
      exp()
    # ratio of jacobian
    q1_jacobian = (exp(log_shape_proposed[1])*exp(log_shape_proposed[2])) / (outcomes[iter-1,1]*outcomes[iter-1,2])
    q1 = (q1_post * q1_jacobian) %>% as.numeric()
  }
  # update shape
  if (runif(1) < q1) {
    outcomes[iter,1] = exp(log_shape_proposed[1])
    outcomes[iter,2] = exp(log_shape_proposed[2])
    a1 = a1 + 1
  }else{
    outcomes[iter,1] = outcomes[iter-1,1]
    outcomes[iter,2] = outcomes[iter-1,2]
  }
  
  #### ---------------------- MH for scale/coefs ---------------------- ####
  # propose new coefs
  coefs_proposed = coefs_sampler(coefs = outcomes[iter-1,3:ncol(outcomes)], 
                                 multiplier = 0.1)
  loglike = log_likelihood(data = is_dummy,
                           alpha = outcomes[iter,1],
                           theta = outcomes[iter,2],
                           coefs = coefs_proposed)
  if(loglike == 0 || is.nan(loglike)){
    q2 = 0
  }else{
    # ratio of posterior 
    q2 = 
      (
        loglike + 
          log(prior(alpha = outcomes[iter,1],
                    theta = outcomes[iter,2],
                    coefs = coefs_proposed)) -
          log_likelihood(data = is_dummy,
                         alpha = outcomes[iter,1],
                         theta = outcomes[iter,2],
                         coefs = outcomes[iter-1,3:ncol(outcomes)]) -
          log(prior(alpha = outcomes[iter,1],
                    theta = outcomes[iter,2],
                    coefs = outcomes[iter-1,3:ncol(outcomes)]))
      ) %>% 
      exp()
  }
  # update scale
  if (runif(1) < q2) {
    outcomes[iter,3:ncol(outcomes)] = coefs_proposed
    a2 = a2 + 1
  }else{
    outcomes[iter,3:ncol(outcomes)] = outcomes[iter-1,3:ncol(outcomes)]
  }
  
  if (iter %% 10000 == 0) {cat(paste0("iteration: ", iter, "\n"))}
}

#### step 3: check acceptance rate and convergence ####

cat(paste0("acceptance rate for shape:  ", a1/B, "\n"))
cat(paste0("acceptance rate for coefs:  ", a2/B, "\n"))
as.data.frame(outcomes) %>% 
  mutate(iteration = c(1:nrow(outcomes))) %>% 
  reshape2::melt(id.vars = "iteration") %>%
  ggplot(aes(x = iteration, y = value, color = variable)) +
  geom_line(size = 0.1) +
  geom_hline(yintercept = 0, color = "brown3", size = 0.5, linetype = "dashed") +
  facet_wrap(~variable, ncol = 3, scales = "free_y") +
  labs(y = "coefficients") +
  theme_minimal() +
  theme(legend.position = "none")


#### step 4: save posterior dist and acceptance rate ####

outcomes = rbind(c(rep(a1/B, 2), rep(a2/B, ncol(outcomes)-2)), outcomes)
# write.csv(outcomes, file = "./data/posterior_aco.csv", row.names = FALSE)



#### CO ============================ ####

#### step 1: filter out CO and dummy categorical features ####

is_dummy <- 
  is_data %>% 
  filter(Population == "CO") %>% 
  mutate(
    intercept = 1,
    infected = ifelse(Treatment == "Infected", 1, 0),
    female = ifelse(Sex == "Female", 1, 0),
    age_14 = ifelse(Age == 14, 1, 0),
    age_28 = ifelse(Age == 28, 1, 0),
    age_42 = ifelse(Age == 42, 1, 0),
    age_56 = ifelse(Age == 56, 1, 0),
    inf_female = infected*female, 
    inf_14 = infected*age_14,
    inf_28 = infected*age_28,
    inf_42 = infected*age_42,
    inf_56 = infected*age_56,
    female_14 = female*age_14,
    female_28 = female*age_28,
    female_42 = female*age_42,
    female_56 = female*age_56
  ) %>% 
  dplyr::select(c(Day, intercept, infected, female, age_14, age_28, age_42, age_56,
                  inf_female, inf_14, inf_28, inf_42, inf_56, female_14, female_28, female_42, female_56))

#### step 2: initial values and draws from MH algorithm ####

B <- 50000
a1 <- a2 <- 0

outcomes <- matrix(0, nrow = B + 1, ncol = ncol(is_dummy)-1+2)
colnames(outcomes) <- c("alpha", "theta", colnames(is_dummy)[-1])
# outcomes[1,] = c(1,5,1.1,rep(0, ncol(outcomes)-3))
# outcomes[1,] = c(3.2,0.4,2.6,rep(0, ncol(outcomes)-3))
# outcomes[1,] = c(3.68270763,0.33296154,2.62699638,-0.69282490,0.27977484,1.36118492,  
#                  1.17942181,0.82778030,0.52509799,-0.08029493,0.44532987,0.17147463,
#                  -0.49091124,0.33804166,-0.15852870,-0.33884779,-0.23518983,-0.27212150)

outcomes[1,] = c(rep(0.5,2),rep(0, ncol(outcomes)-2))

#### step 3: for loop of MH algorithm

for(iter in 2:(B+1)){
  
  #### ---------------------- MH for shape ---------------------- ####
  # propose new shape
  log_shape_proposed = log_shape_sampler2(log_alpha = log(outcomes[iter-1,1]),
                                          log_theta = log(outcomes[iter-1,2]),
                                          multiplier = 1)
  loglike = log_likelihood(data = is_dummy,
                           alpha = exp(log_shape_proposed[1]),
                           theta = exp(log_shape_proposed[2]),
                           coefs = outcomes[iter-1,3:ncol(outcomes)])
  if(loglike == 0){
    q1 = 0
  }else{
    # ratio of posterior 
    q1_post = 
      (
        loglike + 
          log(prior(alpha = exp(log_shape_proposed[1]),
                    theta = exp(log_shape_proposed[2]),
                    coefs = outcomes[iter-1,3:ncol(outcomes)])) -
          log_likelihood(data = is_dummy,
                         alpha = outcomes[iter-1,1],
                         theta = outcomes[iter-1,2],
                         coefs = outcomes[iter-1,3:ncol(outcomes)]) -
          log(prior(alpha = exp(log_shape_proposed[1]),
                    theta = exp(log_shape_proposed[2]),
                    coefs = outcomes[iter-1,3:ncol(outcomes)]))
      ) %>% 
      exp()
    # ratio of jacobian
    q1_jacobian = (exp(log_shape_proposed[1])*exp(log_shape_proposed[2])) / (outcomes[iter-1,1]*outcomes[iter-1,2])
    q1 = (q1_post * q1_jacobian) %>% as.numeric()
  }
  # update shape
  if (runif(1) < q1) {
    outcomes[iter,1] = exp(log_shape_proposed[1])
    outcomes[iter,2] = exp(log_shape_proposed[2])
    a1 = a1 + 1
  }else{
    outcomes[iter,1] = outcomes[iter-1,1]
    outcomes[iter,2] = outcomes[iter-1,2]
  }
  
  #### ---------------------- MH for scale/coefs ---------------------- ####
  # propose new coefs
  coefs_proposed = coefs_sampler2(coefs = outcomes[iter-1,3:ncol(outcomes)], 
                                  multiplier = 1)
  loglike = log_likelihood(data = is_dummy,
                           alpha = outcomes[iter,1],
                           theta = outcomes[iter,2],
                           coefs = coefs_proposed)
  if(loglike == 0 || is.nan(loglike)){
    q2 = 0
  }else{
    # ratio of posterior 
    q2 = 
      (
        loglike + 
          log(prior(alpha = outcomes[iter,1],
                    theta = outcomes[iter,2],
                    coefs = coefs_proposed)) -
          log_likelihood(data = is_dummy,
                         alpha = outcomes[iter,1],
                         theta = outcomes[iter,2],
                         coefs = outcomes[iter-1,3:ncol(outcomes)]) -
          log(prior(alpha = outcomes[iter,1],
                    theta = outcomes[iter,2],
                    coefs = outcomes[iter-1,3:ncol(outcomes)]))
      ) %>% 
      exp()
  }
  # update scale
  if (runif(1) < q2) {
    outcomes[iter,3:ncol(outcomes)] = coefs_proposed
    a2 = a2 + 1
  }else{
    outcomes[iter,3:ncol(outcomes)] = outcomes[iter-1,3:ncol(outcomes)]
  }
  
  if (iter %% 10000 == 0) {cat(paste0("iteration: ", iter, "\n"))}
}

#### step 3: check acceptance rate and convergence ####

cat(paste0("acceptance rate for shape:  ", a1/B, "\n"))
cat(paste0("acceptance rate for coefs:  ", a2/B, "\n"))

# outcomes %>% burn_n_thin_draws(burn_at = 1, jump=1) %>% apply(2, function(c) acf(c))
outcomes %>%
  burn_n_thin_draws(burn_at = 1, jump = 1) %>%
  as.data.frame() %>%
  rownames_to_column("iteration") %>%
  reshape2::melt(id.vars = "iteration") %>%
  ggplot(aes(x = as.numeric(iteration), y = value, color = variable)) +
  geom_line(size = 0.1) +
  facet_wrap(~variable, ncol = 3, scales = "free_y") +
  labs(y = "coefficients") +
  theme_minimal() +
  theme(legend.position = "none")


#### step 4: save posterior dist and acceptance rate ####

outcomes = rbind(c(rep(a1/B, 2), rep(a2/B, ncol(outcomes)-2)), outcomes)
# write.csv(outcomes, file = "./data/posterior_co_random.csv", row.names = FALSE)

