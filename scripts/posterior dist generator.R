
library(tidyverse)
library(data.table)
library(magrittr)
library(ggsci)
library(kableExtra)

source("./scripts/functions.R")

is_data <- readxl::read_excel("data/IS_data.xlsx") %>% clean_data(cols_to_numeric = c("Day after spray","Deaths"))


#### ACO (parameters good 550k needed) ============================= ####

#### step 1: filter out ACO and dummy categorical features ####

is_dummy <- is_data %>% 
  filter(Population == "ACO" & Treatment != "Control") %>% 
  convert_to_death_status() %>% dummy_aco(.) 

#### step 2: initial values and draws from MH algorithm ####

B <- 10000
a1 <- a2 <- 0

outcomes <- matrix(0, nrow = B + 1, ncol = ncol(is_dummy)-1+2)
colnames(outcomes) <- c("alpha", "theta", colnames(is_dummy)[-1])
outcomes[1,] = c(1.890390276,  1.846996258,  1.823753465, -0.615436519,
                 0.030102303,  0.936181470,  0.329697524, -0.003964608,
                 -0.613027522, -0.194533670, 0.108735851, -0.064016412)

#### step 3: for loop of MH algorithm
start_time <- Sys.time()
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
                                 multiplier = 0.05)
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

as.data.frame(outcomes)[-1,] %>% 
  mutate(iteration = c(1:(nrow(outcomes)-1))) %>% 
  reshape2::melt(id.vars = "iteration") %>%
  ggplot(aes(x = iteration, y = value, color = variable)) +
  geom_line(size = 0.1) +
  # geom_hline(yintercept = 0, color = "brown3", size = 0.5, linetype = "dashed") +
  facet_wrap(~variable, ncol = 3, scales = "free_y") +
  labs(y = "coefficients") +
  theme_minimal() +
  theme(legend.position = "none")

coda::effectiveSize(outcomes)

#### step 4: save posterior dist and acceptance rate ####

outcomes = rbind(c(rep(a1/B, 2), rep(a2/B, ncol(outcomes)-2)), outcomes) %>% as_tibble()
write_csv(outcomes, "./data/posterior_aco_2020-04-14.csv")
end_time <- Sys.time()
end_time - start_time






#### CO (parameters good 840-850k needed) ============================================ ####

#### step 1: filter out CO and dummy categorical features ####

is_dummy <- is_data %>% 
  filter(Population == "CO" & Treatment != "Control") %>% 
  convert_to_death_status() %>% 
  dummy_co(.)

#### step 2: initial values and draws from MH algorithm ####

B <- 50000
a1 <- a2 <- 0

outcomes <- matrix(0, nrow = B + 1, ncol = ncol(is_dummy)-1+2)
colnames(outcomes) <- c("alpha", "theta", colnames(is_dummy)[-1])


# extract from posterior_co.csv
# outcomes[1,] <- c(
#   1.097296367,   4.953922497,   1.209283349,  -0.482228059,   0.025895856,   1.832793248,   1.453513721,
#   0.735576928,   0.250055192,  -0.015691754,  -0.884285239,  -0.906604992,   -0.521759576,  -0.309321503,
#   0.158424317,   0.096360510,   0.325413932,   0.152977188,   0.064844156,   0.006385447,  -0.182813001,
#   0.042295215 )

outcomes[1,] <-
  c(
    1.114988, 4.732720,  1.224525, -0.4636519,  0.06126246, 1.839870, 1.456864, 0.7476735, 0.2777832, -0.03781173, -0.9081203,
    -0.9012455, -0.5302593, -0.3442507, 0.13778395,  0.06706877, 0.2957999, 0.10999838,  0.09359469, 0.04177699, -0.18162112, 0.08890605
  )


# 2020_04_12 changed
# outcomes[1,] <- c(
#   3.95744600,    0.31253430,    2.67353700,   -0.88902345,    0.18291240,    
#   1.33515900,    1.04234250,    0.81375680,    0.45575790,    0.11387325 ,   
#   0.59973100,    0.61522040,   -0.31968705,   0.31373570,   -0.10320105,   
#   -0.07326436,   -0.15598775,   -0.09338692,   -0.13881100,   -0.79989680,   -0.10959390,    0.02975302
# )
# inital values started from last draw of recent draws
# outcomes[1, ] <- read_csv('data/posterior_co_trials_recent_run.csv') %>% 
#   last() %>% 
#   unlist()

#### step 3: for loop of MH algorithm
start_time <- Sys.time()
for(iter in 2:(B+1)){
  
  #### ---------------------- MH for shape ---------------------- ####
  # propose new shape
  log_shape_proposed = log_shape_sampler2(log_alpha = log(outcomes[iter-1,1]),
                                          log_theta = log(outcomes[iter-1,2]),
                                          multiplier = 0.07)
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
                                  multiplier = 0.05)
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
outcomes[-1,] %>%
  # burn_n_thin_draws(burn_at = 1, jump = 260) %>%
  as.data.frame() %>%
  rownames_to_column("iteration") %>%
  reshape2::melt(id.vars = "iteration") %>%
  ggplot(aes(x = as.numeric(iteration), y = value, color = variable)) +
  geom_line(size = 0.1) +
  facet_wrap(~variable, ncol = 3, scales = "free_y") +
  labs(y = "coefficients") +
  theme_minimal() +
  theme(legend.position = "none")

end_time <- Sys.time()
end_time - start_time

coda::effectiveSize(outcomes)

#### step 4: save posterior dist and acceptance rate ####

outcomes = rbind(c(rep(a1/B, 2), rep(a2/B, ncol(outcomes)-2)), outcomes)
write.csv(outcomes, file = "./data/posterior_co_2020-04-14.csv", row.names = FALSE)








#### Control vs Uninfected (parameters good 250k needed) ================================ ####

#### step 1: dummy categorical features ####

is_dummy <- is_data %>% filter(Age == 14 & Treatment != "Infected") %>% convert_to_death_status() %>% 
  mutate(intercept = 1,
         control = ifelse(Treatment=="Control", 1, 0),
         female = ifelse(Sex=="Female", 1, 0),
         aco = ifelse(Population == "ACO", 1, 0),
         control_female = control*female,
         control_aco = control*aco,
         female_aco = female*aco,
         control_f_aco = control*female*aco) %>% 
  select(Day,intercept,control,female,aco,control_female,control_aco,female_aco,control_f_aco)
  

#### step 2: initial values and draws from MH algorithm ####

B <- 20000
a1 <- a2 <- 0

outcomes <- matrix(0, nrow = B + 1, ncol = ncol(is_dummy)-1+2)
colnames(outcomes) <- c("alpha", "theta", colnames(is_dummy)[-1])
outcomes[1,] = c(3.88887942,0.71959404,3.89842895,-0.03050519,0.15341679,-0.75850896,0.06685480,-0.03322455,-0.03874913,-0.07996730)


#### step 3: for loop of MH algorithm ####
for(iter in 2:(B+1)){
  
  #### ---------------------- MH for shape ---------------------- ####
  # propose new shape
  log_shape_proposed = log_shape_sampler3(log_alpha = log(outcomes[iter-1,1]),
                                         log_theta = log(outcomes[iter-1,2]),
                                         multiplier = 0.1)
  
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
  coefs_proposed = coefs_sampler3(coefs = outcomes[iter-1,3:ncol(outcomes)], 
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
  rowid_to_column('iteration') %>% 
  reshape2::melt(id.vars = "iteration") %>%
  ggplot(aes(x = iteration, y = value, color = variable)) +
  geom_line(size = 0.1) +
  facet_wrap(~variable, ncol = 3, scales = "free_y") +
  labs(y = "coefficients") +
  theme_minimal() +
  theme(legend.position = "none")

coda::effectiveSize(outcomes)

#### step 4: save posterior dist and acceptance rate ####

outcomes = rbind(c(rep(a1/B, 2), rep(a2/B, ncol(outcomes)-2)), outcomes)
write.csv(outcomes, file = "./data/posterior_control_2020-04-17.csv", row.names = FALSE)



