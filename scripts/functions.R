
# First Edited: 03/19/2020

#### data cleaning --------------------------------------------- ####

## remove day 1 and convert data tpyes
clean_data <- function(data, cols_to_numeric){
  data %>% 
    mutate_all(function(x) as.factor(x)) %>% 
    mutate_at(cols_to_numeric, function(x) as.numeric(as.character(x))) %>% 
    mutate(Sex = recode(Sex, "F" = "Female", "M" = "Male")) %>% 
    filter(`Day after spray` != 1 & Treatment != "Control") %>% 
    droplevels() %>% 
    setnames(old = c("Age at spray", "Day after spray"),
             new = c("Age","Day"))
}

## convert death counts to death status
convert_to_death_status <- function(data){
  data[rep(1:nrow(data), data$Deaths),] %>% 
    mutate(Death = 1) %>% 
    dplyr::select(-Deaths)
}

#### MH ratio of posterior  ------------------------------------- ####

## define prior 
prior <- function(alpha, theta, coefs){
  term_1 = dexp(alpha, rate = 1/3)
  term_2 = dexp(theta, rate = 1/2)
  term_3 = sapply(coefs[2:length(coefs)], function(x) dnorm(x, 0, 100)) %>% prod()
  (term_1 * term_2 * dnorm(coefs[1], 1, 100) * term_3) %>% as.numeric()
}

## define density for each row
density_calculator <- function(t, alpha, theta, sigma){
  term_1 = alpha * theta / sigma
  term_exp = exp(-(t/sigma)^alpha)
  term_3 = (t/sigma)^(alpha-1)
  term_1 * (1-term_exp)^(theta-1) * term_exp * term_3
}

## pack rows to log-likelihood
log_likelihood <- function(data, alpha, theta, coefs){
  sigma = exp(as.matrix(data[,-1]) %*% as.matrix(coefs, ncol=1))
  condition_1 = theta < 1
  condition_2 = sum(exp(-(data[,1]/sigma)^alpha) == 1) %>% as.logical()
  if(condition_1 && condition_2){
    return(0)
  }else{
    density_calculator(data[,1], alpha, theta, sigma) %>% log() %>% sum() %>% return()
  }
}

#### random sampler for log_alpha and log_theta -------------------- ####

log_shape_sampler <- function(log_alpha, log_theta, multiplier){
  m = matrix(c(0.0006053001,-0.001027386,-0.001027386,0.002148208), ncol=2)
  MASS::mvrnorm(1, c(log_alpha, log_theta), multiplier*m)
}

#### random sampler for coefficients ------------------------------- ####

coefs_sampler <- function(coefs, multiplier){
  m = diag(x = length(coefs)) * c(0.0001)
  MASS::mvrnorm(1, coefs, multiplier*m)
}

#### thin posterior draws ----------------------------------------- ####

thin_draws <- function(data, B, b){
  data = outcomes
  b = B/50
  dummy_data = matrix(0, nrow = b+1, ncol = ncol(data))
  dummy_data[1,] = data[1,]
  dummy_data[2:(b+1),] = data[(data$iteration %% (B/b)) == 0, ]
  return(as.data.frame(dummy_data))
}

#### cdf of exponentiated weibull dist  ---------------------------- ####

cdf_calculator <- function(t, alpha, theta, sigma) {
  (1 - exp(-((t/sigma)^alpha)))^theta
}


#### get survival function bounds  ---------------------------- ####

get_bounds_survival <- function(data, time_limit){
  time_sequence = seq(2, time_limit, 1)
  sapply(time_sequence, 
         function(time) {
           apply(data, 1, function(row) cdf_calculator)
         } )
}









