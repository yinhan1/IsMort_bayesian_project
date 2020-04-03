
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

#### raw data plot: survival percent  ---------------------------- ####

## calculate survival percent 
get_survival_percent <- function(data, cut_off = 0.005){
  data_percent <- 
    data %>% 
    group_by(Population, Treatment, Sex, Age, Day) %>% 
    summarise(Deaths = sum(Deaths)) %>% 
    group_by(Population, Treatment, Sex, Age) %>% 
    mutate(Total = sum(Deaths),
           Deaths = cumsum(Deaths)) %>% 
    ungroup() %>% 
    mutate(Percent = 1 - Deaths/Total) %>% 
    select(-c(Total,Deaths))

    bind_rows(
      data_percent, 
      data_percent %>% 
        select(-c(Percent,Day)) %>% 
        unique() %>% 
        mutate(Day = 0, Percent = 1)
    ) %>% 
    filter(Percent > cut_off)
} 

## plot survival percent 
plot_survival_percent <- function(data){
  data %>% 
  mutate(tag = paste(Population, Sex)) %>% 
    ggplot(aes(x = Day + as.numeric(as.character(Age)),
               y = Percent*100, 
               group = interaction(Population,Treatment,Sex,Age),
               color = Population, 
               linetype = Treatment)) +
    geom_line() +
    facet_grid(Sex~Population) +
    scale_color_manual(values = c("ACO" = "brown3", "CO" = "grey30")) +
    scale_linetype_manual(values = c("Infected" = "solid", "Uninfected" = "dashed")) +
    theme_minimal() +
    ylim(0, 107) +
    labs(x = "Days after spray", y = "Survival percent (%)")
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

#### random sampler for log_alpha and log_theta ACO -------------------- ####

log_shape_sampler <- function(log_alpha, log_theta, multiplier){
  m = matrix(c(0.0006053001,-0.001027386,-0.001027386,0.002148208), ncol=2)
  MASS::mvrnorm(1, c(log_alpha, log_theta), multiplier*m)
}

#### random sampler for coefficients ACO ------------------------------- ####

coefs_sampler <- function(coefs, multiplier){
  m = diag(x = length(coefs)) * c(0.0001)
  MASS::mvrnorm(1, coefs, multiplier*m)
}

#### random sampler for log_alpha and log_theta CO -------------------- ####

log_shape_sampler2 <- function(log_alpha, log_theta, multiplier){
  m = matrix(c(8.864063e-04,-0.0094593231,-0.0094593231,0.1105368735), ncol=2)
  MASS::mvrnorm(1, c(log_alpha, log_theta), multiplier*m)
}

#### random sampler for coefficients CO ------------------------------- ####

coefs_sampler2 <- function(coefs, multiplier){
  m = matrix(
    c(
      0.0040710324, -8.751589e-04, -1.141065e-03, -1.718342e-03, -1.479397e-03, -1.382290e-03, -1.243257e-03,  2.083276e-04,  1.100848e-03,  5.990355e-04,  0.0005853875,  0.0003508888,  1.257629e-03,  0.0009751685,  0.0009006624,  1.223964e-03,
      -0.0008751589,  1.639023e-03, -2.274207e-05,  5.508726e-04,  4.678327e-04,  4.892571e-04,  5.571443e-04, -3.925505e-04, -1.492931e-03, -1.338323e-03, -0.0013815837, -0.0014580305,  2.030281e-04,  0.0002801420,  0.0002689931,  2.035145e-04,
      -0.0011410650, -2.274207e-05,  1.171505e-03,  8.403559e-04,  8.058213e-04,  7.957070e-04,  7.565865e-04, -1.191741e-04,  4.465597e-05,  6.596132e-05,  0.0001147259,  0.0001906197, -1.112871e-03, -0.0010523018, -0.0010806770, -1.115909e-03,
      -0.0017183419,  5.508726e-04,  8.403559e-04,  1.599665e-03,  1.148139e-03,  1.147607e-03,  1.140571e-03, -1.474958e-05, -8.894544e-04, -4.740215e-04, -0.0004924252, -0.0004250131, -1.087028e-03, -0.0007958929, -0.0007999547, -9.000300e-04,
      -0.0014793966,  4.678327e-04,  8.058213e-04,  1.148139e-03,  1.547099e-03,  1.053528e-03,  1.091250e-03, -6.957815e-06, -5.566094e-04, -7.095671e-04, -0.0004042833, -0.0003622183, -8.079231e-04, -0.0010771283, -0.0007542152, -8.808369e-04,
      -0.0013822900,  4.892571e-04,  7.957070e-04,  1.147607e-03,  1.053528e-03,  1.442052e-03,  1.084849e-03, -6.260946e-06, -5.841735e-04, -4.387456e-04, -0.0006969650, -0.0004154847, -8.255802e-04, -0.0007305863, -0.0010176703, -8.434370e-04,
      -0.0012432570,  5.571443e-04,  7.565865e-04,  1.140571e-03,  1.091250e-03,  1.084849e-03,  1.803576e-03, -1.791519e-06, -6.181014e-04, -5.105923e-04, -0.0005376494, -0.0008920102, -7.527625e-04, -0.0007187303, -0.0007055809, -1.248724e-03,
       0.0002083276, -3.925505e-04, -1.191741e-04, -1.474958e-05, -6.957815e-06, -6.260946e-06, -1.791519e-06,  4.000704e-04,  1.750714e-04,  1.758909e-04,  0.0001502162,  0.0001619473, -1.084067e-04, -0.0001203937, -0.0001158421, -1.151855e-04,
      0.0011008475, -1.492931e-03,  4.465597e-05, -8.894544e-04, -5.566094e-04,-5.841735e-04, -6.181014e-04,  1.750714e-04,  1.914052e-03,  1.293478e-03,  0.0013561539,  0.0013789665, -9.161988e-05, -0.0001668238, -0.0001537390, -7.691657e-05,
      0.0005990355, -1.338323e-03,  6.596132e-05, -4.740215e-04, -7.095671e-04, -4.387456e-04, -5.105923e-04,  1.758909e-04,  1.293478e-03,  1.658192e-03,  0.0012183922,  0.0012756995, -1.534455e-04, -0.0001486325, -0.0001892409, -1.143662e-04,
      0.0005853875, -1.381584e-03,  1.147259e-04, -4.924252e-04, -4.042833e-04,-6.969650e-04, -5.376494e-04,  1.502162e-04,  1.356154e-03,  1.218392e-03,  0.0017412798,  0.0013672090, -1.613126e-04, -0.0002574529, -0.0002064878, -1.678280e-04,
      0.0003508888, -1.458031e-03,  1.906197e-04, -4.250131e-04, -3.622183e-04, -4.154847e-04, -8.920102e-04,  1.619473e-04,  1.378967e-03,  1.275699e-03,  0.0013672090,  0.0021065395, -2.883503e-04, -0.0003347314, -0.0003086805, -2.506416e-04,
      0.0012576289,  2.030281e-04, -1.112871e-03, -1.087028e-03, -8.079231e-04, -8.255802e-04, -7.527625e-04, -1.084067e-04, -9.161988e-05, -1.534455e-04, -0.0001613126, -0.0002883503,  1.604704e-03,  0.0011258079,  0.0011658313,  1.220864e-03,
      0.0009751685,  2.801420e-04, -1.052302e-03, -7.958929e-04, -1.077128e-03, -7.305863e-04, -7.187303e-04, -1.203937e-04, -1.668238e-04, -1.486325e-04, -0.0002574529, -0.0003347314,  1.125808e-03,  0.0016158744,  0.0011003264,  1.162148e-03,
      0.0009006624,  2.689931e-04, -1.080677e-03, -7.999547e-04, -7.542152e-04, -1.017670e-03, -7.055809e-04, -1.158421e-04, -1.537390e-04, -1.892409e-04, -0.0002064878, -0.0003086805,  1.165831e-03,  0.0011003264,  0.0015963546,  1.157664e-03,
      0.0012239640,  2.035145e-04, -1.115909e-03, -9.000300e-04, -8.808369e-04, -8.434370e-04, -1.248724e-03, -1.151855e-04, -7.691657e-05, -1.143662e-04, -0.0001678280, -0.0002506416,  1.220864e-03,  0.0011621484,  0.0011576638,  1.859971e-03
    ),
    ncol = length(coefs)
  )
  MASS::mvrnorm(1, coefs, multiplier*m)
}


#### burn in and thin posterior draws ----------------------------------------- ####

burn_n_thin_draws <- function(data, jump, burn_at = 10000){
  data_burned <- data[-c(1:burn_at), ]
  rows_to_take <- (c(1:nrow(data_burned)) %% jump) == 0
  data_burned[rows_to_take, ]
}

#### cdf of exponentiated weibull dist  ---------------------------- ####

cdf_calculator <- function(t, alpha, theta, sigma) {
  (1 - exp(-((t/sigma)^alpha)))^theta
}




































