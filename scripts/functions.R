
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
  # m = matrix(c(8.864063e-04,-0.0094593231,-0.0094593231,0.1105368735), ncol=2)
  # new proposed 2020-04-12
  m = matrix(c(-2.576053, -0.005930974, -0.005930974, -7.193297), ncol = 2)
  MASS::mvrnorm(1, c(log_alpha, log_theta), multiplier*m)
}

#### random sampler for coefficients CO ------------------------------- ####

coefs_sampler2 <- function(coefs, multiplier){
  m = matrix(
    c(
       0.005672872, -0.002669458, -0.002691113, -0.002878614, -0.003271481, -0.003010084, -0.002511487,  0.002381041,  0.002543716,  0.003070963,  0.002465176,  0.002210576,  0.002527249,  0.002956663,  0.002911526,  0.002481151,  -0.001537563,  -0.003033601,  -0.002475190,  -0.002038347,
      -0.002669458,  0.005316181,  0.002281056,  0.001947195,  0.002424420,  0.002512071,  0.001929093, -0.005105257, -0.004962621, -0.005487076, -0.005496556, -0.004979651, -0.001952308, -0.002384429, -0.002641449, -0.001850093,   0.004678553,   0.005496978,   0.005421307,   0.004640388,
      -0.002691113,  0.002281056,  0.003131074,  0.002429113,  0.002645602,  0.002733310,  0.002338012, -0.002701463, -0.002202043, -0.002474463, -0.002412784, -0.002204857, -0.002977257, -0.003316759, -0.003419056, -0.002806065,   0.002453809,   0.003127546,   0.002817517,   0.002458429,
      -0.002878614,  0.001947195,  0.002429113,  0.003006421,  0.002567143,  0.002670483,  0.002368381, -0.001847035, -0.002457747, -0.002158871, -0.002052244, -0.001971318, -0.002914302, -0.002543383, -0.002753763, -0.002128606,   0.002158013,   0.002311202,   0.002065475,   0.001649250,
      -0.003271481,  0.002424420,  0.002645602,  0.002567143,  0.003425107,  0.002885716,  0.002487734, -0.002276897, -0.002319177, -0.003281496, -0.002471562, -0.002336591, -0.002475771, -0.003372978, -0.002953902, -0.002329716,   0.001988507,   0.003300746,   0.002388285,   0.002033732,
      -0.003010084,  0.002512071,  0.002733310,  0.002670483,  0.002885716,  0.003497085,  0.002648466, -0.002414070, -0.002403610, -0.002780940, -0.003127605, -0.002523005, -0.002575681, -0.002865820, -0.003590417, -0.002452247,   0.002126090,   0.002883627,   0.003134849,   0.002175550,
      -0.002511487,  0.001929093,  0.002338012,  0.002368381,  0.002487734,  0.002648466,  0.003358532, -0.001766779, -0.001840760, -0.002207491, -0.002160660, -0.002945170, -0.002208050, -0.002508299, -0.002642995, -0.003110995,   0.001573319,   0.002375846,   0.002027300,   0.002484386,
       0.002381041, -0.005105257, -0.002701463, -0.001847035, -0.002276897, -0.002414070, -0.001766779,  0.006001062,  0.004858528,  0.005240927,  0.005349085,  0.004750296,  0.002424078,  0.002898175,  0.003114761,  0.002224544,  -0.005669470,  -0.006397770,  -0.006295366,  -0.005466571,
       0.002543716, -0.004962621, -0.002202043, -0.002457747, -0.002319177, -0.002403610, -0.001840760,  0.004858528,  0.005664104,  0.005189908,  0.005230166,  0.004739592,  0.002566534,  0.002387046,  0.002677152,  0.001731221,  -0.005509289,  -0.005365271,  -0.005343287,  -0.004466786,
       0.003070963, -0.005487076, -0.002474463, -0.002158871, -0.003281496, -0.002780940, -0.002207491,  0.005240927,  0.005189908,  0.006716127,  0.005671574,  0.005218324,  0.002084786,  0.003137643,  0.002905619,  0.002033507,  -0.004846037,  -0.006643600,  -0.005608494,  -0.004771376,
       0.002465176, -0.005496556, -0.002412784, -0.002052244, -0.002471562, -0.003127605, -0.002160660,  0.005349085,  0.005230166,  0.005671574,  0.006696344,  0.005369149,  0.002147722,  0.002489050,  0.003331143,  0.002129007,  -0.005019242,  -0.005704943,  -0.006713798,  -0.005033419,
       0.002210576, -0.004979651, -0.002204857, -0.001971318, -0.002336591, -0.002523005, -0.002945170,  0.004750296,  0.004739592,  0.005218324,  0.005369149,  0.006492015,  0.001881183,  0.002270886,  0.002611430,  0.002791607,  -0.004480240,  -0.005272799,  -0.005245970,  -0.006049108,
       0.002527249, -0.001952308, -0.002977257, -0.002914302, -0.002475771, -0.002575681, -0.002208050,  0.002424078,  0.002566534,  0.002084786,  0.002147722,  0.001881183,  0.003883010,  0.003211362,  0.003362058,  0.002665316,  -0.003246177,  -0.002895465,  -0.002717173,  -0.002192086,
       0.002956663, -0.002384429, -0.003316759, -0.002543383, -0.003372978, -0.002865820, -0.002508299,  0.002898175,  0.002387046,  0.003137643,  0.002489050,  0.002270886,  0.003211362,  0.004536735,  0.003651844,  0.003073851,  -0.002806558,  -0.004285694,  -0.002983509,  -0.002672160,
       0.002911526, -0.002641449, -0.003419056, -0.002753763, -0.002953902, -0.003590417, -0.002642995,  0.003114761,  0.002677152,  0.002905619,  0.003331143,  0.002611430,  0.003362058,  0.003651844,  0.004724280,  0.003085641,  -0.003010241,  -0.003668241,  -0.004231749,  -0.002848987,
       0.002481151, -0.001850093, -0.002806065, -0.002128606, -0.002329716, -0.002452247, -0.003110995,  0.002224544,  0.001731221,  0.002033507,  0.002129007,  0.002791607,  0.002665316,  0.003073851,  0.003085641,  0.004040772,  -0.001982212,  -0.002831457,  -0.002470968,  -0.003451005,
      -0.001537563,  0.004678553,  0.002453809,  0.002158013,  0.001988507,  0.002126090,  0.001573319, -0.005669470, -0.005509289, -0.004846037, -0.005019242, -0.004480240, -0.003246177, -0.002806558, -0.003010241, -0.001982212,   0.007221781,   0.006221453,   0.006174744,   0.005215791,
      -0.003033601,  0.005496978,  0.003127546,  0.002311202,  0.003300746,  0.002883627,  0.002375846, -0.006397770, -0.005365271, -0.006643600, -0.005704943, -0.005272799, -0.002895465, -0.004285694, -0.003668241, -0.002831457,   0.006221453,   0.008604713,   0.006754588,   0.005992424,
      -0.002475190,  0.005421307,  0.002817517,  0.002065475,  0.002388285,  0.003134849,  0.002027300, -0.006295366, -0.005343287, -0.005608494, -0.006713798, -0.005245970, -0.002717173, -0.002983509, -0.004231749, -0.002470968,   0.006174744,   0.006754588,   0.008440304,   0.005903923,
      -0.002038347,  0.004640388,  0.002458429,  0.001649250,  0.002033732,  0.002175550,  0.002484386, -0.005466571, -0.004466786, -0.004771376, -0.005033419, -0.006049108, -0.002192086, -0.002672160, -0.002848987, -0.003451005,   0.005215791,   0.005992424,   0.005903923,   0.007631515
    ),
    ncol = length(coefs)
  )
  MASS::mvrnorm(1, coefs, multiplier*m)
}

#### burn in and thin posterior draws ------------------------------ ####

burn_n_thin_draws <- function(data, jump, burn_at = 10000){
  data_burned <- data[-c(1:burn_at), ]
  rows_to_take <- (c(1:nrow(data_burned)) %% jump) == 0
  data_burned[rows_to_take, ]
}

#### cdf of exponentiated weibull dist  ---------------------------- ####

cdf_calculator <- function(t, alpha, theta, sigma) {
  (1 - exp(-((t/sigma)^alpha)))^theta
}

#### median life function ------------------------------------------ #### 

median_life <- function(u, alpha, theta, sigma){
  sigma * (-log(1-u^(1/theta)))^(1/alpha)
}

#### median residual life function  ------------------------------- #### 

median_residual_life <- function(t, u, alpha, theta, sigma){
  K = (1-u)*(1-cdf_calculator(t, alpha, theta, sigma))
  term_log = log(1-((1-K)^(1/theta)))
  (-sigma^alpha*term_log)^(1/alpha) - t
}


#### dummy categorical covariates --------------------------------- #### 

dummy_aco <- function(data){
  data %>% 
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
}

dummy_co <- function(data){
  data %>% 
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
      female_56 = female*age_56, 
      female_inf_14 = female*infected*age_14,
      female_inf_28 = female*infected*age_28,
      female_inf_42 = female*infected*age_42,
      female_inf_56 = female*infected*age_56
    ) %>% 
    dplyr::select(c(Day, intercept, infected, female, age_14, age_28, age_42, age_56,
                    inf_female, inf_14, inf_28, inf_42, inf_56, female_14, female_28, female_42, female_56,
                    female_inf_14, female_inf_28, female_inf_42, female_inf_56))
}

#### calculate sigma by product of covariates and coefficients --------------- ####

sigma_calculator <- function(row, x_vector){
  exp(sum(row[3:length(row)] * x_vector))
}

#### calculate the survival function ---------------------------------------- ####

survival_function_calculator <- function(data, x_fit, max_day){
  data$sigma <- apply(data, 1, function(r) sigma_calculator(r, x_fit)) 
  seq(2,max_day) %>% 
    sapply(function(day)(
      apply(data, 1, function(row) (
        1 - cdf_calculator(t = day, alpha = row[1], theta = row[2], sigma = row[length(row)])
      ))
    )) %>% 
    as.data.frame()
}

#### calculate the survival intervals ---------------------------------------- ####

survival_interval_calculator <- function(data, x_fit, max_day){
  dummy_data = survival_function_calculator(data, x_fit, max_day)
  apply(dummy_data, 2, function(col) quantile(col, probs = c(0.025, 0.5, 0.975)) ) %>% 
    as.data.frame()
}

#### calculate median residual life -------------------------------------- #### 

median_residual_calculator <- function(data, x_fit, max_day){
  data$sigma= apply(data, 1, function(row) sigma_calculator(row, x_fit))
  seq(2, max_day) %>% 
    sapply(function(day)(
      apply(data, 1, function(row)(
        median_residual_life(t = day, u = 0.5, alpha = row[1], theta = row[2], sigma = row[length(row)])
      ))
    )) %>% 
    as.data.frame()
}

#### calculate the median residual intervals ------------------------- ####

median_residual_interval_calculator <- function(data, x_fit, max_day){
  dummy_data = median_residual_calculator(data, x_fit, max_day)
  apply(dummy_data, 2, function(col) quantile(col, probs = c(0.025, 0.5, 0.975)) ) %>% 
    as.data.frame()
}

#### calculate median life ----------------------------------------- ####

median_calculator <- function(data, x_fit){
  data$sigma <- apply(data, 1, function(row) sigma_calculator(row, x_fit))
  apply(data, 1, function(row)  
    median_residual_life(u = 0.5, alpha = row[1], theta = row[2], sigma = row[length(row)]))
}

#### calculate median intervals ------------------------------------- ####

median_intervals_calculator <- function(data){
  apply(data, 2, function(col) quantile(col, probs = c(0.025, 0.5, 0.975)) ) %>% 
    as.data.frame()
}

#### convert dummied age columns back to one --------------------------------- ####

convert_age_back <- function(age_14, age_28){
  if(age_14 == 1) {return(14)}
  if(age_28 == 1) {return(28)}
  if((age_14==0) && (age_28==0)) {return(42)}
}

convert_age_back2 <- function(age_14, age_28, age_42, age_56){
  if(age_14 == 1) {return(14)}
  if(age_28 == 1) {return(28)}
  if(age_42 == 1) {return(42)}
  if(age_56 == 1) {return(56)}
  if((age_14==0) && (age_28==0) && (age_42==0) && (age_56==0)) {return(72)}
}

#### pack survival function and add initial day survival ----------------- ####

add_survival_initial <- function(data, subgroups, max_day){
  data %>% 
    rbind(Day = rep(c(2:max_day), nrow(subgroups))) %>% 
      t() %>% 
      set_colnames(c("lower","est","upper","Day")) %>% 
      cbind(
        subgroups[rep(1:nrow(subgroups), each = length(c(2:max_day))),]
      ) %>% 
      rbind(
        cbind(
          data.frame(
            lower = rep(1, nrow(subgroups)),
            est = rep(1, nrow(subgroups)),
            upper = rep(1, nrow(subgroups)),
            Day = rep(1, nrow(subgroups))
          ),
          subgroups
        ),
        .
      )
}

#### pack median residual function and add initial median residual life ----------------- ####

add_median_residual_initial <- function(data, subgroups, max_day){
  data %>% 
    rbind(Day = rep(c(2:max_day), nrow(subgroups))) %>% 
    t() %>% 
    set_colnames(c("lower","est","upper","Day")) %>% 
    cbind(
      subgroups[rep(1:nrow(subgroups), each = length(c(2:max_day))),]
    ) %>% 
    filter( (!is.infinite(lower)) & (!is.infinite(upper)) & (!is.infinite(est)))
}


#### plot survival function ---------------------------------------------- ####

plot_survival <- function(data, color){
  data %>% 
    filter(est > 0.007) %>% 
    ggplot(aes(x = Day + Age, 
               linetype = Treatment,
               group = interaction(Treatment,Sex,Age))) +
    geom_line(aes(y = est), color = color, size = 0.5) +
    geom_line(aes(y = lower), color = color, size = 0.5) +
    geom_line(aes(y = upper), color = color, size = 0.5) +
    scale_linetype_manual(values = c("Infected" = "solid", "Uninfected" = "dashed")) +
    facet_grid(Sex~.) +
    theme_minimal() +
    labs(x = "Days after spray", y = "Survival function") +
    xlim(10,100)
}


#### plot median residual function ---------------------------------------------- ####

plot_median_residual <- function(data, type){
  if(type == 1){
    ggplot(data,
           aes(x = Day+Age, y = est, ymin = lower, ymax = upper, color = as.factor(Age),
               group = interaction(Treatment,Sex,Age))) +
      geom_line() +
      geom_errorbar() +
      facet_wrap(Treatment~Sex) +
      scale_color_jco() +
      theme_minimal() +
      labs(x= "Age", y = "Median residual life", color = "Age")
  }else{
    ggplot(data,
           aes(x = Day, y = est, ymin = lower, ymax = upper, color = as.factor(Age),
               group = interaction(Treatment,Sex,Age))) +
      geom_line() +
      geom_errorbar() +
      facet_wrap(Treatment~Sex) +
      scale_color_jco() +
      theme_minimal() +
      labs(y = "Median residual life", color = "Age")
  }
}


#### plot median residual slope  ---------------------------------------------- ####

plot_median_residual_slope <- function(data, type){
  if(type == 1){
    ggplot(data, aes(x = Day+Age, y = abs(est), color = as.factor(Age),
               group = interaction(Treatment,Sex,Age))) +
      geom_line() +
      scale_color_jco() +
      facet_wrap(Treatment~Sex) +
      labs(x= "Age", y = "Slope", color = "Age") + 
      theme_minimal()
  }else{
    ggplot(data, aes(x = Day, y = abs(est), color = as.factor(Age),
               group = interaction(Treatment,Sex,Age))) +
      geom_line() +
      scale_color_jco() +
      facet_wrap(Treatment~Sex) +
      labs(y = "Slope", color = "Age") + 
      theme_minimal()
  }
}

