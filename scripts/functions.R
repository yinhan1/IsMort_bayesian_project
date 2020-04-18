
# First Edited: 03/19/2020

#### data cleaning --------------------------------------------- ####

## remove day 1 and convert data tpyes
clean_data <- function(data, cols_to_numeric){
  data %>% 
    mutate_all(function(x) as.factor(x)) %>% 
    mutate_at(cols_to_numeric, function(x) as.numeric(as.character(x))) %>% 
    mutate(Sex = recode(Sex, "F" = "Female", "M" = "Male")) %>% 
    filter(`Day after spray` != 1) %>% 
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
        mutate(Day = 1, Percent = 1)
    ) %>% 
    filter(Percent > cut_off)
} 

## plot survival percent -- control
plot_survival_percent <- function(data){
  data %>% 
    mutate(tag = paste(Treatment, Population, Sex)) %>% 
    ggplot(aes(x = Day, y = Percent*100, color = tag, linetype = tag)) +
    geom_line(size = 1.2) +
    scale_color_manual(
      values = c(
        "Uninfected ACO Female" = "brown3",
        "Infected ACO Female" = "brown3",
        "Control ACO Female" = "brown3",
        "Uninfected ACO Male" = "brown3",
        "Infected ACO Male" = "brown3",
        "Control ACO Male" = "brown3",
        "Uninfected CO Female" = "grey30",
        "Infected CO Female" = "grey30",
        "Control CO Female" = "grey30",
        "Uninfected CO Male" = "grey30",
        "Infected CO Male" = "grey30",
        "Control CO Male" = "grey30"
      )) +
    scale_linetype_manual(
      values = c(
        "Uninfected ACO Female" = "dashed",
        "Infected ACO Female" = "solid",
        "Control ACO Female" = "dotted",
        "Uninfected ACO Male" = "dashed",
        "Infected ACO Male" = "solid",
        "Control ACO Male" = "dotted",
        "Uninfected CO Female" = "dashed",
        "Infected CO Female" = "solid",
        "Control CO Female" = "dotted",
        "Uninfected CO Male" = "dashed",
        "Infected CO Male" = "solid",
        "Control CO Male" = "dotted"
      )) +
    theme_minimal() +
    labs(x = "Days after spray", y = "Survival percent (%)", color = "", linetype = "")
}

## plot survival percent -- aco
plot_survival_percent2 <- function(data){
  data %>% 
    mutate(tag = paste(Treatment, Population, Sex)) %>% 
    ggplot(aes(x = Day + as.numeric(as.character(Age)), y = Percent*100, color = tag, linetype = tag, 
               group = paste(Treatment, Population, Sex, Age))) +
    geom_line(size = 1.2) +
    scale_color_manual(
      values = c(
        "Uninfected ACO Female" = "brown3",
        "Infected ACO Female" = "brown3",
        "Control ACO Female" = "brown3",
        "Uninfected ACO Male" = "brown3",
        "Infected ACO Male" = "brown3",
        "Control ACO Male" = "brown3",
        "Uninfected CO Female" = "grey30",
        "Infected CO Female" = "grey30",
        "Control CO Female" = "grey30",
        "Uninfected CO Male" = "grey30",
        "Infected CO Male" = "grey30",
        "Control CO Male" = "grey30"
      )) +
    scale_linetype_manual(
      values = c(
        "Uninfected ACO Female" = "dashed",
        "Infected ACO Female" = "solid",
        "Control ACO Female" = "dotted",
        "Uninfected ACO Male" = "dashed",
        "Infected ACO Male" = "solid",
        "Control ACO Male" = "dotted",
        "Uninfected CO Female" = "dashed",
        "Infected CO Female" = "solid",
        "Control CO Female" = "dotted",
        "Uninfected CO Male" = "dashed",
        "Infected CO Male" = "solid",
        "Control CO Male" = "dotted"
      )) +
    theme_minimal() +
    labs(x = "Days after spray", y = "Survival percent (%)", color = "", linetype = "")
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
  m = matrix(c(0.01114262,-0.01573048,-0.01573048,0.02324322), ncol=2)
  MASS::mvrnorm(1, c(log_alpha, log_theta), multiplier*m)
}

#### random sampler for coefficients ACO ------------------------------- ####

coefs_sampler <- function(coefs, multiplier){
  m = matrix(
    c(
      0.0071851664, -0.004962065, -0.0058428651, -0.0056797092, -5.654202e-03,  8.258855e-04,  0.004700024,  0.004150738,  0.0056207377,  0.0056849459,
     -0.0049620650,  0.013112832,  0.0090023144,  0.0069026423,  7.220138e-03, -5.953280e-04, -0.011470043, -0.012786175, -0.0088386152, -0.0093792422,
     -0.0058428651,  0.009002314,  0.0133393838,  0.0080449943,  8.282609e-03, -2.210229e-04, -0.007667616, -0.009040652, -0.0132055722, -0.0136164703,
     -0.0056797092,  0.006902642,  0.0080449943,  0.0076077576,  6.900719e-03,  1.050835e-04, -0.006502139, -0.006730158, -0.0086534516, -0.0081556094,
     -0.0056542019,  0.007220138,  0.0082826094,  0.0069007190,  8.188692e-03,  4.711354e-05, -0.006368861, -0.007743543, -0.0080731500, -0.0094521029,
      0.0008258855, -0.000595328, -0.0002210229,  0.0001050835,  4.711354e-05,  4.749380e-03, -0.001977225, -0.001724045, -0.0008051706, -0.0007537978,
      0.0047000239, -0.011470043, -0.0076676160, -0.0065021392, -6.368861e-03, -1.977225e-03,  0.013363245,  0.012301371,  0.0079438542,  0.0084799558,
      0.0041507376, -0.012786175, -0.0090406525, -0.0067301582, -7.743543e-03, -1.724045e-03,  0.012301371,  0.016950937,  0.0092541111,  0.0099765103,
      0.0056207377, -0.008838615, -0.0132055722, -0.0086534516, -8.073150e-03, -8.051706e-04,  0.007943854,  0.009254111,  0.0146440366,  0.0135491457,
      0.0056849459, -0.009379242, -0.0136164703, -0.0081556094, -9.452103e-03, -7.537978e-04,  0.008479956,  0.009976510,  0.0135491457,  0.0160567796
    ),
    ncol = length(coefs)
  )
  MASS::mvrnorm(1, coefs, multiplier*m)
}

#### random sampler for log_alpha and log_theta CO -------------------- ####

log_shape_sampler2 <- function(log_alpha, log_theta, multiplier){
  # m = matrix(c(8.864063e-04,-0.0094593231,-0.0094593231,0.1105368735), ncol=2)
  # new proposed 2020-04-12
  m = matrix(c(0.0007697761, -0.001871584, -0.001871584, 0.004913074), ncol = 2)
  MASS::mvrnorm(1, c(log_alpha, log_theta), multiplier*m)
}

#### random sampler for coefficients CO ------------------------------- ####

coefs_sampler2 <- function(coefs, multiplier){
  m = matrix(
    c( # cov(co_burned[,-c(1,2)]) posterior_co.csv
      0.005652016, -0.002559801, -0.002776032, -0.002996357, -0.003354646, -0.003084640, -0.002684503,  0.002322912,  0.002469943,  0.002971691,  0.002286614,  0.002031493,  0.002609980,  0.002976414,  0.002931943,  0.002656492,  -0.001475481,  -0.002962345,  -0.002332574,  -0.002016518,
      -0.002559801,  0.005494845,  0.002280648,  0.002018612,  0.002292863,  0.002451683,  0.001962589, -0.005337779, -0.005208700, -0.005496785, -0.005680588, -0.005191097, -0.001978139, -0.002203570, -0.002506969, -0.002000402,   0.004953133,   0.005480354,   0.005589501,   0.005064522,
      -0.002776032,  0.002280648,  0.003135123,  0.002444622,  0.002630967,  0.002661271,  0.002387629, -0.002725250, -0.002165135, -0.002443351, -0.002265278, -0.002114092, -0.002965531, -0.003220677, -0.003293003, -0.002834306,   0.002429708,   0.003085663,   0.002682295,   0.002419584,
      -0.002996357,  0.002018612,  0.002444622,  0.003070282,  0.002642961,  0.002723788,  0.002480706, -0.001900806, -0.002524132, -0.002255357, -0.002064401, -0.001987639, -0.002847756, -0.002446722, -0.002691118, -0.002146468,   0.002129979,   0.002294343,   0.002058065,   0.001647215,
      -0.003354646,  0.002292863,  0.002630967,  0.002642961,  0.003399568,  0.002835663,  0.002511498, -0.002151212, -0.002212778, -0.003079532, -0.002289923, -0.002112533, -0.002479244, -0.003238591, -0.002792185, -0.002320244,   0.001817701,   0.003018408,   0.002186349,   0.001886629,
      -0.003084640,  0.002451683,  0.002661271,  0.002723788,  0.002835663,  0.003460465,  0.002689403, -0.002307309, -0.002403937, -0.002639511, -0.003053661, -0.002357664, -0.002525455, -0.002670162, -0.003404629, -0.002380121,   0.002049861,   0.002625371,   0.002975891,   0.001982815,
      -0.002684503,  0.001962589,  0.002387629,  0.002480706,  0.002511498,  0.002689403,  0.003405908, -0.001802993, -0.001918400, -0.002147106, -0.002139558, -0.002874677, -0.002301237, -0.002427088, -0.002630598, -0.003140322,   0.001596459,   0.002273129,   0.002005465,   0.002493292,
      0.002322912, -0.005337779, -0.002725250, -0.001900806, -0.002151212, -0.002307309, -0.001802993,  0.006270013,  0.005074049,  0.005301141,  0.005490455,  0.004985750,  0.002439715,  0.002753241,  0.002945677,  0.002366814,  -0.005905703,  -0.006458191,  -0.006446063,  -0.005887221,
      0.002469943, -0.005208700, -0.002165135, -0.002524132, -0.002212778, -0.002403937, -0.001918400,  0.005074049,  0.005965899,  0.005329475,  0.005474683,  0.005034119,  0.002484915,  0.002134233,  0.002554533,  0.001829379,  -0.005698194,  -0.005349709,  -0.005537511,  -0.004806354,
      0.002971691, -0.005496785, -0.002443351, -0.002255357, -0.003079532, -0.002639511, -0.002147106,  0.005301141,  0.005329475,  0.006537889,  0.005642654,  0.005191835,  0.002136289,  0.002896811,  0.002687634,  0.002139173,  -0.004952678,  -0.006396931,  -0.005559660,  -0.005010744,
      0.002286614, -0.005680588, -0.002265278, -0.002064401, -0.002289923, -0.003053661, -0.002139558,  0.005490455,  0.005474683,  0.005642654,  0.006959385,  0.005557318,  0.002006578,  0.002179177,  0.003093230,  0.002084871,  -0.005215035,  -0.005572746,  -0.006859371,  -0.005301447,
      0.002031493, -0.005191097, -0.002114092, -0.001987639, -0.002112533, -0.002357664, -0.002874677,  0.004985750,  0.005034119,  0.005191835,  0.005557318,  0.006681979,  0.001885596,  0.002009044,  0.002422554,  0.002839850,  -0.004738902,  -0.005170758,  -0.005429394,  -0.006439851,
      0.002609980, -0.001978139, -0.002965531, -0.002847756, -0.002479244, -0.002525455, -0.002301237,  0.002439715,  0.002484915,  0.002136289,  0.002006578,  0.001885596,  0.003687895,  0.003039454,  0.003222307,  0.002651649,  -0.003081190,  -0.002846157,  -0.002563268,  -0.002134867,
      0.002976414, -0.002203570, -0.003220677, -0.002446722, -0.003238591, -0.002670162, -0.002427088,  0.002753241,  0.002134233,  0.002896811,  0.002179177,  0.002009044,  0.003039454,  0.004335243,  0.003373750,  0.002981896,  -0.002473805,  -0.003996985,  -0.002708386,  -0.002473080,
      0.002931943, -0.002506969, -0.003293003, -0.002691118, -0.002792185, -0.003404629, -0.002630598,  0.002945677,  0.002554533,  0.002687634,  0.003093230,  0.002422554,  0.003222307,  0.003373750,  0.004487640,  0.002966434,  -0.002799120,  -0.003340578,  -0.003981506,  -0.002605861,
      0.002656492, -0.002000402, -0.002834306, -0.002146468, -0.002320244, -0.002380121, -0.003140322,  0.002366814,  0.001829379,  0.002139173,  0.002084871,  0.002839850,  0.002651649,  0.002981896,  0.002966434,  0.004148426,  -0.001957908,  -0.002850689,  -0.002387628,  -0.003596330,
      -0.001475481,  0.004953133,  0.002429708,  0.002129979,  0.001817701,  0.002049861,  0.001596459, -0.005905703, -0.005698194, -0.004952678, -0.005215035, -0.004738902, -0.003081190, -0.002473805, -0.002799120, -0.001957908,   0.007351463,   0.006208475,   0.006302368,   0.005432536,
      -0.002962345,  0.005480354,  0.003085663,  0.002294343,  0.003018408,  0.002625371,  0.002273129, -0.006458191, -0.005349709, -0.006396931, -0.005572746, -0.005170758, -0.002846157, -0.003996985, -0.003340578, -0.002850689,   0.006208475,   0.008364978,   0.006655412,   0.006105928,
      -0.002332574,  0.005589501,  0.002682295,  0.002058065,  0.002186349,  0.002975891,  0.002005465, -0.006446063, -0.005537511, -0.005559660, -0.006859371, -0.005429394, -0.002563268, -0.002708386, -0.003981506, -0.002387628,   0.006302368,   0.006655412,   0.008552421,   0.006137083,
      -0.002016518,  0.005064522,  0.002419584,  0.001647215,  0.001886629,  0.001982815,  0.002493292, -0.005887221, -0.004806354, -0.005010744, -0.005301447, -0.006439851, -0.002134867, -0.002473080, -0.002605861, -0.003596330,   0.005432536,   0.006105928,   0.006137083,   0.008189789
    ),
    ncol = length(coefs)
  )
  MASS::mvrnorm(1, coefs, multiplier*m)
}

#### random sampler for log_alpha and log_theta Control vs Uninfected -------------------- ####

log_shape_sampler3 <- function(log_alpha, log_theta, multiplier){
  m = matrix(c(0.002809302, -0.001473003, -0.001473003, 0.002899220), ncol = 2)
  MASS::mvrnorm(1, c(log_alpha, log_theta), multiplier*m)
}

#### random sampler for coefficients Control vs Uninfected ------------------------------- ####

coefs_sampler3 <- function(coefs, multiplier){
  m = matrix(
    c(
      0.0002875011, -0.0002238994, -0.0002424512, -0.0002853866,  0.0002379977,  0.0002287953,  0.0002735135,  -0.0002389153,
     -0.0002238994,  0.0004814111,  0.0002440348,  0.0002936916, -0.0002731769, -0.0003596304, -0.0004660404,   0.0003524155,
     -0.0002424512,  0.0002440348,  0.0004026862,  0.0002488287, -0.0003800377, -0.0003722545, -0.0002385877,   0.0003437202,
     -0.0002853866,  0.0002936916,  0.0002488287,  0.0004646400, -0.0003982045, -0.0002539445, -0.0004707457,   0.0003990771,
      0.0002379977, -0.0002731769, -0.0003800377, -0.0003982045,  0.0006823800,  0.0003586515,  0.0003935021,  -0.0006441216,
      0.0002287953, -0.0003596304, -0.0003722545, -0.0002539445,  0.0003586515,  0.0005806781,  0.0003546247,  -0.0005491987,
      0.0002735135, -0.0004660404, -0.0002385877, -0.0004707457,  0.0003935021,  0.0003546247,  0.0007410937,  -0.0006045969,
     -0.0002389153,  0.0003524155,  0.0003437202,  0.0003990771, -0.0006441216, -0.0005491987, -0.0006045969,   0.0010752469
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

#### pdf of exponentiated weibull dist  ---------------------------- ####

pdf_calculator <- function(t, alpha, theta, sigma) {
  exp_term = exp(-(t/sigma)^alpha)
  alpha*theta/sigma * (1-exp_term)^{theta-1} * exp_term * (t/sigma)^(alpha-1)
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

#### hazard rate function ---------------------------------------- ####

hazard_rate <- function(t, alpha, theta, sigma){
  pdf_calculator(t, alpha, theta, sigma) / (1 - cdf_calculator(t, alpha, theta, sigma))
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

#### calculate hazard rate -------------------------------------- #### 

hazard_rate_calculator <- function(data, x_fit, max_day){
  data$sigma <- apply(data, 1, function(r) sigma_calculator(r, x_fit)) 
  seq(2,max_day) %>% 
    sapply(function(day)(
      apply(data, 1, function(row) (
        hazard_rate(t = day, alpha = row[1], theta = row[2], sigma = row[length(row)])
      ))
    )) %>% 
    as.data.frame()
}

#### calculate the hazard ratio ---------------------------------------- ####

hazard_ratio_calculator <- function(s1, s2){
  d1 = hazard_function_calculator(posterior_dist, s1, max_day)
  d2 = hazard_function_calculator(posterior_dist, s2, max_day)
  
  d1[sapply(d1, is.infinite)] <- NA
  d2[sapply(d2, is.infinite)] <- NA
  
  cbind(
    Day = c(2:max_day),
    t(apply(d1/d2, 2, function(col) quantile(col, probs = c(0.025,0.5,0.975), na.rm = TRUE)))
  ) %>% 
    set_colnames(c("Day","lower","est","upper")) %>% 
    na.omit()
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

## control
plot_survival_function <- function(data){
  data %>% 
    filter(est > 0.007) %>% 
    mutate(tag = paste(Treatment, Population, Sex)) %>% 
    ggplot(aes(x = Day + Age, y = est, color = tag,  ymin = lower, ymax = upper)) +
    geom_errorbar(size = 0.5) +
    geom_line(aes(linetype = tag), size = 1.2) +
    scale_color_manual(
      values = c(
        "Uninfected ACO Female" = "brown3",
        "Infected ACO Female" = "brown3",
        "Control ACO Female" = "brown3",
        "Uninfected ACO Male" = "brown3",
        "Infected ACO Male" = "brown3",
        "Control ACO Male" = "brown3",
        "Uninfected CO Female" = "grey30",
        "Infected CO Female" = "grey30",
        "Control CO Female" = "grey30",
        "Uninfected CO Male" = "grey30",
        "Infected CO Male" = "grey30",
        "Control CO Male" = "grey30"
      )) +
    scale_linetype_manual(
      values = c(
        "Uninfected ACO Female" = "dashed",
        "Infected ACO Female" = "solid",
        "Control ACO Female" = "dotted",
        "Uninfected ACO Male" = "dashed",
        "Infected ACO Male" = "solid",
        "Control ACO Male" = "dotted",
        "Uninfected CO Female" = "dashed",
        "Infected CO Female" = "solid",
        "Control CO Female" = "dotted",
        "Uninfected CO Male" = "dashed",
        "Infected CO Male" = "solid",
        "Control CO Male" = "dotted"
      )) +
    theme_minimal() +
    labs(x = "Days after spray", y = "Survival function", color = "", linetype = "")
}

## aco
plot_survival_function2 <- function(data){
  data %>% 
    filter(est > 0.007) %>% 
    mutate(tag = paste(Treatment, Population, Sex)) %>% 
    ggplot(aes(x = Day + as.numeric(as.character(Age)), y = est, ymin = lower, ymax = upper, color = tag, linetype = tag,
               group = paste(Treatment, Population, Sex, Age))) +
    geom_errorbar(size = 0.5, alpha = 0.8) +
    geom_line(size = 0.5) +
    scale_color_manual(
      values = c(
        "Uninfected ACO Female" = "brown3",
        "Infected ACO Female" = "brown3",
        "Control ACO Female" = "brown3",
        "Uninfected ACO Male" = "brown3",
        "Infected ACO Male" = "brown3",
        "Control ACO Male" = "brown3",
        "Uninfected CO Female" = "grey30",
        "Infected CO Female" = "grey30",
        "Control CO Female" = "grey30",
        "Uninfected CO Male" = "grey30",
        "Infected CO Male" = "grey30",
        "Control CO Male" = "grey30"
      )) +
    scale_linetype_manual(
      values = c(
        "Uninfected ACO Female" = "dashed",
        "Infected ACO Female" = "solid",
        "Control ACO Female" = "dotted",
        "Uninfected ACO Male" = "dashed",
        "Infected ACO Male" = "solid",
        "Control ACO Male" = "dotted",
        "Uninfected CO Female" = "dashed",
        "Infected CO Female" = "solid",
        "Control CO Female" = "dotted",
        "Uninfected CO Male" = "dashed",
        "Infected CO Male" = "solid",
        "Control CO Male" = "dotted"
      )) +
    theme_minimal() +
    labs(x = "Days after spray", y = "Survival function", color = "", linetype = "")
}


#### plot hazard ratio  ---------------------------------------------- ####

# control
plot_hazard_ratio <- function(data){
  data %>% 
    ggplot(aes(x = Day, y = est, color = tag, linetype = tag, ymin = lower, ymax = upper)) +
    geom_line(size = 1.2) +
    geom_errorbar(size = 0.7, alpha = 0.5) +
    scale_color_manual(
      values = c(
        "ACO Female" = "brown3",
        "ACO Male" = "brown3",
        "CO Female" = "grey30",
        "CO Male" = "grey30"
      )) +
    scale_linetype_manual(
      values = c(
        "ACO Female" = "dashed",
        "ACO Male" = "dashed",
        "CO Female" = "solid",
        "CO Male" = "solid"
      )) +
    theme_minimal() +
    labs(x = "Days after spray", y = "Hazard ratio", color = "", linetype = "")
}

# aco
plot_hazard_ratio <- function(data){
  data %>% 
    ggplot(aes(x = Day, y = est, color = tag, linetype = tag, ymin = lower, ymax = upper)) +
    geom_line(size = 1.2) +
    geom_errorbar(size = 0.7, alpha = 0.5) +
    scale_color_manual(
      values = c(
        "ACO Female" = "brown3",
        "ACO Male" = "brown3",
        "CO Female" = "grey30",
        "CO Male" = "grey30"
      )) +
    scale_linetype_manual(
      values = c(
        "ACO Female" = "dashed",
        "ACO Male" = "dashed",
        "CO Female" = "solid",
        "CO Male" = "solid"
      )) +
    theme_minimal() +
    labs(x = "Days after spray", y = "Hazard ratio", color = "", linetype = "")
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

