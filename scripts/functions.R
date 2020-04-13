
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
  m = matrix(c(0.004765520, -0.005930974, -0.005930974, 0.007679827), ncol = 2)
  MASS::mvrnorm(1, c(log_alpha, log_theta), multiplier*m)
}

#### random sampler for coefficients CO ------------------------------- ####

coefs_sampler2 <- function(coefs, multiplier){
  m = matrix(
    c( 
      0.002379386, -0.002229623, -0.001989456, -0.002001112, -0.002022311, -0.001867147, -0.001931965,  0.002198859,  0.002407254,  0.002433021,  0.002078134,  0.002324412,  0.001980339,  0.002022347,  0.001858222,  0.001886710,  -0.002229246,  -0.002218460,  -0.001999588,  -0.002159403,
     -0.002229623,  0.009461154,  0.002165117,  0.001932411,  0.002025612,  0.001861434,  0.002080501, -0.009816314, -0.009534871, -0.009714387, -0.009583566, -0.009597880, -0.002076397, -0.002336621, -0.002023292, -0.002165693,   0.009901627,   0.010092065,   0.009979032,   0.009774144,
     -0.001989456,  0.002165117,  0.003540308,  0.001993401,  0.002021423,  0.001993054,  0.002013090, -0.003721771, -0.002158945, -0.002281721, -0.002094922, -0.002151621, -0.003561548, -0.003582930, -0.003574370, -0.003494755,   0.003855187,   0.003859802,   0.003604336,   0.003639463,
     -0.002001112,  0.001932411,  0.001993401,  0.002344157,  0.001974453,  0.001949430,  0.001917159, -0.001917510, -0.002332213, -0.002001244, -0.001930394, -0.001880601, -0.002390161, -0.002011322, -0.002011843, -0.001904774,   0.002383728,   0.001994953,   0.001936480,   0.001798980,
     -0.002022311,  0.002025612,  0.002021423,  0.001974453,  0.002486958,  0.001932854,  0.001918111, -0.002128522, -0.002035297, -0.002584179, -0.002055691, -0.001979962, -0.002052116, -0.002526618, -0.001993438, -0.001895594,   0.002200683,   0.002703262,   0.002129507,   0.002029198,
     -0.001867147,  0.001861434,  0.001993054,  0.001949430,  0.001932854,  0.002674803,  0.001889182, -0.001965051, -0.001841576, -0.001878688, -0.002582962, -0.001729556, -0.002032532, -0.002019854, -0.002762743, -0.001892967,   0.002019861,   0.001990437,   0.002701553,   0.001829111,
     -0.001931965,  0.002080501,  0.002013090,  0.001917159,  0.001918111,  0.001889182,  0.002964603, -0.002115664, -0.002049600, -0.002117171, -0.002030819, -0.003108594, -0.002030855, -0.002013338, -0.002006340, -0.002975775,   0.002118606,   0.002168995,   0.002027776,   0.003101394,
      0.002198859, -0.009816314, -0.003721771, -0.001917510, -0.002128522, -0.001965051, -0.002115664,  0.015578818,  0.009808936,  0.010185948,  0.010173289,  0.009855569,  0.003542512,  0.004014998,  0.003700559,  0.003603887,  -0.015584118,  -0.016167212,  -0.015914541,  -0.015268816,
      0.002407254, -0.009534871, -0.002158945, -0.002332213, -0.002035297, -0.001841576, -0.002049600,  0.009808936,  0.010762457,  0.009830076,  0.009553196,  0.009755918,  0.002454702,  0.002345624,  0.001974685,  0.002102229,  -0.011026228,  -0.010065436,  -0.009907881,  -0.009795306,
      0.002433021, -0.009714387, -0.002281721, -0.002001244, -0.002584179, -0.001878688, -0.002117171,  0.010185948,  0.009830076,  0.011711525,  0.009938288,  0.009823622,  0.002213872,  0.002948635,  0.002086057,  0.002286093,  -0.010336566,  -0.012140418,  -0.010479731,  -0.010160698,
      0.002078134, -0.009583566, -0.002094922, -0.001930394, -0.002055691, -0.002582962, -0.002030819,  0.010173289,  0.009553196,  0.009938288,  0.013384946,  0.009511851,  0.001993752,  0.002346276,  0.002733396,  0.002014571,  -0.010180221,  -0.010457777,  -0.014223671,  -0.009923969,
      0.002324412, -0.009597880, -0.002151621, -0.001880601, -0.001979962, -0.001729556, -0.003108594,  0.009855569,  0.009755918,  0.009823622,  0.009511851,  0.013662161,  0.002039434,  0.002382791,  0.001850949,  0.003162372,  -0.009857380,  -0.010160151,  -0.009681874,  -0.013589132,
      0.001980339, -0.002076397, -0.003561548, -0.002390161, -0.002052116, -0.002032532, -0.002030855,  0.003542512,  0.002454702,  0.002213872,  0.001993752,  0.002039434,  0.004314761,  0.003624080,  0.003634575,  0.003531149,  -0.004432179,  -0.003735342,  -0.003426474,  -0.003410062,
      0.002022347, -0.002336621, -0.003582930, -0.002011322, -0.002526618, -0.002019854, -0.002013338,  0.004014998,  0.002345624,  0.002948635,  0.002346276,  0.002382791,  0.003624080,  0.004602512,  0.003616106,  0.003502585,  -0.004178295,  -0.005167502,  -0.003946820,  -0.003986141,
      0.001858222, -0.002023292, -0.003574370, -0.002011843, -0.001993438, -0.002762743, -0.002006340,  0.003700559,  0.001974685,  0.002086057,  0.002733396,  0.001850949,  0.003634575,  0.003616106,  0.005011021,  0.003520571,  -0.003838942,  -0.003751048,  -0.004943673,  -0.003505764,
      0.001886710, -0.002165693, -0.003494755, -0.001904774, -0.001895594, -0.001892967, -0.002975775,  0.003603887,  0.002102229,  0.002286093,  0.002014571,  0.003162372,  0.003531149,  0.003502585,  0.003520571,  0.005455708,  -0.003686820,  -0.003782049,  -0.003243772,  -0.005460419,
     -0.002229246,  0.009901627,  0.003855187,  0.002383728,  0.002200683,  0.002019861,  0.002118606, -0.015584118, -0.011026228, -0.010336566, -0.010180221, -0.009857380, -0.004432179, -0.004178295, -0.003838942, -0.003686820,   0.017808439,   0.016250477,   0.015981048,   0.015150585,
     -0.002218460,  0.010092065,  0.003859802,  0.001994953,  0.002703262,  0.001990437,  0.002168995, -0.016167212, -0.010065436, -0.012140418, -0.010457777, -0.010160151, -0.003735342, -0.005167502, -0.003751048, -0.003782049,   0.016250477,   0.020559826,   0.016412648,   0.015974689,
     -0.001999588,  0.009979032,  0.003604336,  0.001936480,  0.002129507,  0.002701553,  0.002027776, -0.015914541, -0.009907881, -0.010479731, -0.014223671, -0.009681874, -0.003426474, -0.003946820, -0.004943673, -0.003243772,   0.015981048,   0.016412648,   0.023018410,   0.015320688,
     -0.002159403,  0.009774144,  0.003639463,  0.001798980,  0.002029198,  0.001829111,  0.003101394, -0.015268816, -0.009795306, -0.010160698, -0.009923969, -0.013589132, -0.003410062, -0.003986141, -0.003505764, -0.005460419,   0.015150585,   0.015974689,   0.015320688,   0.021530513
    ),
    ncol = length(coefs)
  )
  MASS::mvrnorm(1, coefs, multiplier*m)
}

#### random sampler for log_alpha and log_theta Control vs Uninfected -------------------- ####

log_shape_sampler3 <- function(log_alpha, log_theta, multiplier){
  m = matrix(c(0.001085796, -0.001337716, -0.001337716, 0.001934802), ncol = 2)
  MASS::mvrnorm(1, c(log_alpha, log_theta), multiplier*m)
}

#### random sampler for coefficients Control vs Uninfected ------------------------------- ####

coefs_sampler3 <- function(coefs, multiplier){
  m = matrix(
    c(
      2.352395e-04, -4.678037e-05, -1.840872e-04, -1.416648e-04,  1.435204e-04,  4.439736e-05,  2.535909e-05,
     -4.678037e-05,  2.344465e-04,  8.308017e-06,  1.098870e-05,  8.918215e-05, -1.370062e-04, -1.393824e-04,
     -1.840872e-04,  8.308017e-06,  3.748156e-04,  4.387675e-05, -2.443901e-04, -1.786456e-04,  1.208028e-04,
     -1.416648e-04,  1.098870e-05,  4.387675e-05,  2.129162e-04, -1.260438e-04,  8.088946e-05, -1.146378e-04,
      1.435204e-04,  8.918215e-05, -2.443901e-04, -1.260438e-04,  3.678938e-04, -2.726152e-05, -1.001764e-04,
      4.439736e-05, -1.370062e-04, -1.786456e-04,  8.088946e-05, -2.726152e-05,  3.181270e-04, -6.304845e-05,
      2.535909e-05, -1.393824e-04,  1.208028e-04, -1.146378e-04, -1.001764e-04, -6.304845e-05,  2.778154e-04
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

