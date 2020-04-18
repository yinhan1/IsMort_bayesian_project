
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
    geom_line(size = 0.7) +
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
  m = matrix(c(0.001960230,-0.002874102,-0.0028741025,0.004999855), ncol=2)
  MASS::mvrnorm(1, c(log_alpha, log_theta), multiplier*m)
}

#### random sampler for coefficients ACO ------------------------------- ####

coefs_sampler <- function(coefs, multiplier){
  m = matrix(
    c(
      0.0016331242, -7.546217e-04, -3.756889e-04, -0.0013879285, -0.0011487794, -1.478646e-04,  1.007800e-03,  7.730165e-04,  4.618921e-04,  3.755600e-04,
      -0.0007546217,  2.306518e-03, -1.765128e-05, -0.0006699123, -0.0001506769, -4.931073e-04, -4.392248e-04, -1.304168e-03,  1.048303e-04,  1.510364e-04,
      -0.0003756889, -1.765128e-05,  1.481559e-03,  0.0001205303,  0.0001324929, -4.211933e-04,  4.236033e-04,  2.699192e-04, -1.312530e-03, -1.246737e-03,
      -0.0013879285, -6.699123e-04,  1.205303e-04,  0.0055117460,  0.0039422929,  1.403599e-03, -3.883897e-03, -1.868431e-03, -8.909010e-04, -4.303256e-04,
      -0.0011487794, -1.506769e-04,  1.324929e-04,  0.0039422929,  0.0032579977,  1.083500e-03, -2.783775e-03, -1.669844e-03, -6.955772e-04, -5.674560e-04,
      -0.0001478646, -4.931073e-04, -4.211933e-04,  0.0014035989,  0.0010834997,  1.027390e-03, -9.371255e-04, -4.306868e-04, -1.084219e-04, -2.259865e-05,
      0.0010077997, -4.392248e-04,  4.236033e-04, -0.0038838973, -0.0027837750, -9.371255e-04,  3.938644e-03,  2.361863e-03, -6.162519e-05, -2.301085e-04,
      0.0007730165, -1.304168e-03,  2.699192e-04, -0.0018684307, -0.0016698437, -4.306868e-04,  2.361863e-03,  2.322829e-03, -1.192243e-05, -1.629997e-04,
      0.0004618921,  1.048303e-04, -1.312530e-03, -0.0008909010, -0.0006955772, -1.084219e-04, -6.162519e-05, -1.192243e-05,  1.759685e-03,  1.338649e-03,
      0.0003755600,  1.510364e-04, -1.246737e-03, -0.0004303256, -0.0005674560, -2.259865e-05, -2.301085e-04, -1.629997e-04,  1.338649e-03,  1.554111e-03
      ),
    ncol = length(coefs)
  )
  MASS::mvrnorm(1, coefs, multiplier*m)
}

#### random sampler for log_alpha and log_theta CO -------------------- ####

log_shape_sampler2 <- function(log_alpha, log_theta, multiplier){
  m = matrix(c(0.0004298026,-0.001039285,-0.001039285,0.004004406), ncol=2)
  # new proposed 2020-04-12
  # m = matrix(c(0.004765520, -0.005930974, -0.005930974, 0.007679827), ncol = 2)
  MASS::mvrnorm(1, c(log_alpha, log_theta), multiplier*m)
}

#### random sampler for coefficients CO ------------------------------- ####

coefs_sampler2 <- function(coefs, multiplier){
  m = matrix(
    c( 
       0.005699091, -0.002631979, -0.002782551, -0.003077558, -0.003360253, -0.003090080, -0.002514716,  0.002364193,  0.002570392,  0.002994810,  0.002309234,  0.001852255,  0.002696786,  0.003031882,  0.002957056,  0.002442440,  -0.001535293,  -0.003053749,  -0.002418183,  -0.001744473,
      -0.002631979,  0.005401202,  0.002223637,  0.001964786,  0.002324856,  0.002491203,  0.001914776, -0.005264697, -0.005099403, -0.005469853, -0.005572284, -0.004918411, -0.002025840, -0.002224568, -0.002587876, -0.001851319,   0.004927044,   0.005468460,   0.005667965,   0.004723761,
      -0.002782551,  0.002223637,  0.003158728,  0.002507776,  0.002654173,  0.002689964,  0.002337880, -0.002628378, -0.002137178, -0.002426010, -0.002183843, -0.002059953, -0.003007695, -0.003265178, -0.003305048, -0.002759593,   0.002306060,   0.003073867,   0.002605110,   0.002288300,
      -0.003077558,  0.001964786,  0.002507776,  0.003167965,  0.002710925,  0.002717880,  0.002329440, -0.001887202, -0.002522224, -0.002133822, -0.001917766, -0.001737309, -0.003045548, -0.002606856, -0.002720394, -0.001999414,   0.002166075,   0.002297423,   0.001962871,   0.001423694,
      -0.003360253,  0.002324856,  0.002654173,  0.002710925,  0.003473508,  0.002850420,  0.002456272, -0.002170950, -0.002266112, -0.003146166, -0.002229766, -0.002054657, -0.002574697, -0.003382807, -0.002846395, -0.002226378,   0.001861335,   0.003178775,   0.002192379,   0.001735878,
      -0.003090080,  0.002491203,  0.002689964,  0.002717880,  0.002850420,  0.003501578,  0.002609450, -0.002360295, -0.002382203, -0.002663706, -0.003073346, -0.002385575, -0.002560429, -0.002745616, -0.003468329, -0.002311093,   0.002053160,   0.002756795,   0.003086204,   0.002043395,
      -0.002514716,  0.001914776,  0.002337880,  0.002329440,  0.002456272,  0.002609450,  0.003161222, -0.001824658, -0.001754394, -0.002088580, -0.002052233, -0.002756146, -0.002209102, -0.002470232, -0.002584996, -0.002965662,   0.001561988,   0.002354778,   0.001999860,   0.002520000,
       0.002364193, -0.005264697, -0.002628378, -0.001887202, -0.002170950, -0.002360295, -0.001824658,  0.006130938,  0.005091093,  0.005357606,  0.005460782,  0.004820917,  0.002462211,  0.002727003,  0.002987808,  0.002273074,  -0.005879357,  -0.006434583,  -0.006448277,  -0.005589061,
       0.002570392, -0.005099403, -0.002137178, -0.002522224, -0.002266112, -0.002382203, -0.001754394,  0.005091093,  0.005931711,  0.005205139,  0.005300206,  0.004595046,  0.002703320,  0.002233994,  0.002620211,  0.001612458,  -0.005895538,  -0.005371603,  -0.005594649,  -0.004430475,
       0.002994810, -0.005469853, -0.002426010, -0.002133822, -0.003146166, -0.002663706, -0.002088580,  0.005357606,  0.005205139,  0.006614835,  0.005565063,  0.004945641,  0.002174246,  0.003069763,  0.002812970,  0.002024447,  -0.004990024,  -0.006620122,  -0.005742387,  -0.004702583,
       0.002309234, -0.005572284, -0.002183843, -0.001917766, -0.002229766, -0.003073346, -0.002052233,  0.005460782,  0.005300206,  0.005565063,  0.006831365,  0.005320507,  0.001995426,  0.002119809,  0.003171744,  0.001910341,  -0.005214528,  -0.005551613,  -0.006974356,  -0.005079428,
       0.001852255, -0.004918411, -0.002059953, -0.001737309, -0.002054657, -0.002385575, -0.002756146,  0.004820917,  0.004595046,  0.004945641,  0.005320507,  0.006401710,  0.001827261,  0.002051182,  0.002556525,  0.002737176,  -0.004598002,  -0.005034643,  -0.005390234,  -0.006274907,
       0.002696786, -0.002025840, -0.003007695, -0.003045548, -0.002574697, -0.002560429, -0.002209102,  0.002462211,  0.002703320,  0.002174246,  0.001995426,  0.001827261,  0.003919412,  0.003171430,  0.003255405,  0.002555015,  -0.003240102,  -0.002924006,  -0.002539424,  -0.002031349,
       0.003031882, -0.002224568, -0.003265178, -0.002606856, -0.003382807, -0.002745616, -0.002470232,  0.002727003,  0.002233994,  0.003069763,  0.002119809,  0.002051182,  0.003171430,  0.004483041,  0.003428013,  0.002953833,  -0.002496395,  -0.004171920,  -0.002608533,  -0.002383967,
       0.002957056, -0.002587876, -0.003305048, -0.002720394, -0.002846395, -0.003468329, -0.002584996,  0.002987808,  0.002620211,  0.002812970,  0.003171744,  0.002556525,  0.003255405,  0.003428013,  0.004524961,  0.002926238,  -0.002829247,  -0.003491146,  -0.004076888,  -0.002731424,
       0.002442440, -0.001851319, -0.002759593, -0.001999414, -0.002226378, -0.002311093, -0.002965662,  0.002273074,  0.001612458,  0.002024447,  0.001910341,  0.002737176,  0.002555015,  0.002953833,  0.002926238,  0.004052987,  -0.001871091,  -0.002855858,  -0.002275758,  -0.003628235,
      -0.001535293,  0.004927044,  0.002306060,  0.002166075,  0.001861335,  0.002053160,  0.001561988, -0.005879357, -0.005895538, -0.004990024, -0.005214528, -0.004598002, -0.003240102, -0.002496395, -0.002829247, -0.001871091,   0.007579362,   0.006242812,   0.006360717,   0.005292980,
      -0.003053749,  0.005468460,  0.003073867,  0.002297423,  0.003178775,  0.002756795,  0.002354778, -0.006434583, -0.005371603, -0.006620122, -0.005551613, -0.005034643, -0.002924006, -0.004171920, -0.003491146, -0.002855858,   0.006242812,   0.008583797,   0.006713186,   0.005852787,
      -0.002418183,  0.005667965,  0.002605110,  0.001962871,  0.002192379,  0.003086204,  0.001999860, -0.006448277, -0.005594649, -0.005742387, -0.006974356, -0.005390234, -0.002539424, -0.002608533, -0.004076888, -0.002275758,   0.006360717,   0.006713186,   0.008694456,   0.006015239,
      -0.001744473,  0.004723761,  0.002288300,  0.001423694,  0.001735878,  0.002043395,  0.002520000, -0.005589061, -0.004430475, -0.004702583, -0.005079428, -0.006274907, -0.002031349, -0.002383967, -0.002731424, -0.003628235,   0.005292980,   0.005852787,   0.006015239,   0.008098664
     ),
    ncol = length(coefs)
  )
  MASS::mvrnorm(1, coefs, multiplier*m)
}

#### random sampler for log_alpha and log_theta Control vs Uninfected -------------------- ####

log_shape_sampler3 <- function(log_alpha, log_theta, multiplier){
  m = diag(2) * 0.001
  # m = matrix(c(0.001008922, -0.001898425, -0.001898425, 0.003), ncol = 2)
  MASS::mvrnorm(1, c(log_alpha, log_theta), multiplier*m)
}

#### random sampler for coefficients Control vs Uninfected ------------------------------- ####

coefs_sampler3 <- function(coefs, multiplier){
  m = matrix(
    c(
       0.0003341708, -0.0002924757, -0.0002053593, -0.0002254700,  0.0002652508,  0.0003636899,  0.0002021378, -0.0003303293,
      -0.0002924757,  0.0004877895,  0.0002203237,  0.0002693134, -0.0004130929, -0.0005574818, -0.0002324381,  0.0004829751,
      -0.0002053593,  0.0002203237,  0.0003335476,  0.0002443345, -0.0003272519, -0.0002757294, -0.0003467267,  0.0003615778,
      -0.0002254700,  0.0002693134,  0.0002443345,  0.0004004281, -0.0002766106, -0.0004276332, -0.0003838086,  0.0004263236,
       0.0002652508, -0.0004130929, -0.0003272519, -0.0002766106,  0.0006431313,  0.0004984967,  0.0003512696, -0.0006970985,
       0.0003636899, -0.0005574818, -0.0002757294, -0.0004276332,  0.0004984967,  0.0008934480,  0.0004112550, -0.0008317232,
       0.0002021378, -0.0002324381, -0.0003467267, -0.0003838086,  0.0003512696,  0.0004112550,  0.0005990230, -0.0006049894,
      -0.0003303293,  0.0004829751,  0.0003615778,  0.0004263236, -0.0006970985, -0.0008317232, -0.0006049894,  0.0012039680
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

hazard_ratio_calculator <- function(data, s1, s2, max_day){
  d1 = hazard_rate_calculator(data, s1, max_day)
  d2 = hazard_rate_calculator(data, s2, max_day)
  
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
    geom_errorbar(size = 0.5, alpha = 0.5) +
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


#### calculate immunity measurement ------------------ ####

get_change_on_sigma <- function(data, s1, s2){
  apply(posterior_dist, 1, function(row) 
     {
      sigma_inf = sigma_calculator(row,s1)
      sigma_unf = sigma_calculator(row,s2)
      (sigma_inf / sigma_unf) / sigma_unf 
     }
    ) %>% 
    quantile(probs = c(0.025, 0.5, 0.975)) %>% 
    t() %>% 
    as.data.frame()
}


