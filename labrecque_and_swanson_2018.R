# Added by RJCB
require(pracma) ; require(magrittr) ; require(dplyr)


# GbyE_shape <- "FTO" ; exp_window_shape <- "recent"


MR_longitudinal_sim2 <- function(GbyE_shape, exp_window_shape) {
  
  # Create age vector
  age <- seq(0, 50, length.out = 1001)
  
  # Set lifetime effect at age 50
  lifetime_effect_at_50=2
  
  # Function to create vector of G-A relation for one of four prespecified relationships
  GbyE_fn <- function(age, GbyE_shape) {
    if (GbyE_shape=="unif") {
      # Difference betwee homozygous alleles always 0.5
      return(rep(0.5, times = length(age)))
    } else if (GbyE_shape=="FTO") {
      # Difference between homozygous alleles follows a pattern roughly similar
      # to FTO (i.e. rough bell-shaped curve centred at age 25)
      return(dnorm(x = age, mean = 25, sd = 10)*13.1 + 0.2538)  # >>>>>>>>>>>>>> From appendix:       return(dnorm(x = age, mean = 25, sd = 10)*15 + 0.2538)
    } else if (GbyE_shape=="incr") {
      # Difference between homozygous alleles increases with age
      return(0.1+0.1*age/50)                                    # >>>>>>>>>>>>>> From appendix:      return(0.01+0.9*age/50)
    } else if (GbyE_shape=="decr") {
      # Difference between homozygous alleles decreases with age
      return(0.2-0.1*age/50)                                    # >>>>>>>>>>>>>> From appendix:      return(1-0.9*age/50)
    } else {
      stop("Shape argument for GbyE_fn must be either unif, FTO, incr or decr")
    }
  }
  
  lapply(c("unif","FTO","incr","decr"), FUN=function(x) {
    return(c(trapz(age, GbyE_fn(age = age, GbyE_shape = x)),
             mean(GbyE_fn(age = age, GbyE_shape = x))))
  })
  
  exp_window_fn <- function(age, exp_window_shape) {
    
    if (exp_window_shape=="unif") {
      # The exposure window is constant over time. This is a cumulative model
      # where exposure has the same effect at any time
      return(rep(1, times = length(age)))
    } else if (exp_window_shape=="recent") {
      # Recent exposure has more effect (Note that the curve at 30 is shifted
      # later in the code.)
      return(dnorm(x = age, mean = 50, sd = 10))
    } else if (exp_window_shape=="critical") {
      # Importance of exposure increases with time
      return(dnorm(x = age, mean = 14, sd = 2))  # >>>>>>>>>> From appendix: return(dnorm(x = age, mean = 25, sd = 2))
    } else if (exp_window_shape=="incr") {
      # Importance of exposure decreases with time
      return(age/50)
    } else {
      stop("Shape argument for exp_window_fn must be either unif, recent, critical or incr")
    }
    
  }
  
  # Create BMI difference by genetic variant
  BMI_diff <- GbyE_fn(age = age, GbyE_shape = GbyE_shape)
  
  
  # Create exposure window
  exp_window_wts_unscaled <- exp_window_fn(age = age, exp_window_shape = exp_window_shape)
  
  # Scale exposure window to equal lifetime effect at 50
  exp_window_wts <- lifetime_effect_at_50*exp_window_wts_unscaled/trapz(age, exp_window_wts_unscaled)
  
  
  # Calculate lifetime effect of changing BMI trajectory by 1 unit at age 30 and 
  # age 50. (Note that the if statement shifts the exposure window to where it 
  # should be at age 30.)
  est_lifetime_effect30 <- trapz(age[age<=30], (BMI_diff[age<=30]+1)*exp_window_wts[age<=30]) - trapz(age[age<=30], (BMI_diff[age<=30])*exp_window_wts[age<=30])
  if (exp_window_shape=="recent") est_lifetime_effect30 <- trapz(age[age<=30], (BMI_diff[age<=30]+1)*exp_window_wts[age>=20]) - trapz(age[age<=30], (BMI_diff[age<=30])*exp_window_wts[age>=20])
  est_lifetime_effect50 <- trapz(age, (BMI_diff+1)*exp_window_wts) - trapz(age, (BMI_diff)*exp_window_wts)
  
  # Calculate the reduced form estimate at age 30 and age 50 (Note that the line 
  # with the if statement shift the exposures to where it should be for age 30.)
  reduced_form_est30 <- trapz(age[age<=30], BMI_diff[age<=30]*exp_window_wts[age<=30])
  if (exp_window_shape=="recent") reduced_form_est30 <- trapz(age[age<=30], BMI_diff[age<=30]*exp_window_wts[age>=20])
  reduced_form_est50 <- trapz(age, BMI_diff*exp_window_wts)
  
  # MR estimates at age 30 and age 50
  MR_est30 <- round(reduced_form_est30/BMI_diff[age==30], 2)
  MR_est50 <- round(reduced_form_est50/BMI_diff[age==50], 2)
  
  # Absoulate bias at age 30 and age 50
  abs_bias30 <- MR_est30 - est_lifetime_effect30
  abs_bias50 <- MR_est50 - est_lifetime_effect50
  
  # Relative bias at age 30 and age 50
  rel_bias30 <- ((MR_est30/est_lifetime_effect30)-1)*100
  rel_bias50 <- ((MR_est50/est_lifetime_effect50)-1)*100
  
  # Arrange results
  res <- round(c(est_lifetime_effect30, MR_est30, abs_bias30, rel_bias30, 
                 est_lifetime_effect50, MR_est50, abs_bias50, rel_bias50),2)
  names(res) <- c(paste0(c("true", "MR", "abs", "rel"), 30),
                  paste0(c("true", "MR", "abs", "rel"), 50))
  
  return(res)
  
}


# 

grid <- expand.grid(GbyE_shape=c("unif","incr","decr","FTO"),exp_window_shape=c("unif","recent","critical","incr"))

ds <- do.call(rbind, mapply(grid[ ,1], grid[ ,2],FUN=function(x,y) MR_longitudinal_sim2(GbyE_shape = x, exp_window_shape = y), SIMPLIFY = FALSE)) %>% as.data.frame                
ds$GbyE <- as.character(grid$GbyE_shape)
ds$exp_win <- as.character(grid$exp_window_shape)


ds %<>% select(GbyE, exp_win, everything())
ds$GbyE <- rep(c("Constant", "Increasing", "Decreasing", "FTO"),4)
ds$exp_win <- rep(c("Uniform", "Recent", "Critical", "Increasing"), each=4)
ds <- ds[order(match(ds$GbyE,c("Constant", "Increasing", "Decreasing", "FTO"))),]
row.names(ds) <- NULL
names(ds) <- c("Genetic scenario", "Exposure window", 
               "True effect", "MR estimate", "Absolute bias", "Relative bias (%)",
               "True effect", "MR estimate", "Absolute bias", "Relative bias (%)")

ds[, c(3,4,5,7,8,9)] <- format(round(ds[, c(3,4,5,7,8,9)],1), nsmall = 1)
ds[, c(6,10)] <- format(round(ds[, c(6,10)],0), nsmall = 0)

# Added by RJCB
write.csv(ds, "table_1.csv")
