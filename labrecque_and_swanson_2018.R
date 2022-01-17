require(pracma) ; require(magrittr) ; require(dplyr)

MR_longitudinal_sim2 <- function(GA_shape, exp_window_shape) {
  # Create age vector
  age <- seq(0, 50, length.out = 1001)
  # Set lifetime effect at age 50
  lifetime_effect_at_50=2
  # Function to create vector of G-A relation for one of four prespecified relationships
  GA_fn <- function(age, GA_shape) {
    if (GA_shape=="unif") {
      # Difference between homozygous alleles always 0.5
      return(rep(0.5, times = length(age)))
    } else if (GA_shape=="FTO") {
      # Difference between homozygous alleles follows a pattern roughly similar
      # to FTO (i.e. rough bell-shaped curve centred at age 25)
      return(dnorm(x = age, mean = 25, sd = 10)*15 + 0.2538)
    } else if (GA_shape=="incr") {
      # Difference between homozygous alleles increases with age
      return(0.01+0.9*age/50)
    } else if (GA_shape=="decr") {
      # Difference between homozygous alleles decreases with age
      return(1-0.9*age/50)
    } else {
      stop("Shape argument for GA_fn must be either unif, FTO, incr or decr")
    }
  }
  # Define the the exposure window
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
      return(dnorm(x = age, mean = 25, sd = 2))
    } else if (exp_window_shape=="incr") {
      # Importance of exposure decreases with time
      return(age/50)
    } else {
      stop("Shape argument for exp_window_fn must be either unif, recent, critical or incr")
    }
    5
  }
  # Create BMI difference by genetic variant
  BMI_diff <- GA_fn(age = age, GA_shape = GA_shape)
  # Create exposure window
  exp_window_wts_unscaled <- exp_window_fn(age = age, exp_window_shape = exp_window_shape)
  # Scale exposure window to equal lifetime effect at 50
  exp_window_wts <- lifetime_effect_at_50*exp_window_wts_unscaled/
    trapz(age, exp_window_wts_unscaled)
  # Calculate lifetime effect of changing BMI trajectory by 1 unit at age 30 and
  # age 50. (Note that the if statement shifts the exposure window to where it
  # should be at age 30.)
  est_lifetime_effect30 <- trapz(age[age<=30],
                                 (BMI_diff[age<=30]+1)*exp_window_wts[age<=30]) -
    trapz(age[age<=30],
          (BMI_diff[age<=30])*exp_window_wts[age<=30])
  if (exp_window_shape=="recent") {
    est_lifetime_effect30 <- trapz(age[age<=30],
                                   (BMI_diff[age<=30]+1)*exp_window_wts[age>=20]) -
      trapz(age[age<=30], (BMI_diff[age<=30])*exp_window_wts[age>=20])
  }
  est_lifetime_effect50 <- trapz(age, (BMI_diff+1)*exp_window_wts) -
    trapz(age, (BMI_diff)*exp_window_wts)
  # Calculate the reduced form estimate at age 30 and age 50 (Note that the line
  # with the if statement shift the exposures to where it should be for age 30.)
  reduced_form_est30 <- trapz(age[age<=30], BMI_diff[age<=30]*exp_window_wts[age<=30])
  if (exp_window_shape=="recent") {
    reduced_form_est30 <- trapz(age[age<=30],
                                BMI_diff[age<=30]*exp_window_wts[age>=20])
  }
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
# Create data.frame of all combinations of G-A relationships and exposre windows
grid <- expand.grid(GA_shape=c("unif","FTO","incr","decr"),
                    exp_window_shape=c("unif","recent","critical","incr"))
ds <- do.call(rbind, mapply(grid[ ,1], grid[ ,2],
                            FUN=function(x,y) MR_longitudinal_sim2(GA_shape = x,
                                                                   exp_window_shape = y),
                            SIMPLIFY = FALSE)) %>% as.data.frame
ds$GA <- as.character(grid$GA_shape)
ds$exp_win <- as.character(grid$exp_window_shape)
ds %<>% dplyr::select(., GA, exp_win, dplyr::everything())
ds$GA <- rep(c("Time-fixed", "FTO", "Increasing", "Decreasing"),4)
ds$exp_win <- rep(c("Uniform", "Recent", "Critical", "Increasing"), each=4)
ds <- ds[order(match(ds$GA,c("Time-fixed", "FTO", "Increasing", "Decreasing"))),]
row.names(ds) <- NULL
names(ds)[1:2] <- c("GA shape", "Exposure window")
names(ds)[names(ds) %in% c("rel30","rel50")] <- c("rel30 (%)", "rel50 (%)")
ds
knitr::kable(ds)

# |GA shape   |Exposure window | true30| MR30| abs30| rel30 (%)| true50|  MR50| abs50| rel50 (%)|
# |:----------|:---------------|------:|----:|-----:|---------:|------:|-----:|-----:|---------:|
# |Time-fixed |Uniform         |   1.20| 1.20|  0.00|      0.00|      2|  2.00|  0.00|       0.0|
# |Time-fixed |Recent          |   1.99| 1.99|  0.00|     -0.23|      2|  2.00|  0.00|       0.0|
# |Time-fixed |Critical        |   1.99| 1.99|  0.00|      0.12|      2|  2.00|  0.00|       0.0|
# |Time-fixed |Increasing      |   0.72| 0.72|  0.00|      0.00|      2|  2.00|  0.00|       0.0|
# |FTO        |Uniform         |   1.20| 0.92| -0.28|    -23.33|      2|  3.93|  1.93|      96.5|
# |FTO        |Recent          |   1.99| 1.95| -0.04|     -2.24|      2|  3.03|  1.03|      51.5|
# |FTO        |Critical        |   1.99| 2.14|  0.15|      7.67|      2|  6.00|  4.00|     200.0|
# |FTO        |Increasing      |   0.72| 0.66| -0.06|     -8.33|      2|  3.93|  1.93|      96.5|
# |Increasing |Uniform         |   1.20| 0.61| -0.59|    -49.17|      2|  1.01| -0.99|     -49.5|
# |Increasing |Recent          |   1.99| 1.48| -0.51|    -25.80|      2|  1.68| -0.32|     -16.0|
# |Increasing |Critical        |   1.99| 1.66| -0.33|    -16.48|      2|  1.01| -0.99|     -49.5|
# |Increasing |Increasing      |   0.72| 0.48| -0.24|    -33.33|      2|  1.34| -0.66|     -33.0|
# |Decreasing |Uniform         |   1.20| 1.90|  0.70|     58.33|      2| 11.00|  9.00|     450.0|
# |Decreasing |Recent          |   1.99| 2.61|  0.62|     30.85|      2|  4.87|  2.87|     143.5|
# |Decreasing |Critical        |   1.99| 2.38|  0.39|     19.74|      2| 11.00|  9.00|     450.0|
# |Decreasing |Increasing      |   0.72| 1.00|  0.28|     38.89|      2|  8.00|  6.00|     300.0|
