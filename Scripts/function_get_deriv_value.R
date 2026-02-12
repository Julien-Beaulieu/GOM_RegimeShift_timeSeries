# This script is to make a function that takes a first or second derivative 
# and pulls out the derivatives values during regime shifts. The input is the 
# output of the Deriv function with columns as: year, y_deriv_estim, sd, upper_CI,
# and lower_CI.
# Contributor: Julien Beaulieu

values_deriv2 <- function(data = deriv2_cover_MDS1, rs1 = c(1987,1992.7),
                          rs2 = c(1999.4,2002), rs3 = c(2009.5,2011.7),
                          name = "Sessile MDS1"){
  # create a df for future output
  
  df <- data.frame(var = c(name,name,name),
                   regime_shift = c(1,2,3),
                   upper = c(NA,NA,NA),
                   lower = c(NA,NA,NA))
  
  # make 1 df for data of each regime shift
  data1 <- filter(data, rs1[1] <= data$year)
  data1 <- filter(data1, data1$year <= rs1[2])
  
  data2 <- filter(data, rs2[1] <= data$year)
  data2 <- filter(data2, data2$year <= rs2[2])
  
  data3 <- filter(data, rs3[1] <= data$year)
  data3 <- filter(data3, data3$year <= rs3[2])
  
  # make a column with upper*lower so the max is where both are + or -
  
  data1$vector = data1$upper_CI * data1$lower_CI
  data2$vector = data2$upper_CI * data2$lower_CI
  data3$vector = data3$upper_CI * data3$lower_CI
  
  # finding where vector is max
  
  ind1 <- which.max(data1$vector)
  ind2 <- which.max(data2$vector)
  ind3 <- which.max(data3$vector)
  
  if (data1$vector[ind1] > 0){
    df$upper[1] <- data1$upper_CI[ind1]
    df$lower[1] <- data1$lower_CI[ind1]
  } else{
    ind <- which.max(abs(data1$y_deriv_estim))
    df$upper[1] <- data1$upper_CI[ind]
    df$lower[1] <- data1$lower_CI[ind]
  }
  
  if (data2$vector[ind2] > 0){
    df$upper[2] <- data2$upper_CI[ind2]
    df$lower[2] <- data2$lower_CI[ind2]
  } else{
    ind <- which.max(abs(data2$y_deriv_estim))
    df$upper[2] <- data2$upper_CI[ind]
    df$lower[2] <- data2$lower_CI[ind]
  }
  
  if (data3$vector[ind3] > 0){
    df$upper[3] <- data3$upper_CI[ind3]
    df$lower[3] <- data3$lower_CI[ind3]
  } else{
    ind <- which.max(abs(data3$y_deriv_estim))
    df$upper[3] <- data3$upper_CI[ind]
    df$lower[3] <- data3$lower_CI[ind]
  }
  return(df)  
  
}

