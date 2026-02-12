# This script is a function to calculate First derivative with forward
# difference the input is a dataframe with the year as first colums and the
# draws posterior prediction as column. This is also the output of the extract 
# function. 

# contributor: Julien Beaulieu

Deriv <- function(data, model_type = "other"){
  if(model_type == "other"){
    mty1 <- data[,-1]
    
    #calculate numerateur 
    mty0 <- mty1[-nrow(mty1),]
    mty1 <- mty1[-1,]
    
    num = mty1 - mty0
    
    #calculate denominateur
    x0 <- data$year[-nrow(data)]
    x1 <- data$year[-1]
    
    deno <- (x1-x0)
    
    #calculate derivative
    mt_deriv <- as.data.frame(matrix(ncol = ncol(mty1), nrow = nrow(num)))
    for (i in c(1:ncol(num))) {
      mt_deriv[,i] <- num[,i]/deno
    }  
    
    deriv <-  as.data.frame(matrix(nrow = nrow(mt_deriv), ncol = 5)) 
    names(deriv) <- c("year","y_deriv_estim","sd", "upper_CI", "lower_CI")
    deriv$year <- x1
    
    for(i in c(1:nrow(deriv))) {
      deriv$y_deriv_estim[i] <- mean(as.numeric(mt_deriv[i,]))
      deriv$sd[i] <- sd(as.numeric(mt_deriv[i,]))
      deriv$upper_CI[i] <- ci(as.numeric(mt_deriv[i,]))$CI_high
      deriv$lower_CI[i] <- ci(as.numeric(mt_deriv[i,]))$CI_low
    }
    return(deriv)
  }
  #=============================================================================
  if(model_type == "GI"){
    mty1 <- data[,-1]
    
    #calculate numerateur 
    mty0 <- mty1[-c((nrow(mty1)-3): nrow(mty1)),]
    mty2 <- mty1[-c(1:4),]
    mty1 <- mty1[-c(1,2,nrow(data), nrow(data)-1),]
    
    num = mty2 - 2*mty1 + mty0
    
    #calculate denominateur
    x0 <- data$year[-c((nrow(mty1)-3): nrow(mty1))]
    x1 <- data$year[-c(1,2,nrow(data),(nrow(data)-1))]
    x2 <- data$year[-c(1:4)]
    
    deno <- (x1-x0)*(x2-x1)
    
    #calculate derivative
    mt_deriv <- as.data.frame(matrix(ncol = ncol(mty1), nrow = nrow(num)))
    for (i in c(1:ncol(num))) {
      mt_deriv[,i] <- num[,i]/deno
    }  
    
    deriv <-  as.data.frame(matrix(nrow = nrow(mt_deriv), ncol = 5)) 
    names(deriv) <- c("year","y_deriv_estim","sd", "upper_CI", "lower_CI")
    deriv$year <- x1
    
    for(i in c(1:nrow(deriv))) {
      deriv$y_deriv_estim[i] <- mean(as.numeric(mt_deriv[i,]))
      deriv$sd[i] <- sd(as.numeric(mt_deriv[i,]))
      deriv$upper_CI[i] <- ci(as.numeric(mt_deriv[i,]))$CI_high
      deriv$lower_CI[i] <- ci(as.numeric(mt_deriv[i,]))$CI_low
    }
    return(deriv)
  }
}