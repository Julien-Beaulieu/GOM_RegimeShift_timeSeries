# This script is a funcition to extract spaghetti posterior prediction from 
# brms models. It is used mainly in the script second_derivative_central-diff.R


extract <- function(data, model_type) {
  data <- conditional_smooths(data, spaghetti = T, ndraws = 1000)
  
  if (model_type == "GS") {
    f <- attributes(data$`mu: s(year,bs="tp")`)$spaghetti
    f <- f[, c("year", "estimate__", "sample__")]
    
    f_wide <-
      reshape2::dcast(f, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide[is.nan(f_wide)] <- 0
    return(f_wide)
  }
  
  if (model_type == "I") {
    f <- attributes(data$`mu: s(year,by=site,bs="tp",m=2)`)$spaghetti
    f_SW <- filter(f, cond__ == "site = SW Appledore")
    f_NW <- filter(f, cond__ == "site = NW Appledore")
    f_NE <- filter(f, cond__ == "site = NE Appledore")
    f_SE <- filter(f, cond__ == "site = SE Appledore")
    #f_MC <- filter(f, cond__ == "site = Malaga Cut")
    f_BC <- filter(f, cond__ == "site = Babb's Cove")
    
    f_SW <- f_SW[, c("year", "estimate__", "sample__")]
    f_NW <- f_NW[, c("year", "estimate__", "sample__")]
    f_NE <- f_NE[, c("year", "estimate__", "sample__")]
    f_SE <- f_SE[, c("year", "estimate__", "sample__")]
    #f_MC <- f_MC[,c("year", "estimate__", "sample__")]
    f_BC <- f_BC[, c("year", "estimate__", "sample__")]
    
    #NW
    f_wide_NW <-
      reshape2::dcast(f_NW, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_NW[is.nan(f_wide_NW)] <- 0
    
    #SW
    f_wide_SW <-
      reshape2::dcast(f_SW, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_SW[is.nan(f_wide_SW)] <- 0
    
    #NE
    f_wide_NE <-
      reshape2::dcast(f_NE, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_NE[is.nan(f_wide_NE)] <- 0
    
    #SE
    f_wide_SE <-
      reshape2::dcast(f_SE, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_SE[is.nan(f_wide_SE)] <- 0
    
    #MC
    #f_wide_MC <- reshape2::dcast(f_MC, year ~ sample__, value.var = "estimate__", mean) %>%
    #  as_tibble()
    #is.nan.data.frame <- function(x)
    #  do.call(cbind, lapply(x, is.nan))
    #f_wide_MC[is.nan(f_wide_MC)] <- 0
    
    #BC
    f_wide_BC <-
      reshape2::dcast(f_BC, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_BC[is.nan(f_wide_BC)] <- 0
    
    #return(list(f_wide_BC,f_wide_MC,f_wide_NE,f_wide_NW,f_wide_SE,f_wide_SW))
    return(list(f_wide_BC, f_wide_NE, f_wide_NW, f_wide_SE, f_wide_SW))
  }
  if (model_type == "I_zi") {
    f <- attributes(data$`mu: s(year,by=site,bs="tp",m=2)`)$spaghetti
    f_SW <- filter(f, cond__ == "site = SW Appledore")
    f_NW <- filter(f, cond__ == "site = NW Appledore")
    f_NE <- filter(f, cond__ == "site = NE Appledore")
    f_SE <- filter(f, cond__ == "site = SE Appledore")
    #f_MC <- filter(f, cond__ == "site = Malaga Cut")
    f_BC <- filter(f, cond__ == "site = Babb's Cove")
    
    f_SW <- f_SW[, c("year", "estimate__", "sample__")]
    f_NW <- f_NW[, c("year", "estimate__", "sample__")]
    f_NE <- f_NE[, c("year", "estimate__", "sample__")]
    f_SE <- f_SE[, c("year", "estimate__", "sample__")]
    #f_MC <- f_MC[,c("year", "estimate__", "sample__")]
    f_BC <- f_BC[, c("year", "estimate__", "sample__")]
    
    #NW
    f_wide_NW <-
      reshape2::dcast(f_NW, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_NW[is.nan(f_wide_NW)] <- 0
    
    #SW
    f_wide_SW <-
      reshape2::dcast(f_SW, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_SW[is.nan(f_wide_SW)] <- 0
    
    #NE
    f_wide_NE <-
      reshape2::dcast(f_NE, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_NE[is.nan(f_wide_NE)] <- 0
    
    #SE
    f_wide_SE <-
      reshape2::dcast(f_SE, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_SE[is.nan(f_wide_SE)] <- 0
    
    #MC
    #f_wide_MC <- reshape2::dcast(f_MC, year ~ sample__, value.var = "estimate__", mean) %>%
    #  as_tibble()
    #is.nan.data.frame <- function(x)
    #  do.call(cbind, lapply(x, is.nan))
    #f_wide_MC[is.nan(f_wide_MC)] <- 0
    
    #BC
    f_wide_BC <-
      reshape2::dcast(f_BC, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_BC[is.nan(f_wide_BC)] <- 0
    
    ### Now the zi part
    
    f <-
      attributes(data$`zi: s(year,by=site,bs="tp",m=2)`)$spaghetti
    f_SW_zi <- filter(f, cond__ == "site = SW Appledore")
    f_NW_zi <- filter(f, cond__ == "site = NW Appledore")
    f_NE_zi <- filter(f, cond__ == "site = NE Appledore")
    f_SE_zi <- filter(f, cond__ == "site = SE Appledore")
    f_BC_zi <- filter(f, cond__ == "site = Babb's Cove")
    
    f_SW_zi <- f_SW_zi[, c("year", "estimate__", "sample__")]
    f_NW_zi <- f_NW_zi[, c("year", "estimate__", "sample__")]
    f_NE_zi <- f_NE_zi[, c("year", "estimate__", "sample__")]
    f_SE_zi <- f_SE_zi[, c("year", "estimate__", "sample__")]
    f_BC_zi <- f_BC_zi[, c("year", "estimate__", "sample__")]
    
    #NW
    f_wide_NW_zi <-
      reshape2::dcast(f_NW_zi, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_NW_zi[is.nan(f_wide_NW_zi)] <- 0
    
    #SW
    f_wide_SW_zi <-
      reshape2::dcast(f_SW_zi, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_SW_zi[is.nan(f_wide_SW_zi)] <- 0
    
    #NE
    f_wide_NE_zi <-
      reshape2::dcast(f_NE_zi, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_NE_zi[is.nan(f_wide_NE_zi)] <- 0
    
    #SE
    f_wide_SE_zi <-
      reshape2::dcast(f_SE_zi, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_SE_zi[is.nan(f_wide_SE_zi)] <- 0
    
    #BC
    f_wide_BC_zi <-
      reshape2::dcast(f_BC_zi, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_BC_zi[is.nan(f_wide_BC_zi)] <- 0
    
    #return(list(f_wide_BC,f_wide_MC,f_wide_NE,f_wide_NW,f_wide_SE,f_wide_SW))
    return(
      list(
        f_wide_BC,
        f_wide_NE,
        f_wide_NW,
        f_wide_SE,
        f_wide_SW,
        f_wide_BC_zi,
        f_wide_NE_zi,
        f_wide_NW_zi,
        f_wide_SE_zi,
        f_wide_SW_zi
      )
    )
    
  }
  if (model_type == "I_hu") {
    f <- attributes(data$`mu: s(year,by=site,bs="tp",m=2)`)$spaghetti
    f_SW <- filter(f, cond__ == "site = SW Appledore")
    f_NW <- filter(f, cond__ == "site = NW Appledore")
    f_NE <- filter(f, cond__ == "site = NE Appledore")
    f_SE <- filter(f, cond__ == "site = SE Appledore")
    #f_MC <- filter(f, cond__ == "site = Malaga Cut")
    f_BC <- filter(f, cond__ == "site = Babb's Cove")
    
    f_SW <- f_SW[, c("year", "estimate__", "sample__")]
    f_NW <- f_NW[, c("year", "estimate__", "sample__")]
    f_NE <- f_NE[, c("year", "estimate__", "sample__")]
    f_SE <- f_SE[, c("year", "estimate__", "sample__")]
    #f_MC <- f_MC[,c("year", "estimate__", "sample__")]
    f_BC <- f_BC[, c("year", "estimate__", "sample__")]
    
    #NW
    f_wide_NW <-
      reshape2::dcast(f_NW, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_NW[is.nan(f_wide_NW)] <- 0
    
    #SW
    f_wide_SW <-
      reshape2::dcast(f_SW, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_SW[is.nan(f_wide_SW)] <- 0
    
    #NE
    f_wide_NE <-
      reshape2::dcast(f_NE, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_NE[is.nan(f_wide_NE)] <- 0
    
    #SE
    f_wide_SE <-
      reshape2::dcast(f_SE, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_SE[is.nan(f_wide_SE)] <- 0
    
    #MC
    #f_wide_MC <- reshape2::dcast(f_MC, year ~ sample__, value.var = "estimate__", mean) %>%
    #  as_tibble()
    #is.nan.data.frame <- function(x)
    #  do.call(cbind, lapply(x, is.nan))
    #f_wide_MC[is.nan(f_wide_MC)] <- 0
    
    #BC
    f_wide_BC <-
      reshape2::dcast(f_BC, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_BC[is.nan(f_wide_BC)] <- 0
    
    ### Now the zi part
    
    f <-
      attributes(data$`hu: s(year,by=site,bs="tp",m=2)`)$spaghetti
    f_SW_zi <- filter(f, cond__ == "site = SW Appledore")
    f_NW_zi <- filter(f, cond__ == "site = NW Appledore")
    f_NE_zi <- filter(f, cond__ == "site = NE Appledore")
    f_SE_zi <- filter(f, cond__ == "site = SE Appledore")
    f_BC_zi <- filter(f, cond__ == "site = Babb's Cove")
    
    f_SW_zi <- f_SW_zi[, c("year", "estimate__", "sample__")]
    f_NW_zi <- f_NW_zi[, c("year", "estimate__", "sample__")]
    f_NE_zi <- f_NE_zi[, c("year", "estimate__", "sample__")]
    f_SE_zi <- f_SE_zi[, c("year", "estimate__", "sample__")]
    f_BC_zi <- f_BC_zi[, c("year", "estimate__", "sample__")]
    
    #NW
    f_wide_NW_zi <-
      reshape2::dcast(f_NW_zi, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_NW_zi[is.nan(f_wide_NW_zi)] <- 0
    
    #SW
    f_wide_SW_zi <-
      reshape2::dcast(f_SW_zi, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_SW_zi[is.nan(f_wide_SW_zi)] <- 0
    
    #NE
    f_wide_NE_zi <-
      reshape2::dcast(f_NE_zi, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_NE_zi[is.nan(f_wide_NE_zi)] <- 0
    
    #SE
    f_wide_SE_zi <-
      reshape2::dcast(f_SE_zi, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_SE_zi[is.nan(f_wide_SE_zi)] <- 0
    
    #BC
    f_wide_BC_zi <-
      reshape2::dcast(f_BC_zi, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_BC_zi[is.nan(f_wide_BC_zi)] <- 0
    
    #return(list(f_wide_BC,f_wide_MC,f_wide_NE,f_wide_NW,f_wide_SE,f_wide_SW))
    return(
      list(
        f_wide_BC,
        f_wide_NE,
        f_wide_NW,
        f_wide_SE,
        f_wide_SW,
        f_wide_BC_zi,
        f_wide_NE_zi,
        f_wide_NW_zi,
        f_wide_SE_zi,
        f_wide_SW_zi
      )
    )
  }
  if (model_type == "GI") {
    f <- attributes(data$`mu: s(year,by=site,bs="tp")`)$spaghetti
    f_SW <- filter(f, cond__ == "site = SW Appledore")
    f_NW <- filter(f, cond__ == "site = NW Appledore")
    f_NE <- filter(f, cond__ == "site = NE Appledore")
    f_SE <- filter(f, cond__ == "site = SE Appledore")
    #f_MC <- filter(f, cond__ == "site = Malaga Cut")
    f_BC <- filter(f, cond__ == "site = Babb's Cove")
    
    f_SW <- f_SW[, c("year", "estimate__", "sample__")]
    f_NW <- f_NW[, c("year", "estimate__", "sample__")]
    f_NE <- f_NE[, c("year", "estimate__", "sample__")]
    f_SE <- f_SE[, c("year", "estimate__", "sample__")]
    #f_MC <- f_MC[,c("year", "estimate__", "sample__")]
    f_BC <- f_BC[, c("year", "estimate__", "sample__")]
    
    #NW
    f_wide_NW <-
      reshape2::dcast(f_NW, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_NW[is.nan(f_wide_NW)] <- 0
    
    #SW
    f_wide_SW <-
      reshape2::dcast(f_SW, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_SW[is.nan(f_wide_SW)] <- 0
    
    #NE
    f_wide_NE <-
      reshape2::dcast(f_NE, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_NE[is.nan(f_wide_NE)] <- 0
    
    #SE
    f_wide_SE <-
      reshape2::dcast(f_SE, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_SE[is.nan(f_wide_SE)] <- 0
    
    #MC
    #f_wide_MC <- reshape2::dcast(f_MC, year ~ sample__, value.var = "estimate__", mean) %>%
    #  as_tibble()
    #is.nan.data.frame <- function(x)
    #  do.call(cbind, lapply(x, is.nan))
    #f_wide_MC[is.nan(f_wide_MC)] <- 0
    
    #BC
    f_wide_BC <-
      reshape2::dcast(f_BC, year ~ sample__, value.var = "estimate__", mean) %>%
      as_tibble()
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    f_wide_BC[is.nan(f_wide_BC)] <- 0
    
    return(list(f_wide_BC, f_wide_NE, f_wide_NW, f_wide_SE, f_wide_SW))
  }
}




extract_site_MV <- function(spag = site_cover_MDS1){
  f <- spag %>% rename("year" = YEAR) %>% select(c(year,SITE,cond__,estimate__,sample__))
  f_SW <- filter(f, cond__ == "SITE = SW Appledore")
  f_NW <- filter(f, cond__ == "SITE = NW Appledore")
  f_NE <- filter(f, cond__ == "SITE = NE Appledore")
  f_SE <- filter(f, cond__ == "SITE = SE Appledore")
  f_BC <- filter(f, cond__ == "SITE = Babb's Cove")
  
  f_SW <- f_SW[, c("year", "estimate__", "sample__")]
  f_NW <- f_NW[, c("year", "estimate__", "sample__")]
  f_NE <- f_NE[, c("year", "estimate__", "sample__")]
  f_SE <- f_SE[, c("year", "estimate__", "sample__")]
  f_BC <- f_BC[, c("year", "estimate__", "sample__")]
  
  #NW
  f_wide_NW <-
    reshape2::dcast(f_NW, year ~ sample__, value.var = "estimate__", mean) %>%
    as_tibble()
  f_wide_NW[is.na(f_wide_NW)] <- 0
  
  #SW
  f_wide_SW <-
    reshape2::dcast(f_SW, year ~ sample__, value.var = "estimate__", mean) %>%
    as_tibble()
  f_wide_SW[is.na(f_wide_SW)] <- 0
  
  #NE
  f_wide_NE <-
    reshape2::dcast(f_NE, year ~ sample__, value.var = "estimate__", mean) %>%
    as_tibble()
  f_wide_NE[is.na(f_wide_NE)] <- 0
  
  #SE
  f_wide_SE <-
    reshape2::dcast(f_SE, year ~ sample__, value.var = "estimate__", mean) %>%
    as_tibble()
  f_wide_SE[is.na(f_wide_SE)] <- 0
  
  #BC
  f_wide_BC <-
    reshape2::dcast(f_BC, year ~ sample__, value.var = "estimate__", mean) %>%
    as_tibble()
  f_wide_BC[is.na(f_wide_BC)] <- 0
  
  #return(list(f_wide_BC,f_wide_MC,f_wide_NE,f_wide_NW,f_wide_SE,f_wide_SW))
  return(list(f_wide_BC, f_wide_NE, f_wide_NW, f_wide_SE, f_wide_SW))
}
