### This script makes a function to shift from genus specie to G. specie
### Contributor: Julien Beaulieu


####____________________________________________________________________________
#### Genus specie to G. specie #################################################
####____________________________________________________________________________



genusSp2gSp <- function(data, column){
  data <- data %>% rename("col" = column)
  
  sp_index <- which(str_detect(data$col, " "))
  data_no_sp <- data[-sp_index,]
  
  data_sp <- data[sp_index,] %>% separate(col, c("gen","sp"), sep = " ") %>% 
    mutate(temp = substring(gen,1,1)) %>% 
    unite("col",c(temp,sp), sep = ". ") %>% 
    dplyr::select(-gen)
  
  data <- rbind(data_no_sp,data_sp) %>% rename({{column}} := col)
  
  return(data)
}
