# This script is to make a function that change sp. name format
# from genus_spcie to Genus specie. 


fix_sp_name <- function(data, column = 1) {
  # Capitalize first letter if it starts with lowercase
  colnames(data)[column] <- "specie"
  data$specie <- ifelse(grepl("^[a-z]", data$specie),
                           paste0(toupper(substr(data$specie, 1, 1)), 
                                  substr(data$specie, 2, nchar(data$specie))),
                           data$specie)
  
  data2 <- data %>% separate(specie, c("genus","sp"), sep = "_") %>% 
    unite(specie,c("genus","sp"), sep = " ") %>% 
    mutate(specie = sub(" sp$", "", specie)) %>% 
    mutate(specie = sub(" NA$", "", specie))
  return(data2)
}
