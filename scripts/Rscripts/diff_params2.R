#Just get the different parameters in an ensemble
library(tidyverse)
a<-list.files() %>% .[grep(".csv",.)] %>% read.table(header=T, sep="\t")
b <- a %>% select_if(function(x)length(unique(x))>1) %>% select(-c("X", "name"))%>% lapply(unique)


