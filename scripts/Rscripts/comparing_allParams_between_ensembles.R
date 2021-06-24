library(tidyversse)
s3 <- read.table("small3/param_files/small3_allConds.csv", header=T, sep="\t")
s4 <- read.table("small4/param_files/small4_allConds.csv", header=T, sep="\t")

s4b <- s4 %>%  select_if(names(.) %in% names(s3))
s3b <- s3 %>%  select_if(names(.) %in% names(s4))

int <- intersect(s4b[s4b$X == 196,c(-1, -61)], s3b[c(-1, -61)]) #doesn't work
names <- names(s3b[c(-1, -61)])
s3b <- s3b %>% mutate(diff = sapply(1:nrow(s3b),FUN=function(i)(s3b[i,names] ==  s4b[s4b$X == 196,names]) %>% which %>% length ))
row <- s3b %>% filter(diff == max(diff))
differences <- sapply(1:nrow(row),FUN=function(i)(row[i,names] !=  s4b[s4b$X == 196,names]) %>% which  )
difnames <- names[unique(differences)]

s3b %>% filter(name %in% row$name) %>% select(difnames) #128 and 132
s4b %>% filter(X == 196) %>% select(difnames) #t1_transition_critical_distance == 0.01
