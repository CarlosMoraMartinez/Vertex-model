library(tidyverse)

NAMES <- c("extra", "ind", 
           "e_ind_0", "e_tension_0", "e_len_0", "e_len0_0", "e_energy_0",
           "e_ind_1", "e_tension_1", "e_len_1", "e_len0_1", "e_energy_1",
           "e_ind_2", "e_tension_2", "e_len_2", "e_len0_2", "e_energy_2",
           "e_ind_3", "e_tension_3", "e_len_3", "e_len0_3", "e_energy_3",
           "e_ind_4", "e_tension_4", "e_len_4", "e_len0_4", "e_energy_4", 
           "spr_tension", "spr_len",
           "term2", "term2final") 
setwd("/home/carmoma/vertex/Vertex-model/cut2layer9/outfiles")
 f<- list.files() %>% subset(!grepl("-", .))
 
 
 
 
 
 fi <- f[1]
 term2f <- function(ten, len, len0){ten*(len - len0)^2}
 
 
 for(fi in f){
 a<-read.table(fi, header=F, stringsAsFactors = F) %>% set_names(NAMES)
 for(i in as.character(0:4)){
   vars <- paste(c("e_tension_", "e_len_", "e_len0_", "e_energy_", "calc_energy_", "e_diff_"), i, sep="")
   a[, vars[5]] <- term2f(a[, vars[1]], a[, vars[2]], a[, vars[3]])
   a[, vars[6]] <- a[, vars[5]] - a[, vars[4]]
   cat(i, ": ", max(a[ a[,vars[6]] != 999, vars[6]]), ", ")
 }
 cat("*****\n")
 }
 
 #####
 a <- read.table(f[1], header=F, stringsAsFactors = F) %>% set_names(NAMES)
 b <- read.table(f[10], header=F, stringsAsFactors = F) %>% set_names(NAMES)
 c <- read.table(f[2], header=F, stringsAsFactors = F) %>% set_names(NAMES)
 
 res <- a %>% select(e_ind_0, e_len0_0)
 res$bl0 <- b$e_len0_0[match(res$e_ind_0, b$e_ind_0)]
 res$bl1 <- c$e_len0_0[match(res$e_ind_0, c$e_ind_0)]
 
 
 b <- a %>% filter(e_diff_1 != 0)
 
 
 b2 <- b %>% mutate(move = 1:nrow(.)) %>% gather("varname", "val",e_ind_0:spr_len) %>% 
   mutate(varname = gsub("e_", "E", varname)) %>% 
   separate(varname, c("varname", "num_term"), sep="_") %>% 
   spread(varname, val) %>% 
   mutate(diff = Elen - Elen0)


hist(b2$diff)
hist(b2$Eenergy)


 