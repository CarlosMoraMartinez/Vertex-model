library(tidyverse)
setwd("/home/carmoma/vertex/Vertex-model/dpygrad_mode41_intervein/dpygrad_mode41_intervein_0/etournay1_intervein1/")
setwd("etournay1_unmoveX7/")
setwd("/home/carmoma/vertex/Vertex-model/dpygrad_mode120/dpygrad_mode120_1/etournay1_nosprings3/")

a <- read_tsv("etournay1_nosprings3_moved_40.ptab")
a <- a %>% mutate(ratio = moves_accepted/(moves_accepted+moves_rejected))

a$cols =  colorRampPalette(c('red', 'blue'))(length(a$ratio))[rank(a$ratio)]
plot(a$x, a$y, col=a$cols, type="p", pch=16, cex=1)#, xlim=c(220, 330), ylim=c(100, 200)
plot(a$x, a$y, col=a$cols, type="p", pch=16, cex=1, xlim=c(320, 400), ylim=c(40, 90))

b <- a[a$ratio < 0.25 & !is.na(a$ratio), ]
text(b$x, b$y, as.character(b$ind))

ggplot(a, aes(x=x, y=y, fill=ratio, color=ratio, alpha=0.2))+
  geom_point(size=0.75) +
  scale_fill_continuous(cols)
  scale_color_gradient(low="red", high="blue") +
  scale_fill_gradient(low="red", high="blue") +
  theme(element_blank())

#23620 seems to be the problem

c <- a %>% filter(grepl("-999", cells)) %>% filter(!grepl("^-999|-999,$", cells))
a$have2cells <- ifelse(a$ind %in% c$ind, "red", "blue")
plot(a$x, a$y, col=a$have2cells, type="p", pch=16, cex=ifelse(a$ind %in% c$ind, 2, 0.25))
plot(a$moves_accepted, a$moves_rejected, col=a$have2cells, type="p", pch=16, cex=ifelse(a$ind %in% c$ind, 2, 0.25))
plot(a$x, a$y, col=cols, type="p", pch=16, cex=ifelse(a$ind %in% c$ind, 2, 0.25))

#Other points with 2 cells like 23620 look normal

cells <- c(17,11702)
edges<-c(35104,35105,35103)
neig<-c(130,54,23151)

ct <- read_tsv("etournay1_intervein1_moved_40.celltab") %>% mutate(interest=ifelse(ind %in% cells, "red", "blue"))
plot(ct$centroid_x, ct$centroid_y, col=ct$interest)
#Something is wrong: both cells are very far away and also their edges are:
#67,35105,35103,72,76,80,34402,84,88,                   
#56,35104,60,64,68,35105,52,48,180,176,70,175,49042,181,

v1 <- c(54,23619,58,61,55,1,23151,64,23620)
v2<-c(49,32911,54,23620,130,129,418,132,37,52,36,40,43,46)
a$newcol <- ifelse(a$ind %in% v1 & a$ind %in% v2, "red", ifelse(a$ind %in% v1, "blue", ifelse(a$ind %in% v2, "green", "black")) )
plot(a$x, a$y, col=a$newcol, type="p", pch=16, 
     cex=ifelse(a$newcol != "black", 1, 0.25),
     xlim=c(200,400), ylim=c(0,170))

#35104 and 105 should be next to each other in cell 11702
numf <- 40
basename <- "etournay1_unmoveX7_moved_"
center_vert <- 24005
for(i in 3:numf){
 name <- paste(basename, as.character(i), sep="", collapse="")

  p0 <- read_tsv(paste(name, ".ptab", sep="", collapse=""))
  a<-p0
  cells <- p0[p0$ind==center_vert, "cells"] %>% gsub(",$", "", .) %>% strsplit(",") %>% unlist %>% as.integer()

  c0 <- read_tsv(paste(name, ".celltab", sep="", collapse=""))%>% mutate(interest=ifelse(ind %in% cells, "red", "blue"))
  tovec<-function(x){
   return(gsub(",$", "", x) %>% lapply(function(c){strsplit(c,",")[[1]]}))
  }
  e0 <- read_tsv(paste(name, ".edges", sep="", collapse=""))
  e0$vertices <- tovec(pull(e0, vertices))
  e0$cells <- tovec(pull(e0, cells))

  v1 <- c0$vertices[c0$ind == 18] %>% gsub(",$", "", .)%>% strsplit(",") %>% unlist %>% as.numeric
  v2 <- c0$vertices[c0$ind == 11895] %>% gsub(",$", "", .)%>% strsplit(",") %>% unlist %>% as.numeric


  a$newcol <- ifelse(a$ind %in% v1 & a$ind %in% v2, "red", ifelse(a$ind %in% v1, "blue", ifelse(a$ind %in% v2, "green", "black")) )
  v1p <- a[match(v1, a$ind),]
  v1p <- bind_rows(v1p, v1p[1,])
  v2p <- a[match(v2, a$ind),]
  v2p <- bind_rows(v2p, v2p[1,])
  both <- bind_rows(v1p, v2p)

  #pdf(paste(name, ".pdf", sep="", collapse=""))
    plot(a$x, a$y, col=a$newcol, type="p", pch=16, 
      cex=ifelse(a$newcol != "black", 1, 0.25),
      xlim=c(min(both$x)-10,max(both$x)+10), ylim=c(min(both$y)-10,max(both$y)+10))
    segments(v1p$x[1:(nrow(v1p)-1)], v1p$y[1:(nrow(v1p)-1)], v1p$x[2:(nrow(v1p))], v1p$y[2:(nrow(v1p))], col="blue")
    segments(v2p$x[1:(nrow(v2p)-1)], v2p$y[1:(nrow(v2p)-1)], v2p$x[2:(nrow(v2p))], v2p$y[2:(nrow(v2p))], col="red")

    v1p <- a[match(v1, a$ind),]
    v2p <- a[match(v2, a$ind),]
text(v1p$x+runif(1, -4, 4), v1p$y+runif(1, -2, 2), as.character(v1p$ind))
text(v2p$x[!v2p$ind %in% v1p$ind]+runif(1, -2, 2), v2p$y[!v2p$ind %in% v1p$ind]+runif(1, -2, 2), as.character(v2p$ind))
  #dev.off()
}


plot(c0$centroid_x, c0$centroid_y, col=c0$interest, pch=16)
c0 %>% filter(interest=="red")
#vertices of cell 17: 
vs<-c(49,51,53,54)
#edges of cell 17: 
es <- c(67,63,69,70)
p0 %>% filter(ind %in% vs)


pall <- data.frame()
for(i in 1:numf){
  name <- paste(basename, as.character(i), sep="", collapse="")
  
  p0 <- read_tsv(paste(name, ".ptab", sep="", collapse=""))
  p0$time <- i
  pall <- bind_rows(pall, p0)
}

px <- pall %>% filter(ind == 49)
px$move_trials <- px$moves_accepted + px$moves_rejected
px$ratio <- px$moves_accepted/px$move_trials

ggplot(px, aes(x=time, y=move_trials)) +
  geom_line() +
  geom_line(aes(y=moves_accepted, col="blue"))+
  geom_line(aes(y=moves_rejected, col="red"))

ggplot(px, aes(x=time, y=ratio)) + geom_line()
ggplot(px, aes(x=time, y=energy)) + geom_line()
