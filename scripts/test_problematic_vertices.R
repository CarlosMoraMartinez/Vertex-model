library(tidyverse)

setwd("/home/carmoma/vertex/Vertex-model/test_border")
setwd("/home/carmoma/vertex/Vertex-model/dpygrad_mode120/dpygrad_mode120_1/etournay1_nosprings3/")


a <- read_tsv("v23693.tsv")
table(a$Accepted)
9021/(37594+9021)

b <- a %>% filter(Accepted == "Accepted")
b$energy_bigger <- ifelse(b$new_energy < b$old_energy, "should_accept", "no_accept")
table(b$energy_bigger)
#no_accept should_accept 
#2237          6784 
2237/(2237+6784) #24.7% of movements are not favourable



plot(b$moves_accepted, b$new_energy)
plot(a$old_energy, a$new_energy, col=ifelse(a$Accepted=="Accepted", "blue", "red"))
a$energy_diff <- a$old_energy - a$new_energy
hist(a$energy_diff, xlab = "old energy - new energy (negative->unfavoured mov.)")

f <- list.files("etournay1_unmoveX7/", full.names = T) %>% subset(grepl("_40",.))
c <- read_tsv(f[2])
a <- read_tsv(f[6])


par(mfrow=c(1,2))
plot(c$centroid_x, c$centroid_y, col=ifelse(c$type==1, "red", "blue"),
     xlim=c(300,380), ylim=c(50,90), pch=16)
sux <- c %>% filter(centroid_y < 30)
text(sux$centroid_x, sux$centroid_y, as.character(sux$ind))
#11739 is out

plot(a$x, a$y, type="p", pch=16, cex=0.25, xlim=c(300,380), ylim=c(50,90))
plot(a$x, a$y, type="p", pch=16, cex=1,
     ylim=c(60,70), xlim=c(340,350))
text(a$x, a$y-runif(nrow(a), 0, 3), as.character(a$ind))



### Plot differences with time
setwd("/home/carmoma/vertex/Vertex-model/dpygrad_mode66/dpygrad_mode66_12/etournay1_intervein4/")
f<-list.files() %>% subset(grepl(".ptab", .))
res <- data.frame()
for(i in 1:(length(f)-1)){
  af <- f[grepl(paste("_", as.character(i-1), ".ptab", sep="", collapse=""), f)]
  bf <- f[grepl(paste("_", as.character(i), ".ptab", sep="", collapse=""), f)]
  a <- read_tsv(af) %>% select(ind, x, y, energy, moves_accepted, moves_rejected)
  b <- read_tsv(bf)%>% select(ind, x, y, energy, moves_accepted, moves_rejected)
  c <- merge(a, b, by=c("ind"), all.x=F, all.y=T)
  c$time <- i
  res <- bind_rows(res, c)
  
}

res$der_x <- res$x.x - res$x.y
res$der_y <- res$y.x - res$y.y
res$distance <- sqrt(res$der_x^2 + res$der_y^2)
res$accepted <- (res$moves_accepted.y - res$moves_accepted.x)
res$rejected <- (res$moves_rejected.y - res$moves_rejected.x)
res$ratio <- res$accepted/res$rejected
ggplot(res, aes(x=der_x, y=der_y, alpha = 0.5)) +
  geom_point() +
  facet_wrap(~time, ncol=4)

ggplot(res, aes(x=ratio)) +
  geom_histogram() +
  facet_wrap(~time, ncol=4)

ggplot(res, aes(x=der_y)) +
  geom_histogram() +
  facet_wrap(~time, ncol=4)

ggplot(res, aes(x=distance)) +
  geom_histogram() +
  facet_wrap(~time, ncol=4)


ggplot(res %>% filter(ratio < 0.3), aes(x=ratio)) +
  geom_histogram() +
  facet_wrap(~time, ncol=4)

res8<-res %>% filter(time == 4)
ggplot(res8, aes(x=x.x, y.x, fill = ratio, color=ratio)) +
  geom_point(size=ifelse(res8$ratio < 0.25, 2, 0.2)) +
  scale_fill_gradient2(low="darkred",mid="yellow", high="blue", midpoint=0.5) +
  scale_color_gradient2(low="darkred", mid="yellow",high="blue", midpoint=0.5) +
  theme(element_blank())


res8<-res %>% filter(time == 32)
ggplot(res8, aes(x=x.x, y.x, fill = distance, color=distance, alpha=0.5)) +
  geom_point(size=0.5) +
  scale_fill_gradient2(low="blue",mid="black", high="red", midpoint=1) +
  scale_color_gradient2(low="blue", mid="black",high="red", midpoint=1) +
  theme(element_blank())

ggplot(res, aes(x=x.x, y.x, fill = distance, color=distance, alpha=0.5)) +
  facet_wrap(~time, ncol=4)+
  geom_point(size=0.5) +
  scale_fill_gradient2(low="blue",mid="black", high="red", midpoint=5) +
  scale_color_gradient2(low="blue", mid="black",high="red", midpoint=5) +
  theme(element_blank())

