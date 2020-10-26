library(tidyverse)

setwd("/home/carmoma/vertex/Vertex-model/test_border")

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
