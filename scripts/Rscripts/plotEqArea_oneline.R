
#Plot in several folders


n<-paste("etournay1_3cpv", c("", "_X","_Y", "_XY", "_Ynotime" ), sep="")
colors<-c("blue", "black", "red", "green")

i<-1;f<-paste(n[i],"/",n[i], "_moved_12.celltab", sep="", collapse="");a<-read.table(f, sep="\t", header=T);par(mfrow=c(1,2));plot(a$centroid_x, a$preferred_area, main=n[i], col=colors[a$type+1]);plot(a$centroid_y, a$preferred_area, main=n[i], col=colors[a$type+1])
