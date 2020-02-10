
library(magrittr)


plotAll <- function(folder){
	current <- getwd()
	setwd(folder)
	files <- list.files() 
	files <- files[grep(pattern = ".out", files)]

	res <- data.frame()
	for(f in files){
   	 d <- read.table(f, stringsAsFactors=F, header=T)
  	  name <- strsplit(f, "_")[[1]] 
   	 name <- name[length(name)] %>% gsub(".out", "", .) %>% as.numeric
 	   d$time <- name
 	   res <- rbind(res, d)
	    cat(nrow(res), " ",nrow(d), "\n")
	}

	res2 <- res[res$type == 1, ]

	pdf("timeExp_2.5_coordExp_0.5_1-exp.pdf")
	plot(res2$centroid_x, res2$preferred_area, col = res2$time, pch = 19, cex = 0.2, xlab="Cell centroid horizontal position", ylab = "Equilibrium area", main = "Exp. time = 2.5, exp. coord = 0.5")

	for(t in unique(res2$time)){
 	   aux <- res2[res2$time == t, ]
	    aux <- aux[order(aux$centroid_x) , ]
	    lines(aux$centroid_x, aux$preferred_area, col = t)
	}

	res3 <- res[res$type == 0, ]
	res4 <- res[res$type == 1, ]
	res5 <- res[res$type == 2, ]
	res6 <- res[res$type == 3, ]
	par(mfrow=c(2,2))
	plot(res3$time, res3$preferred_area, col = res3$num_divisions + 1, xlab = "time", ylab = "preferred area", main = "Blade")
	plot(res4$time, res4$preferred_area, col = res4$num_divisions + 1, xlab = "time", ylab = "preferred area", main = "Hinge")
	plot(res5$time, res5$preferred_area, col = res5$num_divisions + 1, xlab = "time", ylab = "preferred area", main = "Blade Vein")
	plot(res6$time, res6$preferred_area, col = res6$num_divisions + 1, xlab = "time", ylab = "preferred area", main = "Hinge Vein")

	dev.off()
	setwd(current)
	cat("\n\n***\n\n")

}


dirs <- paste( "/home/carmoma/vertex/Vertex-model/ensemble17/ensemble17_", as.character(0:15), "/wing2E", sep="")
sapply(dirs, plotAll)


rr <- res5
rr<- rr[rr$num_divisions >= 10,]
dividers<-rr$ind %>% unique
r13 <- res[res$ind %in% dividers, ]
rows= as.integer(round(sqrt(length(dividers))))
par(mfrow=c(rows, rows))
for(d in dividers){
	plot(r13$time[r13$ind == d], r13$preferred_area[r13$ind == d], col = r13$num_divisions[r13$ind == d])
}

