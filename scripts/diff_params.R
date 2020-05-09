


ensemble <- "aprilbud2"
sims <- c(10,11,22,58)



fname <- paste(ensemble, "_allConds.csv", sep="", collapse="")
a<- read.table(fname, header=T, stringsAsFactors=F, sep="\t")

names <- paste(ensemble, as.character(sims), sep="_")
b <- a[a$name %in% names, ]
c <- b[, sapply(b, FUN=function(x)length(unique(x)))>1]


print(c)

write.table(c, paste(ensemble,"_Differences_between_good_sims.csv", sep="", collapse=""), sep="\t", row.names=F, quote=F)
