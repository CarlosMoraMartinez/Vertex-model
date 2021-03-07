

rr <- res2 %>% filter(time== 0)
INITIAL_PERIM <- 0.5
expf <- function(x, exp){
  return(INITIAL_PERIM*exp(-exp*x))
}
expf2 <- function(x, exp){
  return(INITIAL_PERIM*(1 - exp(-exp*x)))
}
rr <- rr[order(rr$centroid_x),]
rr$scaled_c <- (rr$centroid_x - min(rr$centroid_x ))/(max(rr$centroid_x ) - min(rr$centroid_x ))


exponents <- seq(0.01,10,1)
#exponents <- data.frame(exp = exponents,
#                        exp_char=as.character(exponents), 
#                        xpos = 0.5,
#                      values = expf(0.5, exponents))
g1 <- ggplot(rr, aes(x=scaled_c, y = base_eq_perimcontr)) + 
  geom_point() + 
  xlim(0, 1)

colors <-c("red", "blue", "green", "yellow", "black", "gray", "red", "blue", "green", "yellow", "black", "gray","red", "blue", "green", "yellow", "black", "gray")
j <- 1
for(i in exponents){
  g1 <- g1 +  stat_function(x=a, fun=expf, args=c(i), col=colors[j])
  j = j+1
}
g1
