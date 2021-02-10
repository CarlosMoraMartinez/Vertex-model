

library(tidyverse)
library(wesanderson)

dirs <- c("/home/carmoma/vertex/Vertex-model/dpygrad_mode232/dpygrad_mode232_0/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode229/dpygrad_mode229_1/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode229/dpygrad_mode229_2/etournay1_strings8b"
)


d <- dirs[1]

setwd(d)
## READ CELL FILES
f<-list.files() %>%subset(grepl(".ptab", .)) %>% subset(grepl("moved",.))
time <- gsub("etournay1_strings8b_moved_", "", f) %>% gsub(".ptab", "", .) %>% as.numeric
f <- f[order(time)]
time <- time[order(time)]
res<-data.frame();
for(i in 1:length(f)){
  aux<-read_tsv(f[i]);
  aux$time <- time[i]
  aux$name <- f[i]
  if(i>1){
    indmap <- match(aux$ind, aux_previous$ind)
    aux$previous_x <- aux_previous$x[indmap]
    aux$previous_y <- aux_previous$y[indmap]
    aux$previous_name <- aux_previous$name[1]
    aux$previous_time <- aux_previous$time[1]
    aux$displacement_x <- aux$x - aux_previous$x[indmap]
    aux$displacement_y <- aux$y - aux_previous$y[indmap]
    aux$displacement <- sqrt(aux$displacement_x^2 + aux$displacement_y^2)
    aux$moves_accepted_interval <- aux$moves_accepted - aux_previous$moves_accepted[indmap]
    aux$moves_rejected_interval <- aux$moves_rejected - aux_previous$moves_rejected[indmap]
    
    res<-rbind(res, aux)
  }
  aux_previous <- aux
};

#res<-res %>% mutate(type = recode(type+1, "blade", "hinge", "vein blade", "vein hinge"))

res2 <- res %>% filter(!is.na(displacement) & movable ==1)
res2$prop_movements_accepted <- res2$moves_accepted_interval/(res2$moves_accepted_interval+res2$moves_rejected_interval)

#plot(res2$time, res2$displacement)
res3 <- res2 %>% group_by(time) %>% summarise(n = n(),
                                                    mean = mean(displacement), 
                                                    median = median(displacement),
                                                    sd=sd(displacement), 
                                                    var=var(displacement),
                                                    q05 = quantile(displacement, 0.05),
                                                    q25 = quantile(displacement, 0.25),
                                                    q75 = quantile(displacement, 0.75),
                                                    q95 = quantile(displacement, 0.95),
                                              mean_moves = mean(prop_movements_accepted), 
                                              median_moves = median(prop_movements_accepted),
                                              sd_moves=sd(prop_movements_accepted), 
                                              var_moves=var(prop_movements_accepted),
                                              q05_moves = quantile(prop_movements_accepted, 0.05),
                                              q25_moves = quantile(prop_movements_accepted, 0.25),
                                              q75_moves = quantile(prop_movements_accepted, 0.75),
                                              q95_moves = quantile(prop_movements_accepted, 0.95)
) %>% 
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = mean + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = mean - qt((1-0.95)/2, n - 1) * sem,
         MEAN_plus_SD = mean + sd, 
         MEAN_minus_SD = mean - sd)


g2 <- ggplot(res2, aes(x=factor(time), fill=time, y=displacement))+
  geom_boxplot(outlier.shape = NA)+
  ylim(0,15)+
  ggtitle("Vertex movement with time") +
  ylab("Displacement") +
  xlab("time") +
  geom_vline(xintercept=40, linetype = 2, col="red") +
  geom_hline(yintercept= res3$median[res3$time==40], linetype = 2, col="red") +
  theme_light() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") 
  ) 
g2

ggplot(res3, aes(x=time, y=median)) +
  geom_line(aes(x=time, y=median), color="black",  size=1) +
  geom_ribbon(aes(ymin=q25,ymax=q75), fill="grey70", color="grey70",alpha=0.55) +
  geom_ribbon(aes(ymin=q05,ymax=q95), fill="grey70", color="grey70", alpha=0.15, linetype=2) +
  geom_vline(xintercept=40, linetype = 2, col = "red") +
  #geom_vline(xintercept=20, linetype = 2, col = "blue") +
  geom_hline(yintercept= res3$median[res3$time==40], linetype = 2, col="red") +
  #geom_hline(yintercept= res3$median[res3$time==20], linetype = 2, col="blue") +
  ggtitle("Median vertex displacement") +
  ylab("Vertex displacement") +
  #scale_fill_manual(values=wes_palette(n=4, name="GrandBudapest2")) +
  #scale_color_manual(values=wes_palette(n=4, name="GrandBudapest2")) +
  theme_light() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") 
  ) 

##Proportion of movesments accepted
ggplot(res3, aes(x=time, y=median_moves)) +
  geom_line(aes(x=time, y=median_moves), color="black",  size=1) +
  geom_ribbon(aes(ymin=q25_moves,ymax=q75_moves), fill="grey70", color="grey70",alpha=0.55) +
  geom_ribbon(aes(ymin=q05_moves,ymax=q95_moves), fill="grey70", color="grey70", alpha=0.15, linetype=2) +
  ggtitle("Median proportion of accepted moves per vertex") +
  ylab("Proportion of accepted moves") +
  #scale_fill_manual(values=wes_palette(n=4, name="GrandBudapest2")) +
  #scale_color_manual(values=wes_palette(n=4, name="GrandBudapest2")) +
  theme_light() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") 
  ) 



res4 <- res2 %>%gather( "move_type","number", moves_accepted_interval, moves_rejected_interval) 
res4 <- res4 %>%   mutate(move_type = recode(move_type, moves_accepted_interval="moves accepted", moves_rejected_interval="moves rejected"))
g2 <- ggplot(res4, aes(x=factor(time), y=number, fill=move_type, color=move_type))+
  geom_boxplot(outlier.shape = NA)+
  ggtitle("Movements per vertex with time") +
  ylab("Number of movements") +
  xlab("time") +
  #geom_vline(xintercept=40, linetype = 2, col="red") +
  #geom_hline(yintercept= res3$median[res3$time==40], linetype = 2, col="red") +
  theme_light() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") 
  ) +
  scale_fill_manual(values=wes_palette(n=2, name="GrandBudapest1")) +
  scale_color_manual(values=wes_palette(n=2, name="GrandBudapest1")) 
g2

res5 <- res4 %>% group_by(move_type, time) %>% summarise(n = n(),
                                              mean = mean(number), 
                                              median = median(number),
                                              sd=sd(number), 
                                              var=var(number),
                                              q05 = quantile(number, 0.05),
                                              q25 = quantile(number, 0.25),
                                              q75 = quantile(number, 0.75),
                                              q95 = quantile(number, 0.95)
) 

ggplot(res5, aes(x=time, y=median, color=move_type, fill=move_type)) +
  geom_line(aes(x=time, y=median),  size=1) +
  geom_ribbon(aes(ymin=q25,ymax=q75), alpha=0.55) +
  geom_ribbon(aes(ymin=q05,ymax=q95), alpha=0.15, linetype=2) +
  #geom_vline(xintercept=40, linetype = 2, col = "red") +
  #geom_vline(xintercept=20, linetype = 2, col = "blue") +
  #geom_hline(yintercept= res3$median[res3$time==40], linetype = 2, col="red") +
  #geom_hline(yintercept= res3$median[res3$time==20], linetype = 2, col="blue") +
  ggtitle("Movements per vertex with time") +
  ylab("Number of movements") +
  #scale_fill_manual(values=wes_palette(n=4, name="GrandBudapest2")) +
  #scale_color_manual(values=wes_palette(n=4, name="GrandBudapest2")) +
  theme_light() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") 
  ) +
  scale_fill_manual(values=wes_palette(n=2, name="GrandBudapest1")) +
  scale_color_manual(values=wes_palette(n=2, name="GrandBudapest1")) 



#plot displacement by position
xmin <- min(res2$x)
xmax <- max(res2$x)
ymin <- min(res2$y)
ymax <- max(res2$y)
mindesp <- quantile(res2$displacement, 0.05)
maxdesp <- quantile(res2$displacement, 0.95)


if(!dir.exists("displacement_with_time"))dir.create("displacement_with_time")

t<-10
for(t in 1:max(res2$time)){
 res_aux <- res2 %>% filter(time==t)
title <- paste("displacement at t = ", as.character(t), sep="", collapse="")
 gg<- ggplot(res_aux, aes(x=x, y=y, fill=displacement, color=displacement)) +
   geom_point(alpha=1, size=1)+
   xlim(xmin, xmax) +
   ylim(ymin, ymax) +
   theme_bw() +
   theme(
     plot.title = element_text(size=16, face="bold.italic"),
     axis.title.x = element_text( size=16, face="bold"),
     axis.title.y = element_text(size=16, face="bold"),
     strip.text = element_text(size=16, face="bold") 
   ) + 
   ggtitle(title) +
   scale_color_gradient(low="blue", high="red",limits = c(mindesp, maxdesp),oob = scales::squish) +
   scale_fill_gradient(low="blue", high="red",limits = c(mindesp, maxdesp),oob = scales::squish)

fname <- paste("displacement_with_time/displacement_", as.character(t), ".png", sep="", collapse="")
ggsave(fname, gg, "png", width = 25, height = 10)
}
