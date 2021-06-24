
library(tidyverse)
library(wesanderson)

dirs <- c("/home/carmoma/vertex/Vertex-model/dpygrad_mode232/dpygrad_mode232_0/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode229/dpygrad_mode229_1/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode229/dpygrad_mode229_2/etournay1_strings8b"
)


d <- dirs[1]

setwd(d)
## READ CELL FILES
f<-list.files() %>%subset(grepl(".celltab", .)) %>% subset(grepl("moved",.))
time <- gsub("etournay1_strings8b_moved_", "", f) %>% gsub(".celltab", "", .) %>% as.numeric
f <- f[order(time)]
time <- time[order(time)]
res<-data.frame();
for(i in 1:length(f)){
  aux<-read_tsv(f[i]);
  aux$time <- time[i]
  aux$name <- f[i]
  if(i>1){
    indmap <- match(aux$ind, aux_previous$ind)
    aux$previous_x <- aux_previous$centroid_x[indmap]
    aux$previous_y <- aux_previous$centroid_y[indmap]
    aux$previous_name <- aux_previous$name[1]
    aux$previous_time <- aux_previous$time[1]
    aux$displacement_x <- aux$centroid_x - aux_previous$centroid_x[indmap]
    aux$displacement_y <- aux$centroid_y - aux_previous$centroid_y[indmap]
    aux$displacement <- sqrt(aux$displacement_x^2 + aux$displacement_y^2)
    res<-rbind(res, aux)
  }
  aux_previous <- aux
};

res<-res %>% mutate(type = recode(type+1, "blade", "hinge", "vein blade", "vein hinge"))

res2 <- res %>% filter(!is.na(displacement) & previous_x > 10 & previous_y > 8)
plot(res2$time, res2$displacement)
res3 <- res2 %>% group_by(time, type) %>% summarise(n = n(),
                                              mean = mean(displacement), 
                                              median = median(displacement),
                                              sd=sd(displacement), 
                                              var=var(displacement),
                                              q05 = quantile(displacement, 0.05),
                                              q25 = quantile(displacement, 0.25),
                                              q75 = quantile(displacement, 0.75),
                                              q95 = quantile(displacement, 0.95)
                                              ) %>% 
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = mean + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = mean - qt((1-0.95)/2, n - 1) * sem,
         MEAN_plus_SD = mean + sd, 
         MEAN_minus_SD = mean - sd)
  


plot(res3$time, res3$mean)
plot(res3$time, res3$median)
plot(res3$time, res3$sd)

g2 <- ggplot(res2, aes(x=factor(time), fill=time, y=displacement))+
  geom_boxplot(outlier.shape = NA)+
  ylim(0, 10) +
  #ggtitle(names[which(dirs == d)]) +
  #ylab("Perimeter contractility") +
  theme_light() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") 
  ) #+ 
  #scale_color_gradient(low="blue", high="red") +
  #scale_fill_gradient(low="blue", high="red")

ggplot(res3, aes(x=time, y=median, color = type)) +
  geom_line(aes(x=time, y=median, color=type), size=1) +
  #facet_grid(type~., scales="free_y")+
  geom_ribbon(aes(ymin=q25,ymax=q75,fill=type, color=type),alpha=0.25) +
  geom_ribbon(aes(ymin=q05,ymax=q95,fill=type, color=type),alpha=0.05, linetype=2) +
  scale_fill_manual(values=wes_palette(n=4, name="GrandBudapest2")) +
  scale_color_manual(values=wes_palette(n=4, name="GrandBudapest2")) +
  theme_light() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") 
  )
ggplot(res3, aes(x=time, y=median, fill="grey70", color = "grey70")) +
  geom_line(aes(x=time, y=median, color="grey70"), size=1) +
  facet_grid(type~., scales="free_y")+
  geom_ribbon(aes(ymin=q25,ymax=q75,fill="grey70", color="grey70"),alpha=0.25) +
  geom_ribbon(aes(ymin=q05,ymax=q95,fill="grey70", color="grey70"),alpha=0.05, linetype=2) +
  geom_vline(xintercept=40, linetype = 2) +
  geom_vline(xintercept=20, linetype = 2) +
  ggtitle("Median vertex displacement") +
  ylab("Median displacement") +
  #scale_fill_manual(values=wes_palette(n=4, name="GrandBudapest2")) +
  #scale_color_manual(values=wes_palette(n=4, name="GrandBudapest2")) +
  theme_light() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") 
  ) 

ggplot(res3, aes(x=time, y=mean, fill="grey70", color = "grey70")) +
  geom_line(aes(x=time, y=mean, color="grey70"), size=1) +
  facet_grid(type~., scales="free_y")+
  geom_ribbon(aes(ymin=MEAN_plus_SD,ymax=MEAN_minus_SD,fill="grey70", color="grey70"),alpha=0.25) +
  #geom_ribbon(aes(ymin=q05,ymax=q95,fill="grey70", color="grey70"),alpha=0.05, linetype=2) +
  geom_vline(xintercept=40, linetype = 2) +
  geom_vline(xintercept=20, linetype = 2) +
  ggtitle("Mean vertex displacement +/- sd") +
  ylab("Mean displacement") +
  #scale_fill_manual(values=wes_palette(n=4, name="GrandBudapest2")) +
  #scale_color_manual(values=wes_palette(n=4, name="GrandBudapest2")) +
  theme_light() +
  theme(
    plot.title = element_text(size=16, face="bold.italic"),
    axis.title.x = element_text( size=16, face="bold"),
    axis.title.y = element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold") 
  ) 


#Plot displacement with position
b <- res2 %>% filter(time==80)
ggplot(b, aes(x=centroid_x, y=centroid_y, 
              color = displacement, fill=displacement))+
  geom_point(size=b$displacement/max(b$displacement)*2 + 0.1) +
  scale_color_gradient(low="blue", high="red") +
  scale_fill_gradient(low="blue", high="red")

c<-b %>% filter(displacement>50)

d <- res2 %>% filter(time >= 39) %>% filter(time==39 | ind %in% c$ind)
d$moved_color <- ifelse(d$ind %in% c$ind, ifelse(d$time==40, "red", "blue"),"black" )
d$moved_size <- ifelse(d$ind %in% c$ind, ifelse(d$time==40, 2, 1),0.2 )
ggplot(d, aes(x=centroid_x, y=centroid_y, 
              color = moved_color, fill=moved_color))+
  geom_point(size=d$moved_size) 



## READ VERTEX FILES

