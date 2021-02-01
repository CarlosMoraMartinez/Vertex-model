
library(tidyverse)
library(wesanderson)


dirs <- c("/home/carmoma/vertex/Vertex-model/dpygrad_mode227/dpygrad_mode227_0/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode227/dpygrad_mode227_1/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode227/dpygrad_mode227_2/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode227/dpygrad_mode227_3/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode227/dpygrad_mode227_4/etournay1_strings8b")

dirs <- c("/home/carmoma/vertex/Vertex-model/dpygrad_mode229/dpygrad_mode229_0/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode229/dpygrad_mode229_1/etournay1_strings8b",
          "/home/carmoma/vertex/Vertex-model/dpygrad_mode229/dpygrad_mode229_2/etournay1_strings8b"
          )


names <- c("Only temporal gradient", "temporal + PD gradient", "temporal + AP gradient", 
           "temporal + PD + AP gradients", "temporal + PD Isaac's way")

names <- paste("perim_grad", as.character(c(0.5, 1.0, 2.0)), sep="_")


setwd("/home/carmoma/vertex/Vertex-model/dpygrad_mode229/")
for(d in dirs){
  
  setwd(d)
  f<-list.files() %>%subset(grepl(".celltab", .)) %>% subset(grepl("moved",.))
  time <- gsub("etournay1_strings8b_moved_", "", f) %>% gsub(".celltab", "", .)
  res<-data.frame();
  for(i in 1:length(f)){
    aux<-read_tsv(f[i]);
   aux$time <- time[i]
   res<-rbind(res, aux)
  };
  res$time<-as.character(res$time)
  res<-res %>% mutate(type = recode(type+1, "blade", "hinge", "vein blade", "vein hinge"))
  res2 <- res %>%filter(time %in% c("1", "5", "10", "20") & preferred_area > 0 & type %in% c("blade", "hinge") )

  res3 <- res2 %>% gather("coordinate", "centroid_position", centroid_x, centroid_y)
  res3$time <- as.numeric(res3$time)
  g1 <- ggplot(res3, aes(x=centroid_position, y=(preferred_area*(num_divisions+1)), 
                       fill=time, color=time, shape=type))+
    geom_point(size=0.5, stroke=0.5)+
    facet_grid(coordinate~.)+
    ggtitle(names[which(dirs == d)]) +
    ylab("Equilibrium Area") +
    theme_light() +
    theme(
      plot.title = element_text(size=16, face="bold.italic"),
      axis.title.x = element_text( size=16, face="bold"),
      axis.title.y = element_text(size=16, face="bold"),
      strip.text = element_text(size=16, face="bold") 
    ) 

  g2 <- ggplot(res3, aes(x=centroid_position, y=(perim_contract), 
                       fill=time, color=time, shape=type))+
    geom_point(size=0.5, stroke=0.5)+
    facet_grid(coordinate~.)+
    ggtitle(names[which(dirs == d)]) +
    ylab("Perimeter contractility") +
    theme_light() +
    theme(
      plot.title = element_text(size=16, face="bold.italic"),
      axis.title.x = element_text( size=16, face="bold"),
      axis.title.y = element_text(size=16, face="bold"),
      strip.text = element_text(size=16, face="bold") 
    ) + 
    scale_color_gradient(low="blue", high="red") +
    scale_fill_gradient(low="blue", high="red")
  g3 <- ggplot(res3, aes(x=centroid_position, y=K, 
                       fill=time, color=time, shape=type))+
    geom_point(size=0.5, stroke=0.5, alpha=0.5)+
    facet_grid(coordinate~.)+
    ggtitle(names[which(dirs == d)]) +
    ylab("K") +
    theme_light() +
    theme(
      plot.title = element_text(size=16, face="bold.italic"),
      axis.title.x = element_text( size=16, face="bold"),
      axis.title.y = element_text(size=16, face="bold"),
      strip.text = element_text(size=16, face="bold") 
    ) + 
    scale_color_gradient(low="green", high="purple") +
    scale_fill_gradient(low="green", high="purple")
  #scale_fill_grey()
  #scale_color_manual(values=wes_palette(n=4, name="GrandBudapest")) +
  #scale_color_manual(values=wes_palette(n=4, name="GrandBudapest"))
  
  res3$ratio_area <- res3$area / res3$preferred_area
  g4 <- ggplot(res3, aes(x=centroid_position, y=ratio_area, 
                 fill=time, color=time, shape=type))+
    geom_point(size=0.5, stroke=0.5)+
    geom_hline(yintercept=1, linetype = 2) +
    facet_grid(coordinate~.)+
    ggtitle(names[which(dirs == d)]) +
    xlab("Centroid position") +
    ylab("Area / Eq. Area") +
    theme_light() +
    theme(
      plot.title = element_text(size=16, face="bold.italic"),
      axis.title.x = element_text( size=16, face="bold"),
      axis.title.y = element_text(size=16, face="bold"),
      strip.text = element_text(size=16, face="bold") 
    ) 

  setwd("../../")
  pdf(paste(gsub(" ", "_", names[which(dirs == d)]), ".pdf", sep="", collapse = ""))
  print(g1)
  print(g2)
  print(g3)
  print(g4)
  dev.off()
}##for directories


