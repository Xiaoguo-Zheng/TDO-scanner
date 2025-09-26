setwd("C:/Users/zxg/Desktop/rdir/qyb")

library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)

colors2 <- c("#96C3D8","#F5B375","#C0937E","#67A59B","#A5D38F", "#8D75AF","#F19294","#E45D61","#BDA7CB" )
#色板
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4',
               '#AB3282', '#23452F', '#BD956A','#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', 
               '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#CCE0F5', '#CCC9E6', '#625D9E',
               '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#68A180', '#3A6963', '#968175' )#颜色设

#读取数据
data1 <- read.table("Figure1.txt",head=T)
data2 <- read.table("Figure2.txt",head=T,sep="\t")
markers <- read.table("markers.txt",head=T)
marker1 <- read.table("maker1.txt",head=T)
show1 <- read.table("show1.txt",head=T)
marker2 <- read.table("maker2.txt",head=T)
show2 <- read.table("show2.txt",head=T)


head(data1)
head(data2)
head(markers)

#合并数据
dat0 <- merge(data1, markers, by = "Gene", all = TRUE)
dat0 <- subset(dat0,logP<=60)

dat1 <- merge(data1, marker1, by = "Gene", all = TRUE)
dat1 <- merge(dat1, show1, by = "Gene", all = TRUE)
dat1 <- subset(dat1,logP<=60)

dat2 <- merge(data1, marker2, by = "Gene", all = TRUE)
dat2 <- merge(dat2, show2, by = "Gene", all = TRUE)
dat2 <- subset(dat2,logP<=60)

dat3 <- merge(data2, markers, by = "Gene", all = TRUE)

marker_subsets1 <- unique(marker1$Marker)
marker_subsets2 <- unique(marker2$Marker)

##############################################
pdf("sgmut_figures1.pdf",height=7,width=9)
for(marker in marker_subsets1){
p <- ggplot(dat1, aes(x = LFC, y = logP)) +
  geom_point(color = "grey", alpha = 0.7, size = 2) +
  geom_point(data = subset(dat1, Marker == marker), aes(x = LFC, y = logP), color =  "navyblue", alpha=0.7,size = 2) +
  ggrepel::geom_text_repel(subset(dat1,Marker==marker & show==marker),
                           mapping=aes(LFC, logP,label=((subset(dat1,Marker==marker & show==marker))$Gene)),color="black",
                           force = 30,size=8, max.overlaps = 20,seed = 233, min.segment.length = 0,
                           force_pull = 2,box.padding = 0.5,segment.linetype = 3, segment.color = 'black', 
                           segment.alpha = 1, direction = "both",  hjust = 0.5,show.legend = F)+
  labs(x = "", y = "", color = "Marker",title = marker) +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.text = element_text(size = 26),
        axis.title = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(angle = 90, vjust = 0.5),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(colour = "black", size = 0.1)) +
  scale_x_continuous(breaks = seq(-8, 4, by = 2)) +  # 设置X轴刻度间隔为2
  guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_y_continuous(breaks = seq(0, 60, by = 10))   # 设置Y轴刻度间隔为10
print(p)
}
dev.off()


pdf("sgmut_figures2.pdf",height=7,width=7)
for(marker in marker_subsets2){
  p <- ggplot(dat2, aes(x = LFC, y = logP)) +
    geom_point(color = "grey", alpha = 0.7, size = 2) +
    geom_point(data = subset(dat2, Marker == marker), aes(x = LFC, y = logP), color =  "navyblue", alpha=0.7,size = 2) +
    ggrepel::geom_text_repel(subset(dat2,Marker==marker & show==marker),
                             mapping=aes(LFC, logP,label=((subset(dat2,Marker==marker & show==marker))$Gene)),color="black",
                             force = 30,size=8, max.overlaps = 20,seed = 233, min.segment.length = 0,
                             force_pull = 2,box.padding = 0.5,segment.linetype = 3, segment.color = 'black', 
                             segment.alpha = 1, direction = "both",  hjust = 0.5,show.legend = F)+
    labs(x = "", y = "", color = "Marker",title = marker) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
          axis.text = element_text(size = 26),
          axis.title = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(angle = 90, vjust = 0.5),
          axis.ticks.length = unit(0.2, "cm"),
          axis.ticks = element_line(colour = "black", size = 0.1)) +
    scale_x_continuous(breaks = seq(-8, 4, by = 2)) +  # 设置X轴刻度间隔为2
    guides(color = guide_legend(override.aes = list(size = 3))) +
    scale_y_continuous(breaks = seq(0, 60, by = 10))   # 设置Y轴刻度间隔为10
  print(p)
}
dev.off()
##########################################################



















x1 <- ggplot(dat0, aes(x = LFC, y = logP)) +
  geom_point(aes(color = "grey"), alpha = 0.7, size = 2) +
  geom_point(data = subset(dat0, marker == "m1"), aes(x = LFC, y = logP, color = "m1"), size = 2) +
  geom_point(data = subset(dat0, marker == "m2"), aes(x = LFC, y = logP, color = "m2"), size = 2) +
  geom_point(data = subset(dat0, marker == "m3"), aes(x = LFC, y = logP, color = "m3"), size = 2) +
  geom_point(data = subset(dat0, marker == "m4"), aes(x = LFC, y = logP, color = "m4"), size = 2) +
  geom_point(data = subset(dat0, marker == "m5"), aes(x = LFC, y = logP, color = "m5"), size = 2) +
  labs(x = "Log2(FoldChange)", y = "-log10(pvalue)", color = "Marker") +
  scale_color_manual(name = "Marker",
                     values = c("m1" = "#96C3D8", "m2" = "#F5B375", "m3" = "#8D75AF", 
                                "m4" = "#F19294", "m5" = "#E45D61"),
                     breaks = c("m1", "m2", "m3", "m4", "m5"),
                     labels = c("Marker 1", "Marker 2", "Marker 3", "Marker 4", "Marker 5")) +
  labs(x = "Log2(FoldChange)", y = "-log10(pvalue)", color = "Marker",title = "") +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(angle = 90, vjust = 0.5),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(colour = "black", size = 0.1)) +
  scale_x_continuous(breaks = seq(-8, 4, by = 2)) +  # 设置X轴刻度间隔为2
  guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_y_continuous(breaks = seq(0, 60, by = 10))   # 设置Y轴刻度间隔为10








x2 <- ggplot(dat0, aes(x = LFC, y = logP)) +
  geom_point(aes(color = "grey"), alpha = 0.7, size = 2) +
  geom_point(data = subset(dat0, marker == "m1"), aes(x = LFC, y = logP, color = "m1"), size = 2) +
  geom_point(data = subset(dat0, marker == "m2"), aes(x = LFC, y = logP, color = "m2"), size = 2) +
  geom_point(data = subset(dat0, marker == "m3"), aes(x = LFC, y = logP, color = "m3"), size = 2) +
  geom_point(data = subset(dat0, marker == "m4"), aes(x = LFC, y = logP, color = "m4"), size = 2) +
  geom_point(data = subset(dat0, marker == "m5"), aes(x = LFC, y = logP, color = "m5"), size = 2) +
  ggrepel::geom_text_repel(subset(dat0,marker=="m1"),mapping=aes(LFC, logP,label=((subset(dat0,marker=="m1"))$Gene)),color="#96C3D8",
                           force = 80,size=4, max.overlaps = 20,seed = 233, min.segment.length = 0,
                           force_pull = 2,box.padding = 0.1,segment.linetype = 3, segment.color = 'black', 
                           segment.alpha = 1, direction = "both",  hjust = 0.5,show.legend = F)+
  ggrepel::geom_text_repel(subset(dat0,marker=="m2"),mapping=aes(LFC, logP,label=((subset(dat0,marker=="m2"))$Gene)),color="#F5B375",
                           force = 80,size=4, max.overlaps = 20,seed = 233, min.segment.length = 0,
                           force_pull = 2,box.padding = 0.1,segment.linetype = 3, segment.color = 'black', 
                           segment.alpha = 1, direction = "both",  hjust = 0.5,show.legend = F)+
  ggrepel::geom_text_repel(subset(dat0,marker=="m3"),mapping=aes(LFC, logP,label=((subset(dat0,marker=="m3"))$Gene)),color="#8D75AF",
                           force = 80,size=4, max.overlaps = 20,seed = 233, min.segment.length = 0,
                           force_pull = 2,box.padding = 0.1,segment.linetype = 3, segment.color = 'black', 
                           segment.alpha = 1, direction = "both",  hjust = 0.5,show.legend = F)+
  ggrepel::geom_text_repel(subset(dat0,marker=="m4"),mapping=aes(LFC, logP,label=((subset(dat0,marker=="m4"))$Gene)),color="#F19294",
                           force = 80,size=4, max.overlaps = 20,seed = 233, min.segment.length = 0,
                           force_pull = 2,box.padding = 0.1,segment.linetype = 3, segment.color = 'black', 
                           segment.alpha = 1, direction = "both",  hjust = 0.5,show.legend = F)+
  ggrepel::geom_text_repel(subset(dat0,marker=="m5"),mapping=aes(LFC, logP,label=((subset(dat0,marker=="m5"))$Gene)),color="#E45D61",
                           force = 80,size=4, max.overlaps = 20,seed = 233, min.segment.length = 0,
                           force_pull = 2,box.padding = 0.1,segment.linetype = 3, segment.color = 'black', 
                           segment.alpha = 1, direction = "both",  hjust = 0.5,show.legend = F)+
  labs(x = "Log2(FoldChange)", y = "-log10(pvalue)", color = "Marker") +
  scale_color_manual(name = "Marker",
                     values = c("m1" = "#96C3D8", "m2" = "#F5B375", "m3" = "#8D75AF", 
                                "m4" = "#F19294", "m5" = "#E45D61"),
                     breaks = c("m1", "m2", "m3", "m4", "m5"),
                     labels = c("Marker 1", "Marker 2", "Marker 3", "Marker 4", "Marker 5")) +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(angle = 90, vjust = 0.5),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(colour = "black", size = 0.1)) +
  scale_x_continuous(breaks = seq(-8, 4, by = 2)) +  # 设置X轴刻度间隔为2
  scale_y_continuous(breaks = seq(0, 60, by = 10)) +  # 设置Y轴刻度间隔为10
guides(color = guide_legend(override.aes = list(size = 3))) 

  
  
  
  
x3 <-  ggplot(data2, aes(x = LFC_relative_to_W, y = LFC_relative_to_Negative)) +
  geom_point(aes(color = "grey"), alpha = 0.7, size = 2) +
  geom_smooth(se=T,method="lm")+
  geom_point(data = subset(dat3, marker == "m1"), aes(x = LFC_relative_to_W, y = LFC_relative_to_Negative, color = "m1"), size = 2) +
  geom_point(data = subset(dat3, marker == "m2"), aes(x = LFC_relative_to_W, y = LFC_relative_to_Negative, color = "m2"), size = 2) +
  geom_point(data = subset(dat3, marker == "m3"), aes(x = LFC_relative_to_W, y = LFC_relative_to_Negative, color = "m3"), size = 2) +
  geom_point(data = subset(dat3, marker == "m4"), aes(x = LFC_relative_to_W, y = LFC_relative_to_Negative, color = "m4"), size = 2) +
  geom_point(data = subset(dat3, marker == "m5"), aes(x = LFC_relative_to_W, y = LFC_relative_to_Negative, color = "m5"), size = 2) +
  labs(x = "LFC_relative_to_W", y = "LFC_relative_to_Negative", color = "Marker") +
  scale_color_manual(name = "Marker",
                     values = c("m1" = "#96C3D8", "m2" = "#F5B375", "m3" = "#8D75AF", 
                                "m4" = "#F19294", "m5" = "#E45D61"),
                     breaks = c("m1", "m2", "m3", "m4", "m5"),
                     labels = c("Marker 1", "Marker 2", "Marker 3", "Marker 4", "Marker 5")) +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(angle = 90, vjust = 0.5),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(colour = "black", size = 0.1)) +
  scale_x_continuous(breaks = seq(-8, 8, by = 2)) +  # 设置X轴刻度间隔为2
  scale_y_continuous(breaks = seq(-8, 8, by = 2)) +  # 设置Y轴刻度间隔为10
  guides(color = guide_legend(override.aes = list(size = 3))) 
  








ggplot(data2, aes(x = LFC_relative_to_W, y = LFC_relative_to_Negative)) +
  geom_point(aes(color = "grey"), alpha = 0.7, size =2) +
  geom_smooth(se=T,method="lm")
#计算Rsuqare
round(summary(lm(LFC_relative_to_W ~ LFC_relative_to_Negative, data = data2))$r.squared, 2)




x4 <- ggplot(data2, aes(x = LFC_relative_to_W, y = LFC_relative_to_Negative)) +
  geom_point(aes(color = "grey"), alpha = 0.7, size =2) +
  geom_point(data = subset(dat3, marker == "m1"), aes(x = LFC_relative_to_W, y = LFC_relative_to_Negative, color = "m1"), size = 2) +
  geom_point(data = subset(dat3, marker == "m2"), aes(x = LFC_relative_to_W, y = LFC_relative_to_Negative, color = "m2"), size = 2) +
  geom_point(data = subset(dat3, marker == "m3"), aes(x = LFC_relative_to_W, y = LFC_relative_to_Negative, color = "m3"), size = 2) +
  geom_point(data = subset(dat3, marker == "m4"), aes(x = LFC_relative_to_W, y = LFC_relative_to_Negative, color = "m4"), size = 2) +
  geom_point(data = subset(dat3, marker == "m5"), aes(x = LFC_relative_to_W, y = LFC_relative_to_Negative, color = "m5"), size = 2) +
  ggrepel::geom_text_repel(subset(dat3,marker=="m1"),mapping=aes(x = LFC_relative_to_W, y = LFC_relative_to_Negative,label=((subset(dat3,marker=="m1"))$Gene)),color="#96C3D8",
                           force = 80,size=4, max.overlaps = 20,seed = 233, min.segment.length = 5,
                           force_pull = 2,box.padding = 0.1,segment.linetype = 3, segment.color = 'black', 
                           segment.alpha = 1, direction = "both",  hjust = 0.5,show.legend = F)+
  ggrepel::geom_text_repel(subset(dat3,marker=="m2"),mapping=aes(x = LFC_relative_to_W, y = LFC_relative_to_Negative,label=((subset(dat3,marker=="m2"))$Gene)),color="#F5B375",
                           force = 80,size=4, max.overlaps = 20,seed = 233, min.segment.length = 0,
                           force_pull = 2,box.padding = 0.1,segment.linetype = 3, segment.color = 'black', 
                           segment.alpha = 1, direction = "both",  hjust = 0.5,show.legend = F)+
  ggrepel::geom_text_repel(subset(dat3,marker=="m3"),mapping=aes(x = LFC_relative_to_W, y = LFC_relative_to_Negative,label=((subset(dat3,marker=="m3"))$Gene)),color="#8D75AF",
                           force = 80,size=4, max.overlaps = 20,seed = 233, min.segment.length = 0,
                           force_pull = 2,box.padding = 0.1,segment.linetype = 3, segment.color = 'black', 
                           segment.alpha = 1, direction = "both",  hjust = 0.5,show.legend = F)+
  ggrepel::geom_text_repel(subset(dat3,marker=="m4"),mapping=aes(x = LFC_relative_to_W, y = LFC_relative_to_Negative,label=((subset(dat3,marker=="m4"))$Gene)),color="#F19294",
                           force = 80,size=4, max.overlaps = 20,seed = 233, min.segment.length = 0,
                           force_pull = 2,box.padding = 0.1,segment.linetype = 3, segment.color = 'black', 
                           segment.alpha = 1, direction = "both",  hjust = 0.5,show.legend = F)+
  ggrepel::geom_text_repel(subset(dat3,marker=="m5"),mapping=aes(x = LFC_relative_to_W, y = LFC_relative_to_Negative,label=((subset(dat3,marker=="m5"))$Gene)),color="#E45D61",
                           force = 80,size=4, max.overlaps = 20,seed = 233, min.segment.length = 0,
                           force_pull = 2,box.padding = 0.1,segment.linetype = 3, segment.color = 'black', 
                           segment.alpha = 1, direction = "both",  hjust = 0.5,show.legend = F)+
  labs(x = "LFC_relative_to_W", y = "LFC_relative_to_Negative", color = "Marker") +
  scale_color_manual(name = "Marker",
                     values = c("m1" = "#96C3D8", "m2" = "#F5B375", "m3" = "#8D75AF", 
                                "m4" = "#F19294", "m5" = "#E45D61"),
                     breaks = c("m1", "m2", "m3", "m4", "m5"),
                     labels = c("Marker 1", "Marker 2", "Marker 3", "Marker 4", "Marker 5")) +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(angle = 90, vjust = 0.5),
        axis.ticks.length = unit(0.2, "cm"),
        axis.ticks = element_line(colour = "black", size = 0.1)) +
  scale_x_continuous(breaks = seq(-8, 8, by = 2)) +  # 设置X轴刻度间隔为2
  scale_y_continuous(breaks = seq(-8, 8, by = 2)) +  # 设置Y轴刻度间隔为10
  guides(color = guide_legend(override.aes = list(size = 3))) 



  
pdf("sgmut_h7w9.pdf",height=7,width=9)
print(x1)
print(x2)
      print(x3)
            print(x4)
dev.off()
  
  
  