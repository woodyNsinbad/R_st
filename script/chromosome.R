
### PLOT COVERAGE ###
library("ggplot2")
library("scales")
source("~/Desktop/Work/R_script/FormatSI.R")

########
chromosome <- read.delim("~/Desktop/Work/6Y002.bed",header = F)
names(chromosome) <- c("index","chr","Start","End")
chromosome2 <- read.delim("~/Desktop/Work/6Y004.bed",header = F)
names(chromosome2) <- c("index","chr","Start","End")
chromosome3 <- read.delim("~/Desktop/Work/6Y006.bed",header = F)
names(chromosome3) <- c("index","chr","Start","End")
TelCen <- read.delim("~/Desktop/Work/TelCen.r.bed",header = F)
names(TelCen) <- c("index","chr","Start","End")
gap <- read.delim("~/Desktop/Work/gap.txt",header = F)
names(gap) <- c("index","chr","Start","End")
FiveColorSet <- c("#F25F4F","#FFB12E","#1091A4","#9BBD73","#295570")
#FiveColorSet <- RColorBrewer::brewer.pal(5,"Dark2")
######
hg19 <- read.delim("~/Desktop/Work/hg19.chromosome.full.length.txt",header = F)
names(hg19) <- c("index","chr","Start","End")

###### to combine 2 dataset into a plot ####
leg <- data.frame(Names = c("6Y006","6Y004","6Y002","Telomere/centromere","N"),
                  start = c(0,0,0,0,0),
                  end = c(0,0,0,0,0)
                  )


p <- ggplot(data=chromosome) +
    geom_rect(aes(xmin = Start ,ymin= index , xmax = End,ymax = index+0.2),
              fill = FiveColorSet[1],color=FiveColorSet[1]) +
    scale_x_continuous(labels = FormatSI) +
    scale_y_continuous(breaks = c(24:1) , labels = hg19$chr) + 
  
    geom_point(data=leg,aes(x= start ,y = end, color = Names,shape = NA)) +
  
    geom_rect(data=chromosome2,aes(xmin = Start ,ymin=index+0.2,xmax = End,ymax = index+0.4),
             fill = FiveColorSet[2],color=FiveColorSet[2]) +
    geom_rect(data=chromosome3,aes(xmin = Start ,ymin=index+0.4,xmax = End,ymax = index+0.6),
            fill = FiveColorSet[3],color=FiveColorSet[3]) +
    geom_rect(data=gap,aes(xmin = Start ,ymin=index+0,xmax = End,ymax = index+0.6),
            fill = "grey" , color= "grey") +
    geom_rect(data=TelCen,aes(xmin = Start ,ymin=index+0,xmax = End,ymax = index+0.6),
            fill = FiveColorSet[4],color=FiveColorSet[4]) + 
    geom_rect(data=hg19,aes(xmin = Start ,ymin=index-0.05,xmax = End,ymax = index),
          fill = "black",color="black") + 
    theme( panel.grid = element_line(colour = "white"),
          panel.background = element_rect(fill="white")
          ) +
    scale_color_manual(name="", values = c(FiveColorSet[1:3],"grey","lightgreen")) +
    xlab("") +ylab("")

print(p)

