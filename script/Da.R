BNX.data <- read.delim("head -2000 ~/Desktop/Work/temp.data |",header = T)
color_set1 <-c("#CDDB6A","#FA7942","#E95235","#575586","#4D8F48","#81CE86")
png("~/Desktop/Work/Mol.Dis.png",width = 800,height = 600)
library(ggplot2)
p <- ggplot(BNX.data,aes(x = Length / 100000 ))
p <- p + geom_density( alpha = .6 ,
                       fill=color_set1[3])
p <- p + theme_bw()
p <- p + theme(axis.text=element_text(size=12),
	        axis.title=element_text(size=14,face="bold"),
		plot.title = element_text(size=20,face="bold")
	)
p <- p + xlim(c(0,10))
p <- p + ggtitle("Molecule Size Distribution")
p <- p + xlab("Molecule Length(100kb)")
p <- p + ylab("Density")
p
print(p)
#ggsave(p,"~/Desktop/test.jpg")
