#! Rscript ##
library(Gviz)
library(Biostrings)

###############################

outputpng  <- "~/Desktop/Work/R_script/Gviz.test"
FASTAfile <- "~/Desktop/Work/R_script/Tair.chr1.fa"
GFFfile <- "~/Desktop/Work/R_script/tair9_Assembly_gaps.gff"
Mycolor <- c('#4682B4','#FF8C00','#A0522D',
	      	'#87CEEB','#6B8E23','#6A5ACD',
		'#778899','#DAA520','#B22222',
		'#FF0000','#00FF00','#0000FF',
		'#00FFFF','#FF00FF','#FFFF00')



Mycolor2 <- c("#F5B1A8","#F6DA91","#ACD9D2",
		"#B379A2","#E66591","#F47B6D",
		"#446373","#BBC686")
seqobject <- Biostrings::readDNAStringSet(FASTAfile)
axistrack <- Gviz::GenomeAxisTrack(range = seqobject[1]@ranges)

atrack <- Gviz::AnnotationTrack(range = GFFfile,
			  name="Text",
			  genome = "Tair10",
			  background.title = "brown", ### Change the Title background ###
#			  fill = "salmon"
		      )







for(i in 1:length(Mycolor2)){

	png(filename = paste(outputpng,i,".png",sep="") ,width = 1200 ,heigh = 200 )
	###############################
	#seqobject <- Biostrings::readDNAStringSet(FASTAfile)

	Gviz::displayPars(atrack) <- list(
#			  background.panel = "#FFFEDB" , ### change the plot background ###
			  Gene = Mycolor2[i],
			  Gap = Mycolor2[i+1]
	)

	#Gviz::plotTracks(list(axistrack,atrack))
	Gviz::plotTracks(atrack)
}

q()
