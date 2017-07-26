library(vegan)
library(plyr)
library(dplyr)
library(cluster)

### PCA ###

rda(benz.dere2.y ~ benz.dere2)

#### heatmap ####

my.clust <- function(x) {hclust(x,method= "ward.D2")}
my.dist <- function(x) {dist(x,method="binary")}
heatmap(as.matrix(benz.dere2.infilter),
        distfun = my.dist ,
        hclustfun=my.clust ,
        col = RColorBrewer::brewer.pal("RdYlBu",n=11) ,
        scale = "col")

#### hcluster ####

benz.dere2.dist <- dist(benz.dere2.infilter,method= "binary")
benz.dere2.hcl <- hclust(benz.dere2.dist ,method = "ward")
benz.dere2.cluster  <- cutree(benz.dere2.hcl ,k = 5)   ## use cutree ## 
plot( benz.dere2.hcl ,
      fill = brewer.pal("Set1",n = 5)[kmeans(benz.dere2.y[-884],centers =5 )$cluster])

#### kmeans #####

kmeans(benz.dere2.y[-884],centers =5 )
pam(benz.dere2.y[-884], k =5 ) ## Partitioning Around Medoids more powerful than kmeans ##

#### silhouette plot ####

benz.sample.cutree <- cutree(benz.sample.cluster.T,k = 5) ## this is for hcluster 
benz.hclust.si <- silhouette(benz.sample.cutree,
                             benz.dere2.dist)
plot(benz.hclust.si,col = RColorBrewer::brewer.pal("Set1",n=8))

#### GAP ####

clusplot(cutree(benz.sample.cluster,k = 6))
