FormatSI <- function(x, ...) {
  scale.frac <- 1000
  scale.unit <- c("k", "M", "G", "T", "P", "E", "Z", "Y")
  p <- rep(" ", length(x))
  for(i in 1:length(scale.unit)) {
    p[x >= scale.frac] <- scale.unit[i]
    x[x >= scale.frac] <- x[x >= scale.frac] / scale.frac
  }
  return(paste(format(round(x,1), trim=TRUE, scientific=FALSE, ...), p))
}

Index2Chr <- function(x){
  chrmap <-list("24" = "chr1","23" = "chr2","22" = "chr3","21" = "chr4","20" = "chr5",
             "19" = "chr6","18" = "chr7","17" = "chr8","16" = "chr9","15" = "chr10",
             "14" = "chr11","13" = "chr12","12" = "chr13","11" = "chr14","10" = "chr15",
             "9" = "chr16","8" = "chr17","7" = "chr18","6" = "chr19","5" = "chr20",
             "4" = "chr21","3" = "chr22","2" = "chrX","1" = "chrY")
  chrmap2 <-list("1" = "chr1","2" = "chr2","3" = "chr3","4" = "chr4","5" = "chr5",
                "6" = "chr6","7" = "chr7","8" = "chr8","9" = "chr9","10" = "chr10",
                "11" = "chr11","12" = "chr12","13" = "chr13","14" = "chr14","15" = "chr15",
                "16" = "chr16","17" = "chr17","18" = "chr18","19" = "chr19","20" = "chr20",
                "21" = "chr21","22" = "chr22","23" = "chrX","24" = "chrY")
  result <- c()
  for(i in 1:length(x)){
    result[i] = chrmap[[x[i]]]
  }
  return(result)
}

