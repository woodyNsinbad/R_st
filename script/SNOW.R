df <- data.frame(val=1:10, ind=c(rep(2, 5), rep(3, 5))) ## make a example dataset ##

library(doSNOW)
registerDoSNOW(makeCluster(2, type = "SOCK"))
system.time(print(ddply(df, .(ind), function(x) { Sys.sleep(2); sum(x) }, .parallel=FALSE)))
system.time(print(ddply(df, .(ind), function(x) { Sys.sleep(2); sum(x) }, .parallel=TRUE)))

FiveColorSet <- c("#F25F4F","#FFB12E","#1091A4","#9BBD73","#295570")
