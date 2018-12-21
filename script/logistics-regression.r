#mydata <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")
#xtabs(~admit + rank, data = mydata)
mydata$rank <- factor(mydata$rank)
mylogit <- glm(admit ~ gre + gpa + rank, data = mydata, family = "binomial")