#LDA

require(MASS)
require(ggplot2)
require(scales)
require(gridExtra)


dat <- read.table("SD_top_50percent.txt", header=TRUE, row.names=1)
group <- c("B","B","B","B","B","D20","D20","D20","D20","D20","D20","D20","D20","D20","D20","D60","D60","D60","D60","D60","D60","D60","D60","D60","D60")
dat <- dat[,1:25]
sdat <- t(dat[1:11653,])
sdatg <- data.frame(sdat,group=group)
lda <- lda(group ~ ., sdatg, prior=c(5,10,10)/25)
prop.lda = lda$svd^2/sum(lda$svd^2)
prop.lda

plda <- predict(object=lda,newdata=sdatg)
dataset=data.frame(group=sdatg[,"group"],lda=plda$x)
p1 <- ggplot(dataset) + geom_point(aes(lda.LD1, lda.LD2, colour = group, shape = group), size = 2.5) +
  labs(x = paste("LD1 (", percent(prop.lda[1]), ")", sep=""),
       y = paste("LD2 (", percent(prop.lda[2]), ")", sep=""))


