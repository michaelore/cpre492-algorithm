require(plotrix)
circledata <- read.table('00128results.txt')
plot(circledata$V7, circledata$V8, cex=0.2)
par(new=T)
plot(circledata$V1, circledata$V2, pch=19, cex=0.5, col="red")
plot(circledata$V3, circledata$V4, pch=19, cex=0.5, col="red")
plot(circledata$V5, circledata$V6, pch=19, cex=0.5, col="red")
for (i in 1:length(circledata$V9)) {
    draw.circle(circledata$V7[i], circledata$V8[i], circledata$V9[i])
}
