#!/usr/bin/env Rscript
require(plotrix)
circledata <- read.table('00032results.txt')
ylow = min(c(circledata$V2, circledata$V4, circledata$V6))
yhig = max(c(circledata$V2, circledata$V4, circledata$V6))
xlow = min(c(circledata$V1, circledata$V3, circledata$V5))
#xhig = max(c(circledata$V1, circledata$V3, circledata$V5))
xhig = xlow + yhig-ylow
plot(circledata$V7, circledata$V8, xlim=c(xlow, xhig), ylim=c(ylow, yhig),  pch=19, cex=0.4, col="blue")
par(new=T)
plot(circledata$V1, circledata$V2, xlim=c(xlow, xhig), ylim=c(ylow, yhig),  pch=19, cex=0.5, col="red")
par(new=T)
plot(circledata$V3, circledata$V4, xlim=c(xlow, xhig), ylim=c(ylow, yhig),  pch=19, cex=0.5, col="red")
par(new=T)
plot(circledata$V5, circledata$V6, xlim=c(xlow, xhig), ylim=c(ylow, yhig),  pch=19, cex=0.5, col="red")
par(new=T)
for (i in 1:length(circledata$V9)) {
    draw.circle(circledata$V7[i], circledata$V8[i], circledata$V9[i])
}
