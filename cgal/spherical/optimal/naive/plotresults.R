#!/usr/bin/env Rscript
require(plotrix)

args = commandArgs(trailingOnly = TRUE)

tau = 6.283185

stereo = function(lat, lon, clat, clon) {
    phi = lat*tau/360
    lam = lon*tau/360
    phi1 = clat*tau/360
    lam0 = clon*tau/360
    r = 6371009
    k = 2*r/(1 + sin(phi1)*sin(phi)+cos(phi1)*cos(phi)*cos(lam-lam0))
    x = k*cos(phi)*sin(lam-lam0)
    y = k*(cos(phi1)*sin(phi)-sin(phi1)*cos(phi)*cos(lam-lam0))
    return(data.frame(x, y))
}

pointdata = read.table(args[1])
center = tail(pointdata, 1)
stereodata = stereo(pointdata$V1, pointdata$V2, center$V1, center$V2)
circledata = read.table(args[2])
#xmin = min(c(circledata$V1, circledata$V3, stereodata[[1]]))
#xmax = max(c(circledata$V1, circledata$V3, stereodata[[1]]))
#ymin = min(c(circledata$V2, circledata$V4, stereodata[[2]]))
#ymax = max(c(circledata$V2, circledata$V4, stereodata[[2]]))
xmin = min(c(circledata$V3, stereodata[[1]]))
xmax = max(c(circledata$V3, stereodata[[1]]))
ymin = min(c(circledata$V4, stereodata[[2]]))
ymax = max(c(circledata$V4, stereodata[[2]]))
xlow = xmin
xhig = xmax
ylow = ymin
yhig = ymax
plot(stereodata[[1]], stereodata[[2]], xlim=c(xlow, xhig), ylim=c(ylow, yhig),  pch=19, cex=0.5, col="red", asp=1)
par(new=T)
plot(circledata$V1, circledata$V2, xlim=c(xlow, xhig), ylim=c(ylow, yhig),  pch=19, cex=0.4, col="green", asp=1)
par(new=T)
plot(circledata$V3, circledata$V4, xlim=c(xlow, xhig), ylim=c(ylow, yhig),  pch=19, cex=0.4, col="blue", asp=1)
par(new=T)
for (i in 1:length(circledata$V5)) {
    draw.circle(circledata$V3[i], circledata$V4[i], circledata$V5[i])
}

summary(stereodata)
summary(circledata)
