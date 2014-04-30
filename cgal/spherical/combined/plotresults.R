#!/usr/bin/env Rscript
require(plotrix)

x11()

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

stereo_rad = function(lat, lon, clat, clon, crad) {
    r = 6371009
    d_lat = crad/r*360/tau
    center = stereo(lat, lon, clat, clon)
    boundary = stereo(lat+d_lat, lon, clat, clon)
    return(dist(rbind(center, boundary)))
}

pointdata = read.table(args[1])
center = tail(pointdata, 1)
stereopointdata = stereo(pointdata$V1, pointdata$V2, center$V1, center$V2)
circledata = read.table(args[2])
stereocircledata = stereo(circledata$V1, circledata$V2, center$V1, center$V2)
k1_rads = rep(0, length(circledata$V1))
#err_rads = rep(0, length(circledata$V1))
for (i in 1:length(circledata$V1)) {
    k1_rads[i] = stereo_rad(circledata$V1[i], circledata$V2[i], center$V1, center$V2, circledata$V3[i])
    #err_rads[i] = stereo_rad(circledata$V1[i], circledata$V2[i], center$V1, center$V2, circledata$V4[i])
}
xmin = min(c(stereopointdata[[1]]))
xmax = max(c(stereopointdata[[1]]))
ymin = min(c(stereopointdata[[2]]))
ymax = max(c(stereopointdata[[2]]))
xlow = xmin
xhig = xmax
ylow = ymin
yhig = ymax
plot(stereopointdata[[1]], stereopointdata[[2]], xlim=c(xlow, xhig), ylim=c(ylow, yhig),  pch=19, cex=0.2, col="red", asp=1)
par(new=T)
plot(stereocircledata[[1]], stereocircledata[[2]], xlim=c(xlow, xhig), ylim=c(ylow, yhig),  pch=19, cex=0.7, col="blue", asp=1)
par(new=T)
plot(stereocircledata[[1]][1], stereocircledata[[2]][1], xlim=c(xlow, xhig), ylim=c(ylow, yhig),  pch=19, cex=0.8, col="green", asp=1)
par(new=T)

#Sys.sleep(10)

#for (i in 1:length(stereocircledata[[1]])) {
    i = 1
    x = stereocircledata[[1]][i]
    y = stereocircledata[[2]][i]
    if (xlow <= x & x <= xhig & ylow <= y & y <= yhig) {
        #Sys.sleep(1)
        draw.circle(x, y, k1_rads[i])
        #draw.circle(stereocircledata[[1]][i], stereocircledata[[2]][i], err_rads[i])
    }
#}

summary(stereopointdata)
summary(stereocircledata)
#summary(k1_rads)
#summary(err_rads)

Sys.sleep(999999)
