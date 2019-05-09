# plot 4A from Hampson et al. 2009
# https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1000053

# Outbreak size data from Serengeti and Ngorongoro
outbreaks = read.csv("data/Serengeti_outbreaks_2002_2007.csv")
outbreaks2 = read.csv("data/Ngorongoro_outbreaks_2002_2007.csv")

# Outbreak size simulation data
osizes = read.csv("outputs/outbreaksizes.csv", header=T) # use 10,000 runs for smoothing
# osizes=cbind(osizes2, osizes3, osizes4, osizes5) # RAN MULTIPLE TIMES
osizes = osizes[-100,] # Remove the runs that just didn't take off
meansize = apply(osizes, 1, mean)
size97.5 = apply(osizes, 1, quantile, 0.975, na.rm=TRUE)
size95 = apply(osizes, 1, quantile, 0.95, na.rm=TRUE)
size90 = apply(osizes, 1, quantile, 0.90, na.rm=TRUE)
size75 = apply(osizes, 1, quantile, 0.75, na.rm=TRUE)
size50 = apply(osizes, 1, quantile, 0.5, na.rm=TRUE)
size25 = apply(osizes, 1, quantile, 0.25, na.rm=TRUE)
size0.25 = apply(osizes, 1, quantile, 0.025, na.rm=TRUE)

maxcases = 40
covmax = 99
runs = ncol(osizes)
prop = prop2 = pds = matrix(nrow=maxcases, ncol=covmax)
for (i in 1:maxcases){
  prop2[i,] = apply (osizes, 1, function(x) length(which(x<i))/runs)
  print(i)
}

# Plot Fig 4A:
postscript(file = "figs/vacc_image.eps", width = 8, height = 4)
par(mfrow = c(1,1), cex = 0.7, lwd = 0.4, tck = -0.005, mgp = c(1,0.2,0),
    plt = c(0.1, 0.6, 0.1, 0.97), yaxs = "r", xaxs = "i")

# Background shading
image(x=0:98, y=2:40, t(prop2)[1:99,2:40],
      col=gray(0:100/100)[21:100], bty="n",
      xlab="Population-level immunity from vaccination (%)", ylab="Outbreak size",
      ylim=c(1.5,35))
axis(1,lwd=0.5); axis(2,lwd=0.5)
# THINK I SHORTENED THE FIG TO GO UP TO 90% COVERAGE! But have not fixed the aesthetics here yet!

# Contours
contour(x = 0:98, y = 2:35, t(prop2)[1:99,2:35], levels = c(0.50,0.75,0.90,0.95,0.99),
        add = TRUE, method="edge", labcex = 0.4, lwd = 0.5)

# Data on outbreaks in Serengeti and Ngorongoro
points(outbreaks$coverage*100, jitter(outbreaks$ncases), pch = 20, col="blue") # Serengeti
points(outbreaks2$coverage*100, jitter(outbreaks2$ncases), pch = 20, col="red") # Ngorongoro
dev.off()

########################################################################

# plot an example village and the timing of outbreaks e.g. Mosongo
# I've not cleaned this up yet because the raw data are not de-identified
# 9/May/19 - I'll extract de-identified data and work through the rest of this over the next few days!

NYpop = calcsus2(VPDsd$Village, dogs, vday, vaccnos, VPDsd$Vill[21])
NYcov = 100*(NYpop[,3]/apply(NYpop[,2:3], 1, sum))
NYcases = rabid$Sstart[rabid$Village=="Mosongo"]
NYmonthly = hist(NYcases, breaks=0:217, plot=FALSE)$counts[145:204]
xlabels = rep("",11)
xlabels[seq(2,10,2)] = 2002:2006

par(new=T, yaxs="i", xaxs="i", plt=c(0.38, 0.58, 0.69, 0.96), cex=0.45, mgp=c(0.7,0.1,0))
plot(NYpop[,1], NYcov,
     type="l", axes=FALSE, ylim=c(0,60),
     ylab="Population-level immunity \n from vaccination (%)", xlab="")
axis(1, labels=xlabels, at=1.2*seq(0,5,0.5), lwd=0.5)
#mtext(xlabels, 1, line=0, cex=0.45, at=1.2*seq(0,5,0.5))
axis(2, mgp=c(0.5,0.1,0), lwd=0.5)

par(new=T); barplot(NYmonthly, xlim=c(0,60*1.2), ylim=c(0,8), axes=FALSE)
axis(4, lwd=0.5)
mtext("Cases per month", 4, line=1, cex=0.5)



# Plot secondary cases/dog & P(preventing an outbreak) vs vaccination coverage:
par(new=T, yaxs="r", xaxs="i", 	plt=c(0.65, 0.95, 0.55, 0.98), cex=0.7, mgp=c(1.5, 0.2, 0))
plot(popcovSD[,1]*100, jitter(RinformSD), xlim=c(0,100), col="blue",
     ylim=c(0,max(RinformSD, RinformND, na.rm=T)), pch=20, axes=FALSE,
     xlab="", ylab="Secondary cases per rabid dog")
axis(1, lwd=0.5); axis(2, lwd=0.5)
points(popcovND[,1]*100, jitter(RinformND), pch=20, col="red")

par(new=T, yaxs="r", xaxs="i", 	plt=c(0.65, 0.95, 0.1, 0.5), mgp=c(0.8,0.05,0), cex=0.7)
input=0:98
outcome=as.numeric(pOB[5,])
fit = predict(loess(outcome~input))
plot(input, fit, type="l", ylim=c(0,0.5), xlab="Population-level immunity from vaccination (%)",
     ylab="Probability of an outbreak", axes=FALSE, lwd=2)
axis(1, lwd=0.5); axis(2, lwd=0.5)

outcome=as.numeric(pOB[10,])
fit = predict(loess(outcome~input))
lines(input, fit, type="l", ylim=c(0,1), col="orange", lwd=2)

outcome=as.numeric(pOB[20,])
fit = predict(loess(outcome~input))
lines(input, fit, type="l", ylim=c(0,1), col="forest green", lwd=2)

legend(x=50, y=0.4, legend=c("5+ cases","10+ cases", "20+ cases"), bty="n",
       col=c("black", "orange", "forest green"), lty=c(1,1,1), lwd=c(2,2,2))
#lines(0:100, rep(0.95, 101), col="gray")
# dev.off()


# postscript(file = "figs/vacc.eps", width=3, height=5, horizontal=F)
#
# par(mfrow = c(1,1), cex = 0.7, lwd = 0.4, plt = c(0.1, 0.96, 0.07, 0.49), tck = -0.005, mgp=c(1,0.2,0), yaxs = "r", xaxs = "i")
# plot(popcovSD[,1]*100, jitter(RinformSD),
#      xlim = c(0,100), ylim = c(0, max(RinformSD, RinformND, na.rm=T)),
#      xlab = "Vaccination coverage (%)", ylab = "Secondary cases per rabid dog",
#      pch = 20, axes = F)
# axis(1); axis(2)
# points(popcovND[,1]*100, jitter(RinformND),
#        pch = 20, col = "red")
#
# par(new = T, plt = c(0.1, 0.96, 0.54, 0.96), yaxs = "i")
# plot(coverage, jitter(outbreaks$ncases), pch=20, xlab="", ylab="Outbreak size",
#      xlim = c(0,100), ylim = c(0, max(outbreaks$ncases)+2), axes = F);
# axis(1); axis(2)
# points(coverage2, outbreaks2$ncases, pch = 20, col = "red")
# lines(0:98, size95, col="gray")
# lines(0:98, meansize, col="gray", lwd=2)
#
# dev.off()
