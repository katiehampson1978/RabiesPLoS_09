#MODEL IN CONTINUOUS TIME # revised by KH 8 May 2019
source("R/IBMfunctions.R")
# FUNCTIONS: assign number of bites and if they get rabies (biteIDs)
#	    assign serial intervals that get rabies
#	    run model given a defined no of events, seeded cases, a population and an enddate


#PARAMETERS
SIshape = 1.372851; SIrate = 0.05797735 # Gamma Serial interval
plot(hist(rgamma(5000, shape=SIshape, rate=SIrate), breaks = 0:200)) # plot SI
mean(rgamma(5000, shape=SIshape, rate=SIrate)) # check mean is ~23 days
SI=dgamma(0:500, shape=SIshape, rate=SIrate) # actually get new estimates from ParamEst.R

prabies=0.4746377 			# Prob(rabies|bite)

bitemean = 2.24649; bitevar = 6.298767;	#biting distribution (neg binomial)
bitesize = (bitemean^2)/(bitevar-bitemean)

# ALGORITHM (infections plus demography):
# Seed model with rabid dogs at dates in a queue
# determine earliest event (litter, death, infection)
# At earliest case assign rules for rabid animal
	# (bites, p(rabies), serial intervals) and put new cases in queue
# at earliest death remove animal and its events
# at earliest litter create new individuals with all their events
# move to next event in queue
# IN THIS MODEL THERE ARE NO WASTED BITES (SUSCEPTIBLES THAT DIE ARE JUST REMOVED..FREQ DEPENDENT)

#INITIALISE MODEL
population = makepopulation(500, 0)
seeds = makeseeds(population, 2) # Use 2 seeds to examine size of outbreaks....
Nevents = 2000 # No. events to monitor
end = 500 # end date to stop simulation
nruns = 10
nmonths = length(seq(0, end+30, 30.5))	#No. epidemics to run and no. months
eruns = pruns = matrix(ncol = nmonths, nrow = nruns) # matrices to store info on epidemics
# test = runpopmodel(Nevents, seeds, population, end)

#__________________________________________________________________________________
#Now run simulation and calc size of outbreaks for diff vacc coverages
#or run simulation to see probability of outbreaks at different vacc coverages
	# using Negative binomial versus poission biters!!!

# I'VE NOT TIMED THIS BUT TAKES A MIN OR AT EACH COVERAGE!
nruns = 5000	# Set up for lots of epidemics
eruns = pruns = matrix(ncol=nmonths, nrow=nruns)

vaccrange=seq(0,.99,0.01) 			#set up range of vaccination coverages to loop across
cases=20					# Set up matrix of Probability of outbreak from 1-20 cases
pOB = matrix(nrow=cases, ncol=length(vaccrange))
Osizes = matrix(nrow=length(vaccrange), ncol=nruns)


for (i in 1:length(vaccrange)){
	for (j in 1:nruns){
		population=makepopulation(500, vaccrange[i])
		seeds=makeseeds(population, 2)
#		seeds=makeseeds(population, 1) # change to one seed to get the probability of an incursion!
		erun = runpopmodel(Nevents, seeds, population, end)
#		erun = runpopmodelP(Nevents, seeds, population, end) #also run with poisson biting!
		pruns[j,] = erun[,1]	# population size
		eruns[j,] = erun[,2]	# epidemic size
		print(c(i,j))
	}

# If looking at the probability of an outbreak put in next line - otherwise lose it!
#	for (k in 1:cases){pOB[k,i]=length(which(apply(eruns,1,sum)>k))/nruns}
	Osizes[i,] = apply(eruns,1,sum)
}

write.csv(Osizes, file="output/outbreaksizes.csv", row.names = FALSE)
# pOB=pOB[,-c(98, 99, 100)] #get rid of NAs
#write.csv(pOB, file="pOutbreak2.csv", row.names = FALSE)#NEGATIVE BINOMIAL
#write.csv(pOB, file="pOutbreakPoisson.csv", row.names = FALSE)#POISSON


#__________________________________________________________________________________
#Run epidemic repeatedly - with no vaccination coverage - store runs in matrix
# Then run curve fit to test the effect of sampling
Nevents = 2000; end=500		#No. events to monitor; endate to stop simulation
nruns=1000; nmonths=length(seq(0, end+30, 30.5))	#No. epidemics to run and no. months
eruns=pruns=matrix(ncol=nmonths, nrow=nruns)

for (i in 1:nruns){
	population=makepopulation(500, 0)
	seeds=makeseeds(population, 2)
#	erun = runpopmodel(Nevents, seeds, population, end)	#Negative binomial biting
	erun = runpopmodelP(Nevents, seeds, population, end)	#Poisson biting
	pruns[i,] = erun[,1]	#population size
	eruns[i,] = erun[,2]	#epidemic size
	print(i)
}

write.csv(eruns, file="erunsP.csv", row.names = FALSE)
#write.csv(eruns, file="eruns.csv", row.names = F)

#plot final outbreak size
Osizes=apply(eruns, 1, sum)
Obins=hist(Osizes, breaks=c(seq(0,50,5), 1000), plot=F)$counts
barplot(Obins)

#plot the probability of an outbreak of different sizes
poutbreak=numeric(length=100)
for (i in 3:100){poutbreak[i]=length(which(Osizes>i))/nruns}
plot(3:100, poutbreak[3:100], type="l", xlab="Outbreak size", ylab="Probability", ylim=c(0,1))
#poutbreak=length(which(Osizes>2))/1000 #40% chance of an outbreak of less than 10

plot(1:18, eruns[1,], type="l", ylim=c(0,100))
for (i in 1:1000){lines(1:18, eruns[i,])}




# first see if epidemic occurred (>30 cases) and if R0 can be fit....
#Function to estimate lambda
SIshape=1.43920570; SIrate=0.05894706
SI=dgamma(0:500, shape=SIshape, rate=SIrate)

library("MASS")
source("Model\\R0fit_functions.R")

R0s=fit=aics1=aics2=numeric(length=1000)
for(i in 1:1000){
	peakmonth=which(eruns[i,]==max(eruns[i,]))[1]
	timeseq=(1:peakmonth)*30.5
	fit[i]=ifelse(length(timeseq)>3, 1, 0)
	R0s[i]=ifelse(fit[i]==1 & eruns[i,peakmonth]>20, Roest(eruns[i,1:peakmonth], timeseq, SI)[3], 0)
	aics1[i]=ifelse(fit[i]==1 & eruns[i,peakmonth]>20, Roest(eruns[i,1:peakmonth], timeseq, SI)[4], 0)
	aics2[i]=ifelse(fit[i]==1 & eruns[i,peakmonth]>20, Roest2(eruns[i,1:peakmonth], timeseq, SI)[4], 0)
}

runs=which(R0s>0); sampleR0s=matrix(nrow=length(runs), ncol=10)
for (i in 1:length(runs)){
	cases=rep(1:18, eruns[runs[i],])
	for (j in 1:10){
		samplecases=sample(cases, round(length(cases)*.8))
		sampleseq=hist(samplecases, 1:18, plot=F)$counts
		peakmonth=which(sampleseq==max(sampleseq))[1]
		timeseq=(1:peakmonth)*30.5
		sampleR0s[i,j]=ifelse(length(timeseq)>3, Roest(sampleseq[1:peakmonth], timeseq, SI)[3], NA)
	}
}

#This shows that there is negligible potential for bias induced by sampling
mean80= mean(sampleR0s, na.rm=T)
hist(mean80-R0s[runs], breaks=seq(-0.6, 0.2, 0.01))
#but this sampling is actually somewhat biased because if you randomly observe only 80% of cases it means you
#are observing proportionally less during the larger part of the epidemic. Instead could just see what effect
#is if you select 80% of cases for every month of the epidemic.


#__________________________________________________________________________________
#now fit epidemics with negative binomial and poisson errors and check which is best
# SO NEED TO CONDITION ON EPIDEMIC OCCURRING

Nevents = 2000; end=500		#No. events to monitor; endate to stop simulation
nruns=1000; nmonths=length(seq(0, end+30, 30.5))	#No. epidemics to run and no. months
eruns=pruns=matrix(ncol=nmonths, nrow=nruns)

for (i in 1:nruns){ # first make 1000 epidemic runs - NEGATIVE BINOMIAL!!!
	erun=matrix(0, nrow=2,ncol=2)
	while(sum(erun[,2])<30){
		population=makepopulation(500, 0)
		seeds=makeseeds(population, 2)
		erun = runpopmodel(Nevents, seeds, population, end)
	}
	eruns[i,] = erun[,2]	#epidemic size
	print(i)
}
write.csv(eruns, file="conditionederuns.csv", row.names = FALSE)

for (i in 1:nruns){#first make 1000 epidemic runs - POISSON!!!
	erun=matrix(0, nrow=2,ncol=2)
	while(sum(erun[,2])<30){
		population=makepopulation(500, 0)
		seeds=makeseeds(population, 2)
		erun = runpopmodelP(Nevents, seeds, population, end)
	}
	eruns[i,] = erun[,2]	#epidemic size
	print(i)
}
write.csv(eruns, file="conditionederunsP.csv", row.names = FALSE)


R0s=fit=aics1=aics2=numeric(length=1000)
for(i in 1:1000){
	peakmonth=which(eruns[i,]==max(eruns[i,]))[1]
	timeseq=(1:peakmonth)*30.5
	fit[i]=ifelse(length(timeseq)>3, 1, 0)
	R0s[i]=ifelse(fit[i]==1 & eruns[i,peakmonth]>20, Roest(eruns[i,1:peakmonth], timeseq, SI)[3], 0)
	aics1[i]=ifelse(fit[i]==1 & eruns[i,peakmonth]>20, Roest(eruns[i,1:peakmonth], timeseq, SI)[4], 0)
	aics2[i]=ifelse(fit[i]==1 & eruns[i,peakmonth]>20, Roest2(eruns[i,1:peakmonth], timeseq, SI)[4], 0)
}
#work out how many out of 1000 were best described by negative binomial errors
which(aics1<aics2)
sum(aics1-aics2) #which performed better overall?


#########################################################################
# now run epidemics with negative binomial and poisson distributions and cf 95%CIs
erunsP=read.csv("conditionederunsP.csv")
eruns=read.csv("conditionederuns.csv")

R0s=R0sP=fit=fitp=numeric(length=1000)
for(i in 1:1000){
	peakmonth=which(eruns[i,]==max(eruns[i,]))[1]
	timeseq=(1:peakmonth)*30.5
	fit[i]=ifelse(length(timeseq)>3, 1, 0)
	R0s[i]=ifelse(fit[i]==1 & eruns[i,peakmonth]>20, Roest(eruns[i,1:peakmonth], timeseq, SI)[3], 0)
}

for(i in 1:1000){
	peakmonth=which(erunsP[i,]==max(erunsP[i,]))[1]
	timeseq=(1:peakmonth)*30.5
	fit[i]=ifelse(length(timeseq)>3, 1, 0)
	R0sP[i]=ifelse(fit[i]==1 & erunsP[i,peakmonth]>20, Roest(erunsP[i,1:peakmonth], timeseq, SI)[3], 0)
}

for (i in 1:nruns){
	erun= runpopmodel(Nevents, seeds, population, end)
	pruns[i,] = erun[,1]	#population size
	eruns[i,] = erun[,2]	#epidemic size
	print(i)
}

#OUTPUT RESULTS
par(mfrow=c(2,1))
plot(1:nmonths, eruns[1,], type="l", ylim=c(0,100))
for (i in 1:nruns) {lines(1:nmonths, eruns[i,])}
plot(1:nmonths, pruns[1,], type="l", ylim=c(0,2000))
for (i in 1:nruns) {lines(1:nmonths, pruns[i,])}

#Plot the range of epidemics observed...
epidemic=incidence=vector(mode="numeric", nrow(eruns))
plot(1:nmonths, rep(0, nmonths), type="l", ylim=c(0,100), col="white") #Set up frame
for (i in 1:nrow(eruns)) {
	ncases=sum(eruns[i,]) 			#for each epidemic sum the number of cases
	if (sum(eruns[i, 18:nmonths])>0) { 	#if there were still cases after 18 months:
		points(1:nmonths, eruns[i,], col="gray", pch=20)#plot the monthly no. cases
		epidemic[i]=1 			#assign there to be an epidemic
	}#pc of cases per head of the population every year (12 mo)
	incidence[i]=((ncases*12/ncol(eruns))*100)/mean(t(pruns[i,]))
}

pepidemic=sum(epidemic)/nrow(eruns)
lines(1:nmonths, mean(eruns[epidemic==1, ]), lwd=2)
lines(1:nmonths, apply(eruns[epidemic==1,], 2, quantile, prob=0.25))
lines(1:nmonths, apply(eruns[epidemic==1,], 2, quantile, prob=0.75))

#If a dog gets rabies, it bites on its first day of symptoms, then disappears for rest of time
#If a dog infects others and they develop symptoms after we stop keeping track, ignore



