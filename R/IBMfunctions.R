#FUNCTIONS for IBM with demography
#each rabid dog makes predetermined number of bites - some go on to get rabies
#but if dogs that died of rabies are still in the pool then susceptible depletion still plays a role...

#make a dataframe with all the important popualtion demographics in it
makepopulation=function(popsize, vaccprop){
	pop=1:popsize			#IDs
	sex=rbinom(popsize, 1, 0.5)	#sexes
	death = rexp(popsize, 0.45)*365 #death times
	RFs=floor(sum(sex)*0.638)	# Use stable age distn to assign No. ads & juvs (RFs = reproductive females)
	ldRFs= runif(RFs, 0, 365) # Litter dates of RFs
	jFs=sum(sex)-RFs # juvenile females
	ldjFs= runif(jFs, 365, 640) # litter dates of juvenile females
	ld = c(ldRFs, ldjFs)		# females litter date (up to 2 years on - could be young dog)
	ldall=ifelse(sex==0, NA, ld)	# put litter dates in one long vector
	vacc=rbinom(popsize, 1, vaccprop)
	as.data.frame(cbind(pop, sex, death, ldall, vacc))
}


makeseeds = function(pop, N){ # set up the seeds of infection (making sure they are not vaccinated!)
	seeddates = runif(N, 0, 30) # seeds occur in the same month
	seedIDs = sample(pop$pop[pop$vacc==0], N, replace=FALSE) #IDs of the seeds
	seeds=as.data.frame(cbind(seedIDs, seeddates))
	names(seeds)=c("ID", "date")
	seeds
}


#function to assign the bites made by the rabid dogs and whether they result in a case
biteIDs = function(N, caseID){ 		# Assign IDs to dogs bitten by a rabid dog
	Nposs = N[-(which(N==caseID))] 		# but make sure dog cannot bite itself!
	bites = rnbinom(1, size=bitesize, mu=bitemean)

	# IT IS NOT POSS TO TAKE SAMPLE LARGER THAN POPULATION!:
	if(bites<length(Nposs)){
	  bitten = sample(Nposs, bites, replace=FALSE)
	}else{ bitten=Nposs }

	cases = runif(length(bitten), 0, 1) < prabies
	as.data.frame(cbind(bitten, cases))
}


#Assign the bites made by rabid dogs and whether they result in a case
# BUT WITH POISSON DISTRIBUTION NOT NEG. BINOM.
biteIDsP = function(N, caseID){ 		#Assign IDs to dogs bitten by a rabid dog
	Nposs=N[-(which(N==caseID))] 		#but dog cannot bite itself!
#	bites=rnbinom(1, size=bitesize, mu=bitemean)
	bites=rpois(1, lambda=bitemean)
	#SORT OUT SO ITS NOT POSS TO TAKE SAMPLE LARGER THAN POPULATION
	if(bites<length(Nposs)){bitten=sample(Nposs, bites, replace=FALSE)}else{bitten=Nposs}
	cases=runif(length(bitten), 0, 1) <prabies
	as.data.frame(cbind(bitten, cases))
}

#function to assign the serial interval for each bitten dog that goes onto get rabies
symptoms = function(caseIDlist){rgamma(length(caseIDlist), shape = SIshape, rate = SIrate)}

#INFECTION EVENT (negative binomial)
infect=function(cases, cd, cID, populationDF){
	bitten = biteIDs(populationDF$pop, cID)	# assign those bitten by the case
	getrabies = bitten$bitten[bitten$cases>0]	# and which get rabies
	vaccstat = populationDF$vacc[match(getrabies, populationDF$pop)]
	ID = getrabies[vaccstat<1]

	if (length(ID)>0) {			# IF INFECTION PRODUCES NEW CASES
		date = symptoms(ID)+cd 		# give them onset dates
		newcases = cbind(ID, date)	# make matrix of new cases
		cases = rbind(cases, newcases)	# add future cases to queue
	}
	cases = cases[-which(cases$ID == cID),]	# remove cases that occurred (and that would have occurred)
	cases					# report all the cases
}

#INFECTION EVENT - but POISSON BITING!
infectP = function(cases, cd, cID, populationDF){
	bitten = biteIDsP(populationDF$pop, cID)	# assign those bitten by the case
	getrabies = bitten$bitten[bitten$cases>0]	# and which get rabies
	vaccstat = populationDF$vacc[match(getrabies, populationDF$pop)]
	ID = getrabies[vaccstat<1]

	if (length(ID)>0) {			#IF INFECTION PRODUCES NEW CASES
		date = symptoms(ID)+cd 		#give them onset dates
		newcases = cbind(ID, date)	#make matrix of new cases
		cases = rbind(cases, newcases)	#add future cases to queue
	}
	cases = cases[-which(cases$ID == cID),]	#remove cases that occurred (and that would have occurred)
	cases					# report all the cases
}


#LITTER EVENT
givebirth = function(pmax, l1){
	litter = rbinom(1, 8, 0.25)
	if(litter > 0){			# IF LITTER PRODUCES OFFSPRING:
		pop = (pmax+1):(pmax+litter)		# give newborn IDs
		sex = rbinom(litter, 1, 0.5)		#...sexes
		death = l1+(rexp(litter, 0.45)*365)	#...death dates
		ld = l1+runif(sum(sex),230,640)		#...date of 1st litter (if female!)
		ldall = ifelse(sex == 0, NA, ld)		# GIVE ALL newborns 1st litter date!
		vacc = rep(0, length(litter))
		nrows = cbind(pop, sex, death, ldall, vacc) # report litter with details
	} else {
	  nrows = matrix(ncol = 4, nrow = 0)
	  }
	nrows
}


# function to run model with demography and infections
# Need to enter more events than you expect before the enddate
# If you run out of cases the loop stops....

runpopmodel=function(n.events, caseDF, populationDF, enddate) {
	cdates = cIDs = vector(mode="numeric")	#SET EVERYTHING
	event = edates = 0
	pmax = max(populationDF$pop)
	psize= nrow(populationDF)
	cases = caseDF

	for (i in 1:n.events){	#FOR A SET NO. OF ORDERED EVENTS
#		if(psize[length(psize)]<2){break}			# if population size gets really small break
		cd = ifelse(nrow(cases)<1, enddate+40, min(cases$date)) 	#select next case date
		if (cd > enddate){break}				# or time so loop will break
		l1 = min(populationDF$ldall, na.rm=T)		#next litter date
		d1 = min(populationDF$death) 				# next death date
 		edate = min(d1, cd, l1)					#SELECT NEXT EVENT
		if (edate > enddate) {break}		#STOP if next event is after ENDDATE

		if (cd==edate) {			# IF NEXT EVENT IS AN INFECTION
			cID = cases$ID[which(cases$date == cd)]	# determine ID of infection
			cases = infect(cases, cd, cID, populationDF) # generate the new list of cases and dates
			cdates = c(cdates, cd)			#list dates when cases occurred
			cIDs = c(cIDs, populationDF$pop[cID])	#list IDs of cases (just assume they die!)
			populationDF = populationDF[-which(populationDF$pop==cID),]   	#UPDATE POPULATION AND DEMOGRAPHICS
		}
#		print(c(nrow(cases),i,cID))
		if (l1 == edate) {			#IF NEXT EVENT IS A LITTER
			#SET FEMALES NEXT LITTER DATE!
			populationDF$ldall[which(populationDF$ldall==l1)]=l1+(runif(1, 270, 540))
			populationDF=rbind(populationDF, givebirth(pmax, l1))#UPDATE POPULATION AND DEMOGRAPHICS
			pmax=max(populationDF$pop)
		}

		if (d1 == edate) {			#IF NEXT EVENT IS A DEATH
			deathID = which(populationDF$death==d1)
			ndeath = populationDF$pop[deathID] #Remove rabies cases if they die naturally first!
			if (length(which(cases$ID==ndeath))>0) {cases=cases[-which(cases$ID==ndeath),]}
			populationDF=populationDF[-deathID,]	#UPDATE POPULATION AND DEMOGRAPHICS
		}
		psize = c(psize, nrow(populationDF)) #UPDATE
		edates = c(edates, edate)
		event=event+1
		if(event==n.events) {break} #END WHEN REACH LAST EVENT
	}

	#Output timeseries of cases and populations size every 30 days
	caseTS = hist(cdates, breaks=seq(0, enddate+60, 30.5), plot=FALSE)$counts
	mpopi = sapply(seq(1, enddate+30, 30.5), function (x){ max(which(edates<x))})
	mpsize = psize[mpopi]
	cbind(mpsize, caseTS)

	#################################################################
#	pophistory = cbind(edates, psize)
#	mresults=list(pophistory, caseTS)
#	mresults #OUTPUT THE NUMBER OF RABIES CASES AND THE POPULATION HISTORY
}


#################################################################
# use poisson biting not Neg binom!
runpopmodelP=function(n.events, caseDF, populationDF, enddate) {
	cdates = cIDs = vector(mode="numeric")	#SET EVERYTHING
	event = edates = 0
	pmax=max(populationDF$pop)
	psize=nrow(populationDF)
	cases=caseDF

	for (i in 1:n.events){	#FOR A SET NO. OF ORDERED EVENTS
#		if(psize[length(psize)]<2){break}			#if population size gets really small break
		cd=ifelse(nrow(cases)<1, enddate+40, min(cases$date)) 	#select next case date
		if(cd>enddate){break}				#or time so loop will break
		l1 = min(populationDF$ldall, na.rm=T)			#next litter date
		d1 = min(populationDF$death) 				#next death date
 		edate = min(d1, cd, l1)					#SELECT NEXT EVENT
		if (edate > enddate) {break}		#STOP if next event is after ENDDATE

		if (cd==edate) {			#IF NEXT EVENT IS An INFECTION
			cID=cases$ID[which(cases$date==cd)]	#determine ID of infection
			cases=infectP(cases, cd, cID, populationDF)#generate the new list of cases and dates
			cdates=c(cdates, cd)			#list dates when cases occurred
			cIDs=c(cIDs, populationDF$pop[cID])	#list IDs of cases (just assume they die!)
			populationDF=populationDF[-which(populationDF$pop==cID),]   	#UPDATE POPULATION AND DEMOGRAPHICS
		}
#		print(c(nrow(cases),i,cID))
		if (l1 == edate) {			#IF NEXT EVENT IS A LITTER
			#SET FEMALES NEXT LITTER DATE!
			populationDF$ldall[which(populationDF$ldall==l1)]=l1+(runif(1, 270, 540))
			populationDF=rbind(populationDF, givebirth(pmax, l1))#UPDATE POPULATION AND DEMOGRAPHICS
			pmax=max(populationDF$pop)
		}

		if (d1 == edate) {			#IF NEXT EVENT IS A DEATH
			deathID = which(populationDF$death==d1)
			ndeath = populationDF$pop[deathID] #Remove rabies cases if they die naturally first!
			if (length(which(cases$ID==ndeath))>0) {cases=cases[-which(cases$ID==ndeath),]}
			populationDF=populationDF[-deathID,]	#UPDATE POPULATION AND DEMOGRAPHICS
		}
		psize = c(psize, nrow(populationDF)) #UPDATE
		edates = c(edates, edate)
		event=event+1
		if(event==n.events) {break} #END WHEN REACH LAST EVENT
	}

	#Output timeseries of cases and populations size every 30 days
	caseTS = hist(cdates, breaks=seq(0, enddate+60, 30.5), plot=FALSE)$counts
	mpopi = sapply(seq(1, enddate+30, 30.5), function (x){ max(which(edates<x))})
	mpsize = psize[mpopi]
	cbind(mpsize, caseTS)
#################################################################
#	pophistory = cbind(edates, psize)
#	mresults=list(pophistory, caseTS)
#	mresults #OUTPUT THE NUMBER OF RABIES CASES AND THE POPULATION HISTORY
}


