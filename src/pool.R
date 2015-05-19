### Pool polls one party at a time, with logit-scale random-walk to respect boundaries
### From README:
###	* Input (1): priorout.RData
###	* Input (2): data contained in working/polling
###	* output (1): .RData file containing an 
###   object named poolout which is a list of four matrices, named 
###   mean, hi, lo, and var, where each matrix has 
###   dimension nParties (R) * nDays (C)
###	* output (2): .RData file containing an 
###   object named houseout which is a list of four matrices, named 
###   mean, hi, lo, and var, where each matrix has 
###   dimension nParties * nHouses
###	* output (3): .RData file containing an 
###   object named poolposterior which is 
###   an array of nParties * nDays * nIter

### Libraries
library(rjags)
library(MASS)
### Config
nOut <- 2000
nThin <- 5
nIter <- nOut * nThin
nBurn <- nIter / 2

### Get priors for end of walk
load(paste0("working/historical/",forecastYear,"priorout.RData"))

### Get priors for relationship of polling and outcome
load(paste0("working/historical/", forecastYear,"calibrout.RData"))

### Get polling data for current election
load(paste0("working/polling/",forecastYear,"pollcountdata.RData"))

### Jags script

model <- '
	model
	{
		for (i in 1:nPolls) {
			for (j in 1:nParties) {
			    mu[i,j] <- exp(house[org[i],j] + alphastar[date[i],j])
			}
			y[i,1:nParties] ~ dmulti(mu[i,1:nParties], size[i])
		}

		### Walk starts one year before the election
		alphastar[1,1] <- 0
		for (j in 2:nParties) {
			alphastar[1,j] ~ dunif(-10,10)
		}	

		### Walk goes forwards in time
		for (i in 2:nPeriods) {
			alphastar[i,1] <- 0
			for (j in 2:nParties) {
				alphastar[i,j] ~ dnorm(alphastar[i-1,j], tau)
			}
		}	
		
		## Construct pooled poll estimate	
		for (i in 1:nPeriods) {
			for (j in 1:nParties) {
				ealphastar[i,j] <- exp(alphastar[i,j])
				alpha[i,j] <- ealphastar[i,j]/sum(ealphastar[i,])
			}
		}		
		
		### Use calibration model to predict election results at each day
		### Uses beta distributed errors to avoid illegal values
		for (i in 1:nPeriods) {
			for (j in 1:nParties) {
				election[i,j] <- election.unnorm[i,j]/sum(election.unnorm[i,1:nParties])
				election.unnorm[i,j] ~ dbeta(election.alpha[i,j],election.beta[i,j])
				
				election.alpha[i,j] <- (election.fitted[i,j])*(election.fitted[i,j] - election.fitted[i,j]*election.fitted[i,j] - (field[j] * priorsigmasq[i]*(election.fitted[i,j]*(1-election.fitted[i,j]))/0.21))/((field[j] * priorsigmasq[i]*(election.fitted[i,j]*(1-election.fitted[i,j]))/0.21))
				election.beta[i,j] <- (1-election.fitted[i,j])*(election.fitted[i,j] - election.fitted[i,j]*election.fitted[i,j] - (field[j] * priorsigmasq[i]*(election.fitted[i,j]*(1-election.fitted[i,j]))/0.21))/((field[j] * priorsigmasq[i]*(election.fitted[i,j]*(1-election.fitted[i,j]))/0.21))
						
				election.fitted[i,j] <- priorbeta[i]*alpha[i,j]+(1-priorbeta[i])*lag[j]
			}
			priorsigmasq[i] <- pow(priorsigma[i],2)
		}
			
		### Step size
		tau <- pow(sigma,-2)
		sigma ~ dunif(0,0.2)
			
		
		### Centered house effects
		for (i in 1:nOrgs) {
			for (j in 1:nParties) {
				house.star[i,j] ~ dunif(-0.5,0.5)
				house[i,j] <- house.star[i,j] - inprod(house.star[,j],orgweights)
			}
		}
		
	}
'

# what are the parties? 
party.names <- c("Conservatives","Labour","Liberal Democrats",
	"SNP","Plaid Cymru","Greens","UKIP","Other")
	
# what fraction of seats does each party contest?	
field <- c(632,632,632,59,40,632,632,632)/632	

size <- rowSums(polls[,party.names])

# what absolute weights do we put on each house in constructing the polling average?

if (forecastYear == 2010){
	
	HouseWeights <- rep(1,length(levels(polls$Company)))
	HouseWeights <- HouseWeights/sum(HouseWeights)
	
} else {
		
	HouseWeights <- rep(1,length(levels(polls$Company)))
	
	# zero out weight on old methodologies and flaky pollsters
	HouseWeights <- HouseWeights*(!is.element(levels(polls$Company),c("YouGov","Ashcroft","Populus","Panelbase")))
		
	HouseWeights <- HouseWeights/sum(HouseWeights)	
		
}		
		


### Prep the data for Jags
###

forJags <- list(y = polls[,party.names],
	nParties = 8,
	nPeriods = 365,
	nPolls = nrow(polls),
	nOrgs = length(levels(polls$Company)),
	date = polls$date,
	org = as.numeric(polls$Company),
	orgweights = HouseWeights,
	size = size,
	lag = priorout$lag,
	priorbeta = calibrout$coef,
	priorsigma = calibrout$sigma,
	field = field)

jags.mod <- jags.model(file = textConnection(model),
	data=forJags,
	n.chains=1,
	n.adapt=100)

jags.out <- coda.samples(jags.mod,n.burnin=nBurn,
	n.iter=nIter,thin=nThin,
	variable.names=c("alpha","house","sigma","election"))

### Create objects

poolout <- list(mean = matrix(NA,nrow = 8,ncol = 365),
	hi = matrix(NA,nrow = 8,ncol = 365),
	lo = matrix(NA,nrow = 8,ncol = 365),
	var = matrix(NA,nrow = 8,ncol = 365))
	
electionout <- list(mean = matrix(NA,nrow = 8,ncol = 1),
	hi = matrix(NA,nrow = 8,ncol = 365),
	lo = matrix(NA,nrow = 8,ncol = 365),
	var = matrix(NA,nrow = 8,ncol = 365))	

houseout <- list(mean = matrix(NA,nrow = 8,ncol = forJags$nOrgs),
	hi = matrix(NA,nrow = 8,ncol = forJags$nOrgs),
	lo = matrix(NA,nrow = 8,ncol = forJags$nOrgs),
	var = matrix(NA,nrow = 8,ncol = forJags$nOrgs))

poolposterior <- array(NA,dim = c(8,365,nOut))
electionposterior <- array(NA,dim = c(8,365,nOut))

### Populate poolout

# First get the mean

alpha.pos <- grep("alpha",colnames(jags.out[[1]]))
alpha.bar <- colMeans(jags.out[[1]][,alpha.pos])

poolout$mean <- matrix(alpha.bar,
	nrow = 8, ncol = 365,
	byrow = TRUE)

# Now get the remaining quantities

alpha.var <- apply(jags.out[[1]][,alpha.pos],2,var)
alpha.lo <- apply(jags.out[[1]][,alpha.pos],2,quantile,probs=.025)
alpha.hi <- apply(jags.out[[1]][,alpha.pos],2,quantile,probs=.975)

poolout$var <- matrix(alpha.var,
	nrow = 8, ncol = 365,
	byrow = TRUE)
poolout$hi <- matrix(alpha.hi,
	nrow = 8, ncol = 365,
	byrow = TRUE)
poolout$lo <- matrix(alpha.lo,
	nrow = 8, ncol = 365,
	byrow = TRUE)
	
### Populate electionout	

election.pos <- grep("election",colnames(jags.out[[1]]))
election.bar <- colMeans(jags.out[[1]][, election.pos])	

electionout$mean <- matrix(election.bar,
	nrow = 8, ncol = 365,
	byrow = TRUE)
	
# Now get the remaining quantities	

election.var <- apply(jags.out[[1]][,election.pos],2,var)
election.hi <- apply(jags.out[[1]][,election.pos],2,quantile,probs=.975)
election.lo <- apply(jags.out[[1]][,election.pos],2,quantile,probs=.025)

electionout$var <- matrix(election.var,
	nrow = 8, ncol = 365,
	byrow = TRUE)
electionout$hi <- matrix(election.hi,
	nrow = 8, ncol = 365,
	byrow = TRUE)
electionout$lo <- matrix(election.lo,
	nrow = 8, ncol = 365,
	byrow = TRUE)

### Convert house effects to vote scale and populate houseout

house.pos <- grep("house",colnames(jags.out[[1]]))
house.posterior <- array(jags.out[[1]][,house.pos],dim=c(nOut,forJags$nOrg,forJags$nParties))

today <- poolout$mean[,max(polls$date)]
for (sim in 1:nOut){
	for (i in 1:ncol(houseout$mean)){
		tmp <- log(today) + house.posterior[sim,i,]
		tmp <- exp(tmp)/sum(exp(tmp))
		house.posterior[sim,i,] <- tmp - today
	}
}

house.bar <- apply(house.posterior,c(2,3),mean)
house.var <- apply(house.posterior,c(2,3),var)
house.hi <- apply(house.posterior,c(2,3),quantile,probs=.975)
house.lo <- apply(house.posterior,c(2,3),quantile,probs=.025)

houseout$mean <- matrix(house.bar,
	nrow = 8, ncol = forJags$nOrgs,
	byrow = TRUE,
	dimnames=list(party.names ,levels(polls$Company)))
houseout$var <- matrix(house.var,
	nrow = 8, ncol = forJags$nOrgs,
	byrow = TRUE,
	dimnames=list(party.names ,levels(polls$Company)))
houseout$hi <- matrix(house.hi,
	nrow = 8, ncol = forJags$nOrgs,
	byrow = TRUE,
	dimnames=list(party.names ,levels(polls$Company)))
houseout$lo <- matrix(house.lo,
	nrow = 8, ncol = forJags$nOrgs,
	byrow = TRUE,
	dimnames=list(party.names ,levels(polls$Company)))
	
### Calculate pollster "forecast" ranges

house.posterior <- array(jags.out[[1]][,house.pos],dim=c(nOut,forJags$nOrg,forJags$nParties))

averageSizeByPollster <- aggregate(forJags$size,list(polls$Company),median)
averageSizeByPollster <- round(averageSizeByPollster[,2]/10)*10
names(averageSizeByPollster) <- levels(polls$Company)
todaysPoolPosterior <- jags.out[[1]][,grep(paste0("alpha[",max(polls$date),","),colnames(jags.out[[1]]),fixed=TRUE)]	
pollsterforecasts <- array(NA,dim=c(forJags$nParties,3,forJags$nOrgs))	
for (k in 1:dim(house.posterior)[2]){
	# apply house effects to the pool posterior
	pollsterPosterior <- apply(log(todaysPoolPosterior) + house.posterior[,k,],1,function(x) exp(x)/sum(exp(x)))
	# add multinomial noise, given average size polls
	pollSims <- apply(pollsterPosterior,2,function(x) rmultinom(1,averageSizeByPollster[k],x)/averageSizeByPollster[k])

	pollsterforecasts[,,k] <- t(apply(pollSims,1,quantile,c(0.05,0.50,0.95)))
	
	}	

dimnames(pollsterforecasts)[[1]] <- party.names
dimnames(pollsterforecasts)[[2]] <- c(0.05,0.50,0.95)
dimnames(pollsterforecasts)[[3]] <- levels(polls$Company)	

### Sort out the posterior
for (i in 1:nOut) {
	poolposterior[,,i] <- matrix(jags.out[[1]][i,alpha.pos],
		nrow = 8,
		ncol = 365,
		byrow = TRUE) 	
}

### Sort out the posterior
for (i in 1:nOut) {
	electionposterior[,,i] <- matrix(jags.out[[1]][i,election.pos],
		nrow = 8,
		ncol = 365,
		byrow = TRUE) 	
}

# sigma.bar <- colMeans(jags.out[[1]][,grep("sigma",colnames(jags.out[[1]]))])

lastpolldate <- max(forJags$date)

###	* output (1): .RData file containing an 
###   object named poolout which is a list of four matrices, named 
###   mean, hi, lo, and var, where each matrix has 
###   dimension nParties (R) * nDays (C)
###	* output (2): .RData file containing an 
###   object named electionout which is a list of four matrices, named 
###   mean, hi, lo, and var, where each matrix has 
###   dimension nParties (R) * nDays (C)
###	* output (3): .RData file containing an 
###   object named houseout which is a list of four matrices, named 
###   mean, hi, lo, and var, where each matrix has 
###   dimension nParties * nHouses
###	* output (4): .RData file containing an 
###   object named poolposterior which is 
###   an array of nParties * nDays * nIter
###	* output (5): .RData file containing an 
###   object named electionposterior which is 
###   an array of nParties * nDays * nIter

###

company.names <- levels(polls$Company)
save(poolout,file = paste0("working/polling/",forecastYear,"poolout.RData"))
save(electionout,lastpolldate,file = paste0("working/polling/",forecastYear,"electionout.RData"))
save(houseout,company.names,file = paste0("working/polling/",forecastYear,"houseout.RData"))
save(pollsterforecasts, averageSizeByPollster,file = paste0("working/polling/",forecastYear,"pollsterforecasts.RData"))
save(poolposterior,file = paste0("working/polling/",forecastYear,"poolposterior.RData"))
save(electionposterior,file = paste0("working/polling/",forecastYear,"electionposterior.RData"))

