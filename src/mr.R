###  - MR (mr_*.R); 
### 	* input (1): data contained in working/individual
### 	* input (2): data contained in working/constituency
### 	* input (3): calibration contained in working/historical
### 	* output (1): .RData file containing an object named 
###			devout which is a list of four matrices named 
###			mean, hi, lo and var, where each matrix has 
###			dimension nParties * nConstituencies
### 	* output (2): .RData file containing an object named 
###			devposterior which is an array of
###			nParties * nConstituencies * nIter

elex.date <- as.Date("2015-05-07")
origin.date <- elex.date - 365
currentdate <- as.numeric(Sys.Date() - origin.date)
if (currentdate > 365) currentdate <- 365

### Libraries
library(rjags)
load.module("glm")
### Config
nOut <- 2000
nThin <- 5
nIter <- nOut * nThin
nBurn <- nIter / 2

### 	* input (0a): constituency-level specific to generic polling adjustment factors "GenericPollPropAdjFactor"

load(file=paste0("working/constituency/",forecastYear,"genericpollpropadjfactor.RData"))

### 	* input (0a): rmse of UNS estimate

load(file=paste0("working/constituency/",forecastYear,"RMSEUNS.RData"))

### 	* input (1): constituency-level polling data

load("working/constituency/2015constituencypollingdata.RData")

###
### 	* input (2): constituency covariate data
###

constituency <- read.csv("data/constituency/canonical_seatvars_2011census.csv",
	header=TRUE)
eudisapproval <- read.csv("data/constituency/eudisapp-estimates.csv",
	header=TRUE)	
yesvote <- read.csv("data/constituency/imputedyesvotebyconstituency.csv",
	header=TRUE)		
standingdown <- read.csv("data/constituency/standing_down.csv",
	header = TRUE)	

constituency <- merge(constituency,eudisapproval[,c("refno","eudisapp.hat")])	
constituency <- merge(constituency, standingdown[,c("refno","StandingDown","FirstTerm")])
constituency <- merge(constituency, yesvote[,c("refno","yes")],all.x=TRUE)
	
	epvote <- read.csv("data/constituency/2014_by_wmin.csv",
		header = TRUE)
	epvote <- epvote[,c("refno","Conservatives","Labour","Liberal.Democrats","SNP","Plaid.Cymru","Greens","UKIP")]
	names(epvote) <- c("refno","con14","lab14","ld14","snp14","pc14","grn14","ukip14")
	constituency$yes <- replace((constituency$yes - mean(constituency$yes,na.rm=TRUE))/sd(constituency$yes,na.rm=TRUE),is.na(constituency$yes),0)

	lauthvote <- read.csv("data/constituency/latest_elex_by_wmin_voteshare.csv",
		header = TRUE)
	lauthvote <- lauthvote[,c("refno","Conservatives","Labour","Liberal.Democrats","SNP","Plaid.Cymru","Greens","UKIP")]
	names(lauthvote) <- c("refno","conla","labla","ldla","snpla","pcla","grnla","ukipla")


constituency <- merge(constituency,epvote)
constituency <- merge(constituency,lauthvote)
	
### Merge with the full range of refnos
constituency <- constituency[order(constituency$refno),]

### Clean up past-vote share variables

	partieslag <- c("con10","lab10","ld10","snp10","pc10","grn10","ukip10")
	eplag <- c("con14","lab14","ld14","snp14","pc14","grn14","ukip14")
	lauthlag <- c("conla","labla","ldla","snpla","pcla","grnla","ukipla")

	incumbency.matrix <- as.matrix(constituency[,c("inc15con","inc15lab",
	"inc15ld","inc15snp","inc15pc",
	"inc15grn","inc15ukip","inc15oth")])

for (p in partieslag) {
	tmp <- constituency[,p]
	tmp [is.na(tmp)] <- 0
	constituency[,p] <- tmp
}
### Convert past vote shares to [0,1]
for (p in partieslag) {
	constituency[,p] <- constituency[,p] / 100
}

for (p in eplag) {
	constituency[,p] <- constituency[,p] / 100
}	

### Not needed for local authority results

### Deal with candidates standing down
	### first, get `standing down' by constituency
	standingdown <- matrix(constituency$StandingDown, 
		ncol = ncol(incumbency.matrix),
		nrow = 632,
		byrow = FALSE)
	### now multiply element wise by incumbency
	standingdown <- standingdown * incumbency.matrix 

	### Second, get sophomore by constituency
	sophomore <- matrix(constituency$FirstTerm, 
	#sophomore <- matrix(0,
		ncol = ncol(incumbency.matrix),
		nrow = 632,
		byrow = FALSE)
	### now multiply element wise by incumbency
	sophomore <- sophomore * incumbency.matrix 


##  tiny vote share for all parties in all constituencies
eps <- 0.00001

prevvote <- as.matrix(constituency[,partieslag])
rownames(prevvote) <- constituency$refno
save(prevvote,file = paste0("working/constituency/",forecastYear,"prevvote.RData"))

prevvote <- cbind(prevvote,1-rowSums(prevvote))  # add other category
prevvote <- replace(prevvote, prevvote <= eps,eps)  # fix negative and zero totals
prevvote <- prevvote/rowSums(prevvote)

# KLUDGEs in lagged votes to deal with speakers seat transitions and nigel farage's '10
	# speaker's seat
	tmp <- which(constituency$YouGovName == "Buckingham")
	prevvote[tmp,] <- c(0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.99993)
	incumbency.matrix[tmp,] <- c(0,0,0,0,0,0,0,1)

	# Bob Spink
	tmp <- which(constituency$YouGovName == "Castle Point")
	prevvote[tmp,7] <- 0.27

	# nigel-specific hack, no longer needed:
	#prevvote[tmp,7] <- prevvote[tmp,7] - 0.15
	#prevvote[tmp,] <- prevvote[tmp,] / sum(prevvote[tmp,])

	
prevwin <- 	incumbency.matrix
save(prevwin,file = paste0("working/constituency/",forecastYear,"prevwin.RData"))	

logitprevvote <- matrix(0,nrow(prevvote),ncol(prevvote))
for (j in 2:ncol(prevvote)){
	logitprevvote[,j] <- log(prevvote[,j]/prevvote[,1])
	}


epvote <- as.matrix(constituency[,eplag])	
epvote <- cbind(epvote,1-rowSums(epvote)) 
epvote <- replace(epvote, epvote < eps,eps) 
epvote <- epvote/rowSums(epvote)

logitepvote <- matrix(0,nrow(epvote),ncol(epvote))
for (j in 2:ncol(epvote)){
	logitepvote[,j] <- log(epvote[,j]/epvote[,1])
}

lauthvote <- as.matrix(constituency[,lauthlag])	
lauthvote <- cbind(lauthvote,1-rowSums(lauthvote)) 
lauthvote <- replace(lauthvote, lauthvote < eps, eps) 
lauthvote <- lauthvote/rowSums(lauthvote)

logitlauthvote <- matrix(0,nrow(lauthvote),ncol(lauthvote))
for (j in 2:ncol(lauthvote)){
	logitlauthvote[,j] <- log(lauthvote[,j]/lauthvote[,1])
}


### Scale some continuous variables
cont.vars <- c("nonwhite","log.earn","relig.christian",
	"relig.none","relig.refused","relig.other","private",
	"owns","female","married","education","socgrd","age","eudisapp.hat")
for (v in cont.vars) {
	tmp <- constituency[,v]
	tmp <- as.vector(scale(tmp))
	constituency[,v] <- tmp
}

yrefnos <- constituencypollingdata$refno
ypolldates <- round((constituencypollingdata$startdate + constituencypollingdata$enddate)/2)
ypollhouses <- constituencypollingdata$house
ygeneric <- as.numeric(constituencypollingdata$generic) + 1
ygenericadj <- rbind(rep(1,8), GenericPollPropAdjFactor)
ygenericladj <- array(c(rep(0,8), GenericPollLinearAdjFactor[,1],rep(1,8), GenericPollLinearAdjFactor[,2]),c(8,2,2))

y <- constituencypollingdata[,c("Conservatives",
	"Labour",
	"Liberal.Democrats",
	"SNP",
	"Plaid.Cymru",
	"Greens",
	"UKIP",
	"Other")]

for (i in 1:ncol(y)) {
	tmp <- y[,i]
	tmp[is.na(tmp)] <- 0
	y[,i] <- tmp
}
nRespondents <- rowSums(y)

for (i in 1:nrow(y)) {
	if (is.na(ypolldates[i])) ypolldates[i] <- ypolldates[i-1]
}


# Fudge Polling Data for Speaker's Seat
refnospeaker <- constituency$refno[which(constituency$YouGovName == "Buckingham")]
speakerseatpolls <- which(yrefnos == refnospeaker)
for (i in 1:length(speakerseatpolls)){	
	y[speakerseatpolls[i],] <- c(0,0,0,0,0,0,0,sum(y[speakerseatpolls[i],]))
	}


### Prepare to use non-NI refnos to link data
### yref = in range 1 .. 632
### yrefnos = original PA number
yref <- rep(NA,length(yrefnos))
for (i in 1:length(yrefnos)){
	yref[i] <- which(yrefnos[i] == constituency$refno)
	}
refnos <- constituency$refno


### Now create the matrix X which will hold all this stuff
X <- as.matrix(constituency[,c(cont.vars)])


### And the region dummies and interactions with log.density
constituency$region <- factor(constituency$region)
constituency$region <- relevel(constituency$region,"South East")
regdum <- model.matrix(~region-1,data = constituency)

Region <- as.numeric(constituency$region)
nRegions <- max(Region)

### 	* input (3): data contained in working/individual
### Get priors for relationship of polling and outcome
load(paste0("working/historical/", forecastYear,"calibrout.RData"))

# Set lagged vote vs poll weights
pollweights <- rep(1,8)

###
### 	* input (3): pooled poll estimates
###

load(file = paste0("working/polling/",forecastYear,"poolposterior.RData"))
pool.est <- apply(poolposterior,c(1,2),mean)

## calculate mean national poll levels for the time period of each constituency-level poll

logitgbvote <- matrix(0,nrow(y),ncol(y))
for (i in 1:nrow(y)){
	pool.tmp <- apply(matrix(pool.est[,constituencypollingdata$startdate[i]:constituencypollingdata$enddate[i]],nrow=nrow(pool.est)),1,mean)
	logitgbvote[i,] <- log(pool.tmp/pool.tmp[1])
	}

## Compute relative to today	

todaystarget <- pool.est[,currentdate]
todaystargetlogit <- log(todaystarget/todaystarget[1])	
logitgbvote <- logitgbvote - matrix(todaystargetlogit,nrow(logitgbvote),ncol(logitgbvote),byrow=TRUE)

### Write out the model

model <- '
	model {
		
		for (k in 1:nConstituencyPolls) {
			y[k,1:nParties] ~ dmulti(pi[k,1:nParties],nRespondents[k])
			for (j in 1:nParties) {
				pi[k,j] <- ygenericladj[j,ygeneric[k],1] + (pi.uncalib[k,j] / sum(pi.uncalib[k,1:nParties]))*ygenericladj[j,ygeneric[k],2]
				pi.uncalib[k,j] <- exp(logitgbvote[k,j] + house[org[k],j] + mu[yref[k],j])
			}	
		}
		
		for (i in 1:nConstituencies) {
			for (j in 1:nParties) {
				lp[i,j] <- inprod(alphaprev[j,1:nParties],logitprevvote[i,1:nParties]) + alphaep[j]*logitepvote[i,j] + alphalauth[j]*logitlauthvote[i,j] + inprod(X[i,1:nVars],beta[j,1:nVars]) + delta[j,Region[i]] + gamma*incumbency[i,j] + gamma2 * standingdown[i,j] + gamma3 * sophomore[i,j]
				mu[i,j] ~ dnorm(lp[i,j],tau[j])
				emu[i,j] <- exp(mu[i,j])
				pstar[i,j] <- emu[i,j] / sum(emu[i,1:nParties])
			}
		}
		
		for (i in 1:nConstituencies) {
			for (j in 1:nParties) {
				p[i,j] <- pu[i,j] / sum(pu[i,1:nParties])
				## (1-priorbeta) = 0 for 2015
				pu[i,j] <- pstar[i,j]*priorbeta[j] + uns[i,j]*(1-priorbeta[j])
				uns[i,j] ~ dbeta(uns.alpha[i,j],uns.beta[i,j])		
				uns.alpha[i,j] <- (prevvote[i,j])*((prevvote[i,j] - prevvote[i,j]*prevvote[i,j] - (sigmasquns*prevvote[i,j]*(1-prevvote[i,j])/0.21))/(sigmasquns*prevvote[i,j]*(1-prevvote[i,j])/0.21))
				uns.beta[i,j] <- (1-prevvote[i,j])*((prevvote[i,j] - prevvote[i,j]*prevvote[i,j] - (sigmasquns*prevvote[i,j]*(1-prevvote[i,j])/0.21))/(sigmasquns*prevvote[i,j]*(1-prevvote[i,j])/0.21))
			}	
		}

		for (j in 1:nParties) {
			beta[j,1:nVars] ~ dmnorm(b0,B0)		
			for (l in 1:nRegions) {
				delta[j,l] ~ dnorm(0,1)	
			}	
			for (jj in 1:nParties) {
				alphaprev[j,jj] ~ dnorm(alphaprevmean[j,jj], alphaprevprec)	
			}		
			alphaep[j] ~ dnorm(0,eppriorprec[j])
			alphalauth[j] ~ dnorm(0,lauthpriorprec[j])

			tau[j] <- pow(sigma[j],-2)
			sigma[j] ~ dunif(0,1)
		}		
		gamma ~ dnorm(0,1)
		gamma2 ~ dnorm(0,1)
		gamma3 ~ dnorm(0,1)
						
		### House effects for all but base house level in prepcon script
				
		for (j in 1:nParties) {
			house[1,j] <- 0
			for (i in 2:nOrgs) {			
				house[i,j] ~ dnorm(0, tauhouse)
			}
		}
		
		tauhouse <- pow(sigmahouse,-2)
		sigmahouse ~ dunif(0,0.1)
		
	}

'

### Select whether to use EP election data
### Specify tight prior for all non-UKIP parties in 2015
	eppriorprec <- rep(1e8,8)
	eppriorprec[7] <- 1

### For prior on local authority vote, take double prior on GE vote
	lauthpriorprec <- rep(32,8)

### Create the list holding the data

forJags <- list(
	y = y,
	yref = yref,
	ygeneric = ygeneric,
	ygenericladj = ygenericladj,	
	nRespondents = nRespondents,
	nConstituencies = nrow(X),
	nConstituencyPolls = nrow(y),
	nParties = ncol(logitprevvote),
	nVars = ncol(X),
	nOrgs = length(levels(ypollhouses)),
	nRegions = nRegions,	
	org = as.numeric(ypollhouses),
	X = X,
	Region = Region,
	alphaprevmean = diag(ncol(logitprevvote)),
	alphaprevprec = 16,
	logitprevvote = logitprevvote,
	logitepvote = logitepvote,	
	logitlauthvote = logitlauthvote, 
	logitgbvote = logitgbvote,
	incumbency = incumbency.matrix,
	sophomore = sophomore,
	standingdown = standingdown,
	b0 = rep(0,ncol(X)),
	B0 = diag(16,ncol(X)),
	priorbeta = pollweights,
	sigmasquns = (RMSEUNS^2)/4,  ## div by 4 for approx 4 parties, vars add, kludge
	prevvote = prevvote,
	eppriorprec = eppriorprec,
	lauthpriorprec = lauthpriorprec
)

 ### Try and create some inits
initfunc <- function() {
	beta <- matrix(NA,
		nrow = forJags$nParties, 
		ncol = ncol(X))
	for (i in 2:nrow(beta)) {
		beta[i,] <- rnorm(ncol(X),0,.25)
	}
	return(list(beta=beta))
}


jags.mod <- jags.model(file = textConnection(model),
	data=forJags,
	n.chains=1,
	# inits = initfunc,
	n.adapt=2000)


jags.out <- coda.samples(jags.mod,n.burnin=nBurn,
	n.iter=nIter,thin=nThin,
	variable.names=c("p","house","alphaprev","alphaep","beta","gamma","gamma2","gamma3","delta","sigma"))

p.pos <- grep("^p",colnames(jags.out[[1]]))

p.bar <- colMeans(jags.out[[1]][,p.pos])
p.var <- apply(jags.out[[1]][,p.pos],2,var)
p.hi <- apply(jags.out[[1]][,p.pos],2,quantile,probs=.975)
p.lo <- apply(jags.out[[1]][,p.pos],2,quantile,probs=.025)

### Create objects

devout <- list(mean = matrix(NA,nrow = 8,ncol = forJags$nConstituencies),
	hi = matrix(NA,nrow = 8,ncol = forJags$nConstituencies),
	lo = matrix(NA,nrow = 8,ncol = forJags$nConstituencies),
	var = matrix(NA,nrow = 8,ncol = forJags$nConstituencies))

devout$mean <- matrix(p.bar,
	nrow = 8, ncol = forJags$nConstituencies,
	byrow = TRUE)

devout$var <- matrix(p.var,
	nrow = 8, ncol = forJags$nConstituencies,
	byrow = TRUE)

devout$hi <- matrix(p.hi,
	nrow = 8, ncol = forJags$nConstituencies,
	byrow = TRUE)

devout$lo <- matrix(p.lo,
	nrow = 8, ncol = forJags$nConstituencies,
	byrow = TRUE)
	
houseout <- list(mean = matrix(NA,nrow = 8,ncol = forJags$nOrgs),
	hi = matrix(NA,nrow = 8,ncol = forJags$nOrgs),
	lo = matrix(NA,nrow = 8,ncol = forJags$nOrgs),
	var = matrix(NA,nrow = 8,ncol = forJags$nOrgs))
	
house.pos <- grep("house",colnames(jags.out[[1]]))
house.bar <- colMeans(jags.out[[1]][,house.pos])
house.var <- apply(jags.out[[1]][,house.pos],2,var)
house.hi <- apply(jags.out[[1]][,house.pos],2,quantile,probs=.975)
house.lo <- apply(jags.out[[1]][,house.pos],2,quantile,probs=.025)

houseout$mean <- matrix(house.bar,
	nrow = 8, ncol = forJags$nOrgs,
	byrow = TRUE)
houseout$var <- matrix(house.var,
	nrow = 8, ncol = forJags$nOrgs,
	byrow = TRUE)
houseout$hi <- matrix(house.hi,
	nrow = 8, ncol = forJags$nOrgs,
	byrow = TRUE)
houseout$lo <- matrix(house.lo,
	nrow = 8, ncol = forJags$nOrgs,
	byrow = TRUE)	
	
colnames(houseout$mean)	<- colnames(houseout$var) <- colnames(houseout$hi)	<- colnames(houseout$lo)	<-	levels(ypollhouses)

	
beta.pos <- grep("^beta",colnames(jags.out[[1]]))	
beta.mat <- matrix(colMeans(jags.out[[1]][,beta.pos]),ncol=8, byrow = TRUE)
beta.hi <- matrix(apply(jags.out[[1]][,beta.pos],2,quantile,probs=0.95),
	ncol = 8, byrow = TRUE)
beta.lo <- matrix(apply(jags.out[[1]][,beta.pos],2,quantile,probs= 0.05),
	ncol = 8, byrow = TRUE)

rownames(beta.mat) <- rownames(beta.hi) <- rownames(beta.lo) <- colnames(X)

alphaprev.pos <- grep("^alphaprev",colnames(jags.out[[1]]))	
alphaprev.est <- matrix(colMeans(jags.out[[1]][, alphaprev.pos]),8,8)

alphaep.pos <- grep("^alphaep",colnames(jags.out[[1]]))	
alphaep.est <- colMeans(jags.out[[1]][, alphaep.pos])

gamma.pos <- grep("^gamma",colnames(jags.out[[1]]))	
gamma.mat <- matrix(colMeans(jags.out[[1]][,gamma.pos]),ncol=1, byrow = TRUE)
gamma.hi <- matrix(apply(jags.out[[1]][,gamma.pos],2,quantile,probs=0.95),
	ncol = 1, byrow = TRUE)
gamma.lo <- matrix(apply(jags.out[[1]][,gamma.pos],2,quantile,probs= 0.05),
	ncol = 1, byrow = TRUE)
	
delta.pos <- grep("^delta",colnames(jags.out[[1]]))	
delta.mat <- matrix(colMeans(jags.out[[1]][,delta.pos]),nrow=8, byrow = TRUE)
delta.hi <- matrix(apply(jags.out[[1]][,delta.pos],2,quantile,probs=0.95),
	ncol = 8, byrow = TRUE)
delta.lo <- matrix(apply(jags.out[[1]][,delta.pos],2,quantile,probs= 0.05),
	ncol = 8, byrow = TRUE)	

rownames(beta.mat) <- rownames(beta.hi) <- rownames(beta.lo) <- colnames(X)

sigma.pos <- grep("^sigma",colnames(jags.out[[1]]))	
sigma.est <- colMeans(jags.out[[1]][, sigma.pos])

### Posterior

### 	* output (2): .RData file containing an object named 
###			devposterior which is an array of
###			nParties * nConstituencies * nIter

devposterior <- array(NA,dim=c(8,forJags$nConstituencies,nOut))
for (i in 1:nOut) {
	devposterior[,,i] <- matrix(jags.out[[1]][i,p.pos],
		nrow = 8,
		ncol = forJags$nConstituencies,
		byrow = TRUE) 	
}


save(beta.mat,beta.hi,beta.lo, gamma.mat,gamma.hi,gamma.lo, alphaprev.est,alphaep.est, 
	file = paste0("working/individual/",forecastYear,"-mrbetas.RData"))

save(devout,
	refnos,file = paste0("working/constituency/",forecastYear,"devout.RData"))

save(devposterior,refnos,
	file = paste0("working/constituency/",forecastYear,"devposterior.RData"))


