#############################################
### HISTORICAL NATIONAL OUTCOME VARIATION ###
############################################

historicaldata <- read.csv("data/historical/historicalvotecounts.csv",check.names=FALSE)
historicalgovs <- read.csv("data/historical/historicalgovernments.csv",check.names=FALSE)

# calculate historical vote shares
historicalcounts <- historicaldata[,c(-1,-9)]
historicalcounts <- cbind(historicalcounts,historicaldata[9]-rowSums(historicalcounts,na.rm=TRUE))
colnames(historicalcounts)[8] <- "Other"
historicalcounts <- replace(historicalcounts,is.na(historicalcounts),0)
historicalshares <- historicalcounts/rowSums(historicalcounts)
rownames(historicalshares) <- historicaldata[,1]
party <- colnames(historicalshares)
#priormean <- as.numeric(historicalshares[1,])
priormean <- as.numeric(colMeans(historicalshares[1:7,]))

# save lagged vote shares
lag <- as.numeric(t(historicalshares[1,]))

# calculate historical election to election changes
historicaldiffs <- apply(historicalshares[nrow(historicalshares):1,],2,diff)
historicaldiffs <- replace(historicaldiffs ,historicaldiffs == 0,NA)
zerocenteredvar <- function(x) mean(x^2,na.rm=TRUE)
priorvar <- as.numeric(apply(historicaldiffs,2, zerocenteredvar))

# regress historical changes for conservatives and labour on government status
conlabchanges <- c(historicaldiffs[,1],historicaldiffs[,2])
conlaggovs <- c(rev(historicalgovs[-1, 2]),rev(historicalgovs[-1, 3]))
govpenalty <- coef(lm(conlabchanges~ conlaggovs))[2]

# adjust priormean to penalize current government parties and redistribute proportionately to other parties
#latestgov <- c(historicalgovs[1,-1],0)
#priormean[latestgov == 1] <- priormean[latestgov == 1] + govpenalty
#priormean[latestgov != 1] <- priormean[latestgov != 1] - govpenalty*priormean[latestgov != 1]/sum(priormean[latestgov != 1])

# construct output: .RData file containing an object named priorout with nParties rows, and three columns: party, priormean, priorvar

priorout <- data.frame(party,lag,priormean,priorvar)

save(priorout,file=paste0("working/historical/",forecastYear,"priorout.RData"))



####################################
### ESTIMATE UNS ERROR MAGNITUDE ###
####################################

### Read in the data
dat <- read.csv("data/constituency/dorling_data.csv", 
	header = TRUE)
rownames(dat) <- dat$name1997
### Subset to post-war elections
dat <- dat[,c(1,2,165:ncol(dat))]

### Patch on the norris data
norris <- read.csv("data/constituency/norris_1992to2005.csv",header=TRUE)
### And the lookup
lookup <- read.csv("data/constituency/norris2dorling.csv",header=TRUE)
norris <- merge(norris,lookup,all=TRUE)
rm(lookup)
### Convert to numeric
vt.vars <- grep("vt01|vt05",names(norris))
for (v in vt.vars) {
	norris[,v] <- as.numeric(gsub("\\D","",norris[,v]))
}
norris$natvt05 <- rowSums(norris[,c("pcvt05","snpvt05")],na.rm = TRUE)
norris$natvt01 <- rowSums(norris[,c("pcvt01","snpvt01")],na.rm = TRUE)

norris$othvt05 <- rowSums(norris[,c("bnpvt05","grnvt05","ukipvt05","othvt05")],na.rm = TRUE)
norris$othvt01 <- rowSums(norris[,c("bnpvt01","grnvt01","ukipvt01","othvt01")],na.rm = TRUE)
norris$bnpvt05 <- norris$grnvt05 <- norris$ukipvt05 <- norris$totvt05 <- norris$majvt05 <- norris$snpvt05 <- norris$pcvt05 <- NULL
norris$bnpvt01 <- norris$grnvt01 <- norris$ukipvt01 <- norris$totvt01 <- norris$majvt01 <- norris$snpvt01 <- norris$pcvt01 <- NULL

vt.vars <- grep("vt01|vt05",names(norris),value = TRUE)
norris <- norris[,c(vt.vars,"DorlingCode")]
names(norris) <- gsub("vt","20",names(norris))
names(norris) <- gsub("ld","lib",names(norris))

dat <- merge(dat,norris,
	all.x=T,
	all.y=FALSE,
	by.x = "PCA1997",
	by.y = "DorlingCode")

### Delete some vars
constituencies <- dat$name1997
dat$name1997 <- NULL
dat$PCA1997 <- NULL

## Remove the information on the size of the electorate
dat <- dat[,!grepl("^ele",names(dat),ignore.case = TRUE)]
nConstituencies <- length(constituencies)
parties <- c("con","lab","lib","nat","oth")
nParties <- length(parties)

###
elections <- unique(substr(names(dat)[3:ncol(dat)],4,100))
nElections <- length(elections)

### Create array to hold votes of dimensions nConstituencies by nParties by nElections
### Create blank matrix
vote.mat <- array(NA,
	dim = c(nConstituencies,nParties,nElections),
	dimnames = list(1:nConstituencies,parties,elections))

### Populate it
for (e in elections) {
	for (p in parties) {
		tmp.votes <- dat[,(grepl(e,names(dat)) & grepl(p,names(dat)))]
		### Assume vote count of zero means the party didn't stand
		tmp.votes[which(tmp.votes==0)] <- NA
		vote.mat[,p,e] <- tmp.votes
	}
}

dim(vote.mat)

UNS.rmse <- matrix(NA,dim(vote.mat)[3],dim(vote.mat)[2])
rownames(UNS.rmse) <- dimnames(vote.mat)[[3]]
colnames(UNS.rmse) <- dimnames(vote.mat)[[2]]
for (e in 2:length(elections)) {
	for (p in 1:length(parties)) {
		constswing <- vote.mat[,p,e]/rowSums(vote.mat[,,e],na.rm=TRUE) - vote.mat[,p,e-1]/rowSums(vote.mat[,,e-1],na.rm=TRUE)
		natswing <- sum(vote.mat[,p,e],na.rm=TRUE)/sum(vote.mat[,,e],na.rm=TRUE) - sum(vote.mat[,p,e-1],na.rm=TRUE)/sum(vote.mat[,,e-1],na.rm=TRUE)
		UNS.rmse[e,p] <- sqrt(mean((constswing - natswing)^2,na.rm=TRUE))
	}
}	

RMSEUNS <- mean(UNS.rmse[11:15,1:3])
save(vote.mat,file=paste0("working/constituency/",forecastYear,"_historical_votemat.RData"))
save(RMSEUNS,file=paste0("working/constituency/",forecastYear,"RMSEUNS.RData"))

################################
### CALIBRATE NATIONAL POLLS ###
################################

### Libraries
library(rjags)

ElectionYears <- c(1979,1983,1987,1992,1997,2001,2005,2010)
ElectionYearCount <- length(ElectionYears)
CalibrationResults <- data.frame(Year=rep(0, ElectionYearCount*4),Party=rep("", ElectionYearCount*4),CampaignDay=rep(0, ElectionYearCount*4),ElectionDay=rep(0, ElectionYearCount*4),ElectionDaySE=rep(0, ElectionYearCount*4),ElectionOutcome=rep(0, ElectionYearCount*4),LaggedElectionOutcome=rep(0, ElectionYearCount*4),stringsAsFactors=FALSE)
HistoricalPools <- array(NA,dim=c(365,4, ElectionYearCount))
HistoricalSDs <- array(NA,dim=c(365,4, ElectionYearCount))
counter <- 0

for (ElectionYear in ElectionYears){

### Load validation data

ValiData <- read.csv("data/historical/historicalvotecounts.csv")
temp <- ValiData[which(ValiData$Year == ElectionYear),]
temp <- as.numeric(c(temp[2],temp[3],temp[4],temp[9]-temp[2]-temp[3]-temp[4]))
electionoutcome <- temp/sum(temp)
temp <- ValiData[which(ValiData$Year == ElectionYear)+1,]
temp <- as.numeric(c(temp[2],temp[3],temp[4],temp[9]-temp[2]-temp[3]-temp[4]))
laggedelectionoutcome <- temp/sum(temp)

### Load polling data

Data <- read.csv(paste0("data/polling/polling-report-",ElectionYear,".csv"))
Data <- Data[rowSums(is.na(Data)) == 0,]
if (ElectionYear == 2010) Data$NumDate <- as.numeric(as.Date(Data$Date)) - min(as.numeric(as.Date(Data$Date))) + 1
if (ElectionYear == 2005 | ElectionYear == 1997) Data$NumDate <- as.numeric(as.Date(Data$Date,format="%d/%m/%Y")) - min(as.numeric(as.Date(Data$Date,format="%d/%m/%Y"))) + 1
if (ElectionYear == 2001) {
	Data$Date[Data$Date == "??/06/01"] <- "05/06/2001"
	Data$NumDate <- as.numeric(as.Date(Data$Date,format="%d/%m/%Y")) - min(as.numeric(as.Date(Data$Date,format="%d/%m/%Y"))) + 1
	}
if (ElectionYear == 1992 | ElectionYear == 1987 | ElectionYear == 1983 | ElectionYear == 1979) {
	wrongformat <- grep("*",Data$Date,fixed=TRUE)
	fixedformat <- as.Date(sub("*","",Data$Date[wrongformat],fixed=TRUE))
	Data$Date <- as.Date(Data$Date,format="%d/%m/%Y")
	Data$Date[wrongformat] <- fixedformat
	Data$NumDate <- as.numeric(Data$Date) - min(as.numeric(Data$Date)) + 1
	}	
### Limit to Last Year of Data Before Election

Data <- subset(Data,(max(	Data$NumDate)-	Data$NumDate) < 365)	
Data$NumDate <- Data$NumDate - max(	Data$NumDate) + 365

Data$House <- factor(sub("/.*","",Data$House))	

party.names <- c("Con","Lab","LD")
Data[,party.names] <- Data[,party.names]/100
Data$Other <- 1 - rowSums(Data[,party.names])
party.names <- c("Con","Lab","LD","Other")


### Jags script

model <- '
	model
	{
		for (i in 1:nPolls) {
			for (j in 1:nParties) {
			    mu[i,j] <- house[org[i],j] + alpha[date[i],j]
			    y[i,j] ~ dnorm(mu[i,j], 1000)
			}
		}

		### Walk starts one year before the election
		for (j in 1:nParties) {
			alphastar[1,j] ~ dnorm(0,0.0001)
			ealphastar[1,j] <- exp(alphastar[1,j])
			alpha[1,j] <- ealphastar[1,j]/sum(ealphastar[1,])
		}

		### Walk goes forwards in time
		for (i in 2:nPeriods) {
			for (j in 1:nParties) {
				alphastar[i,j] ~ dnorm(alphastar[i-1,j], tau[j])
				ealphastar[i,j] <- exp(alphastar[i,j])
				alpha[i,j] <- ealphastar[i,j]/sum(ealphastar[i,])
			}
		}
		
		### Step size
		for (j in 1:nParties) {
			tau[j] <- 1/pow(sigma[j], 2)
			sigma[j] ~ dunif(0, 0.1)
		}
		
		### Double centered house effects
		for (i in 1:nOrgs) {
			for (j in 1:nParties) {
				house.star[i,j] ~ dnorm(0, 10000)
				house[i,j] <- house.star[i,j] - mean(house.star[,j]) - mean(house.star[i,]) + mean(house.star[,])
			}
		}
		
	}
'



### Prep the data for Jags

forJags <- list(y = as.matrix(Data[,party.names]),
	nParties = 4,
	nPeriods = 365,
	nPolls = nrow(Data),
	nOrgs = length(unique(Data$House)),
	date = Data$NumDate,
	org = as.numeric(Data$House))

### Adaptation

model.jags <- jags.model(file=textConnection(model)
				,data=forJags,n.adapt= 100)	
				
### Burn-in		
			
update(model.jags, 1000)	

### Sample

posterior <- jags.samples(model.jags, 
				c("alpha","house","sigma"), 
				n.iter= 1000, thin = 1, progress.bar="text")	

### Examine output

HistoricalPools[,,counter/4+1] <- apply(posterior$alpha,c(1,2),mean)
HistoricalSDs[,,counter/4+1] <- apply(posterior$alpha,c(1,2),sd)

CalibrationResults$Year[counter+1:4] <- ElectionYear
CalibrationResults$Party[counter+1:4] <- as.character(party.names)
CalibrationResults$ElectionOutcome[counter+1:4] <- electionoutcome
CalibrationResults$LaggedElectionOutcome[counter+1:4] <- laggedelectionoutcome
counter <- counter + 4

} # end loop over elections

### 
save(HistoricalPools,file = paste0("working/historical/",forecastYear,"_historical_pools.RData"))

### Add Government Data

CalibrationResults$Gov <- rep(NA,nrow(CalibrationResults))
GovtData <- read.csv("data/historical/historicalgovernments.csv")[,1:5]
colnames(GovtData) <- c("Year","Con","Lab","LD","Other")
for (i in 1:nrow(CalibrationResults)){
	CalibrationResults$Gov[i] <- GovtData[which(CalibrationResults$Year[i] == GovtData$Year)+1,2+((i-1) %% 4)]
}

### Compute Calibration Results

CalibrationResults$ElectionDiff <- CalibrationResults$ElectionOutcome - CalibrationResults$LaggedElectionOutcome
calibrout <- list(coef=rep(NA,365),prec=rep(NA,365),sigma=rep(NA,365))
for (i in 1:365){
	PollDiff <- as.numeric(HistoricalPools[i,,]) - CalibrationResults$LaggedElectionOutcome
	lmtmp <- lm(ElectionDiff ~ PollDiff - 1,data=CalibrationResults,subset=(Year < forecastYear) & Party != "Other")
	
	calibrout$coef[i] <- as.numeric(coef(lmtmp))
	calibrout$prec[i] <- solve(vcov(lmtmp))
	calibrout$sigma[i] <- summary(lmtmp)$sigma

}

save(CalibrationResults,file = paste0("working/historical/",forecastYear,"CalibrationResults.RData"))
save(calibrout,file = paste0("working/historical/",forecastYear,"calibrout_unsmoothed.RData"))

# Local linear smooth with kernel running from 1.5x-10 time to election to 0.5x+10 time to election
tmp <- rep(NA,365)
for (i in 365:1){
	w <- 1*(abs(365:1 - i) < (10+0.5*i))
	x <- 1:365
	tmp[366-i] <- fitted(lm(calibrout$coef~x,weights=w))[366-i]	
	}
calibrout$coef <- tmp	
calibrout$coef[1:which.min(calibrout$coef)] <- min(calibrout$coef)

### Correct Sigma for Pooled Polls Uncertainty (another hacky kludge)
calibrout$sigma <- sqrt(calibrout$sigma^2 - apply(HistoricalSDs[,1:3,],c(1),mean)^2)

### Smooth Sigma Calibration Results with Hacky Kludge

sigmamodel <- lm(calibrout$sigma~ 1)
fittedsigma <- fitted(sigmamodel)
calibrout$sigma[1:rev(which(calibrout$sigma  > fittedsigma))[1]] <- fitted(sigmamodel)[1:rev(which(calibrout$sigma  > fittedsigma))[1]]



### Plot Calibration Results

if (FALSE){
	sigmabyday <- rep(0,365)
	poolcoefbyday <- rep(0,365)
	lagcoefbyday <- rep(0,365)
	for (i in 1:365){	
		poolcoefbyday[i] <- calibrout$coef[i]
		sigmabyday[i] <- calibrout$sigma[i]
	}
	
	plot(1:365,sigmabyday)
	plot(1:365, poolcoefbyday,type="l",col="blue",ylim=c(0,1))
#	lines(1:365, lagcoefbyday,col="red")
}	

### Save Calibration Model

save(calibrout,file = paste0("working/historical/",forecastYear,"calibrout.RData"))

