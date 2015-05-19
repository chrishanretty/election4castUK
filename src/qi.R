### Set alpha level for prediction intervals (e.g. 0.1 for 90% intervals)
predalpha <- 0.1

### * input (1): .RData file named devposterior
### * input (2): .RData file named poolposterior
### * input (3): .RData file named houseout
### * input (4): .RData file named poolout
### * input (5): .RData file named electionout
### * input (6): .csv file with historical vote counts
### * output (1): .csv which has size nConstituencies by nParties, where each cell has predicted vote share
### * output (2): .csv which has size nConstituencies by nParties, where each cell has predicted probability of victory
### * output (3): .csv which has number of rows equal to nParties, with columns named mean, hi, lo, and Pr(x is largest); Pr(x>325); Pr(x>316); Pr(maj. w/o Scotland); 
### * output (4): merge of output (1) and (2) with information in data/constituency/constituencyinfo.csv
### * output (5): map of UK with constituencies coloured by party, and alpha level set to predicted probability of winning
### * output (6): trajectory of party support
### * output (7): plot of house effects
### * output (8): maps of predicted winners and party support
### * output (9): bar charts of party seat distributions


### Libraries
try(detach("package:car", unload=TRUE),silent=TRUE)
library(reshape)
library(ggplot2)
library(ggmap)
require(fGarch)
## library(gridSVG)
library(proto)
library(plyr)
library(scales)
library(munsell)
library(grid)
library(rgeos)
library(SortableHTMLTables)
library(knitr)
library(pander)
library(gridExtra)

watermark <- function(grob) {
	grob <- arrangeGrob(grob,
		 sub = textGrob("electionforecast.co.uk", x = 1, hjust = 1, vjust=0.1,
		                            gp = gpar(fontface = "italic", fontsize = 18)))
	return(grob)
}

### Date parsing

elex.date <- as.Date("2015-05-07")

origin.date <- elex.date - 365

party.names <- c("Conservatives","Labour","Liberal Democrats",
	"SNP","Plaid Cymru","Greens","UKIP","Other")

short.party.names <- c("Con","Lab","LD",
	"SNP","PC","GRN","UKIP","Oth")

### ###########################################################################
### * input (0): TRUE/FALSE variable named "poststratification"
### ###########################################################################
if (!exists("poststratification")) {
	poststratification <- FALSE
}


### ###########################################################################
### * input (1): .RData file named devposterior
### * input (2): .RData file named poolposterior
### * input (3): .RData file named houseout
### * input (4): .RData file named poolout
### ###########################################################################


load(paste0("working/constituency/",forecastYear,"devposterior.RData"))

### Zero out pred.s for parties not standing in a seat
standingCandidates <- read.csv("data/constituency/standing_candidates.csv", 
	header = T)

for (j in 1:632) {
	tmp <- devposterior[,j,]
	tozero <- which(standingCandidates[j,-1] == 0)
	tmp[tozero,] <- 0.0001
	devposterior[,j,] <- tmp
}
## Once again, kludge for the speaker's seat
devposterior[,101,] <- matrix(c(0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.9993),nrow=dim(devposterior)[1],ncol=dim(devposterior)[3])
## 

load(paste0("working/polling/",forecastYear,"poolposterior.RData"))

load(paste0("working/polling/",forecastYear,"electionposterior.RData"))

load(paste0("working/polling/",forecastYear,"poolout.RData"))

load(paste0("working/polling/",forecastYear,"electionout.RData"))

load(paste0("working/polling/",forecastYear,"houseout.RData"))

load(paste0("working/polling/",forecastYear,"pollsterforecasts.RData"))

load(paste0("working/polling/",forecastYear,"pollcountdata.RData"))

load(paste0("working/constituency/",forecastYear,"prevvote.RData"))

load(paste0("working/constituency/",forecastYear,"prevwin.RData"))

# load(paste0("working/constituency/",forecastYear,"niposterior.RData"))

canonicalseatvars <- read.csv("data/constituency/canonical_seatvars_2011census.csv",
	header=TRUE)
### Order canonical seat variables by the refnos
canonicalseatvars <- canonicalseatvars[charmatch(refnos, canonicalseatvars$refno),]

historicaldata <- read.csv("data/historical/historicalvotecounts.csv",check.names=FALSE)
historicaldata[,2:8] <- historicaldata[,2:8]  / historicaldata[,9]  

### ###########################################################################
### Start creating stuff
### Specifically, add the slices from the devposterior
### To the slices from the pooling posterior for the last day
### ###########################################################################

### Estimates prediction for today, unless today is after elex.date, 
### in which case estimates prediction on election day
if (Sys.Date() > elex.date){
	choosepredictday <- 365
	} else {
	choosepredictday <- as.numeric(Sys.Date() - origin.date)
	}	

### Re-create electionday as an array of 8, 632, posterior draws
predictday <- array(NA,dim=c(8,632,dim(electionposterior)[3]))
for (i in 1:632) {
	predictday[,i,] <- electionposterior[,choosepredictday,]
}


### ################################
### Merge the output of the mr 
### model with the pooled polls
### ################################

turnoutweights <- canonicalseatvars$turn10
turnoutweights[which(is.na(canonicalseatvars$turn10))] <- canonicalseatvars$turn05[which(is.na(canonicalseatvars$turn10))]	
turnoutweights <- 632*turnoutweights/sum(turnoutweights)

swingmodel <- "GNS"

if (swingmodel == "UNS"){ # Uniform National Swing 
	
	poolmr <- aperm(array(apply(predictday,c(1,3),mean) - apply(devposterior,c(1,3),mean),dim(devposterior)[c(1,3,2)]),c(1,3,2))
	
	comb.posterior <- devposterior + poolmr
	
	comb.posterior <- replace(comb.posterior, comb.posterior < 0,0)
	
}	

if (swingmodel == "PNS"){ # Proportional National Swing 

	poolmr <- aperm(array(apply(predictday,c(1,3),mean)/apply(devposterior,c(1,3),mean),dim(devposterior)[c(1,3,2)]),c(1,3,2))

	comb.posterior <- devposterior * poolmr
	
}	

	logit <- function(p) log((p)/(1-p))
	ilogit <- function(x) exp(x)/(1+exp(x))
	
if (swingmodel == "LNS"){ # Logistic National Swing
	
	logisticswing <- function(con,target){
		obj <- function(x) return(mean(ilogit(logit(con) + x)) - target)
		return(ilogit(logit(con)+uniroot(obj,c(-10,10))$root))
		}
	
	comb.posterior <- devposterior
	for (sim in 1:dim(devposterior)[3]){
		for (party in 1:dim(devposterior)[1]){
			comb.posterior[party,,sim] <- logisticswing(devposterior[party,,sim],predictday[party,1,sim])	
		}
	}
	
}	

if (swingmodel == "GNS"){ # Generalized Error National Swing with Turnout Weights
	
	library(fGarch)
	
	nu <- 10
	gammainvnu <- gamma(1/nu)
	gamma3invnu <- gamma(3/nu)
	
	# generr <- function(p) qged(p,mean=0,sd=1,nu=nu)
	generr <- function (p) {
   		lambda = sqrt(2^(-2/nu) * (gammainvnu)/(gamma3invnu))
    	q = lambda * (2 * qgamma((abs(2 * p - 1)), 1/nu))^(1/nu)
    	return(q * sign(2 * p - 1))
    	}
    	
	# igenerr <- function(q) pged(q,mean=0,sd=1,nu=nu)
	igenerr <- function (q) {
    	lambda = sqrt(2^(-2/nu) * (gammainvnu)/(gamma3invnu))
    	g = nu/(lambda * (2^(1+1/nu)) * gammainvnu)
    	h = 2^(1/nu) * lambda * g * gammainvnu/nu
    	s = 0.5 * (abs(q)/lambda)^nu
    	return(0.5 + sign(q) * h * pgamma(s, 1/nu))
		}

	generrorswing <- function(con,turnout,target){
		obj <- function(x) return(mean(turnout*igenerr(generr(con) + x)) - target)
		return(igenerr(generr(con)+uniroot(obj,c(-nu,nu))$root))
		}
	
	print("Reconciling National Vote Shares with Constituency Vote Shares...")
	comb.posterior <- devposterior
	for (sim in 1:dim(devposterior)[3]){
		for (party in 1:dim(devposterior)[1]){
			comb.posterior[party,,sim] <- generrorswing(devposterior[party,,sim],turnoutweights,predictday[party,1,sim])	
		}
		comb.posterior[,,sim] <- comb.posterior[,,sim] / matrix(colSums(comb.posterior[,,sim],na.rm = TRUE),dim(devposterior)[1],dim(devposterior)[2],byrow=TRUE)
		if (sim %% 100 == 0) print(paste0(sim," iterations reconciled..."))
	}
	print("Done.")
	
}	

### ############################################################################
### * output (1): .csv which has size nConstituencies by nParties, 
### where each cell has predicted vote share
### ############################################################################

predvote.bar <- apply(comb.posterior,c(1,2),mean)
predvote.bar <- t(predvote.bar)
### Enforce sum-to-one constrain
### predvote.bar <- predvote.bar / rowSums(predvote.bar,na.rm = TRUE)
rownames(predvote.bar) <- as.character(refnos)
colnames(predvote.bar) <- party.names

predvote.int <- apply(comb.posterior,c(1,2),quantile,c(predalpha/2,1-predalpha/2))
predvote.int <- aperm(predvote.int,c(3,2,1))

write.csv(predvote.bar,
	file= paste0("outputs/csv/",forecastYear,"_predicted_vote_share_by_seat.csv"),
	row.names = TRUE)	
	
nowvote.bar <- apply(devposterior,c(1,2),mean)
nowvote.bar <- t(nowvote.bar)
rownames(nowvote.bar) <- as.character(refnos)
colnames(nowvote.bar) <- party.names

nowvote.se <- apply(devposterior,c(1,2),sd)
nowvote.se <- t(nowvote.se)
rownames(nowvote.se) <- as.character(refnos)
colnames(nowvote.se) <- party.names

write.csv(nowvote.bar,
	file= paste0("outputs/csv/",forecastYear,"_current_vote_share_by_seat.csv"),
	row.names = TRUE)	
	
write.csv(nowvote.se,
	file= paste0("outputs/csv/",forecastYear,"_current_vote_share_by_seat.se.csv"),
	row.names = TRUE)		

### ###########################################################################
### * output (2): .csv which has size nConstituencies by nParties, 
### where each cell has predicted probability of victory
### ###########################################################################

predwin <- matrix(NA,
	nrow = dim(comb.posterior)[2],
	ncol = 8)

winner <- apply(comb.posterior,c(2,3),which.max)
for (j in 1:8){
	predwin[,j] <- apply(winner,1,function(x)mean(x==j)) 
}

winnersbyseat <- table(factor(apply(predwin ,1,which.max),1:8))
names(winnersbyseat) <- party.names
colnames(predwin) <- party.names
rownames(predwin) <- refnos

write.csv(predwin,
	file= paste0("outputs/csv/",forecastYear,"_predicted_probability_winning_by_seat.csv"),
	row.names = TRUE)
	
	
nowwin <- matrix(NA,
	nrow = dim(devposterior)[2],
	ncol = 8)

leader <- apply(devposterior,c(2,3),which.max)
for (j in 1:8){
	nowwin[,j] <- apply(leader,1,function(x)mean(x==j)) 
}

leadersbyseat <- table(factor(apply(nowwin ,1,which.max),1:8))
names(leadersbyseat) <- party.names
colnames(nowwin) <- party.names
rownames(nowwin) <- refnos

write.csv(predwin,
	file= paste0("outputs/csv/",forecastYear,"_current_probability_winning_by_seat.csv"),
	row.names = TRUE)	
	
# Seat swings by party pairs

PartySwingTable <- table(factor(party.names[apply(prevwin,1,which.max)],levels=party.names),factor(party.names[apply(predwin,1,which.max)],levels=party.names))
rownames(PartySwingTable) <- paste0("from ",short.party.names)
colnames(PartySwingTable) <- paste0("to ",short.party.names)	

### ############################################################################
### * output (3): .csv which has number of rows equal to nParties, with 
### columns named mean, hi, lo, and Pr(x is largest); 
### Pr(x>325); Pr(x>316); Pr(maj. w/o Scotland); 
### ############################################################################

canonicalseatvars <- canonicalseatvars[,c("refno","YouGovName","region","win05","win10")]
names(canonicalseatvars) <- c("RefNo","Seat","Region","Win05","Win10")
canonicalseatvars <- canonicalseatvars[order(canonicalseatvars$RefNo),]
canonicalseatvars$Win05 <- car:::recode(canonicalseatvars$Win05,
	"'Conservative'='Conservatives';
	'Labour'='Labour';
	'Lib-Dem'='Liberal Democrats';
	'PC'='Plaid Cymru';
	'SNP'='SNP'")
canonicalseatvars$Win10 <- car:::recode(canonicalseatvars$Win10,
	"'Conservative'='Conservatives';
	'Labour'='Labour';
	'Lib-Dem'='Liberal Democrats';
	'PC'='Plaid Cymru';
	'SNP'='SNP'")	
	
exScotland <- canonicalseatvars$Region != "Scotland"

outcomes.df <- data.frame(Party = party.names,
	Mean = rep(NA,length(party.names)),
	Lo = rep(NA,length(party.names)),
	Hi = rep(NA,length(party.names)),
	Maj326 = rep(NA,length(party.names)),
	Maj316 = rep(NA,length(party.names)),
	Largest = rep(NA,length(party.names)))

seats.by.iter <- lapply(1:8,function(x)apply(winner,2,function(y)sum(y==x)))
seats.by.iter <- do.call("cbind",seats.by.iter)

for (j in 1:8){
	outcomes.df$Mean[j] <- mean(seats.by.iter[,j]) 
	outcomes.df$Lo[j] <- round(quantile(seats.by.iter[,j],probs= predalpha/2))
	outcomes.df$Hi[j] <- round(quantile(seats.by.iter[,j],probs= 1 - predalpha/2))
	outcomes.df$Maj326[j] <- mean(seats.by.iter[,j] > 325)
	outcomes.df$Maj316[j] <- mean(seats.by.iter[,j] > 316)
	outcomes.df$Largest[j] <- mean(apply(seats.by.iter,1,which.max) == j)
}

# Apply Webster's method to make total seats equal 632.
D <- 632
tmptotal <- sum(round(632*outcomes.df$Mean/D))
while (tmptotal != 632){
	if (tmptotal < 632) D <- D - 0.0004 
	if (tmptotal > 632) D <- D + 0.0004
	tmptotal <- sum(round(632*outcomes.df$Mean/D))
	print(D)
	}
outcomes.df$Mean <- round(632*outcomes.df$Mean/D)






current.df <- data.frame(Party = party.names,
	Mean = rep(NA,length(party.names)),
	Lo = rep(NA,length(party.names)),
	Hi = rep(NA,length(party.names)),
	Maj326 = rep(NA,length(party.names)),
	Maj316 = rep(NA,length(party.names)),
	Largest = rep(NA,length(party.names)))

current.seats.by.iter <- t(apply(apply(devposterior,c(2,3),which.max),2,tabulate,nbins=8))

for (j in 1:8){
	current.df$Mean[j] <- mean(current.seats.by.iter[,j]) 
	current.df$Lo[j] <- round(quantile(current.seats.by.iter[,j],probs= predalpha/2))
	current.df$Hi[j] <- round(quantile(current.seats.by.iter[,j],probs= 1 - predalpha/2))
	current.df$Maj326[j] <- mean(current.seats.by.iter[,j] > 325)
	current.df$Maj316[j] <- mean(current.seats.by.iter[,j] > 316)
	current.df$Largest[j] <- mean(apply(current.seats.by.iter,1,which.max) == j)
}

# Apply Webster's method to make total seats equal 632.
D <- 632
tmptotal <- sum(round(632* current.df$Mean/D))
while (tmptotal != 632){
	if (tmptotal < 632) D <- D - 0.0004 
	if (tmptotal > 632) D <- D + 0.0004
	tmptotal <- sum(round(632* current.df$Mean/D))
	print(D)
	}
current.df$Mean <- round(632* current.df$Mean/D)







outcomesExS.df <- data.frame(Party = party.names,
	Mean = rep(NA,length(party.names)),
	Lo = rep(NA,length(party.names)),
	Hi = rep(NA,length(party.names)),
	Maj296 = rep(NA,length(party.names)),
	Largest = rep(NA,length(party.names)))
	
seatsExS.by.iter <- lapply(1:8,function(x)apply(winner[exScotland,],2,function(y)sum(y==x)))
seatsExS.by.iter <- do.call("cbind", seatsExS.by.iter)	

for (j in 1:8){
	outcomesExS.df$Mean[j] <- mean(seatsExS.by.iter[,j]) 
	outcomesExS.df$Lo[j] <- quantile(seatsExS.by.iter[,j],probs= predalpha/2)
	outcomesExS.df$Hi[j] <- quantile(seatsExS.by.iter[,j],probs= 1 - predalpha/2)
	outcomesExS.df$Maj296[j] <- mean(seatsExS.by.iter[,j] > 296)
	outcomesExS.df$Largest[j] <- mean(apply(seatsExS.by.iter,1,which.max) == j)
}

# Apply Webster's method to make total seats equal 632.
D <- 573
tmptotal <- sum(round(573* outcomesExS.df$Mean/D))
while (tmptotal != 573){
	if (tmptotal < 573) D <- D - 0.0004 
	if (tmptotal > 573) D <- D + 0.0004
	tmptotal <- sum(round(573* outcomesExS.df$Mean/D))
	print(D)
	}
outcomesExS.df$Mean <- round(573* outcomesExS.df$Mean/D)

### Prepare Vote by Region Tables ###

reg.df <- merge(predvote.bar,canonicalseatvars,
	by.x = "row.names",by.y = "RefNo",
	all.x=T,all.y=F)
reg.df$Turnout <- turnoutweights
regional.predvote <- matrix(NA,nlevels(reg.df$Region),8)
colnames(regional.predvote) <- short.party.names
rownames(regional.predvote) <- levels(reg.df$Region)
for (i in 1:nlevels(reg.df$Region)){
	regional.predvote[i,] <- colSums((matrix(reg.df$Turnout,632,8)*reg.df[,party.names])[as.numeric(reg.df$Region) == i,])/sum(reg.df$Turnout[as.numeric(reg.df$Region) == i])
	
	}
pretty.regional.predvote <- round(regional.predvote*100)
write.csv(pretty.regional.predvote,
	file = paste0("outputs/csv/",forecastYear,"_predicted_outcomes_regional.csv"),
	row.names= TRUE)

reg.df <- merge(nowvote.bar,canonicalseatvars,
	by.x = "row.names",by.y = "RefNo",
	all.x=T,all.y=F)
reg.df$Turnout <- turnoutweights
regional.nowvote <- matrix(NA,nlevels(reg.df$Region),8)
colnames(regional.nowvote) <- short.party.names
rownames(regional.nowvote) <- levels(reg.df$Region)
for (i in 1:nlevels(reg.df$Region)){
	regional.nowvote[i,] <- colSums((matrix(reg.df$Turnout,632,8)*reg.df[,party.names])[as.numeric(reg.df$Region) == i,])/sum(reg.df$Turnout[as.numeric(reg.df$Region) == i])
	
	}
pretty.regional.nowvote <- round(regional.nowvote*100)
write.csv(pretty.regional.nowvote,
	file = paste0("outputs/csv/",forecastYear,"_current_outcomes_regional.csv"),
	row.names= TRUE)



write.csv(outcomes.df,
	file = paste0("outputs/csv/",forecastYear,"_predicted_outcomes_nationwide.csv"),
	row.names= FALSE)
	
three.party.absolute.error.by.iter <- rowSums(abs((seats.by.iter - matrix(outcomes.df$Mean,nrow=nrow(seats.by.iter),ncol=8,byrow=TRUE))[,1:3]))
three.party.absolute.error.mean <- mean(three.party.absolute.error.by.iter)
three.party.absolute.error.lo <- quantile(three.party.absolute.error.by.iter,probs= predalpha/2)
three.party.absolute.error.hi <- quantile(three.party.absolute.error.by.iter,probs= 1 - predalpha/2)

### Prepare Pretty Predicted Vote Tables

predvote.barpretty <- as.data.frame(round(predvote.bar*100,0))
predvote.barpretty$PredictedWinner <- apply(predvote.bar,1,which.max)
predvote.barpretty$PredictedWinner <- party.names[predvote.barpretty$PredictedWinner]

predvote.barpretty$RefNo <- refnos
predvote.barpretty <- merge(predvote.barpretty,canonicalseatvars,
	by = "RefNo",
	all.x=T,all.y=F)
predvote.barpretty$Type <- ifelse(predvote.barpretty$Win10 == predvote.barpretty$PredictedWinner,
	paste0(predvote.barpretty$PredictedWinner," hold"),
	paste0(predvote.barpretty$PredictedWinner," gain"))		

### CSV output
write.csv(predvote.barpretty,
	file = paste0("outputs/csv/",forecastYear,"_predicted_vote_by_seat-pretty.csv"),
	row.names= FALSE)

### Sortable HTML table
sortable.html.table(predvote.barpretty,
	output.directory = "outputs/html/",
	output.file = paste0(forecastYear,"_predicted_vote_by_seat.html"),
	page.title = "Vote share predictions")
	
### Change sort order
changeSortOrder <- function(infile) {
	instring <- readLines(infile)
	instring <- sub("tablesorter();","tablesorter({sortInitialOrder: 'desc'});",instring,fixed=T)
	writeLines(instring,con = infile)
}
addColour <- function(infile) {
	instring <- readLines(infile)
	instring <- sub("style.css","style2.css",instring,fixed=T)
	writeLines(instring,con = infile)
}

changeSortOrder(paste0("outputs/html/",forecastYear,"_predicted_vote_by_seat.html"))
addColour(paste0("outputs/html/",forecastYear,"_predicted_vote_by_seat.html"))

### Simplify for Webpage Table

if (forecastYear == 2015){

	predvote.web <- predvote.barpretty[,c("Conservatives","Labour","Liberal Democrats","SNP","Plaid Cymru","Greens","UKIP","Other","Seat","Region","Win10")]
	for (j in 1:8) predvote.web[,j] <- round(predvote.web[,j] )
	names(predvote.web)[11] <- "2010"
		
	### Sortable HTML table
	sortable.html.table(predvote.web,
		output.directory = "outputs/website/2015/tables/",
		output.file = paste0("predicted_vote_by_seat.html"),
		page.title = "Vote share predictions")	
	changeSortOrder(paste0("outputs/website/2015/tables/","predicted_vote_by_seat.html"))
	addColour(paste0("outputs/website/2015/tables/","predicted_vote_by_seat.html"))
}	


### Prepare Pretty Vote Interval Tables



predvote.intpretty <- as.data.frame(apply(round(100*predvote.int),c(1,2),paste,collapse="-"), stringsAsFactors=FALSE)
rownames(predvote.intpretty) <- rownames(predvote.bar)
colnames(predvote.intpretty) <- colnames(predvote.bar)
predvote.intpretty$PredictedWinner <- apply(predvote.bar,1,which.max)
predvote.intpretty$PredictedWinner <- party.names[predvote.barpretty$PredictedWinner]

predvote.intpretty$RefNo <- refnos
predvote.intpretty <- merge(predvote.intpretty,canonicalseatvars,
                            by = "RefNo",
                            all.x=T,all.y=F)
predvote.intpretty$Type <- ifelse(predvote.barpretty$Win10 == predvote.barpretty$PredictedWinner,
                                                            paste0(predvote.barpretty$PredictedWinner," hold"),
                                                            paste0(predvote.barpretty$PredictedWinner," gain"))		

### CSV output
write.csv(predvote.intpretty,
          file = paste0("outputs/csv/",forecastYear,"_predicted_interval_by_seat-pretty.csv"),
          row.names= FALSE)

### Sortable HTML table
sortable.html.table(predvote.intpretty,
                    output.directory = "outputs/html/",
                    output.file = paste0(forecastYear,"_predicted_interval_by_seat.html"),
                    page.title = "Vote share intervals")


changeSortOrder(paste0("outputs/html/",forecastYear,"_predicted_interval_by_seat.html"))
addColour(paste0("outputs/html/",forecastYear,"_predicted_interval_by_seat.html"))
### Simplify for Webpage Table

if (forecastYear == 2015){
  
  intvote.web <- predvote.intpretty[,c("Conservatives","Labour","Liberal Democrats","SNP","Plaid Cymru","Greens","UKIP","Other","Seat","Region","Win10")]
  names(predvote.web)[11] <- "2010"
  
  ### Sortable HTML table
  sortable.html.table(intvote.web,
                      output.directory = "outputs/website/2015/tables/",
                      output.file = paste0("predicted_interval_by_seat.html"),
                      page.title = "Vote share intervals")	
  changeSortOrder(paste0("outputs/website/2015/tables/","predicted_interval_by_seat.html"))
	addColour(paste0("outputs/website/2015/tables/","predicted_interval_by_seat.html"))
}	


### Prepare Pretty Current Vote Tables

nowvote.barpretty <- as.data.frame(round(nowvote.bar*100,0))
nowvote.barpretty$PredictedWinner <- apply(nowvote.bar,1,which.max)
nowvote.barpretty$PredictedWinner <- party.names[nowvote.barpretty$PredictedWinner]

nowvote.barpretty$RefNo <- refnos
nowvote.barpretty <- merge(nowvote.barpretty,canonicalseatvars,
	by = "RefNo",
	all.x=T,all.y=F)
nowvote.barpretty$Type <- ifelse(nowvote.barpretty$Win10 == nowvote.barpretty$PredictedWinner,
	paste0(nowvote.barpretty$PredictedWinner," hold"),
	paste0(nowvote.barpretty$PredictedWinner," gain"))		

### CSV output
write.csv(nowvote.barpretty,
	file = paste0("outputs/csv/",forecastYear,"_current_vote_by_seat-pretty.csv"),
	row.names= FALSE)

### Sortable HTML table
sortable.html.table(nowvote.barpretty,
	output.directory = "outputs/html/",
	output.file = paste0(forecastYear,"_current_vote_by_seat.html"),
	page.title = "Forecast Vote share")
changeSortOrder(paste0("outputs/html/",forecastYear,"_current_vote_by_seat.html"))
addColour(paste0("outputs/html/",forecastYear,"_current_vote_by_seat.html"))

	
### Simplify for Webpage Table

if (forecastYear == 2015){

	nowvote.web <- nowvote.barpretty[,c("Conservatives","Labour","Liberal Democrats","SNP","Plaid Cymru","Greens","UKIP","Other","Seat","Region","Win10")]
	for (j in 1:8) nowvote.web[,j] <- round(nowvote.web[,j] )
	names(nowvote.web)[11] <- "2010"
		
	### Sortable HTML table
	sortable.html.table(nowvote.web,
		output.directory = "outputs/website/2015/tables/",
		output.file = paste0("current_vote_by_seat.html"),
		page.title = "Current vote share")	
	changeSortOrder(paste0("outputs/website/2015/tables/","current_vote_by_seat.html"))
	addColour(paste0("outputs/website/2015/tables/","current_vote_by_seat.html"))
	
}	

### Prepare Pretty Predicted Win Tables

predwin.pretty <- as.data.frame(round(predwin*100,0))

predwin.pretty$RefNo <- refnos
predwin.pretty <- merge(predwin.pretty,canonicalseatvars,
	by = "RefNo",
	all.x=T,all.y=F)

### CSV output

write.csv(predwin.pretty,
	file = paste0("outputs/csv/",forecastYear,"_predicted_winner_by_seat-pretty.csv"),
	row.names= FALSE)

### Sortable HTML table
sortable.html.table(predwin.pretty,
	output.directory = "outputs/html/",
	output.file = paste0(forecastYear,"_predicted_probability_winning_by_seat.html"),
	page.title = "Predicted probability of winning seat")
changeSortOrder(paste0("outputs/html/",forecastYear,"_predicted_probability_winning_by_seat.html"))
addColour(paste0("outputs/html/",forecastYear,"_predicted_probability_winning_by_seat.html"))
	
### Prepare Pretty Current Win Tables	

nowwin.pretty <- as.data.frame(round(nowwin*100,0))

nowwin.pretty$RefNo <- refnos
nowwin.pretty <- merge(nowwin.pretty,canonicalseatvars,
	by = "RefNo",
	all.x=T,all.y=F)

### CSV output

write.csv(nowwin.pretty,
	file = paste0("outputs/csv/",forecastYear,"_current_winner_by_seat-pretty.csv"),
	row.names= FALSE)

### Sortable HTML table
sortable.html.table(nowwin.pretty,
	output.directory = "outputs/html/",
	output.file = paste0(forecastYear,"_current_probability_winning_by_seat.html"),
	page.title = "Current probability of winning seat")	
changeSortOrder(paste0("outputs/html/",forecastYear,"_current_probability_winning_by_seat.html"))
addColour(paste0("outputs/html/",forecastYear,"_current_probability_winning_by_seat.html"))
	
### Simplify for Webpage Table

if (forecastYear == 2015){

	predwin.web <- predwin.pretty[,c("Conservatives","Labour","Liberal Democrats","SNP","Plaid Cymru","Greens","UKIP","Other","Seat","Region","Win10")]
	names(predwin.web)[11] <- "2010"	
	
	### Sortable HTML table
	sortable.html.table(predwin.web,
		output.directory = "outputs/website/2015/tables/",
		output.file = paste0("predicted_probability_by_seat.html"),
		page.title = "Predicted probability of winning seat")

	changeSortOrder(paste0("outputs/website/2015/tables/","predicted_probability_by_seat.html"))
	addColour(paste0("outputs/website/2015/tables/","predicted_probability_by_seat.html"))

	nowwin.web <- nowwin.pretty[,c("Conservatives","Labour","Liberal Democrats","SNP","Plaid Cymru","Greens","UKIP","Other","Seat","Region","Win10")]
	names(nowwin.web)[11] <- "2010"	
	
	### Sortable HTML table
	sortable.html.table(nowwin.web,
		output.directory = "outputs/website/2015/tables/",
		output.file = paste0("current_probability_by_seat.html"),
		page.title = "Current probability of winning seat")				
	changeSortOrder(paste0("outputs/website/2015/tables/","current_probability_by_seat.html"))
	addColour(paste0("outputs/website/2015/tables/","current_probability_by_seat.html"))

	### Create tables of most likely gains and losses for each party
	
	GainLossTables <- list()
	threshold <- 0.1
	tmp <- 1-rowSums(prevvote)
	prevvote <- cbind(prevvote,replace(tmp,tmp < 0,0))
	
	for (i in 1:8){
		GainLossTables[[short.party.names[i]]] <- list()
		
		lasttmp <- prevwin[,i] == 1
		
		potgainstmp <- which((predwin[,i] > threshold) & (!lasttmp))	
		if (length(potgainstmp) > 0){
			tmptable <- data.frame(Seat=predwin.pretty$Seat[potgainstmp],Region=predwin.pretty$Region[potgainstmp],pGain=predwin[potgainstmp,i],From=short.party.names[apply(matrix(prevwin[potgainstmp,],ncol=8),1,which.max)],Margin10=NA,Margin15=NA, predvote.web[potgainstmp ,1:8])
			for (j in potgainstmp){
				tmptable$Margin10[which(j == potgainstmp)] <- round(-1*(prevvote[j,which.max(prevwin[j,])] - prevvote[j,i])*100)
				tmptable$Margin15[which(j == potgainstmp)] <- round(-1*((predvote.bar[j,which.max(prevwin[j,])] - predvote.bar[j,i]))*100)
				}	
			tmptable <- tmptable[order(tmptable$pGain,decreasing=TRUE),]
			tmptable$pGain <- paste0(round(tmptable$pGain*100),"%")
			tmptable$pGain <- replace(tmptable$pGain,tmptable$pGain == "100%","99%")
			rownames(tmptable) <- NULL
			colnames(tmptable) <- c("Seat","Region","Prob","From","2010 Margin","2015 Margin", short.party.names)
			GainLossTables[[short.party.names[i]]]$Gains <- tmptable
		}
		
		
		### Sortable HTML table
		sortable.html.table(GainLossTables[[short.party.names[i]]]$Gains,
		output.directory = "outputs/website/2015/tables/",
		output.file = paste0(short.party.names[i],"_seat_gains.html"),
		page.title = paste(party.names[i],"Seat Gains"))
		changeSortOrder(paste0("outputs/website/2015/tables/",short.party.names[i],"_seat_gains.html"))
		
	
		atrisktmp <- which((predwin[,i] < (1-threshold)) & (lasttmp))	
		if (length(atrisktmp) > 0){	
			tmptable <- data.frame(Seat=predwin.pretty$Seat[atrisktmp],Region=predwin.pretty$Region[atrisktmp], pLoss =1-predwin[atrisktmp,i],To=(short.party.names[-i])[apply(matrix(predwin[atrisktmp,-i],ncol=7),1,which.max)],Margin10=NA,Margin15=NA,predvote.web[atrisktmp ,1:8])
			for (j in atrisktmp){
				targetter <- which(tmptable$To[which(j == atrisktmp)] == short.party.names)
				tmptable$Margin10[which(j == atrisktmp)] <- round(-1*(prevvote[j, targetter] - prevvote[j,i])*100)
				tmptable$Margin15[which(j == atrisktmp)] <- round(-1*((predvote.bar[j, targetter] - predvote.bar[j,i]))*100)
				}	
			tmptable <- tmptable[order(tmptable$pLoss,decreasing=TRUE),]
			tmptable$pLoss <- paste0(round(tmptable$pLoss*100),"%")	
			tmptable$pLoss <- replace(tmptable$pLoss,tmptable$pLoss == "100%","99%")	
			rownames(tmptable) <- NULL
			colnames(tmptable) <- c("Seat","Region","Prob","To","2010 Margin","2015 Margin",short.party.names)
			GainLossTables[[short.party.names[i]]]$Losses <- tmptable
		}
		
		### Sortable HTML table
		sortable.html.table(GainLossTables[[short.party.names[i]]]$Losses,
		output.directory = "outputs/website/2015/tables/",
		output.file = paste0(short.party.names[i],"_seat_losses.html"),
		page.title = paste(party.names[i],"Seat Losses"))
		changeSortOrder(paste0("outputs/website/2015/tables/",short.party.names[i],"_seat_losses.html"))
		
	}

}

## Compute predicted seats by region

reg.df <- merge(predwin,canonicalseatvars,
	by.x = "row.names",by.y = "RefNo",
	all.x=T,all.y=F)
regional.predwin <- matrix(NA,nlevels(reg.df$Region),8)
colnames(regional.predwin) <- short.party.names
rownames(regional.predwin) <- levels(reg.df$Region)
for (i in 1:nlevels(reg.df$Region)){
	regional.predwin[i,] <- tabulate(apply(reg.df[,party.names],1,which.max)[as.numeric(reg.df$Region) == i],8)
	
	}

national.predwin <- rbind(colSums(regional.predwin[!is.element(rownames(regional.predwin),c("Scotland","Wales")),]),regional.predwin[is.element(rownames(regional.predwin),c("Scotland","Wales")),])
rownames(national.predwin)[1] <- "England"

## Compute predicted seats by 2010 result

retain.df <- merge(predwin,canonicalseatvars,
	by.x = "row.names",by.y = "RefNo",
	all.x=T,all.y=F)
retain.df$Win10 <- factor(retain.df$Win10,levels=party.names)	
retain.df$Win10[which(retain.df$Seat == "Brighton, Pavilion")] <- "Greens"
retain.predwin <- matrix(NA,length(short.party.names),length(short.party.names))
colnames(retain.predwin) <- short.party.names
rownames(retain.predwin) <- paste("2010",short.party.names)
for (i in 1:length(short.party.names)){
	retain.predwin[i,] <- tabulate(apply(retain.df[,party.names],1,which.max)[as.numeric(retain.df$Win10) == i],length(short.party.names))	
	}

### ############################################################################
### ### * output (6): trajectory of party support
### ############################################################################

### Convert the object to something we can use for plotting

### Create vector of party colours and define function for making transparent	
party.colours <- c("#0087DC","#D50000","#FDBB30","#FFFF00","#008142","#99CC33","#70147A","#DDDDDD")

fade.colour <- function(col,alpha){
	components <- col2rgb(col)
	return(rgb(components[1,]/255,components[2,]/255,components[3,]/255,alpha=alpha))
	}	

### Translate party colours to Munsell space, take slice, and convert back to hex
### This helps for choropleths
party.rgb <- col2rgb(party.colours)
party.rgb <- t(party.rgb)
party.rgb <- as.numeric(party.rgb)

party.munsells <- rgb2mnsl(matrix(party.rgb,ncol=3))
party.munsells.major <- gsub(" .*","",party.munsells)

party.colours.lo <- vector("numeric",8)
party.colours.hi <- vector("numeric",8)

for (i in 1:8) {
	tmp <- munsell:::munsell.map[which(munsell:::munsell.map$hue == party.munsells.major[i]),]
	### Exclude 'N'
	tmp <- subset(tmp,!grepl("N ",tmp$name))
	### exclude those without a name
	tmp <- subset(tmp,!is.na(tmp$name))
	### 
	if (nrow(tmp)==0) {
		party.colours.lo[i] <- "#777777"
		party.colours.hi[i] <- "#DDDDDD"

	} else {
		party.colours.lo[i] <- tmp$hex[which.min(tmp$chroma)]
		party.colours.hi[i] <- tmp$hex[which.max(tmp$chroma)]
	}
}



poolout.df <- as.data.frame(do.call("rbind",poolout))
poolout.df$variable <- rep(names(poolout),each = 8)
poolout.df$party <- rep(party.names,times = 4)

### Melt, then cast
poolout.df <- melt(poolout.df)
names(poolout.df)[3] <- "period"
poolout.df$date <- as.numeric(gsub("\\D","",poolout.df$period))
poolout.df$date <- as.Date(poolout.df$date, origin = origin.date)
poolout.df <- cast(poolout.df,party+period+date~variable)

poolout.df$fillcolour <- party.colours[match(poolout.df$party,party.names)]
poolout.df$linecolour <- poolout.df$fillcolour 
poolout.df$linecolour[which(match(poolout.df$party,party.names) == 4)] <- "#00000044"

### Prepare the plot

poolplot <- ggplot(subset(poolout.df,date <= as.Date(choosepredictday,origin = origin.date)),
	aes(x=date,y = mean, ymax = hi, ymin = lo, group = party,fill=fillcolour,colour=linecolour)) + 
	geom_smooth(stat="identity") + 	
	theme_bw() + 
	scale_fill_identity() + 
	scale_colour_identity() + 
	scale_x_date("",limits=c(origin.date,as.Date(choosepredictday,origin = origin.date))) + 
	scale_y_continuous("Support (%)",labels=percent,limits=c(0,0.4)) + 
	theme(axis.text.x=element_text(angle=90)) +
	facet_wrap(~party)

poolplot <- watermark(poolplot)

### Print off in PDF, EPS, and SVG

pdf(file=paste0("outputs/pdf/",forecastYear,"_pool.pdf"))
print(poolplot)
dev.off()

svg(file=paste0("outputs/svg/",forecastYear,"_pool.svg"))
print(poolplot)
dev.off()

### Top three (combined)
top3.df <- subset(poolout.df,party %in% c("Labour","Conservatives","Liberal Democrats") & date <= as.Date(choosepredictday,origin = origin.date))
top3.df$party <- factor(top3.df$party)

top3poolplot <- ggplot(top3.df,
	aes(x=date,y = mean, ymax = hi, ymin = lo, group = party,fill=fillcolour,colour=linecolour)) + 
	geom_smooth(stat="identity") + 	
	theme_bw() + 
	scale_fill_identity("Party",guide="legend",labels = c("Conservatives","Labour","Liberal Democrats")) + 
	scale_colour_identity() + 
	scale_x_date("",limits=c(origin.date,as.Date(choosepredictday,origin = origin.date))) + 
	scale_y_continuous("Support (%)",labels=percent,limits=c(0,0.5)) +
	theme(axis.text.x=element_text(angle=90),legend.position = "bottom") 

top3poolplot <- watermark(top3poolplot)

### Print off in PDF, EPS, and SVG

pdf(file=paste0("outputs/pdf/",forecastYear,"_pool_top_3.pdf"))
print(top3poolplot)
dev.off()

svg(file=paste0("outputs/svg/",forecastYear,"_pool_top_3.svg"))
print(top3poolplot)
dev.off()


### Top four (combined)

top4.df <- subset(poolout.df,party %in% c("Labour","Conservatives","Liberal Democrats","UKIP") & date <= as.Date(choosepredictday,origin = origin.date))
top4.df$party <- factor(top4.df$party)

top4poolplot <- ggplot(top4.df,
	aes(x=date,y = mean, ymax = hi, ymin = lo, group = party,fill=fillcolour,colour=linecolour)) + 
	geom_smooth(stat="identity") + 	
	theme_bw() + 
	scale_fill_identity("Party",guide="legend",labels = c("Conservatives","UKIP","Labour","Liberal Democrats")) + 
	scale_colour_identity() + 
	scale_x_date("",limits=c(origin.date,as.Date(choosepredictday,origin = origin.date))) + 
	scale_y_continuous("Support (%)",labels=percent,limits=c(0,0.5)) +
	theme(axis.text.x=element_text(angle=90),legend.position = "bottom") 

top4poolplot <- watermark(top4poolplot)

### Print off in PDF, EPS, and SVG

pdf(file=paste0("outputs/pdf/",forecastYear,"_pool_top_4.pdf"))
print(top4poolplot)
dev.off()

svg(file=paste0("outputs/svg/",forecastYear,"_pool_top_4.svg"))
print(top4poolplot)
dev.off()


### Top 5

top5.df <- subset(poolout.df,party %in% c("Labour","Conservatives","Liberal Democrats","UKIP","Greens") & date <= as.Date(choosepredictday,origin = origin.date))
top5.df$party <- factor(top5.df$party)

top5poolplot <- ggplot(top5.df,
	aes(x=date,y = mean, ymax = hi, ymin = lo, group = party,fill=fillcolour,colour=linecolour)) + 
	geom_smooth(stat="identity") + 	
	theme_bw() + 
	scale_fill_identity("Party",guide="legend",labels = c("Conservatives","UKIP","Green","Labour","Liberal Democrats")) + 
	scale_colour_identity() + 
	scale_x_date("",limits=c(origin.date,as.Date(choosepredictday,origin = origin.date))) + 
	scale_y_continuous("Support (%)",labels=percent,limits=c(0,0.5)) +
	theme(axis.text.x=element_text(angle=90),legend.position = "bottom") 

top5poolplot <- watermark(top5poolplot)

### Print off in PDF, EPS, and SVG

pdf(file=paste0("outputs/pdf/",forecastYear,"_pool_top_5.pdf"))
print(top5poolplot)
dev.off()

svg(file=paste0("outputs/svg/",forecastYear,"_pool_top_5.svg"))
print(top5poolplot)
dev.off()

### Plots of Election Predictions

electionout.df <- as.data.frame(do.call("rbind",electionout))
electionout.df$variable <- rep(names(electionout),each = 8)
electionout.df$party <- rep(party.names,times = 4)

### Melt, then cast
electionout.df <- melt(electionout.df)
names(electionout.df)[3] <- "period"
electionout.df$date <- as.numeric(gsub("\\D","",electionout.df$period))
electionout.df$date <- as.Date(electionout.df$date, origin = origin.date)
electionout.df <- cast(electionout.df,party+period+date~variable)

	electionout.df$fillcolour <- party.colours[match(electionout.df$party,party.names)]
	electionout.df$linecolour <- electionout.df$fillcolour 
	electionout.df$linecolour[which(match(electionout.df$party,party.names) == 4)] <- "#00000044"
	

### Get the historical data
		lagdata <- historicaldata[1,-1]
		historicaldata <- data.frame(value = rep(NA,8),
		party = party.names,
		fillcolour = rep("#000000",8),
		linecolour = rep("#000000",8))
### Prepare the plot

electionplot <- ggplot(subset(electionout.df,date <= as.Date(choosepredictday,origin = origin.date)),
	aes(x=date,y = mean, ymax = hi, ymin = lo, group = party,fill=fillcolour,colour=linecolour)) + 
	geom_smooth(stat="identity") + 	
	theme_bw() + 
	geom_abline(data = historicaldata, aes(intercept = value, group = party, fill=fillcolour,colour=linecolour), linetype = 2, slope = 0) + 
	scale_fill_identity() + 
	scale_colour_identity() + 
	scale_x_date("",limits=c(origin.date,as.Date(choosepredictday,origin = origin.date))) + 
	scale_y_continuous("UK Vote Share (%)",labels=percent,limits=c(0,0.5)) + 
	theme(axis.text.x=element_text(angle=90)) +
	facet_wrap(~party)

electionplot <- watermark(electionplot)

### Print off in PDF, EPS, and SVG

pdf(file=paste0("outputs/pdf/",forecastYear,"_election.pdf"))
print(electionplot)
dev.off()

svg(file=paste0("outputs/svg/",forecastYear,"_election.svg"))
print(electionplot)
dev.off()


### Top three (combined)

top3.df <- subset(electionout.df,party %in% c("Labour","Conservatives","Liberal Democrats") & date <= as.Date(choosepredictday,origin = origin.date))
top3.df$party <- factor(top3.df$party)
top3hist <- subset(historicaldata,party %in% c("Labour","Conservatives","Liberal Democrats"))
top3hist$party <- factor(top3hist$party)


top3electionplot <- ggplot(top3.df,
	aes(x=date,y = mean, ymax = hi, ymin = lo, group = party,fill=fillcolour)) + 
	geom_smooth(aes(colour= linecolour),stat="identity") + 	
	theme_bw() + 
	geom_abline(data = top3hist, aes(intercept = value, group = party), colour="#000000", linetype = 2, slope = 0) +
	scale_fill_identity("Party",guide="legend",labels = c("Conservatives","UKIP","Labour","Liberal Democrats")) + 
	scale_colour_identity(guide="none") + 
	scale_x_date("",limits=c(origin.date,as.Date(choosepredictday,origin = origin.date))) + 
	scale_y_continuous("UK Vote Forecast (%)",labels=percent,limits=c(0,0.5)) +
	theme(axis.text.x=element_text(angle=90), legend.position = "bottom")  
	
top3electionplot <- watermark(top3electionplot)

### Print off in PDF, EPS, and SVG

pdf(file=paste0("outputs/pdf/",forecastYear,"_election_top_3.pdf"))
print(top3electionplot)
dev.off()

svg(file=paste0("outputs/svg/",forecastYear,"_election_top_3.svg"))
print(top3electionplot)
dev.off()


### Top four (combined)

top4.df <- subset(electionout.df,party %in% c("Labour","Conservatives","Liberal Democrats","UKIP") & date <= as.Date(choosepredictday,origin = origin.date))
top4.df$party <- factor(top4.df$party)
top4hist <- subset(historicaldata,party %in% c("Labour","Conservatives","Liberal Democrats","UKIP"))
top4hist$party <- factor(top4hist$party)

top4electionplot <- ggplot(top4.df,
	aes(x=date,y = mean, ymax = hi, ymin = lo, group = party,fill=fillcolour)) + 
	geom_smooth(aes(colour= linecolour),stat="identity") + 	
	theme_bw() + 
	geom_abline(data = top4hist, aes(intercept = value, group = party), colour="#000000", linetype = 2, slope = 0) +
	scale_fill_identity("Party",guide="legend",labels = c("Conservatives","UKIP","Labour","Liberal Democrats")) + 
	scale_colour_identity(guide="none") + 
	scale_x_date("",limits=c(origin.date,as.Date(choosepredictday,origin = origin.date))) + 
	scale_y_continuous("UK Vote Forecast (%)",labels=percent,limits=c(0,0.5)) +
	theme(axis.text.x=element_text(angle=90),legend.position = "bottom") 

top4electionplot <- watermark(top4electionplot)

### Print off in PDF, EPS, and SVG

pdf(file=paste0("outputs/pdf/",forecastYear,"_election_top_4.pdf"))
print(top4electionplot)
dev.off()

svg(file=paste0("outputs/svg/",forecastYear,"_election_top_4.svg"))
print(top4electionplot)
dev.off()


top5.df <- subset(electionout.df,party %in% c("Labour","Conservatives","Liberal Democrats","UKIP","Greens") & date <= as.Date(choosepredictday,origin = origin.date))
top5.df$party <- factor(top5.df$party)
top5hist <- subset(historicaldata,party %in% c("Labour","Conservatives","Liberal Democrats","UKIP","Greens"))
top5hist$party <- factor(top5hist$party)

top5electionplot <- ggplot(top5.df,
	aes(x=date,y = mean, ymax = hi, ymin = lo, group = party,fill=fillcolour)) + 
	geom_smooth(aes(colour= linecolour),stat="identity") + 	
	theme_bw() + 
	geom_abline(data = top5hist, aes(intercept = value, group = party), colour="#000000", linetype = 2, slope = 0) +
	scale_fill_identity("Party",guide="legend",labels = c("Conservatives","UKIP","Green","Labour","Liberal Democrats")) + 
	scale_colour_identity(guide="none") + 
	scale_x_date("",limits=c(origin.date,as.Date(choosepredictday,origin = origin.date))) + 
	scale_y_continuous("UK Vote Forecast (%)",labels=percent,limits=c(0,0.5)) +
	theme(axis.text.x=element_text(angle=90),legend.position = "bottom") 

top5electionplot <- watermark(top5electionplot)

### Print off in PDF, EPS, and SVG

pdf(file=paste0("outputs/pdf/",forecastYear,"_election_top_5.pdf"))
print(top5electionplot)
dev.off()

svg(file=paste0("outputs/svg/",forecastYear,"_election_top_5.svg"))
print(top5electionplot)
dev.off()




## create seat forecast history plots

if (forecastYear == 2015){

	load(paste0("working/polling/",forecastYear,"seatout.RData"))

	### Update with today's forecast
	
	seatout$mean[, choosepredictday] <- outcomes.df$Mean
	seatout$hi[, choosepredictday] <- outcomes.df$Hi
	seatout$lo[, choosepredictday] <- outcomes.df$Lo
	
	### Fill in NAs from last few days, if found
	
	if (diff(rev(rev(which(!is.na(seatout$mean[1,])))[1:2])) > 1){
		tmpday <- rev(which(!is.na(seatout$mean[1,])))[2] + 1
		} else {
		tmpday <- choosepredictday
		}	
	
	while(tmpday < choosepredictday){
		seatout$mean[, tmpday] <- seatout$mean[, choosepredictday]
		seatout$hi[, tmpday] <- seatout$hi[, choosepredictday] 
		seatout$lo[, tmpday] <- seatout$lo[, choosepredictday]
		tmpday <- tmpday + 1
		}
		
	save(seatout,file=paste0("working/polling/",forecastYear,"seatout.RData"))	
		
	### Plots of Seat Predictions
	
	seatout.df <- as.data.frame(do.call("rbind",seatout))
	seatout.df$variable <- rep(names(seatout),each = 8)
	seatout.df$party <- rep(party.names,times = 4)
	
	### Melt, then cast
	seatout.df <- melt(seatout.df)
	names(seatout.df)[3] <- "period"
	seatout.df$date <- as.numeric(gsub("\\D","",seatout.df$period))
	seatout.df$date <- as.Date(seatout.df$date, origin = origin.date)
	seatout.df <- cast(seatout.df,party+period+date~variable)
	
	seatout.df$fillcolour <- party.colours[match(seatout.df$party,party.names)]
	seatout.df$linecolour <- seatout.df$fillcolour 
	seatout.df$linecolour[which(match(seatout.df$party,party.names) == 4)] <- "#00000044"
	
	### Prepare the plot
	
	seatplot <- ggplot(subset(seatout.df,date <= as.Date(choosepredictday,origin = origin.date)),
		aes(x=date,y = mean, ymax = hi, ymin = lo, group = party,fill=fillcolour,colour=linecolour)) + 
		geom_smooth(stat="identity") + 	
		theme_bw() + 
		scale_fill_identity() + 
		scale_colour_identity() + 
		scale_x_date("",limits=c(origin.date,as.Date(choosepredictday,origin = origin.date))) + 
		scale_y_continuous("GB Seat Forecast",limits=c(0,375)) + 
		theme(axis.text.x=element_text(angle=90)) +
		facet_wrap(~party)
	
	seatplot <- watermark(seatplot)
	
	### Print off in PDF, EPS, and SVG

	pdf(file=paste0("outputs/pdf/2015_seat_forecasts.pdf"))
	print(seatplot)
	dev.off()
	
	svg(file=paste0("outputs/svg/2015_seat_forecasts.svg"))
	print(seatplot)
	dev.off()
	
	
	
	
	### Lab, Con, LD, SNP, UKIP
	
	top5.df <- subset(seatout.df,party %in% c("Labour","Conservatives","Liberal Democrats","SNP","UKIP") & date <= as.Date(choosepredictday,origin = origin.date))
	top5.df$party <- factor(top5.df$party)
	
	
	top5seatplot <- ggplot(top5.df,
		aes(x=date,y = mean, ymax = hi, ymin = lo, group = party,fill=fillcolour)) + 
		geom_smooth(aes(colour= linecolour),stat="identity") + 	
		theme_bw() + 
		scale_fill_identity("Party",guide="legend",labels = c("Conservatives","UKIP","Labour","Liberal Democrats","SNP")) + 
		scale_colour_identity(guide="none") + 
		scale_x_date("",limits=c(origin.date,as.Date(choosepredictday,origin = origin.date))) + 
		scale_y_continuous("GB Seat Forecast",limits=c(0,375)) +
		theme(axis.text.x=element_text(angle=90), legend.position = "bottom")  
		
	top5seatplot <- watermark(top5seatplot)
	
	### Print off in PDF, EPS, and SVG
	
	pdf(file=paste0("outputs/pdf/2015_seat_forecasts_top_5.pdf"))
	print(top5seatplot)
	dev.off()
	
	svg(file=paste0("outputs/svg/2015_seat_forecasts_top_5.svg"))
	print(top5seatplot)
	dev.off()
	
	

}



### ############################################################################
### ### * output (7): plot of house effects
### ############################################################################

houseout.df <- data.frame(mean=as.numeric(houseout$mean)*100,lo=as.numeric(houseout$lo)*100,hi=as.numeric(houseout$hi)*100)
houseout.df$party <- rep(party.names,times = ncol(houseout$mean))
houseout.df$company <- rep(colnames(houseout$mean),each = nrow(houseout$mean))
	
houseplot <- ggplot(houseout.df,aes(x=company,y=mean,ymax=hi,ymin=lo,group=party)) +
	geom_pointrange() + 
	scale_x_discrete("Polling company") + 
	scale_y_continuous("House effects (+/- %)") + 
	theme_bw() + 
	theme(axis.text.x=element_text(angle=90)) + 
	facet_wrap(~party)	

houseplot <- watermark(houseplot)

pdf(file=paste0("outputs/pdf/",forecastYear,"_house_effect_plot.pdf"))
print(houseplot)
dev.off()

svg(file=paste0("outputs/svg/",forecastYear,"_house_effect_plot.svg"))
print(houseplot)
dev.off()

### Pollster-specific forecasts

houseout2.df <- data.frame(mean=as.numeric(pollsterforecasts[,2,])*100,lo=as.numeric(pollsterforecasts[,1,])*100,hi=as.numeric(pollsterforecasts[,3,])*100)
houseout2.df$party <- rep(party.names,times = dim(pollsterforecasts)[3])
houseout2.df$company <- rep(dimnames(pollsterforecasts)[[3]],each = dim(pollsterforecasts)[1])
	
houseplot2 <- ggplot(houseout2.df,aes(x=company,y=mean,ymax=hi,ymin=lo,group=party)) +
	geom_pointrange() + 
	scale_x_discrete("Polling company") + 
	scale_y_continuous("House effects (+/- %)") + 
	theme_bw() + 
	theme(axis.text.x=element_text(angle=90)) + 
	facet_wrap(~party)	

houseplot2 <- watermark(houseplot)

pdf(file=paste0("outputs/pdf/",forecastYear,"_house_pred_plot.pdf"))
print(houseplot2)
dev.off()

svg(file=paste0("outputs/svg/",forecastYear,"_house_pred_plot.svg"))
print(houseplot2)
dev.off()

### Create Table of Ranges

pollsterforecasts.pretty <- round(pollsterforecasts,2)*100
pollsterforecasts.pretty <- apply(pollsterforecasts.pretty,c(1,3),function(x) paste(x[1],"-",x[3]))
rownames(pollsterforecasts.pretty) <- short.party.names
pollsterforecasts.pretty <- rbind(pollsterforecasts.pretty, averageSizeByPollster)
rownames(pollsterforecasts.pretty)[nrow(pollsterforecasts.pretty)] <- "Avg N"
keeppollsters <- setdiff(colnames(pollsterforecasts.pretty),c("YouGov","Ashcroft","Populus"))
pollsterforecasts.pretty <- pollsterforecasts.pretty[, keeppollsters]
colnames(pollsterforecasts.pretty) <- gsub(" 2","",colnames(pollsterforecasts.pretty),fixed=TRUE)


## Plots

## Load recent polls

forecast.date <- Sys.Date()
forecast.date <- forecast.date - origin.date

pollsub <- polls[polls$date%in%(forecast.date-14):forecast.date,]
pollsub[,-c(1:4)] <- t(apply(pollsub[,-c(1:4)],1,function(x)x/sum(x)))

pollsub $Company <- gsub(" 2","", pollsub $Company)

## Remove old house effects
houseout2.df <- houseout2.df[!houseout2.df$company%in%c("YouGov","Ashcroft","Populus"),]
houseout2.df$company <- gsub(" 2","",houseout2.df$company)

houseout2.df$party <- as.factor(houseout2.df$party)

svg(file=paste0("outputs/svg/",forecastYear,"_pollster_forecasts_combined.svg"),16,10)
par(mar=c(8,4,4,2),mfrow=c(2,4),cex=1.2)

for(p in party.names){
tmp <- houseout2.df[houseout2.df$party==p,]
tmp$company.numeric <- as.numeric(as.factor(tmp$company))
plot(tmp$company.numeric,tmp$mean, 
		col="white",
		pch=19,
		bty="n",
		main=p,
		ylim=c(0,40),
		xlab="",
		ylab="Pollster Forecast",
		xaxt="n")
segments(tmp$company.numeric,tmp$lo,y1=tmp$hi, col=party.colours[match(p,party.names)], lwd=3)
axis(1,at=tmp$company.numeric,labels=tmp$company, las=2)
pollsub$Company.numeric <- tmp$company.numeric[match(pollsub$Company,tmp$company)]
points(pollsub $Company.numeric, pollsub[,p]*100,pch="—",cex=0.8,col=alpha(party.colours[match(p,party.names)],0.75))

}
dev.off()


for(p in levels(houseout2.df$party)){

tmp <- houseout2.df[houseout2.df$party==p,]
tmp$company.numeric <- as.numeric(as.factor(tmp$company))

svg(file=paste0("outputs/svg/",forecastYear,"_pollster_forecasts_",short.party.names[match(p,party.names)],".svg"),4,5)
par(mar=c(8,4,4,2),mfrow=c(1,1),cex=1.2)
plot(tmp$company.numeric,tmp$mean, 
		col="white",
		pch=19,
		bty="n",
		main=p,
		ylim=c(0,40),
		xlab="",
		ylab="Pollster Forecast",
		xaxt="n")
segments(tmp$company.numeric,tmp$lo,y1=tmp$hi, col=party.colours[match(p,party.names)], lwd=3)
axis(1,at=tmp$company.numeric,labels=tmp$company, las=2)

pollsub$Company.numeric <- tmp$company.numeric[match(pollsub$Company,tmp$company)]
points(pollsub$Company.numeric, pollsub[,p]*100,pch="—",cex=0.8,col=alpha(party.colours[match(p,party.names)],0.75))
	dev.off()
}




### ############################################################################
### * output (8): maps of predicted winners and party support
### ############################################################################

### THIS SECTION OF CODE REMOVED
### WE CAN'T REDISTRIBUTE THE NECESSARY MAP FILES

### ###########################################################################
### * output (9): bar charts of party seat distributions
### ###########################################################################

seat.tally <- lapply(1:8,function(x){
	res <- apply(winner,2,function(y)sum(y==x))
	return(res)
})
seat.tally <- do.call("cbind",seat.tally)
seat.tally <- as.data.frame(seat.tally)
names(seat.tally) <- party.names
seat.tally$iter <- 1:nrow(seat.tally)


seat.tally <- melt(seat.tally,id.vars="iter")
seat.tally$party.number <- match(seat.tally$variable,party.names)
seat.tally$party.fill <- party.colours[seat.tally$party.number]
seathist <- ggplot(seat.tally,aes(x=value,fill=party.fill)) + 
	facet_wrap(~variable,scales="free") + 
	scale_fill_identity() + 
	geom_bar(aes(y=..density..),binwidth=1,origin=-0.5) + 
	theme_classic()

seathist <- watermark(seathist)

pdf(file=paste0("outputs/pdf/",forecastYear,"_seat_histogram.pdf"))
print(seathist)
dev.off()

svg(file=paste0("outputs/svg/",forecastYear,"_seat_histogram.svg"))
print(seathist)
dev.off()

### ###########################################################################
### * output (9)(b): bar charts of party seat distributions
### ###########################################################################

seat.tally2 <- ddply(seat.tally,.(variable),function(df) {
	data.frame(mean = mean(df$value),
		lo = quantile(df$value,probs = predalpha/2),
		hi = quantile(df$value,probs = 1-predalpha/2),
		party.fill = unique(df$party.fill),
		party.number = unique(df$party.number))
})

watermark(ggplot(seat.tally2,aes(x = variable, y= mean, ymax = hi, ymin = lo, fill = party.fill)) +
	geom_bar(stat="identity") + 
	geom_errorbar(width=0.5) + 
	scale_fill_identity(guide="none") + 
	theme_bw())
	


	

### ###########################################################################
### * output (10): a quick report 
### ###########################################################################

round1 <- function(x) formatC(x,digits=1,format="f")
round2 <- function(x) formatC(x,digits=2,format="f")

if (forecastYear == 2010) lagseats <- c(210,349,62,6,2,0,0,3)
if (forecastYear == 2015) lagseats <- c(306,258,57,6,3,1,0,1)

pretty.seats.df <- outcomes.df[,c(1,3,2,4)]
pretty.seats.df[,2:4]<- round(pretty.seats.df[,2:4])
pretty.seats.df[,5] <- pretty.seats.df[,3] - lagseats
colnames(pretty.seats.df) <- c("Party","Lo","Seats","Hi","Swing")

pretty.votes.df <- pretty.seats.df 
pretty.votes.df[,2] <- paste0(round1(100*apply(predictday[,1,],1,quantile, predalpha/2)),"%")
pretty.votes.df[,3] <- paste0(round1(100*apply(predictday[,1,],1,mean)),"%")
pretty.votes.df[,4] <- paste0(round1(100*apply(predictday[,1,],1,quantile,1-predalpha/2)),"%")
pretty.votes.df[,5] <- paste0(round1(100*apply(predictday[,1,],1,mean) - 100*as.numeric(c(lagdata[1:7],1-sum(lagdata[1:7])))),"%")
colnames(pretty.votes.df) <- c("Party","Lo","Votes","Hi","Swing")

nowcastposterior <- apply(devposterior,3,"%*%", turnoutweights)/632

pretty.current.df <- cbind(current.df[c(1,3,2,4)], current.df[,2] - lagseats,
	paste0(round1(100*apply(nowcastposterior,1,quantile, predalpha/2)),"%"),
	paste0(round1(100*apply(nowcastposterior,1,mean)),"%"),
	paste0(round1(100*apply(nowcastposterior,1,quantile,1-predalpha/2)),"%"),
	paste0(round1(100*apply(nowcastposterior,1,mean) - 100*as.numeric(c(lagdata[1:7],1-sum(lagdata[1:7])))),"%")
	)
colnames(pretty.current.df) <- c("Party","Lo","Seats","Hi","Swing","Lo","Votes","Hi","Swing")	
	

shipping.majority <- c(mean(seats.by.iter[,1] > 325),mean(seats.by.iter[,2] > 325))
shipping.plurality <- c(mean(seats.by.iter[,1] > seats.by.iter[,2]),mean(seats.by.iter[,2] > seats.by.iter[,1]))
names(shipping.majority) <- names(shipping.plurality) <- c("Conservatives","Labour")

