###  - Pre-process polls (preppolls.R);
### 	Takes data from data/polling
###		Cleans in 
###		Spits it out in working/polling
###	 - First three variables are date, Company, SampleSize
###	 - followed by parties in established order 

### Libraries
library(car)


	
	elex.date <- as.Date("2015-05-07")
	origin.date <- elex.date - 365
	
	## Convert YouGov data to counts
	
	ygraw <- read.csv("data/polling/2015/2015YouGovdata.csv")	
	
	ygtmp <- data.frame(
		Conservatives=ygraw$Conservatives*ygraw$samplesize/100,
		Labour=ygraw$Labour*ygraw$samplesize/100,
		Liberal.Democrats=ygraw$Liberal.Democrats*ygraw$samplesize/100,
		SNP=ygraw$NationalPartiesinScotland*ygraw$samplesizeScotland/100,
		Plaid.Cymru=ygraw$NationalPartiesinWalesMidlands*ygraw$samplesizeWalesMidlands/100,
		Greens=ygraw$Greens*ygraw$samplesize/100,
		UKIP=ygraw$UKIP*ygraw$samplesize/100
		)	
	ygtmp$Other=ygraw$samplesize - rowSums(ygtmp)
	ygtmp <- round(ygtmp*0.75) # adjustment for don't know and refusals
		
	ygstart <- as.Date(ygraw$start,format="%d/%m/%Y") - origin.date 
	ygend <- as.Date(ygraw$end,format="%d/%m/%Y") - origin.date 	
	ygdate <- as.integer(ceiling((ygstart+ygend)/2))
	
	ygCompany <- as.character(ygraw$house)
	
	## Merge count data	
	
	counts <- read.csv("data/polling/2015/2015pollcounts.csv")	
	
	tmpstart <- as.Date(counts$start,format="%d/%m/%Y") - origin.date
	tmpend <- as.Date(counts $end,format="%d/%m/%Y") - origin.date
	tmpdate <- as.integer(ceiling((tmpstart + tmpend)/2))
	
	polls <- data.frame(
		startdate=c(ygstart, tmpstart),
		enddate=c(ygend, tmpend),
		date=c(ygdate,tmpdate),
		Company=c(as.character(ygCompany),as.character(counts$house)),
		Conservatives=c(ygtmp$Conservatives,counts$Conservatives),
		Labour=c(ygtmp$Labour,counts$Labour),
		Liberal.Democrats=c(ygtmp$Liberal.Democrats,counts$Liberal.Democrats),
		SNP=c(ygtmp$SNP,counts$SNP),
		Plaid.Cymru=c(ygtmp$Plaid.Cymru,counts$Plaid.Cymru),
		Greens=c(ygtmp$Greens,counts$Greens),
		UKIP=c(ygtmp$UKIP,counts$UKIP),
		Other=c(ygtmp$Other,counts$Other)
		)
		
	polls$Other <- replace(polls$Other,polls$Other < 0,0)	
		
	names(polls) <- car:::recode(names(polls),
		"'Liberal.Democrats'='Liberal Democrats';
		'Plaid.Cymru'='Plaid Cymru'")
		
	polls$Company <- relevel(polls$Company,"YouGov 2")  	
	
	save(polls,file="working/polling/2015pollcountdata.RData")


