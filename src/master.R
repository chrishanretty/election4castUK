forecastYear <- 2015
upload <- F
fullrun <- T
optionname <- ""

### ###########################################################################
### Run through all processes
### ###########################################################################

if (fullrun) source("src/prepconpolls.R",echo=TRUE)
source("src/prepukpolls.R",echo=TRUE)
if (fullrun) source("src/genprior.R",echo=TRUE)
source("src/pool.R",echo=TRUE)
if (fullrun) source("src/mr.R",echo=TRUE)
source("src/qi.R",echo=TRUE)
