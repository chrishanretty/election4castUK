ElectionForecast.co.uk
===

This project contains an amended version of the source code used to generate the forecasts for electionforecast.co.uk

The code has been amended in three respects:

 - code preparing constituency polling, some of which was made available by YouGov, has been omitted
 - code mapping the forecasts has been omitted
 - code dealing with our forecasts for Northern Ireland has been omitted

In the first and second case, we have omitted code because we are unable to redistribute the underlying data (YouGov campaign data for the constituency polling; shapefiles for the mapping).

Although the project contains data on YouGov constituency subsamples (included alongside information from Ashcroft polls, and found in `working/constituency`), these have been reweighted by us to match Census demographics, and cannot be considered "raw data". 

We thank YouGov for giving us access to their data, and would also like to thank Anthony Wells for permission to include his data on historical vote intention, which can be found at ukpollingreport.co.uk

Usage
---

The file `master.R` should do everything. 


