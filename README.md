# CRCGNPortal

Scripts for running interactive web portal for exploring CRCGN data

##Instructions for running the app

First you need to create "datadirs.json", see "datadirs.json\_example",
this file stores the directories of data - do not commit this to git (add to .gitignore).

Then, simply run the shiny server. 

###To run from command line:

```
Rscript source.R
```

See interactive app at url given by shiny output

###To run within R:
start up R
```
library(shiny)
runApp("../CRCGNPortal")
```

This will automatically open a web browser while app is running locally. 


