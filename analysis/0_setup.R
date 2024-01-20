if(!require("tidyverse")) install.packages("tidyverse")
if(!require("devtools")) install.packages("devtools")
if(!require("IVODE")) devtools::install_github("GBarnsley/IVODE")
if(!require("epimixr")) devtools::install_github("sbfnk/epimixr")
if(!require("rmarkdown")) install.packages("rmarkdown")

date_start <- as_date("2000-01-01") #also t = 0

date_crisis_start <- as_date("2023-10-07")

date_projection_start <- as_date("2024-02-07")

date_projection_end <- date_projection_start + ((6/12) * 365)

#add infection and disease
#carmen code review
#send to francesco

library(IVODE)
library(tidyverse)
library(epimixr)
