if(!require("devtools")) install.packages("devtools")
if(!require("tidyverse")) install.packages("tidyverse")
if(!require("IVODE")) devtools::install_github("GBarnsley/IVODE")

date_start <- as_date("2000-01-01") #also t = 0

date_projection_start <- as_date("2023-10-07")

date_projection_end <- as_date("2024-02-01")

library(IVODE)
library(tidyverse)