if(!require("tidyverse")) install.packages("tidyverse")
if(!require("ggpubr")) install.packages("ggpubr")
if(!require("viridis")) install.packages("viridis")
if(!require("devtools")) install.packages("devtools")
if(!require("IVODE")) devtools::install_github("GBarnsley/IVODE")
if(!require("epimixr")) devtools::install_github("sbfnk/epimixr")
if(!require("socialmixr")) devtools::install_github("sbfnk/epimixr")
if(!require("rmarkdown")) install.packages("rmarkdown")

date_start <- as_date("2000-01-01") #also t = 0

date_crisis_start <- as_date("2023-10-07")

date_projection_start <- as_date("2024-02-07")

projection_months <- 6

date_projection_end <- date_projection_start + ((projection_months/12) * 365)

N <- 1000

library(IVODE)
library(tidyverse)
library(epimixr)
library(socialmixr)

#make directories
if(!dir.exists("data/output")) dir.create("data/output")
if(!dir.exists("data/derived")) dir.create("data/derived")
if(!dir.exists("plots")) dir.create("plots")