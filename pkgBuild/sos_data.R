library("data.table")

setwd("~/Documents/School&Work/epaPost/timeScales/pkgBuild")

dat <- fread("../inst/extdata/Daily_SoS_All_Orgd&Cleaned_Data_2015.csv")

sos_data <- dat[,list(Year, Lake, DoY, DateTime, Temp_HYLB, Chla_Conc_HYLB, BGA_Conc_HYLB)]

save(sos_data, file="../data/sos_data.RData")