library("data.table")

setwd("~/Documents/School&Work/epaPost/timeScales/pkgBuild")

# dat <- fread("../inst/extdata/Daily_SoS_All_Orgd&Cleaned_Data_2015.csv")
# sos_data <- dat[,list(Year, Lake, DoY, DateTime, Temp_HYLB, Chla_Conc_HYLB, BGA_Conc_HYLB)]
# save(sos_data, file="../data/sos_data.RData")

load("../inst/extdata/LRT_Squeal2_allYears_HighFreq_SondeData_MARSSinterps_v2.Rdata")
dat <- as.data.table(LRT_Squeal2_allYears_HighFreq_SondeData_MARSSinterps_v2)
sos_data <- dat[,list(Year, Lake, DoY, Temp_HYLB, Chl_HYLB, Chl_logged_HYLB, BGA_HYLB, BGA_logged_HYLB)]
sos_data[,Lake:=list(c("L"="Paul", "R"="Peter", "T"="Tuesday")[Lake])]
setorder(sos_data, Year, Lake, DoY)
save(sos_data, file="../data/sos_data.RData", compress='xz')