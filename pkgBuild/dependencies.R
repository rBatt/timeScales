
update_dependencies <- function(){
	devtools::use_package("data.table", type="Depends") # Basis for handling all data sets
	
	devtools::use_package("rbLib", type="Imports")
	devtools::use_package("stats", type="Imports")
	devtools::use_package("methods", type="Imports")
	devtools::use_package("zoo", type="Imports")
	devtools::use_package("forecast", type="Imports")	
	devtools::use_package("rootSolve", type="Imports")
	devtools::use_package("R2jags", type="Imports")
	devtools::use_package("fields", type="Imports")
	devtools::use_package("deSolve", type="Imports")
	devtools::use_package("phaseR", type="Imports")
	devtools::use_package("foreach", type="Imports")
	devtools::use_package("doParallel", type="Imports")
}