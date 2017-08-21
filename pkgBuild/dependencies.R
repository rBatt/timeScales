
update_dependencies <- function(){
	devtools::use_package("data.table", type="Depends") # Basis for handling all data sets
	
	devtools::use_package("rbLib", type="Imports")
	devtools::use_package("stats", type="Imports")
	devtools::use_package("methods", type="Imports")
	devtools::use_package("zoo", type="Imports")
	devtools::use_package("forecast", type="Imports")	
}