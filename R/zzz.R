#  zzz.R
#
# This function is simply copied from mice package.


# on attach CDM
.onAttach <- function(libname,pkgname){
  d <- packageDescription("BIFIEsurvey")
  packageStartupMessage("|---------------------------------------------------------",
		   "--------\n"  ,
		paste("| " ,d$Package," " , d$Version," (",d$Date,")",sep="") ,
		"                                       " , 
        "\n| Maintainer: Alexander Robitzsch <a.robitzsch at bifie.at >      " ,
		"\n| http://www.bifie.at                                             ",
		"\n|---------------------------------------------------" ,
		"--------------\n" )
	}
version <- function(pkg="BIFIEsurvey"){
  lib <- dirname(system.file(package = pkg))
  d <- packageDescription(pkg)
  return(paste(d$Package,d$Version,d$Date,lib))
}
