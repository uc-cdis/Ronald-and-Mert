# install and soure R packages
install_r_prereqs <- function(){
	pkgTest("RCurl")
	pkgTest("devtools")
	pkgTest("RJSONIO")
	pkgTest("matR")
	if(!require("DESeq"))
	{
		source("https://bioconductor.org/biocLite.R")
		biocLite("DESeq")
	}
    install_github(repo="MG-RAST/matR", dependencies=FALSE, ref="early-release")
    dependencies()
}
pkgTest <- function(x)
{
    if (!require(x,character.only = TRUE))
    {
      	install.packages(x,dep=TRUE)
        if(!require(x,character.only = TRUE)) stop("Package not found")
    }
}
############################################################################################################################
