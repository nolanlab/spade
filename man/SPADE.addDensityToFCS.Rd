\name{SPADE.addDensityToFCS}
\alias{SPADE.addDensityToFCS}
\title{
Annotate FCS file with local density of each observation
}
\description{
Compute the local density of observation and incorporate the result as a new parameter
to the FCS file. The local density is modeled as an integer count of the number of other
observations within a specified distance of the observation.
}
\usage{
SPADE.addDensityToFCS(infilename, outfilename, cols = NULL, arcsinh_cofactor = 5, kernel_mult = 5, apprx_mult = 1.5, med_samples = 2000, comp=TRUE)
}
\arguments{
  \item{infilename}{
Name of the input FCS file
}
  \item{outfilename}{
Name of the output FCS file
}
  \item{cols}{
Usually a vector of strings specifying the columns to be used in the density
calculation, e.g., c("(Cd110)D","(Cs111)D"). Strings will be matched against
the parameter names extracted from the FCS file. The default=NULL will
use all parameters. 
}
  \item{arcsinh_cofactor}{
Cofactor used in arcsinh transform \code{asinh(data/arcsinh_cofactor)} of data
}
  \item{kernel_mult}{
Multiplier of the minimum median distance within which other observations are
counted towards the density
}
  \item{apprx_mult}{
Multiplier of the minimum median distance within which observations are
approximated to have the same density
}
  \item{med_samples}{
Number of observations used to estimate the minimum median distance
}
  \item{comp}{
Apply compensation matrix if present in SPILL or SPILLOVER keywords
}
}
\value{
The name of the written file is returned
}
\author{
Michael Linderman
}
\examples{
	# Load two-parameters sample data included in package
	data_file_path = paste(installed.packages()["spade","LibPath"],"spade","extdata","SimulatedRawData.fcs",sep=.Platform$file.sep)
	
	output_dir <- tempdir()
	
	# Compute and annotate FCS file with density
	density_file_path <- paste(output_dir,.Platform$file.sep,basename(data_file_path),".density.fcs",sep="")
	SPADE.addDensityToFCS(data_file_path, density_file_path, cols=c("marker1","marker2"))
}
