\name{SPADE.downsampleFCS}
\alias{SPADE.downsampleFCS}
\title{
Downsample observcations in a FCS file according to density parameter
}
\description{
Downsample the observations in a FCS file according to a previously computed
density parameter. The goal is to produce a smaller set of observations with
similar density. Downsampling is independent of how the density is modeled. 
}
\usage{
SPADE.downsampleFCS(infilename, outfilename,
		      exclude_pctile = 0.01, target_pctile = 0.05, 
		      desired_samples = NULL)
}
\arguments{
  \item{infilename}{
Name of the input FCS file. Must have a parameter named "density".
}
  \item{outfilename}{
Name of the output FCS file
}
  \item{exclude_pctile}{
Numeric value in [0,1]. Densities below this percentile will be excluded.
}
  \item{target_pctile}{
Numeric value in [0,1]. Densities below this percentile, but above
\option{exclude_pctile} will be retained. Only meaningful if
\option{desired_samples} is \code{NULL}.
}
  \item{desired_samples}{
Desired number of samples. If set to integer value, the target percentile will
be set internally to downsample to approximately the desired number of
samples.
}
}
\value{
The name of the written file is returned
}
\author{
Michael Linderman
}
\seealso{
\code{\link{SPADE.addDensityToFCS}}
}
\note{
Underlying implementations have been parallelized with OpenMP. Set
\code{OMP_NUM_THREADS} in environment to control the number of threads used.
}

\examples{
	# Not run
	## Load two-parameters sample data included in package
	#data_file_path = paste(installed.packages()["spade","LibPath"],"spade","extdata","SimulatedRawData.fcs",sep=.Platform$file.sep)

	#output_dir <- tempdir()
	#
	## Compute and annotate FCS file with density
	#density_file_path <- paste(output_dir,.Platform$file.sep,basename(data_file_path),".density.fcs",sep="")
	#SPADE.addDensityToFCS(data_file_path, density_file_path, cols=c("marker1","marker2"))

	## Downsample FCS file based on density
	#downsample_file_path <- paste(output_dir,.Platform$file.sep,basename(data_file_path),".density.fcs",sep="")
	#SPADE.downsampleFCS(density_file_path, downsample_file_path)
}
