\name{SPADE.driver}
\alias{SPADE.driver}
\title{
SPADE workflow driver
}
\description{
A function to drive the SPADE workflow. Produces graphs annoated with parameter medians and fold change.
}
\usage{
SPADE.driver(files, file_pattern="*.fcs", out_dir=".", cluster_cols=NULL, panels=NULL, comp=TRUE, arcsinh_cofactor=NULL, transforms=flowCore::arcsinhTransform(a=0, b=0.2), downsampling_samples=20000, downsampling_exclude_pctile=0.01, downsampling_target_pctile=0.05, k=200, clustering_samples=50000, layout=igraph0:::layout.kamada.kawai, pctile_color=c(0.02,0.98)) 
}
\arguments{
  \item{files}{
Either a vector of FCS files, or a directory. If a directory, all of the *.fcs files in the directory are processed.
}
  \item{file_pattern}{
Wildcard pattern to match file if \code{files} is a director
}
  \item{out_dir}{
Directory where output files are written. Will be created if it does not exist.
}
	\item{cluster_cols}{
Usually a vector of strings specifying the columns to be used in the clustering, e.g., c("(Cd110)D","(Cs111)D"). Strings will be matched against the parameter names extracted from the FCS file. The default=NULL will use all parameters.

}
  \item{panels}{
List of panels for median and fold change calculations. See details for specific structure. If NULL, medians are computed for all markers in all files.
}
  \item{comp}{
Apply compensation matrix if present in SPILL or SPILLOVER keywords
}
	\item{arcsinh_cofactor}{
DEPRECATED. Cofactor used in arcsinh transform \code{asinh(data/arcsinh_cofactor)} of data.
}
	\item{transforms}{
Transform object to apply to data. A single transform object will be applied to all channels. To apply different transforms to specific channels use a named vector of transform objects (where names are parameter names). 
}
  \item{downsampling_samples}{
Desired number of samples remaining after downsampling files
}
  \item{downsampling_exclude_pctile}{
Numeric value in [0,1]. Densities below this percentile will be excluded during downsampling.
}
  \item{downsampling_target_pctile}{
Numeric value in [0,1]. Densities below this percentile, but above 'exclude_pctile' will be retained during downsampling. Only meaningful if 'downsampling_samples' is 'NULL'.
}
  \item{k}{
Desirec number of clusters. Algorithm might create between [k/2,3k/2] clusters.
}
  \item{clustering_samples}{
Desired number of samples to be used in clustering.
}
  \item{layout}{
Layout function
}
  \item{pctile_color}{
A two element vector specifying lower and upper percentiles that should be used to set the color scale. Values below and above these percentiles will be forced to the 'smallest' and 'largest' color respectively. Not in effect if 'scale' is specified. Relevant for downstream tools that used global value ranges produced by driver.
}
}
\details{
	The \code{panels} argument must be null or a list of panel descriptors, which are themselves lists containing at minimum a vector of panel files and median cols. An example minimum panels argument would be \code{ list( list(panel_files="basal.fcs", median_cols=NULL)) }. \code{panel_files} is a single file name or vector of file names in the experiment. \code{median_cols} is similar to the \code{cluster_cols} argument. Each panel descriptor can optionally specifiy reference_files and columns for fold change analysis. An example full panel descriptor would be \code{ list( list(panel_files=c("basal.fcs", "stim.fcs"), median_cols=NULL, reference_files="basal.fcs", fold_cols=("p1"))) }. \code{fold_cols} is similar to the \code{cluster_cols} argument. \code{reference_files} is a single file name or vector of file names in the experiment and in the \code{panel_files} for this experiment.
	\code{median_cols} and \code{fold_cols} are only interpreted in the context on their panel files, and so partially overlapping panels are possible. However, all the files specified within a panel must have the cluster, median and fold change parameters specified.
}
\value{
NULL
}
\author{
Michael Linderman
}
\examples{
	# Load two-parameters sample data included in package
	data_file_path = system.file(file.path("extdata","SimulatedRawData.fcs"), package = "spade")

	# Run basic SPADE analyses, clustering on two parameters. Annotated graphs will be
	# in output_dir. See SPADE.plot.trees to generate PDFs of annotated graphs.
	output_dir <- tempdir()
	SPADE.driver(data_file_path, out_dir=output_dir, cluster_cols=c("marker1","marker2"))
}
