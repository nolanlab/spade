\name{SPADE.driver}
\alias{SPADE.driver}
\title{
SPADE workflow driver
}
\description{
A function to drive the SPADE workflow. Produces graphs annoated with parameter medians and fold change.
}
\usage{
SPADE.driver(files, file_pattern="*.fcs", out_dir=".", cluster_cols=NULL, arcsinh_cofactor=5.0, layout=SPADE.layout.arch, median_cols=NULL, reference_files=NULL, fold_cols=NULL, downsampling_samples=20000, downsampling_exclude_pctile=0.01, downsampling_target_pctile=0.05, k=200, clustering_samples=50000, comp=TRUE, pctile_color=c(0.02,0.98)) 
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
  \item{arcsinh_cofactor}{
Cofactor used in arcsinh transform \code{asinh(data/arcsinh_cofactor)} of data
}
  \item{layout}{
Layout function
}
  \item{median_cols}{
Similar to \code{cluster_cols}, specifies columns for which to calculate medians. MST graphs are annotated with medians.
}
  \item{reference_files}{
Vector of reference FCS files to use in calculating the fold change. When multiple reference files are used, the observations in those files are concatenated prior to computing the reference medians. 
}
  \item{fold_cols}{
Similar to \code{cluster_cols}, specifies columns for which to calculate fold change. MST graphs are annotated with fold change.
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
  \item{comp}{
Apply compensation matrix if present in SPILL or SPILLOVER keywords
}
 \item{pctile_color}{
A two element vector specifying lower and upper percentiles that should be used to set the color scale. Values below and above these percentiles will be forced to the 'smallest' and 'largest' color respectively. Not in effect if 'scale' is specified. Relevant for downstream tools that used global value ranges produced by driver.
}
}
\details{
}
\value{
}
\references{
}
\author{
Michael Linderman
}
\note{
}
\seealso{
}
\examples{
}
