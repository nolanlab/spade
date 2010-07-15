\name{FlowSPD.downsampleFCS}
\alias{FlowSPD.downsampleFCS}
\title{
Downsample observcations in a FCS file according to density parameter
}
\description{
Downsample the observations in a FCS file according to a previously computed
density parameter. The goal is to produce a smaller set of observations with
similar density. Downsampling is independent of how the density is modeled. 
}
\usage{
FlowSPD.downsampleFCS(infilename, outfilename,
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
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
The name of the written file is returned
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Michael Linderman
}
\seealso{
\code{\link{FlowSPD.addDensityToFCS}}
}
%%\examples{
%%}
