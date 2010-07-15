\name{FlowSPD.addDensityToFCS}
\alias{FlowSPD.addDensityToFCS}
\title{
Annotate FCS file with local density of each observation
}
\description{
Compute the local density of observation and incorporate the result as a new parameter
to the FCS file. The local density is modeled as an integer count of the number of other
observations within a specified distance of the observation.
}
\usage{
FlowSPD.addDensityToFCS(infilename, outfilename, cols = NULL, arcsinh_cofactor = 5, kernel_mult = 5, apprx_mult = 1.5, med_samples = 2000)
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
calculation, e.g., c("115-CD45","110-CD3"). Strings will be matched against
the parameter descriptions extracted from the FCS file. The default=NULL will
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
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The name of the written file is returned
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Michael Linderman
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
%%\examples{
%%}
