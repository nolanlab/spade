\name{SPADE.plot.trees}
\alias{SPADE.plot.trees}
\title{
Plot trees with annotated vertices
}
\description{
Plot trees for each vertex annotation setting vertex size and color based on the particular annotation.
}
\usage{
SPADE.plot.trees(files, file_pattern = "*.gml", out_dir = ".", layout = SPADE.layout.arch, attr_pattern = "median|fold", pctile_keep=c(0.02,0.98), normalize=NULL)
}
\arguments{
  \item{files}{
Either a vector of gml graph files or a directory. If a directory, all of the files matching the \code{pattern} wildcard pattern are processesd.
}
  \item{file_pattern}{
Wildcard pattern to match files if \code{files} is a directory.
}
  \item{out_dir}{
Directory where output files are written. Will be created if it does not exist.
}
  \item{layout}{
Either a function or a numeric matrix specifying how vertices are placed on plot. If it is a matrix, the matrix must have two columns, x and y position, and as many rows as vertices. If \code{layout} is a function, it will be called with an igraph graph as the single parameter.
}
  \item{attr_pattern}{
A regular expression that matches the attributes that should be plotted for each graph.
}
  \item{pctile_keep}{
	  A two element vector specify lower and upper percentiles beyond which attribute values should be discarded.
  }
  \item{normalize}{
	  One of NULL or 'local'. If 'local', attributes are normalized on a per-file, per-attribute basis into the range [-1,1] with zero remaining zero (e.g., [-5,10] will normalize to [-.5,1]).
  }
}
\details{
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
