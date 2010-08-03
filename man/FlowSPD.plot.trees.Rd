\name{FlowSPD.plot.trees}
\alias{FlowSPD.plot.trees}
\title{
Plot trees with annotated vertices
}
\description{
Plot trees for each vertex annotation setting vertex size and color based on the particular annotation.
}
\usage{
FlowSPD.plot.trees(files, file_pattern = "*.gml", out_dir = ".", layout = FlowSPD.layout.arch, attr_pattern = "median|fold")
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
