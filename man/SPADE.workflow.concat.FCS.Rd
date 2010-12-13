\name{SPADE.workflow.concat.FCS}
\alias{SPADE.workflow.concat.FCS}
\title{
SPADE workflow helper for computing medians on concatenated FCS files
}
\description{
A helper function for computing the medians of a set of FCS files concatenated together.
}
\usage{
SPADE.workflow.concat.FCS(files, file_pattern = "*.density.fcs.cluster.fcs", out_dir = ".", cols = NULL, layout = SPADE.layout.arch, arcsinh_cofactor = 5, in_graph_file = "mst.gml", out_graph_file = "concat.gml", cluster_cols=NULL)
}
\arguments{
  \item{files}{
Either a vector of FCS files, or a directory. If a directory, all of the *.fcs files in the directory are processed.
}
  \item{file_pattern}{
Wildcard pattern to match file if \code{files} is a directory
}
  \item{out_dir}{
Directory where output files are written. Will be created if it does not exist.
}
  \item{cols}{
Usually a vector of strings specifying the columns to compute medians for, e.g., c("115-CD45","110-CD3"). Strings will be matched against the parameter descriptions extracted from the FCS file. The default=NULL will use all parameters.
}
  \item{layout}{
Layout function or matrix of node positions
}
  \item{arcsinh_cofactor}{
Cofactor used in arcsinh transform \code{asinh(data/arcsinh_cofactor)} of data
}
  \item{in_graph_file}{
Input GML graph file. If \code{files} is a directory, \code{in_graph_file} will treated as a path relative to \code{files}.
}
  \item{out_graph_file}{
File name for annotated GML file output. Will be placed in \code{out_dir}.
}
  \item{cluster_cols}{
A vector of strings specifying parameters that should be marked as being used in clustering.
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
\code{\link{SPADE.driver}}
}
