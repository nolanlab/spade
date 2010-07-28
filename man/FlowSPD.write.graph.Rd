\name{FlowSPD.write.graph}
\alias{FlowSPD.write.graph}
\title{
Writing the graph to a file in some format
}
\description{
General function for exporting graphs to foreign file formats, however at present only the GML format is implemented.
}
\usage{
FlowSPD.write.graph(graph, file = "", format = c("gml"))
}
\arguments{
  \item{graph}{
The graph to export
}
  \item{file}{
A connection or a string giving the file name to write the graph to.
}
  \item{format}{
Character string giving the file format.
}
}
\details{
GML is general textual format for graphs.

The vertex and edge attributes are written to the file if they are numeric or strings. Currently only the graphics struct is supported, and only for vertices; \code{graphics.x} indicates an \code{x} attribute in the \code{graphics} struct.
}
\value{
A NULL, invisibly
}
\references{
}
\author{
Michael Linderman
}
\note{
}
\seealso{
\code{\link{write.graph}}
}
\examples{
}

