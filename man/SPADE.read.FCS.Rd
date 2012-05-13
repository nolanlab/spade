\name{SPADE.read.FCS}
\alias{SPADE.read.FCS}
\title{
Read an FCS File
}
\description{
Read and FCS file, correcting problematic column names and applying compensation if requested. A wrapper around \link[flowCore]{read.FCS}.
}
\usage{
SPADE.read.FCS(file, comp = TRUE, verbose = FALSE, ...)
}
\arguments{
  \item{file}{
FCS filename
}
  \item{comp}{
boolean indicating whether to apply compensation
}
  \item{verbose}{
boolean indicating whether to not suppress warnings from upstream FCS file reader.
}
  \item{\dots}{
additional arguments are passed through to \link[flowCore]{read.FCS}
}
}
\details{
A wrapper around \link[flowCore]{read.FCS} that copies FCS parameter "short name" to description when none is supplied, and applies compensation matrix stored in the SPILL or SPILLOVER keyword.
}
\value{
	a \link[flowCore]{flowFrame} object.
}
\author{
Michael Linderman <michael.linderman@mssm.edu>
}
\seealso{
\code{\link[flowCore]{read.FCS}}
}
\examples{
	fcs <- SPADE.read.FCS(system.file(file.path("extdata","SimulatedRawData.fcs"), package = "spade"))
}
