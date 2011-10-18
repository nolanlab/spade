\name{SPADE.installPlugin}
\alias{SPADE.installPlugin}
\title{
Install CytoSPADE Cytoscape plugin
}
\description{
Install, CytoSPADE, the Cytoscape plugin for working with SPADE that is distributed with the SPADE R package. CytSPADE provides a GUI for setting-up SPADE analyses and interactively visualizing the results. 
}
\usage{
SPADE.installPlugin(cytoscape_path)
}
\arguments{
  \item{cytoscape_path}{
Path to your Cytoscape install, e.g., on OSX it is typically something like '/Applications/Cytoscape_v2.8.1'
}
}
\details{
Copies the Cytoscape plugin file distributed with the SPADE R package to the Cytoscape plugin directory.
}
\value{
Logical indicating success of the copy operation.
}
\author{
Michael Linderman
}
\examples{
# On OSX:
# SPADE.installPlugin("/Applications/Cytoscape_v2.8.1/")
}
