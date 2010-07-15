#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#
# Run end-to-end FlowSPD for a single FCS file 

# Path to flowSPD package installation
PACKAGE=lib/
# Path to temporary directory
TMPDIR=/tmp/
#List of markers in double quotes surrounded by single quotes 
MARKERS='"115-CD45", "110-CD3", "111-CD3", "112-CD3", "114-CD3", "139-CD45RA", "142-CD19", "144-CD11b", "145-CD4", "146-CD8", "148-CD34", "147-CD20", "CD123", "158-CD33", "CD38", "CD90"'

# Input and output files
INFILE=$1

GRAPHFILE=$1.gml
CLUSTERFILE=$1.cluster
OUTFILE=$1.clustered.fcs

SCRIPT="library(\"flowSPD\",lib.loc=\"${PACKAGE}\");\
cols<-c(${MARKERS});\
FlowSPD.addDensityToFCS(\"${INFILE}\",\"${TMPDIR}${INFILE}.density.fcs\",cols=cols);\
FlowSPD.downsampleFCS(\"${TMPDIR}${INFILE}.density.fcs\",\"${TMPDIR}${INFILE}.downsample.fcs\");\
FlowSPD.FCSToTree(c(\"${TMPDIR}${INFILE}.downsample.fcs\"),\"${GRAPHFILE}\",\"${CLUSTERFILE}\",cols=cols);\
FlowSPD.addClusterToFCS(\"${TMPDIR}${INFILE}.density.fcs\",\"${OUTFILE}\",\"${CLUSTERFILE}\",cols=cols);"

R -e "$SCRIPT"

