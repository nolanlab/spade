#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#
# Add density parameter to FCS file specificied on commandline

# Path to flowSPD package installation
PACKAGE=lib/
# Path to temporary directory
TMPDIR=tmp/
#List of markers in double quotes surrounded by single quotes 
MARKERS='"115-CD45", "110-CD3", "111-CD3", "112-CD3", "114-CD3", "139-CD45RA", "142-CD19", "144-CD11b", "145-CD4", "146-CD8", "148-CD34", "147-CD20", "CD123", "158-CD33", "CD38", "CD90"'


INFILE=$1
OUTFILE=$1.density.fcs

SCRIPT="library(\"flowSPD\",lib.loc=\"${PACKAGE}\");\
cols <- c(${MARKERS});\
FlowSPD.addDensityToFCS(\"${INFILE}\",\"${OUTFILE}\",cols=cols);"

R -e "$SCRIPT"
