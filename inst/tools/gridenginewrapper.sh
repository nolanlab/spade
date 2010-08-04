#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -j y
#
#
# gridenginewrapper.sh
# Wrapper for Grid Engine execution of R job on adamant.stanford.edu
# Erin Simonds - esimonds@stanford.edu
# v1.0 - August 4, 2010
# 
# Usage:
# qsub -m e -M esimonds@stanford.edu -pe threaded 16 gridenginewrapper.sh

# Set R_PATH to the path to R on your system
R_PATH='/share/apps/R/R-2.11.1/bin/R'


# The following line will reserve all available threads on the compute node, thus reserving it for your job
export OMP_NUM_THREADS=$NSLOTS

# cd '/home/esimonds/Rdata/100802 - CyTOF55 PBMC RELEASE2 singlets - stimulated ratio trees'
$R_PATH -f runSPADE.R
