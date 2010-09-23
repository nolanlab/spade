# runSPADE:  R wrapper script for SPADE tree construction
# Erin Simonds - esimonds@stanford.edu
# Version 2.0 - August 2, 2010
#
# Command line instructions:
#   1) Customize the parameters that are CAPITALIZED below 
#   2) In a bash shell, navigate to the folder containing this script and the FCS file to be analyzed
#   3) At the bash command line, run: R -f runSPADE.R
#
# Interactive instructions:
#   1) Customize the parameters that are CAPITALIZED below 
#   2) In an interactive R session, change the working directory to a folder containing this script and the FCS file to be analyzed
#   3) At the R command line, run: source("runSPADE.R")
#
#
#
#

# Set this to the FCS file or folder to be processed.  If a folder is specified, all FCS files will be processed..
FILE_TO_PROCESS="FCS_files/"

# Set this to the markers used for clustering and tree construction.  Usually these are the surface markers.
SURFACE_MARKERS=c("Cell_length", "DNA-1(Ir191)D", "110:114(Merged)", "CD123(Eu151)D", "CD14(Gd160)D", "CD20(Sm147)D", "CD33(Nd148)D", "CD4(Nd145)D", "CD45(In115)D", "CD7(Yb176)D", "HLADR(Yb174)D", "IgM(Yb171)D")

# Set this to the markers used for fold-change calculations.  Usually these are phospho-specific markers.
FUNCTIONAL_MARKERS=c("pAkt(Sm152)D", "pBTK(Er166)D", "pErk1_2(Er168)D", "pLat(Er170)D", "pNFkb(Nd142)D", "pp38(Nd144)D", "pPLCg2(Er167)D", "pS6(Yb172)D", "pShp2(Sm154)D", "pSlp76(Dy164)D", "pStat1(Eu153)D", "pStat3(Gd158)D", "pStat5(Nd150)D", "pZap70(Gd156)D")

# Set this to all markers for which you want the basal medians.
ALL_MARKERS=c("Cell_length", "DNA-1(Ir191)D", "110:114(Merged)", "CD123(Eu151)D", "CD14(Gd160)D", "CD20(Sm147)D", "CD33(Nd148)D", "CD4(Nd145)D", "CD45(In115)D", "CD7(Yb176)D", "HLADR(Yb174)D", "IgM(Yb171)D", "pAkt(Sm152)D", "pBTK(Er166)D", "pErk1_2(Er168)D", "pLat(Er170)D", "pNFkb(Nd142)D", "pp38(Nd144)D", "pPLCg2(Er167)D", "pS6(Yb172)D", "pShp2(Sm154)D", "pSlp76(Dy164)D", "pStat1(Eu153)D", "pStat3(Gd158)D", "pStat5(Nd150)D", "pZap70(Gd156)D")

# Set this to the path to your local flowSPD package installation -- if using a global installation, set to NULL.
LIBRARY_PATH="lib/"

# Set this to the reference file to be used for fold-change calculations.
REFERENCE_FILE="Unstimulated.fcs"

# Set this to the desired number of events remaining after downsampling files
DOWNSAMPLED_EVENTS=50000

# Set this to the target number of clusters.  Algorithm will create clusters within 50% of this value. 
TARGET_CLUSTERS=200

# Path to output directory -- you probably don't need to change this.
OUTPUT_DIR="output/"

# Path to temporary directory -- you probably don't need to change this.
TMPDIR="/tmp/"


###################  No need to modify anything below this point ########################

#Run SPADE analysis workflow

args = commandArgs();

library("flowSPD",lib.loc=LIBRARY_PATH)

if (!exists('FILE_TO_PROCESS')){
	FILE_TO_PROCESS <- args[4]
}

FlowSPD.driver(FILE_TO_PROCESS, file_pattern="*.fcs", out_dir=OUTPUT_DIR, cluster_cols=SURFACE_MARKERS, median_cols=ALL_MARKERS, reference_file=REFERENCE_FILE, fold_cols=FUNCTIONAL_MARKERS, downsampling_samples=DOWNSAMPLED_EVENTS, k=TARGET_CLUSTERS)
FlowSPD.plot.trees(OUTPUT_DIR,file_pattern="*fcs*gml",out_dir=OUTPUT_DIR)
