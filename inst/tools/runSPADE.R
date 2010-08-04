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
SURFACE_MARKERS=c("Cell Length", "191-DNA", "115-CD45", "139-CD45RA", "142-CD19", "144-CD11b", "145-CD4", "146-CD8", "148-CD33", "147-CD20", "164-IgM", "167-CD38", "170-CD56", "110_114-CD3")

# Set this to the markers used for fold-change calculations.  Usually these are phospho-specific markers.
FUNCTIONAL_MARKERS=c("103-Viability", "141-pPLCgamma2", "150-pSTAT5", "152-Ki67", "154-pSHP2", "151-pERK1/2", "153-pMAPKAPK2", "156-pZAP70/Syk", "158-pBtk/Itk", "160-pSLP-76", "159-pSTAT3", "165-pNFkB", "166-IkBalpha", "168-pH3", "169-pP38", "171-pLck", "172-pS6", "174-pSrcFK", "176-pCREB", "175-pCrkL")

# Set this to all markers for which you want the basal medians.
ALL_MARKERS=c("Cell Length", "191-DNA", "103-Viability", "115-CD45", "139-CD45RA", "141-pPLCgamma2", "142-CD19", "144-CD11b", "145-CD4", "146-CD8", "148-CD33", "150-pSTAT5", "147-CD20", "152-Ki67", "154-pSHP2", "151-pERK1/2", "153-pMAPKAPK2", "156-pZAP70/Syk", "158-pBtk/Itk", "160-pSLP-76", "159-pSTAT3", "164-IgM", "165-pNFkB", "166-IkBalpha", "167-CD38", "168-pH3", "170-CD56", "169-pP38", "171-pLck", "172-pS6", "174-pSrcFK", "176-pCREB", "175-pCrkL", "110_114-CD3")

# Set this to the path to your local flowSPD package installation -- if using a global installation, set to NULL.
LIBRARY_PATH="lib/"

# Set this to the reference file to be used for fold-change calculations.
REFERENCE_FILE="CyTOF55_PBMC_Day4_1_Unstim1_Singlets.fcs"

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

FlowSPD.driver(FILE_TO_PROCESS, file_pattern="*.fcs", out_dir=OUTPUT_DIR, cluster_cols=SURFACE_MARKERS, median_cols=ALL_MARKERS, reference_file=REFERENCE_FILE, fold_cols=FUNCTIONAL_MARKERS)
FlowSPD.plot.trees(OUTPUT_DIR,file_pattern="*fcs*gml",out_dir=OUTPUT_DIR)
