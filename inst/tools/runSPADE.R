#!/share/apps/R/R-2.11.1/bin/Rscript
# ^^set this to your Rscript path
#
# runSPADE:  R wrapper script for SPADE tree construction
# Erin Simonds - esimonds@stanford.edu
# Version 2.3 - November 9, 2010
#
# Command line instructions:
#   0) Make sure the first line of this file is your Rscript path
#
#   1) Customize the parameters that are CAPITALIZED below 
#
#   2) In a bash shell, navigate to the folder containing this script and the FCS file(s) to be analyzed
#
#   3a) For normal (non-cluster) use: At the command line, run:
#	$ runSPADE.R -num_threads=X [-file_to_process=Y]
#	where X is the number of threads you wish to use
#	Y is the name of the file to process, if not specified in this script
#
#   3b) For Sun Gridengine: At the command line, run:
#	$ qsub -cwd -j y -m e -M username@domain.ext -pe threaded A runSPADE.R -num_threads=X [-file_to_process=Y]
#	where username@domain.ext is your e-mail address to e-mail when the job is done,
#	A is the number of slots to reserve with Gridengine,
#	X is the number of threads you wish to use in SPADE (usually the same as A).
#	Y is the name of the file to process, if not specified in this script
#
#   3c) For Platform LSF: At the command line, run:
#	$ bsub -n A -R "span[hosts=1]" runSPADE.R -num_threads=X [-file_to_process=Y]
#	A is the number of slots to reserve with LSF,
#	X is the number of threads you wish to use in SPADE (usually the same as A).
#	Y is the name of the file to process, if not specified in this script
#
# Interactive instructions:
#   1) Customize the parameters that are CAPITALIZED below 
#   2) In an interactive R session, change the working directory to a folder containing this script and the FCS file to be analyzed
#   3) At the R command line, run: source("runSPADE.R")
#
#

# Set this to the FCS file or folder to be processed.  If a folder is specified, all FCS files will be processed.
FILE_TO_PROCESS="FCS_files/"

# Set this to the reagent names ($P*S keywords from the FCS file) to be used for clustering and tree construction.  Usually these are the surface markers.
SURFACE_MARKERS=c("Cell_length","CD45","CD19","CD47","CD34","CD33","CD38","CD56","CD20","CD10","DNA","110_114-CD3")

# Set this to the reagent names ($P*S keywords from the FCS file) to be used for fold-change calculations.  Usually these are phospho-specific markers.  If you don't have any, set to NULL
FUNCTIONAL_MARKERS=c("pP38","pStat5","pErk","pSTAT3","pSLP76","IkB","pPLCg2","cCaspase3","pS6","pSrc","pCreb","pCrkL","pc-Cbl","pShp2","pBtk/Itk","pSyk")

# Set this to all markers for which you want the basal medians.
ALL_MARKERS=c(SURFACE_MARKERS, FUNCTIONAL_MARKERS)

# Set this to the path to your local spade package installation -- if using a global installation, set to NULL.
LIBRARY_PATH="lib/"

# Set this to the reference file(s) to be used for fold-change calculations.  If multiple files are specified, the medians will be averaged.  If you don't have one, set to NULL
REFERENCE_FILE=c("Basal1.fcs","Basal2.fcs")

# Set this to the desired arcsinh cofactor to transform cytometry data [transformed_value = arcsinh(original_value/cofactor)]
# CyTOF recommended value:  5
# Fluorescence recommended value:  150
ARCSINH_COFACTOR=5

# Set this to the desired number of events remaining after downsampling files.  Recommended:  50000
DOWNSAMPLED_EVENTS=50000

# Set this to the desired number of events to use for clustering.  Recommended:  50000
CLUSTERING_SAMPLES=50000

# Set this to the desired threshold for outliers (unit is the local density for each cell).  A higher value will discard more events.    Recommended:  0.01
DOWNSAMPLING_EXCLUDE_PCTILE=0.01

# Set this to the target number of clusters.  Algorithm will create clusters within 50% of this value.  Recommended:  200
TARGET_CLUSTERS=200

# Path to output directory.  Recommended:  "output/"
OUTPUT_DIR="output/"

# Path to temporary directory.  Recommended:  "/tmp/"
TMPDIR="/tmp/"

###################  No need to modify anything below this point ########################

#Run SPADE analysis workflow

#Default value:
NUM_THREADS <- 1

for (e in commandArgs()) {
	ta <- strsplit(e,"=",fixed=TRUE)
	if( ta[[1]][1] == "-num_threads") {
		NUM_THREADS <- ta[[1]][2]
	}
	if( ta[[1]][1] == "-file_to_process") {
		FILE_TO_PROCESS <- ta[[1]][2]
	}
}

Sys.setenv("OMP_NUM_THREADS"=NUM_THREADS)

library("spade",lib.loc=LIBRARY_PATH)

# Set this to the desired layout function.  Recommended:  SPADE.layout.arch
# Option 1:  SPADE.layout.arch  ...this is fast, but allows overlapping nodes and edges
# Option 2:  SPADE.layout.arch_layout  ...this is slow, but prevents overlapping nodes and edges
LAYOUT_FUNCTION=SPADE.layout.arch

SPADE.driver(FILE_TO_PROCESS, file_pattern="*.fcs", out_dir=OUTPUT_DIR, cluster_cols=SURFACE_MARKERS, arcsinh_cofactor=ARCSINH_COFACTOR, layout=LAYOUT_FUNCTION, median_cols=ALL_MARKERS, reference_files=REFERENCE_FILE, fold_cols=FUNCTIONAL_MARKERS, downsampling_samples=DOWNSAMPLED_EVENTS, downsampling_exclude_pctile=DOWNSAMPLING_EXCLUDE_PCTILE, k=TARGET_CLUSTERS, clustering_samples=CLUSTERING_SAMPLES)

LAYOUT_TABLE <- read.table(paste(OUTPUT_DIR,"layout.table",sep=""))

SPADE.plot.trees(OUTPUT_DIR,file_pattern="*fcs*gml",layout=as.matrix(LAYOUT_TABLE),out_dir=paste(OUTPUT_DIR,"pdf",sep=""))

Sys.unsetenv("OMP_NUM_THREADS")
