#!/share/apps/R/R-2.11.1/bin/Rscript
# ^^set this to your Rscript path
#
# runSPADE:  R wrapper script for SPADE tree construction
# Erin Simonds - esimonds@stanford.edu
# Version 2.6 - September 24, 2013
#
# Command line instructions:
#   1) Make sure the first line of this file is your Rscript path (found in the same directory as R)
#
#   2) Customize the parameters that are CAPITALIZED below 
#
#   3) In a command shell, navigate to the folder containing this script and the FCS file(s) to be analyzed
#
#   4) Make this script executable.  At the command line, run:
#   $ chmod +x runSPADE.R
#
#   5a) For normal use (not in a load sharing or compute cluster): At the command line, run:
#	$ ./runSPADE.R [-num_threads=X] [-file_to_process=Y]
#	where X is the number of threads you wish to use (default = 1)
#	Y is the name of the file to process (default = use file(s) specified in this script)
#		Note that parameters in brackets may be omitted.
#
#   5b) For Sun Gridengine: At the command line, run:
#	$ qsub -cwd -j y -b y -m e -M username@domain.ext [-pe threaded A] ./runSPADE.R [-num_threads=X] [-file_to_process=Y]
#	where username@domain.ext is your e-mail address to e-mail when the job is done,
#	A is the number of slots to reserve with Gridengine,
#	X is the number of threads to use in SPADE (usually the same as A) (default = 1).
#	Y is the name of the file to process (default = use file(s) specified in this script)
#		Note that parameters in brackets may be omitted.
#
#   5c) For Platform LSF: At the command line, run:
#	$ bsub [-n A -R "span[hosts=1]"] ./runSPADE.R [-num_threads=X] [-file_to_process=Y]
#	A is the number of slots to reserve with LSF,
#	X is the number of threads to use in SPADE (usually the same as A) (default = 1).
#	Y is the name of the file to process (default = use file(s) specified in this script)
#		Note that parameters in brackets may be omitted. (Note that the brackets surrounding "hosts=1" do not indicate an optional parameter.)
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
CLUSTERING_MARKERS=c("CD3(Cd114)Dd", "CD19(Nd142)Dd", "CD11c(Tb159)Dd", "CD14(Dy160)Dd", "CD20(Dy161)Dd")

# Set this to panels. For example, a panel may be one patient with their stimulated and unstimulated samples.
PANELS=list(
list(
panel_files=c( "a_Ungated.fcs","b_Ungated.fcs" ),
median_cols=NULL,
reference_files=c("a_Ungated.fcs"),
fold_cols=c("Time","Cell_length","CD3(Cd110)Dd","CD3(Cd111)Dd","CD3(Cd112)Dd","CD3(Cd114)Dd","CD19(Nd142)Dd","CD33(Nd148)Dd","CD235a(Sm152)Dd","CD11c(Tb159)Dd","DNA(Ir191)Dd","DNA(Ir193)Dd","CD14(Dy160)Dd","CD20(Dy161)Dd","CD8a(Dy162)Dd","CD45(Dy163)Dd","HLA-DR(Dy164)Dd","CD4(Gd158)Dd","Cd(110,111,112,114)")
)
)

# Set this to the path to your local spade package installation -- if using a global installation, set to NULL.
LIBRARY_PATH="lib/"

# Set this to the transforms for each channel
TRANSFORMS=c("Cd(110,111,112,114)"=flowCore::arcsinhTransform(a=0, b=0.2), "CD11c(Tb159)Dd"=flowCore::arcsinhTransform(a=0, b=0.2), "CD14(Dy160)Dd"=flowCore::arcsinhTransform(a=0, b=0.2), "CD19(Nd142)Dd"=flowCore::arcsinhTransform(a=0, b=0.2), "CD20(Dy161)Dd"=flowCore::arcsinhTransform(a=0, b=0.2), "CD235a(Sm152)Dd"=flowCore::arcsinhTransform(a=0, b=0.2), "CD3(Cd110)Dd"=flowCore::arcsinhTransform(a=0, b=0.2), "CD3(Cd111)Dd"=flowCore::arcsinhTransform(a=0, b=0.2), "CD3(Cd112)Dd"=flowCore::arcsinhTransform(a=0, b=0.2), "CD3(Cd114)Dd"=flowCore::arcsinhTransform(a=0, b=0.2), "CD33(Nd148)Dd"=flowCore::arcsinhTransform(a=0, b=0.2), "CD4(Gd158)Dd"=flowCore::arcsinhTransform(a=0, b=0.2), "CD45(Dy163)Dd"=flowCore::arcsinhTransform(a=0, b=0.2), "CD8a(Dy162)Dd"=flowCore::arcsinhTransform(a=0, b=0.2), "Cell_length"=flowCore::linearTransform(a=1.0), "DNA(Ir191)Dd"=flowCore::arcsinhTransform(a=0, b=0.2), "DNA(Ir193)Dd"=flowCore::arcsinhTransform(a=0, b=0.2), "HLA-DR(Dy164)Dd"=flowCore::arcsinhTransform(a=0, b=0.2), "Time"=flowCore::linearTransform(a=1.0))

# Pick one of these and set the rest to NULL.
DOWNSAMPLING_TARGET_NUMBER=NULL
DOWNSAMPLING_TARGET_PCTILE=NULL
DOWNSAMPLING_TARGET_PERCENT=0.1 # 10% -- recommended.

# Set this to the desired number of events to use for clustering.  Recommended:  50000
CLUSTERING_SAMPLES=50000

# Set this to the desired threshold for outliers (unit is the local density for each cell).  A higher value will discard more events.    Recommended:  0.01
DOWNSAMPLING_EXCLUDE_PCTILE=0.01

# Set this to the target number of clusters.  Algorithm will create clusters within 50% of this value.  Recommended:  200
TARGET_CLUSTERS=200

# This is a linear multiplier for node size.
NODE_SIZE_SCALE_FACTOR=1.2

# Path to output directory.  Recommended:  "output/"
OUTPUT_DIR="output/"

# Path to temporary directory.  Recommended:  "/tmp/"
TMPDIR="/tmp/"

# Set this to the desired layout function.  Recommended:  SPADE.layout.arch
# Option 1:  SPADE.layout.arch  ...this is fast, but allows overlapping nodes and edges
# Option 2:  SPADE.layout.arch_layout  ...this is slow, but prevents overlapping nodes and edges
LAYOUT_FUNCTION=SPADE.layout.arch

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

SPADE.driver(FILE_TO_PROCESS, file_pattern="*", out_dir=OUTPUT_DIR, cluster_cols=CLUSTERING_MARKERS, panels=PANELS, transforms=TRANSFORMS, layout=LAYOUT_FUNCTION, downsampling_target_percent=DOWNSAMPLING_TARGET_PERCENT, downsampling_target_number=DOWNSAMPLING_TARGET_NUMBER, downsampling_target_pctile=DOWNSAMPLING_TARGET_PCTILE, downsampling_exclude_pctile=DOWNSAMPLING_EXCLUDE_PCTILE, k=TARGET_CLUSTERS, clustering_samples=CLUSTERING_SAMPLES)

LAYOUT_TABLE <- read.table(paste(OUTPUT_DIR,"layout.table",sep=""))
MST_GRAPH <- read.graph(paste(OUTPUT_DIR,"mst.gml",sep=""),format="gml")
SPADE.plot.trees(MST_GRAPH,OUTPUT_DIR,file_pattern="*fcs*Rsave",layout=as.matrix(LAYOUT_TABLE),out_dir=paste(OUTPUT_DIR,"pdf",sep=""))

Sys.unsetenv("OMP_NUM_THREADS")
