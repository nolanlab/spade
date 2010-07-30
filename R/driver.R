FlowSPD.driver <- function(files, out_dir=".", cluster_cols=NULL, arcsinh_cofactor=5.0, layout=FlowSPD.layout.arch, median_cols=NULL, reference_file = NULL, fold_cols=NULL) {
    if (length(files) == 1 && file.info(files)$isdir) {
	files <- dir(files,pattern="*.fcs")
    }

    out_dir_info <- file.info(out_dir)
    if (is.na(out_dir_info$isdir)) {
	dir.create(out_dir)
    }
    if (!file.info(out_dir)$isdir) {
	stop(paste("out_dir:",out_dir,"is not a directory"))
    }
    out_dir <- paste(out_dir,.Platform$file,sep="")

    density_files <- c()
    sampled_files <- c()
    for (f in files) {
	f_density <- paste(out_dir,f,".density.fcs",sep="")
	f_sampled <- paste(out_dir,f,".downsample.fcs",sep="")
	
	FlowSPD.addDensityToFCS(f, f_density, cols=cluster_cols, arcsinh_cofactor=arcsinh_cofactor)
	FlowSPD.downsampleFCS(f_density, f_sampled)
	
	density_files <- c(density_files, f_density)
	sampled_files <- c(sampled_files, f_sampled)	
    }

    clust_file <- paste(out_dir,"clusters.table",sep="")
    graph_file <- paste(out_dir,"mst.gml",sep="")
    FlowSPD.FCSToTree(sampled_files, graph_file, clust_file, cols=cluster_cols, arcsinh_cofactor=arcsinh_cofactor)

    sampled_files <- c()
    for (f in density_files) {
	f_sampled <- paste(f,".cluster.fcs",sep="")
	FlowSPD.addClusterToFCS(f, f_sampled, clust_file, cols=cluster_cols, arcsinh_cofactor=arcsinh_cofactor)
	sampled_files <- c(sampled_files, f_sampled)
    }

    graph  <- read.graph(graph_file, format="gml")
    layout <- layout(graph)
    
    reference_medians = NULL
    if (!is.null(reference_file)) {
	reference_medians <- FlowSPD.markerMedians(paste(out_dir,reference_file,".density.fcs.cluster.fcs",sep=""), cols=fold_cols, archsinh_cofactor=arcsinh_cofactor)
    }

    for (f in sampled_files) {
	a <- FlowSPD.markerMedians(f, cols=median_cols, arcsinh_cofactor=arcsinh_cofactor)
	g <- FlowSPD.annotateGraph(graph, layout=layout, a)
	FlowSPD.write.graph(g, paste(f,".medians.gml",sep=""), format="gml")
	
	if (!is.null(reference_medians)) {
	    a <- FlowSPD.markerMedians(f, cols=fold_cols, arcsinh_cofactor=arcsinh_cofactor)
	    a <- list(counts=a$counts, fold=(a$medians-reference_medians$medians))
	    g <- FlowSPD.annotateGraph(graph, layout=layout, a)
	    FlowSPD.write.graph(g, paste(f,".fold.gml",sep=""), format="gml")
	}
    }
    
    invisible(NULL)
}
