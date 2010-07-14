# Cluster observations into ~k clusters
FlowSPD.cluster <- function(tbl, k) {
    if (nrow(tbl) > 50000) {
	warning("Potentially too many observations for the clustering step",immediate=TRUE);
    }

    # Transpose table before call into row major order
    clust <- .Call("FSPD_cluster",t(tbl),as.integer(k))
    
    centers = c()
    for (i in c(1:max(clust$assgn))) {
	centers = rbind(centers,colMeans(tbl[clust$assgn == i,]))
    }
    centers
}

FlowSPD.clustersToMST <- function(centers, method="manhattan") {
    adjacency  <- dist(centers, method=method)
    full_graph <- graph.adjacency(as.matrix(adjacency),mode="undirected",weighted=TRUE)
    mst_graph  <- minimum.spanning.tree(full_graph)
    mst_graph
}

FlowSPD.writeGraph <- function(graph, outfilename) {
     write.graph(graph, outfilename, format="gml")
}

FlowSPD.FCSToTree <- function(infilenames, graphfilename, clusterfilename, 
    cols=NULL, k=200, arcsinh_cofactor=5.0, desired_samples=40000) {
    
    data = c()
    for (f in infilenames) {
	# Load in FCS file
	in_fcs  <- read.FCS(f);
	in_data <- exprs(in_fcs);

	params <- parameters(in_fcs);
	pd     <- pData(params);

	# Select out the desired columns
	if (is.null(cols)) { 
	    cols = as.vector(pd$desc) 
	}
	idxs <- match(cols,pd$desc)
	if (any(is.na(idxs))) { 
	    stop("Invalid column specifier") 
	}
	
	data <- rbind(data,in_data[,idxs])
	colnames(data) <- pd$desc[idxs]
    }

    # Downsample data if neccessary 
    if (nrow(data) > desired_samples) {
	data <- data[sample(1:nrow(data),desired_samples),]
    }

    # Compute the cluster centers
    centers <- FlowSPD.cluster(asinh(data/arcsinh_cofactor),k);
    
    # Write out the MST and cluster centers to specified files
    FlowSPD.writeGraph(FlowSPD.clustersToMST(centers),graphfilename);
    write.table(centers,file=clusterfilename,row.names=FALSE,col.names=colnames(data))
}
 
