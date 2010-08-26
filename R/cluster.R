# Cluster observations into ~k clusters
FlowSPD.cluster <- function(tbl, k) {
    if (nrow(tbl) > 60000) {
		warning("Potentially too many observations for the clustering step",immediate=TRUE);
    }

    # Transpose table before call into row major order
    clust <- .Call("FSPD_cluster",t(tbl),as.integer(k))
  
	# Invalid clusters have assgn == 0
	centers = c()
    is.na(clust$assgn) <- which(clust$assgn == 0)
	for (i in c(1:max(clust$assgn, na.rm=TRUE))) {  
		obs <- which(clust$assgn == i)
		if (length(obs) > 1) {
			centers <- rbind(centers,colMeans(tbl[obs,]))
			clust$assgn[obs] <- nrow(centers)
		} else {
			is.na(clust$assgn) <- obs
		}
    }
    return(list(centers=centers,assign=clust$assgn))
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

FlowSPD.FCSToTree <- function(infilenames, outfilename, graphfilename, clusterfilename, 
    cols=NULL, k=200, arcsinh_cofactor=5.0, desired_samples=50000) {
    
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

    # Compute the cluster centers, marking any single observation clusters as NA
    clust <- FlowSPD.cluster(asinh(data/arcsinh_cofactor),k);
    
    # Write out FCS file downsampled data used in clustering, along with assignment
    # Strip out observations in single observation clusters
    ff <- FlowSPD.build.flowFrame(subset(cbind(data, cluster=clust$assign),!is.na(clust$assign)))
    write.FCS(ff, outfilename) 

    # Write out the MST and cluster centers to specified files ignoring single observation clusters
    FlowSPD.writeGraph(FlowSPD.clustersToMST(clust$centers),graphfilename);
    write.table(clust$centers,file=clusterfilename,row.names=FALSE,col.names=colnames(data))
}
 
