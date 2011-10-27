# Cluster observations into ~k clusters
SPADE.cluster <- function(tbl, k) {
    if (nrow(tbl) > 60000) {
		warning("Potentially too many observations for the clustering step",immediate=TRUE);
    }

    # Transpose table before call into row major order
    clust <- .Call("SPADE_cluster",t(tbl),as.integer(1))



	# Invalid clusters have assgn == 0
	is.na(clust$assgn) <- which(clust$assgn == 0)

	# Build out hclust object
	hclust <- list()
	hclust$merge  <- clust$merge
	hclust$height <- clust$height
	hclust$order  <- 1:nrow(tbl)
	hclust$labels <- 1:nrow(tbl)
	class(hclust) <- "hclust"
	
	clust$assgn <- cutree(hclust,k=k);

	centers = c()	
	for (i in c(1:max(clust$assgn, na.rm=TRUE))) {  
		obs <- which(clust$assgn == i)
		if (length(obs) > 1) {
			centers <- rbind(centers,colMeans(tbl[obs,,drop=FALSE]))
			clust$assgn[obs] <- nrow(centers)
		} else {
			is.na(clust$assgn) <- obs
		}
    }


  return(list(centers=centers,assign=clust$assgn,hclust=hclust))
}


SPADE.clustersToMST <- function(centers, method="manhattan") {
    adjacency  <- dist(centers, method=method)
    full_graph <- graph.adjacency(as.matrix(adjacency),mode="undirected",weighted=TRUE)
    mst_graph  <- minimum.spanning.tree(full_graph)
    mst_graph
}

SPADE.writeGraph <- function(graph, outfilename) {
     write.graph(graph, outfilename, format="gml")
}

SPADE.FCSToTree <- function(infilenames, outfilename, graphfilename, clusterfilename, 
    cols=NULL, k=200, arcsinh_cofactor=5.0, desired_samples=50000,comp=TRUE) {
    
    data = c()
    for (f in infilenames) {
		# Load in FCS file
		in_fcs  <- SPADE.read.FCS(f,comp=comp);
		in_data <- exprs(in_fcs);

		params <- parameters(in_fcs);
		pd     <- pData(params);

		# Select out the desired columns
		if (is.null(cols)) { 
			cols <- as.vector(pd$name) 
		}
		idxs <- match(cols,pd$name)
		if (any(is.na(idxs))) { 
			stop("Invalid column specifier") 
		}
	
		data <- rbind(data,in_data[,idxs,drop=FALSE])
		colnames(data) <- pd$name[idxs]
    }

    # Downsample data if neccessary 
    if (nrow(data) > desired_samples) {
		data <- data[sample(1:nrow(data),desired_samples),]
    }

    # Compute the cluster centers, marking any single observation clusters as NA
    clust <- SPADE.cluster(asinh(data/arcsinh_cofactor),k);

	# Generate path for merge order
    mergeOrderPath = paste(dirname(outfilename),"/","merge_order.txt",sep="");

	# Write the DEFAULT merge order
	write.table(clust$hclust$merge,file=mergeOrderPath,sep="\t",quote=F,row.names=F,col.names=F)

    # Write out FCS file downsampled data used in clustering, along with assignment
    # Strip out observations in single observation clusters
    ff <- SPADE.build.flowFrame(subset(cbind(data, cluster=clust$assign),!is.na(clust$assign)))
    write.FCS(ff, outfilename) 

    # Write out the MST and cluster centers to specified files ignoring single observation clusters
    SPADE.writeGraph(SPADE.clustersToMST(clust$centers),graphfilename);
    write.table(clust$centers,file=clusterfilename,row.names=FALSE,col.names=colnames(data))
}


SPADE.cluster.normalizeHclust = function(spadeClusterObject){

	# Find events that have not been assigned to a cluster after clustering
	removedEvents = which(is.na(spadeClusterObject$assign))
	
	# If events removed...
	if (length(removedEvents)>0){
		hclust = list();
		numObs = length(spadeClusterObject$assign)
		merge = spadeClusterObject$hclust$merge
		
		# Shift down (closer to zero) negative singleton indices by one
		# for each missing event
		for (event in sort(removedEvents,decreasing=T)){
			merge[merge < (-event)] = (merge[merge < (-event)]+1)
		}
		
		# Chop off empty merges at end of merge matrix
		hclust$merge = merge[(-((numObs-1):(numObs-length(removedEvents)))),]
		
		# Chop off empty heights at end of height matrix
		hclust$height = spadeClusterObject$hclust$height[ -((numObs-1):(numObs-length(removedEvents))) ]

		# Remove non-clustered events from labels
		hclust$labels = (1:numObs)[-removedEvents]

		hclust$order = 1:length(hclust$labels)
		
		# Add an extra list parameter indicating which events were removed
		hclust$removedEvents = removedEvents;

		class(hclust) = "hclust"
	} else {
		hclust = spadeClusterObject$hclust
		hclust$removedEvents = NULL;
	}
	hclust
}
 
