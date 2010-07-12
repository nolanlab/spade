# Cluster observations into ~k clusters
FlowSPD.cluster <- function(tbl, k) {
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

FlowSPD.FCSToTree <- function(infilename, graphfilename, clusterfilename, cols=NULL, k=2, arcsinh_cofactor=5.0) {
    # Load in FCS file
    in_fcs  <- read.FCS(infilename);
    in_data <- exprs(in_fcs);

    params <- parameters(in_fcs);
    pd     <- pData(params);

    # Select out the desired columns
    idxs <- c();
    if (!is.null(cols)) {
	for (cl in cols) {
	    idxs <- c(idxs,which(pd$desc == cl,arr.ind<-TRUE));
	}
    } else {
	idxs <- c(1:nrow(pd));
    }  

    centers <- FlowSPD.cluster(asinh(in_data[,idxs]/arcsinh_cofactor),k);
    FlowSPD.writeGraph(FlowSPD.clustersToMST(centers),graphfilename);
    write.table(centers,file=clusterfilename,row.names=FALSE,col.names=pd$desc[idxs])
}
 
