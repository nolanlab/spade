SPADE.annotateMarkers <- function(files, cols=NULL, arcsinh_cofactor=5.0) {
	warning("Deprecated: Use SPADE.markerMedians instead")
	SPADE.markerMedians(files, cols=cols, arcsinh_cofactor=arcsinh_cofactor)
}

SPADE.markerMedians <- function(files, num.clusters, cols=NULL, arcsinh_cofactor=NULL, transforms=flowCore::arcsinhTransform(a=0, b=0.2), cluster_cols=NULL, comp=TRUE) {

	if (!is.null(arcsinh_cofactor)) {
		warning("arcsinh_cofactor is deprecated, use transform=flowCore::arcsinhTransform(...) instead")
		transforms <- flowCore::arcsinhTransform(a=0, b=1/arcsinh_cofactor)
	}

	data  <- c()
	
	files <- as.vector(files)
	for (f in files) {
		# Load in FCS file
		in_fcs  <- SPADE.read.FCS(f,comp=comp);
		in_data <- exprs(in_fcs);

		params <- parameters(in_fcs);
		pd     <- pData(params);
	
		# Select out the desired columns
		if (is.null(cols)) {
			cols <- as.vector(pd$name) 
		}
		if (!"cluster" %in% cols) {
			cols <- c(cols,"cluster")
		}

		idxs <- match(cols,pd$name)
		if (any(is.na(idxs))) { 
			stop("Invalid column specifier") 
		}	

		data <- rbind(data, in_data[, idxs,drop=FALSE])

	}
	

	clst <- data[,"cluster"]
	data <- data[,colnames(data)!="cluster",drop=FALSE]
	data_t <- SPADE.transform.matrix(data, transforms) 

	# TODO: Weird things were being done to the naming, and this breaks that so we can do the transforms cleanly...
	colnames(data) <- sapply(colnames(data),function(x) { 
		if (x %in% cluster_cols)
			x <- paste(x,"clust",sep="_")
		x
	})
	colnames(data_t) = colnames(data)

	ids  <- 1:num.clusters
	if (any(is.na(match(unique(clst),ids)))) {
		stop("More clusters in FCS files than indicated")
	}

	count   <- matrix(0,  nrow=num.clusters, ncol=1, dimnames=list(ids, "count"))
	medians <- matrix(NA, nrow=num.clusters, ncol=ncol(data), dimnames=list(ids,colnames(data)))
	raw_medians <- matrix(NA, nrow=num.clusters, ncol=ncol(data), dimnames=list(ids,colnames(data)))
	cvs     <- matrix(NA, nrow=num.clusters, ncol=ncol(data), dimnames=list(ids,colnames(data)))
	for (i in ids) {
		data_s  <- subset(data, clst == i)
		data_s_t <- subset(data_t, clst == i)
		
		count[i,1]  <- nrow(data_s_t)
		medians[i,] <- apply(data_s_t, 2, median)
		raw_medians[i,] <- apply(data_s, 2, median)
		cvs[i,]     <- apply(data_s_t, 2, function(d) { 100*sd(d)/abs(mean(d)) })
	} 
	percenttotal <- matrix((count / sum(count)) * 100.0, nrow=num.clusters, ncol=1, dimnames=list(ids, "percenttotal"))
	list(count=count, medians=medians, raw_medians=raw_medians, cvs=cvs, percenttotal=percenttotal)
}

SPADE.layout.arch <-  function(mst_graph) {
	if (!is.igraph(mst_graph)) {
		stop("Input has to be igraph object")
	}
	if (!is.connected(mst_graph)) {	
		stop("Cannot handle graph that has disjoint components")
	}
		
	# Make sure it is a tree, no circles
	if (girth(mst_graph)$girth > 0) {
		stop("Cannot handle graphs with cycles");
	}
	
	# Find the distance between nodes measured in "hops"
	hops <- c()
	for (v in V(mst_graph)) {
		hops <- rbind(hops,unlist(lapply(get.shortest.paths(mst_graph,v),length))-1)
	}


	# Compute the positions for each vertices
	# --------------------------------------------------------------------------
	v_pos <- array(0,c(vcount(mst_graph),2))

	# The longest path is the backbone arch
	terminals <- which(hops == max(hops), arr.ind=TRUE)[1,] - 1  # igraph vertices are 0-indexed
	back_bone <- unlist(get.shortest.paths(mst_graph, from=terminals["row"], to=terminals["col"]))
	  
	# Layout the backbone arch along an arc
	bb_span <- pi * .55  # span in radians of back bone arch
	bb_unit <- 50.0  # unit distance bewteen nodes
	angles  <- seq(pi/2-bb_span/2, pi/2+bb_span/2, length.out=length(back_bone)) 
	v_pos[back_bone+1,] <- bb_unit*length(back_bone)*cbind(cos(angles),-sin(angles))

	# Layout the side chains as trees normal to the backbone
	for (v in back_bone) {
		# Find subset of vertices that compose side chain by deleting links between current
		# backbone vertex and rest of the backbone and then performing subcomponent
		# Note: E(mst_graph,P=c(mapply(c,v,n))) == E(mst_graph)[v %--% n] but is much faster
		n <- intersect(neighbors(mst_graph, v), back_bone)
		side_v <- sort(subcomponent(delete.edges(mst_graph,E(mst_graph,P=c(mapply(c,v,n)))),v))
		
		# Compute layout for side chains and integrate it into overall layout 
		# Note: Relies on side_v being in sorted order
		if (length(side_v) > 1) {
			# Convert side chains to directed graph, deleting edges with decreasing hop distance
			side_h <- hops[v+1,side_v+1]  # hops between back_bone node and side chain
			side_g <- as.directed(subgraph(mst_graph,side_v),mode="mutual")
			e <- get.edges(side_g,E(side_g))+1  # edges as a matrix, recall igraph is 0-indexed
			side_g <- delete.edges(side_g,subset(E(side_g),side_h[e[,1]] > side_h[e[,2]]))	   

			# Layout side chain
			# -----------------------------------------------------------------
			root <- which.min(side_h)-1
			layout <- layout.reingold.tilford(side_g,root=root)
	
			# rotate tree to be normal to back bone
			polar <- cbind(atan2(layout[,2],layout[,1]), sqrt(rowSums(layout^2)))
			polar[,1] <- polar[,1] + angles[back_bone==v] - pi/2
			layout <- bb_unit*polar[,2]*cbind(cos(polar[,1]),-sin(polar[,1]))	    
	
			# translate layout to back_bone 	    
			layout <- layout + matrix(v_pos[v+1,] - layout[root+1,], nrow=nrow(layout), ncol=2, byrow=TRUE)
			v_pos[side_v+1,] <- layout
		}
	}
	v_pos
}

SPADE.annotateGraph <- function(graph, layout=NULL, anno=NULL) {
	if (!is.igraph(graph)) {
		stop("Not a graph object")
	}

	if (!is.null(layout) && is.matrix(layout)) {
		if (nrow(layout) != vcount(graph) || ncol(layout) != 2) {
			stop("Ill-formated layout matrix, must 2 columns (x,y) and as many rows as vertices")
		}
		# Over write non-struct graphics attributes if they exist
		v_attr <- list.vertex.attributes(graph)
		for (i in grep("^graphics$",v_attr))
			graph <- remove.vertex.attribute(graph, v_attr[i])
		
		graph <- set.vertex.attribute(graph, "graphics.x", value=layout[,1])
		graph <- set.vertex.attribute(graph, "graphics.y", value=layout[,2])
	}

	if (is.null(anno)) {
		return(graph)
	} else if (!is.list(anno)) {
		stop("anno must be a list with named entries");
	}
	
	for (i in seq_len(length(anno))) {
		l <- anno[[i]]
		if (!is.matrix(l)) {
			stop(paste("Argument:",quote(l),"must be a matrix"))
		}
		vt <- V(graph)[match(rownames(l),V(graph)$name)]  # Vertex IDS are 1 indexed
		for (c in colnames(l)) {
			graph <- set.vertex.attribute(graph,ifelse(names(anno)[i] == c,c,paste(names(anno)[i],c,sep="")),index=vt, value=l[,c])
		}
	}
	graph
}

SPADE.write.graph <- function(graph, file="", format = c("gml")) {
	if (!is.igraph(graph)) {
		stop("Not a graph object")
	}
	
	if (file == "") {
		file <- stdout()
	} else if (is.character(file)) {
		file <- file(file, "w")
		on.exit(close(file))
	} else if (!isOpen(file, "w")) {
		open(file, "w")
		on.exit(close(file))
	}
	if (!inherits(file, "connection")) {
		stop("'file' must be a character string or a connection")
	}

	write.gml <- function(graph,file) {
	
		write.attr <- function(name, attr) {
			# Strip out non-alphanumeric characters to avoid Cytoscape parsing errors
			name <- gsub("[^A-Za-z0-9_]","",name)
			if (length(grep("^[0-9]",name))) {
				name <- paste("spade",name,sep="")
			}
			if (is.na(attr) || is.nan(attr))
				stop("Unexpected NA or NaN attribute")
			else if (is.character(attr) && nchar(attr) > 0)
				paste(name," \"",attr,"\"",sep="")
			else if (is.integer(attr))
				paste(name,attr)
			else
				paste(name,formatC(attr,format="f"))
		}
		
		writeLines(c("graph [", paste("directed",ifelse(is.directed(graph),1,0))),con=file)

		# Identify known "structs"	
		v_attr <- list.vertex.attributes(graph)
		v_attr_g <- v_attr[grep("graphics[.]",v_attr)]  # graphics attributes
		v_attr <- setdiff(v_attr, c(v_attr_g, "id"))  
		
		for (v in V(graph)) {
			writeLines("node [",con=file)
			
			writeLines(paste("id",v),con=file)
			for (a in v_attr) {
				val <- get.vertex.attribute(graph,a,index=v)
				if (!is.na(val) && !is.nan(val))
					writeLines(write.attr(a,val),con=file)
			}

			if (length(v_attr_g) > 0) {
				writeLines("graphics [",con=file)
				for (a in v_attr_g) {
					parts <- unlist(strsplit(a,"[.]"))
					val <- get.vertex.attribute(graph,a,index=v)
					if (!is.na(val) && !is.nan(val))
						writeLines(write.attr(parts[2],val),con=file)
				}
				writeLines("]",con=file)
			}
			
			writeLines("]",con=file)
		}
		
		# Identify known "structs"	
		e_attr <- list.edge.attributes(graph)
		if (length(grep("[.]",e_attr)) > 0) {
			stop("Unsupported struct in edge attributes")
		}

		for (e in E(graph)) {
			writeLines("edge [",con=file)
		
			pts <- get.edges(graph,e)
			writeLines(c(paste("source",pts[1]), paste("target",pts[2])),con=file)
			for (a in e_attr) {
				val <- get.edge.attribute(graph,a,index=e)
				if (!is.na(val) && !is.nan(val))
					writeLines(write.attr(a,val),con=file)
			}	

			writeLines("]",con=file)	
		}
		writeLines("]",con=file)
	}

	res <- switch(format,   
		gml = write.gml(graph,file), 
		stop(paste("Unsupported output format:",format)))
	invisible(res)			
}
