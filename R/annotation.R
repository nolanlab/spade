SPADE.annotateMarkers <- function(files, cols=NULL, arcsinh_cofactor=5.0) {
    warning("Deprecated: Use SPADE.markerMedians instead")
    SPADE.markerMedians(files, cols=cols, arcsinh_cofactor=arcsinh_cofactor)
}

SPADE.markerMedians <- function(files, cols=NULL, arcsinh_cofactor=5.0) {

	data  <- c()
	c_ids <- c()

	files <- as.vector(files)
	for (f in files) {
		# Load in FCS file
		in_fcs  <- SPADE.read.FCS(f);
		in_data <- exprs(in_fcs);

		params <- parameters(in_fcs);
		pd     <- pData(params);

		# Find cluster column
		c_idx <- match("cluster",pd$desc)
		if (any(is.na(c_idx))) {
			stop("No cluster parameter in FCS file")
		}
		c_asn <- in_data[,c_idx]
		c_ids <- sort(union(c_ids, unique(c_asn))) # Used cluster indices

		# Select out the desired columns
		if (is.null(cols)) {
			cols <- as.vector(pd$desc) 
		}
		idxs <- match(cols,pd$desc)
		if (any(is.na(idxs))) { 
			stop("Invalid column specifier") 
		}
		idxs <- subset(idxs, idxs != c_idx)  # Strip out cluster column
	
		if (!is.null(colnames(data)) && !setequal(colnames(data),pd$desc[idxs])) {
			stop("Files have different columns")
		}
		data <- rbind(data, in_data[, idxs])
		colnames(data) <- pd$desc[idxs]
	}
	   
	count  <- c()
	medians <- c()

	for (i in c_ids) {
		data_s  <- asinh(subset(data, c_asn == i) / arcsinh_cofactor)
		count   <- rbind(count, nrow(data_s))
		medians <- rbind(medians, apply(data_s, 2, median))
	}

	colnames(count)   <- c("count")
	rownames(count)   <- c_ids
	colnames(medians) <- colnames(data)
	rownames(medians) <- c_ids

    list(count=count, medians=medians)
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
		vt <- V(graph)[match(rownames(l),V(graph)$name)-1]  # Vertex IDS are 0 indexed
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
    
    if (file == "")
	file <- stdout()
    else if (is.character(file)) {
	file <- file(file, "w")
	on.exit(close(file))
    }
    else if (!isOpen(file, "w")) {
	open(file, "w")
	on.exit(close(file))
    }
    if (!inherits(file, "connection")) {
	stop("'file' must be a character string or a connection")
    }

    write.gml <- function(graph,file) {
	
	write.attr <- function(name, attr) {
	    # Strip out non-alphanumeric characters to avoid Cytoscape parsing errors
	    name <- gsub("[^A-Za-z0-9]","",name)
	    if (length(grep("^[0-9]",name))) {
		name <- paste("spade",name,sep="")
	    }
	    if (is.na(attr))
		paste(name,"NaN")
	    else if (is.character(attr))
		paste(name," \"",attr,"\"",sep="")
	    else
		paste(name,attr)
	}
	
	writeLines(c("graph [", paste("directed",ifelse(is.directed(graph),1,0))),con=file)

	# Identify known "structs"	
	v_attr <- list.vertex.attributes(graph)
	v_attr_g <- v_attr[grep("graphics[.]",v_attr)]  # graphics attributes
	v_attr <- setdiff(v_attr, v_attr_g)
	if (length(grep("[.]",v_attr)) > 0) {
	    stop("Unsupported struct in vertex attributes")
	}
	
	for (v in V(graph)) {
	    writeLines("node [",con=file)
	 
	    for (a in v_attr) {
		writeLines(write.attr(a,get.vertex.attribute(graph,a,index=v)),con=file)
	    } 	    

	    if (length(v_attr_g) > 0) {
			writeLines("graphics [",con=file)
			for (a in v_attr_g) {
				parts <- unlist(strsplit(a,"[.]"))
				writeLines(write.attr(parts[2],get.vertex.attribute(graph,a,index=v)),con=file)
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
		writeLines(paste(a,get.edge.attribute(graph,a,index=e)),con=file)
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

SPADE.layout.arch_layout <- function(mst_graph)  
{
	if (!is.igraph(mst_graph))
	{
		cat("Input has to be igraph object")
		return(NA)
	}
	if (!is.connected(mst_graph)) 
	{	
		cat("Cannot handle graph that has disjoint components\n")
		return(NA)
	}
	

	# get number of nodes, and binary adjacency matrix
	num_nodes = length(V(mst_graph));
	tmp = get.adjacency(mst_graph) + t(get.adjacency(mst_graph)); 
	diag(tmp) = 0;
	adj =  array(as.double(tmp!=0), c(num_nodes,num_nodes));

	
	# make sure it is a tree, no circles
	tmp_adj = adj;
	diag(tmp_adj)=0;
	is_node_visited = array(0,c(num_nodes,1)); is_node_visited[1]=1;
	e = array(0,c(num_nodes,1)); e[1]=1;
	for (i in 1:num_nodes-1)
	{
		e_new = tmp_adj %*% e;
		tmp_adj[e!=0,]=0; tmp_adj[,e!=0]=0;
		e = e_new;
		if (any(e!=0 & is_node_visited==1))
		{
			cat("Cannot handle graphs with cycles\n");
			return(NA)
		}
		is_node_visited[e!=0,1]=1;
	}




	# find the shortest distance between each pairs of nodes
	shortest_hop = array(0,c(num_nodes,num_nodes));
	diag(adj) = 1;
	for (i in 1:num_nodes)
	{
		e = array(0,c(num_nodes,1));
		e[i]=1; 
		hop=0;
		while (any(shortest_hop[i,]==0))
		{
			e = array(as.double(adj %*% e > 0), c(num_nodes,1));
			hop = hop + 1;
			shortest_hop[i,t(e)>0 & shortest_hop[i,]==0]=hop;
		}
	}
	diag(adj)=0;
	diag(shortest_hop)=0;


	# find the longest path, which is the backbone arch
	tmp = which(shortest_hop==max(as.vector(shortest_hop)),arr.in=TRUE)[1,];
	backbone_start = min(tmp); backbone_end = max(tmp);
	back_bones = which(shortest_hop[backbone_start,]+shortest_hop[backbone_end,]==shortest_hop[backbone_start,backbone_end]);
	back_bones = back_bones[order(shortest_hop[backbone_start,back_bones])];


	# find sidechains and subsidechains
	side_chains = list();
	side_chain_roots = back_bones;
	counter = 1;
	while (counter<=length(side_chain_roots))
	{
		root_node = side_chain_roots[counter];
		first_neighbors = setdiff(which(adj[root_node,]!=0),side_chain_roots);
		if (length(first_neighbors)==0) {counter = counter + 1;	next;}
		for (i in 1:length(first_neighbors))
		{
			# find the side chain starting from root_node to first_neighbors[i] to as far as it can go
			subtree_nodes_through_this_neighbor = which(shortest_hop[root_node,]>shortest_hop[first_neighbors[i],]);
			end_node = subtree_nodes_through_this_neighbor[which(shortest_hop[root_node,subtree_nodes_through_this_neighbor]==max(shortest_hop[root_node,subtree_nodes_through_this_neighbor]))[1]];
			sub_backbone = which(shortest_hop[root_node,]+shortest_hop[end_node,]==shortest_hop[root_node,end_node]);
			sub_backbone = sub_backbone[order(shortest_hop[root_node,sub_backbone])];			

			side_chain_roots = c(side_chain_roots,sub_backbone[-1]);
			side_chains[[length(side_chains)+1]] = sub_backbone;
		}			
	}

#print(back_bones-1);
#if (length(side_chains)!=0) {for (i in 1:length(side_chains)) {print(side_chains[[i]]-1);}}

	# determine node location of the back_bones
	node_positions = array(0,c(2,num_nodes));
	position_assigned_flag = array(0,c(1,num_nodes));
	backbone_node_angles = 1:length(back_bones);
	backbone_node_angles = backbone_node_angles - mean(backbone_node_angles);
	backbone_node_angles = backbone_node_angles/max(backbone_node_angles)*(pi/90*25);
	node_positions[1,back_bones] = sin(backbone_node_angles);
	node_positions[2,back_bones] = -cos(backbone_node_angles);
	normalization_const = sqrt(as.vector(t(node_positions[,back_bones[1]]-node_positions[,back_bones[2]])%*%(node_positions[,back_bones[1]]-node_positions[,back_bones[2]])));
	node_positions = node_positions/normalization_const;
	position_assigned_flag[back_bones]=1;

	cat("generating layout for ", num_nodes, " nodes ... ",sum(position_assigned_flag), "\n");

if (length(side_chains)!=0)
{
	# determine the order of subbackbones to draw
	centerness = array(-1,c(1,length(side_chains)));
	for (k in 1:length(side_chains))
	{
		if (length(intersect(side_chains[[k]][1],back_bones))==0) {break;}
		centerness[k] = abs(shortest_hop[side_chains[[k]][1],back_bones[1]]-shortest_hop[side_chains[[k]][1],back_bones[length(back_bones)]]);
	}
	centerness = centerness[centerness!=-1];

	# determine node location for each node in the side chains
	side_chain_ordering = order(centerness);
	while(length(side_chain_ordering)<length(side_chains))
	{side_chain_ordering = c(side_chain_ordering,length(side_chain_ordering)+1);}


	for (i in side_chain_ordering)
	{
		for (j in 1:length(side_chains[[i]]))
		{
			if (position_assigned_flag[side_chains[[i]][j]]==1) {next;}
			new_node = side_chains[[i]][j];
			attaching_node = which(adj[new_node,]!=0 & position_assigned_flag==1);
			r = seq(from=0.3, to=0.9, by=0.04);
			theta = seq(from=0, to=2*pi, by=pi/180);
			potential_position_force = array(0,c(length(r),length(theta),2));
			for (m in 1:length(r))
			{
				for (n in 1:length(theta))
				{
					potential_position = node_positions[,attaching_node] + as.vector(r[m])*array(c(cos(theta[n]),sin(theta[n])),c(2,1));
					repel_vector = as.vector(potential_position) - node_positions[,position_assigned_flag==1];
					repel_force = rowSums(repel_vector/(array(1,c(2,1))%*%sqrt(colSums(repel_vector^2))^5));
					cos_alfa = (repel_force)%*%array(c(cos(theta[n]),sin(theta[n])),c(2,1))/sqrt(sum(repel_force^2));
					sin_alfa = sqrt(1-cos_alfa^2);

					potential_position_force[m,n,1]=sqrt(sum(repel_force^2))*sin_alfa;  # force perpendicular to the string
					potential_position_force[m,n,2]=sqrt(sum(repel_force^2))*cos_alfa;  # force along the sting
				}
			}

			best_for_each_layer = array(0,c(0,2));
			best_string_force_each_layer=array(0,c(0,1));
			for (m in 1:length(r))
			{
				n_s = which(potential_position_force[m, ,2]>=0);
				if (length(n_s)==0) {next;}
				I = which(potential_position_force[m,n_s,1]==min(potential_position_force[m,n_s,1]))[1];
				best_for_each_layer = rbind(best_for_each_layer,c(m,n_s[I]));
				best_string_force_each_layer = rbind(best_string_force_each_layer, potential_position_force[m,n_s[I],2]);
			}
			I = which(best_string_force_each_layer==min(best_string_force_each_layer))[1];
			best_r = as.vector(r[best_for_each_layer[I,1]]);
			best_theta = as.vector(theta[best_for_each_layer[I,2]]);

			node_positions[,new_node] = node_positions[,attaching_node] + best_r*array(c(cos(best_theta),sin(best_theta)),c(2,1));
			position_assigned_flag[side_chains[[i]][j]]=1;
	
			cat("generating layout for ", num_nodes, " nodes ... ",sum(position_assigned_flag), "\n");

		}
	}

}


node_positions = node_positions - c(median(node_positions[1,]),median(node_positions[2,]));
node_positions[1,] = node_positions[1,]/max(abs(node_positions[1,]))*50;
node_positions[2,] = node_positions[2,]/max(abs(node_positions[2,]))*50;

t(node_positions)
}


