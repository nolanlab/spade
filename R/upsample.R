# Transpose data before call to put in row major order
SPADE.assignToCluster <- function(tbl, cluster_data, cluster_assign)
	.Call("SPADE_assign",t(tbl),t(cluster_data),as.integer(cluster_assign))

SPADE.addClusterToFCS <- function(
	infilename, 
	outfilename, 
	clusterfilename,
	cols=NULL, 
	arcsinh_cofactor=NULL,
	transforms=flowCore::arcsinhTransform(a=0, b=0.2),	
	comp=TRUE
) {


	if (!is.null(arcsinh_cofactor)) {
		warning("arcsinh_cofactor is deprecated, use transform=flowCore::arcsinhTransform(...) instead")
		transforms <- flowCore::arcsinhTransform(a=0, b=1/arcsinh_cofactor)
	}
	

	# Load in FCS file
	in_fcs  <- SPADE.read.FCS(infilename,comp=comp);
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
	
	# Load in clustered FCS file
	cluster_fcs  <- SPADE.read.FCS(clusterfilename, comp=comp)
	cluster_data <- exprs(cluster_fcs)    
	
	cluster_params <- parameters(cluster_fcs)
	cluster_pd     <- pData(cluster_params)
	
	c_idxs <- match(cols, cluster_pd$name)
	na.fail(c_idxs)
	
	# Assign observations to clusters
	assign <- SPADE.assignToCluster(
		SPADE.transform.matrix(in_data[,idxs], transforms), 
		SPADE.transform.matrix(cluster_data[,c_idxs], transforms), 
		cluster_data[,"cluster"]
	)
	
	# Reload FCS file without transformation, so it can be compactly rewritten...
	in_fcs <- SPADE.read.FCS(infilename,comp=FALSE,transform=FALSE)
	in_data <- exprs(in_fcs);
	
	# Add column named "cluster" to the FCS file
	channel_number <- ncol(in_fcs)+1;
	channel_id     <- paste("$P",channel_number,sep="");
	channel_name   <- "cluster";
	channel_range  <- max(assign)+1;
	
	plist <- matrix(c(channel_name,channel_name,channel_range,0,channel_range-1));
	rownames(plist) <- c("name","desc","range","minRange","maxRange");
	colnames(plist) <- c(channel_id);
	
	pd <- rbind(pd,t(plist));
	pData(params) <- pd;
	
	out_data <- cbind(in_data,"cluster"=assign);
	out_frame <- flowFrame(out_data,params,description=description(in_fcs));
	
	keyval <- list();
	keyval[[paste("$P",channel_number,"B",sep="")]] <- "32";			# Number of bits
	keyval[[paste("$P",channel_number,"R",sep="")]] <- toString(channel_range); # Range
	keyval[[paste("$P",channel_number,"E",sep="")]] <- "0,0";			# Exponent
	keyval[[paste("$P",channel_number,"N",sep="")]] <- channel_name;		# Name
	keyval[[paste("$P",channel_number,"S",sep="")]] <- channel_name;		# Desc	
	keyword(out_frame) <- keyval;
	
	write.FCS(out_frame,outfilename);
}
