# Transpose data before call to put in row major order
FlowSPD.assignToCluster <- function(tbl, clusters)
    .Call("FSPD_assign",t(tbl),t(clusters))

FlowSPD.addClusterToFCS <- function(infilename, outfilename, clusterfilename,
    cols=NULL, arcsinh_cofactor=5.0) {

    # Load in FCS file
    in_fcs  <- read.FCS(infilename);
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

    # Load in clusters
    clusters <- read.table(clusterfilename, header=TRUE, check.names=FALSE)
    if (!all(colnames(clusters) == pd$desc[idxs])) {
	stop("Selected columns don't match cluster columns")
    }

    # Assign observations to clusters
    assign <- FlowSPD.assignToCluster(asinh(in_data[,idxs]/arcsinh_cofactor),clusters)

    # Reload FCS file without transformation, so it can be accurately rewritten...
    in_fcs <- read.FCS(infilename,transform=FALSE)
    in_data <- exprs(in_fcs);

    # Add column named "density" to the FCS file
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
