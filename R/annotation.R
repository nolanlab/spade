FlowSPD.annotateMarkers <- function(infilename, cols=NULL, arcsinh_cofactor=5.0) {

    # Load in FCS file
    in_fcs  <- read.FCS(infilename);
    in_data <- exprs(in_fcs);

    params <- parameters(in_fcs);
    pd     <- pData(params);

    # Find cluster column
    c_idx <- match("cluster",pd$desc)
    if (is.na(c_idx)) {
	stop("No cluster parameter in FCS file")
    }
    c_asn <- in_data[,c_idx]
    c_ids <- sort(unique(c_asn)) # Used cluster indices

    # Select out the desired columns
    if (is.null(cols)) {
	cols <- as.vector(pd$desc) 
    }
    idxs <- match(cols,pd$desc)
    if (any(is.na(idxs))) { 
	stop("Invalid column specifier") 
    }
    idxs <- subset(idxs, idxs != c_idx)  # Strip out cluster column
    data_c <- in_data[,idxs]             # Select out data columns    

    counts  <- c()
    medians <- c()
   
    for (i in c_ids) {
	data_s <- asinh(subset(data_c, c_asn == i) / arcsinh_cofactor)
	
	counts  <- rbind(counts, nrow(data_s))
	medians <- rbind(medians, apply(data_s, 2, median))
    }
 
    colnames(counts)  <- c("count")
    rownames(counts)  <- c_ids
    colnames(medians) <- pd$desc[idxs]
    rownames(medians) <- c_ids

    list(counts=counts, medians=medians)
}
