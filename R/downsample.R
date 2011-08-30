# Transpose table before call to put in row major order
SPADE.density <- function(tbl, kernel_mult=5.0, apprx_mult=1.5, med_samples=2000)
    .Call("SPADE_density",t(tbl),kernel_mult,apprx_mult,med_samples)

SPADE.addDensityToFCS <- function(infilename, outfilename, 
    cols=NULL, arcsinh_cofactor=5.0, kernel_mult=5.0, apprx_mult=1.5, med_samples=2000, comp=TRUE) {

    # Load in FCS file, touching up missing descriptions if needed
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

    # Compute the density
    density <- SPADE.density(	asinh(in_data[,idxs]/arcsinh_cofactor),
				kernel_mult=kernel_mult,apprx_mult=apprx_mult,med_samples=med_samples);
		if (max(density) == 0.0)
			warning(paste(infilename,"has degenerate densities, possibly due to many identical observations",sep=" "))

    # Reload FCS file without transformation, so it can be accurately rewritten...
    in_fcs <- SPADE.read.FCS(infilename,comp=FALSE,transform=FALSE)
    in_data <- exprs(in_fcs);

    # Add column named "density" to the FCS file
    channel_number <- ncol(in_fcs)+1;
    channel_id     <- paste("$P",channel_number,sep="");
    channel_name   <- "density";
    channel_range  <- max(density)+1;


    plist <- matrix(c(channel_name,channel_name,channel_range,0,channel_range-1));
    rownames(plist) <- c("name","desc","range","minRange","maxRange");
    colnames(plist) <- c(channel_id);

    pd <- rbind(pd,t(plist));
    pData(params) <- pd;

    out_data <- cbind(in_data,"density"=density);
    out_frame <- flowFrame(out_data,params,description=description(in_fcs));

    keyval <- list() 
    keyval[[paste("$P",channel_number,"B",sep="")]] <- "32";			# Number of bits
    keyval[[paste("$P",channel_number,"R",sep="")]] <- toString(channel_range); # Range
    keyval[[paste("$P",channel_number,"E",sep="")]] <- "0,0";			# Exponent
    keyval[[paste("$P",channel_number,"N",sep="")]] <- channel_name;		# Name
    keyval[[paste("$P",channel_number,"S",sep="")]] <- channel_name;		# Desc	
    keyword(out_frame) <- keyval;

    write.FCS(out_frame,outfilename);
}

SPADE.downsampleFCS <- function(infilename, outfilename, exclude_pctile=0.01, target_pctile=0.05, desired_samples=NULL) {
    # Load in FCS file
    in_fcs  <- SPADE.read.FCS(infilename,comp=FALSE,transform=FALSE);
    in_data <- exprs(in_fcs);

    params <- parameters(in_fcs);
    pd     <- pData(params);

    d_idx <- match("density",pd$name)
    if (is.na(d_idx)) {
	stop("No density parameter in FCS file")
    }
    
    # boundary[1]: exclusion, boundary[2]: potential target
    boundary <- quantile(in_data[,d_idx],c(exclude_pctile,target_pctile),names=FALSE)
    
    out_data <- subset(in_data, in_data[,d_idx] > boundary[1]) # Exclusion    

    density <- out_data[,d_idx]
    if (is.null(desired_samples)) {
		boundary <- boundary[2]
		out_data <- subset(out_data,boundary/density > runif(nrow(out_data)))
    } else if (desired_samples < nrow(out_data)) {
		# Need to find target density such there are approximately desired_samples
		# remaining after downsampling. To do so we solve for the density such that
		# the sum of samples below that density plus the expected value of
		# samples retained above that density equals approximately the desired
		# number of samples
		density_s <- sort(density)
		cdf       <- rev(cumsum(1.0/rev(density_s)))
		
		# Default solution if target density smaller than any present
		boundary <- desired_samples/cdf[1] 
		if (boundary > density_s[1]) {  # Boundary actually falls amongst densities present
			targets <- (desired_samples-1:length(density_s)) / cdf 
			boundary <- targets[which.min(targets-density_s > 0)]
		}
		out_data  <- subset(out_data,boundary/density > runif(length(density)))
    }

    out_frame <- flowFrame(out_data,params,description=description(in_fcs))
    write.FCS(out_frame,outfilename)
}
