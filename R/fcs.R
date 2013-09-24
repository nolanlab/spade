SPADE.read.FCS <- function(file, comp=TRUE, verbose=FALSE, ...) {
	if (verbose)
		fcs <- read.FCS(file, ...)
	else 
		fcs <- suppressWarnings(read.FCS(file, ...))
	
	params <- parameters(fcs)
	pd     <- pData(params)

	# Replace any null descs with names (for FSC-A, FSC-W, SSC-A)
	bad_col <- grep("^[a-zA-Z0-9]+",pd$desc,invert=TRUE)
	if (length(bad_col) > 0) {
		keyval <- keyword(fcs)
		for (i in bad_col) {
			pd$desc[i] <- pd$name[i]
			keyval[[paste("$P",i,"S",sep="")]] <- pd$name[i]
		}
		pData(params) <- pd;
		fcs <- flowFrame(exprs(fcs),params,description=description(fcs));
		keyword(fcs) <- keyval
	}

	# Compensate data if SPILL or SPILLOVER present, stripping compensation matrix 
	# out of the flowFrame, i.e we should only have to do this once
	apply.comp <- function(in_fcs, keyword) {
		comp_fcs <- compensate(in_fcs, description(in_fcs)[[keyword]])
		flowFrame(exprs(comp_fcs), parameters(comp_fcs), description(comp_fcs)[grep("SPILL",names(description(comp_fcs)),invert=TRUE)])
	}
	
	if (comp && !is.null(description(fcs)$SPILL)) {
		fcs <- apply.comp(fcs, "SPILL")	
	} else if (comp && !is.null(description(fcs)$SPILLOVER)) {
		fcs <- apply.comp(fcs, "SPILLOVER")	
	}
	
	fcs
}

SPADE.transform.matrix <- function(mat, tform=NULL) {
	if (is.null(tform)) {
		mat  # No-op
	} else {
		if (class(tform) == "transform") {
			apply(mat, 2, tform)
		} else {
			for (name in intersect(colnames(mat), names(tform))) {
				mat[,name] <- tform[[name]](mat[,name])
			}
			mat
		}
	} 
}

SPADE.transform.FCS <- function(ff, tform=NULL) {
	if (is.null(tform)) {
		ff  # No-op
	} else {

		if (class(tform) == "transform") {
			new_exprs <- apply(exprs(ff), 2, tform)
			new_range <- apply(range(ff), 2, tform)
		} else {
			new_exprs <- exprs(ff)
			new_range <- range(ff)
			for (name in intersect(colnames(ff), names(tform))) {
				new_exprs[,name] <- tform[[name]](new_exprs[,name])
				new_range[,name] <- tform[[name]](new_range[,name])
			}
			new_range <- as.matrix(new_range)
		}
				
		new_par <- parameters(ff)
		new_par$minRange <- new_range[1,]
		new_par$maxRange <- new_range[2,]
		new_par$range    <- (new_range[2,] - new_range[1,]) + 1

		flowFrame(new_exprs, new_par, description(ff))
	} 
}

SPADE.build.flowFrame <- function(x) {
	if (!is.matrix(x)) {
		stop("Input must be matrix")
  }

  # Build metadata for FCS file
  pd <- c()  # 'params' phenoData
  dl <- list()  # 'description' list
	
  dl[["$DATATYPE"]] <- "F"
  for (c in 1:ncol(x)) {
		c_name <- colnames(x)[c]
		    
		c_min <- min(x[,c])
		c_max <- max(x[,c])
		c_rng <- c_max - c_min + 1

		pl <- matrix(c(c_name, c_name, c_rng, c_min, c_max),nrow=1)
		colnames(pl) <- c("name", "desc", "range", "minRange", "maxRange")
		rownames(pl) <- paste("$P",c,sep="") 
		pd <- rbind(pd, pl)
		
		dl[[paste("$P",c,"B",sep="")]] <- "32";	    # Number of bits
		dl[[paste("$P",c,"R",sep="")]] <- toString(c_rng); # Range
		dl[[paste("$P",c,"E",sep="")]] <- "0,0";	    # Exponent
		dl[[paste("$P",c,"N",sep="")]] <- c_name;	    # Name
		dl[[paste("$P",c,"S",sep="")]] <- c_name;	    # Desc	
	}	
   
  flowFrame(x, as(data.frame(pd), "AnnotatedDataFrame"), description=dl)
}

SPADE.removeExistingDensityAndClusterColumns <- function(file) {
	# Do not comp or transform ... make this step invisible.
	input_file <- suppressWarnings(read.FCS(file))

	input_file_names <- names(input_file)

	if ("<cluster> cluster" %in% input_file_names ||
		"<density> density" %in% input_file_names) {
		# Drop those columns
		cleaned <- input_file[,!(input_file_names %in% c("<cluster> cluster", "<density> density"))]

		# Rename the original file. Increment the suffix if it already exists.
		suffix <- as.integer(regmatches(file, gregexpr("(\\d+)$", file, perl=TRUE)))

		if (is.na(suffix)) {
			suffix <- ".orig1"
			new_file_name <- paste(file, suffix, sep = "")
			# NB: file.rename is a namespaced function, not a method of the argument.
			file.rename(file, new_file_name)
		} else {
			suffix <- paste(".orig", suffix + 1, sep = "")
			new_file_name <- sub("(.orig\\d+)$", suffix, file);
			file.rename(file, new_file_name)
		}

		# Save with the original file name
		write.FCS(cleaned, file)
	}
}