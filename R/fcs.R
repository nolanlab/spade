SPADE.read.FCS <- function(file, comp=TRUE, ...) {
	fcs <- read.FCS(file, ...)
	params <- parameters(fcs)
	pd <- pData(params)

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
