# Helper functions

SPADE.strip.sep <- function(name) {
    ifelse(substr(name,nchar(name),nchar(name))==.Platform$file,substr(name,1,nchar(name)-1),name)
}

SPADE.normalize.out_dir <- function(out_dir) {
	out_dir_info <- file.info(out_dir)
		if (is.na(out_dir_info$isdir)) {
			dir.create(out_dir)
		}
	if (!file.info(out_dir)$isdir) {
		stop(paste("out_dir:",out_dir,"is not a directory"))
	}
	out_dir <- paste(SPADE.strip.sep(out_dir),.Platform$file,sep="")
}

# Driver convenience functions

SPADE.driver <- function(files, file_pattern="*.fcs", out_dir=".", cluster_cols=NULL, arcsinh_cofactor=5.0, layout=SPADE.layout.arch, median_cols=NULL, reference_files=NULL, fold_cols=NULL, downsampling_samples=20000, downsampling_exclude_pctile=0.01, downsampling_target_pctile=0.05, k=200, clustering_samples=50000) {
	if (length(files) == 1 && file.info(files)$isdir) {
		files <- dir(SPADE.strip.sep(files),full.names=TRUE,pattern=glob2rx(file_pattern))
	}
	if (length(files) == 0) {
		stop("No input files found")
	} 
	out_dir_info <- file.info(out_dir)
		if (is.na(out_dir_info$isdir)) {
			dir.create(out_dir)
		}
	if (!file.info(out_dir)$isdir) {
		stop(paste("out_dir:",out_dir,"is not a directory"))
	}
	out_dir <- paste(SPADE.strip.sep(out_dir),.Platform$file,sep="")

	density_files <- c()
	sampled_files <- c()
	for (f in files) {
		cat("Downsampling file:",f,"\n")
		f_density <- paste(out_dir,basename(f),".density.fcs",sep="")
		f_sampled <- paste(out_dir,basename(f),".downsample.fcs",sep="")

		SPADE.addDensityToFCS(f, f_density, cols=cluster_cols, arcsinh_cofactor=arcsinh_cofactor)
		SPADE.downsampleFCS(f_density, f_sampled, 
							exclude_pctile=downsampling_exclude_pctile,
							target_pctile=downsampling_target_pctile,
							desired_samples=downsampling_samples)

		density_files <- c(density_files, f_density)
		sampled_files <- c(sampled_files, f_sampled)	
	}

	cat("Clustering files...\n")
	cells_file <- paste(out_dir,"clusters.fcs",sep="")
	clust_file <- paste(out_dir,"clusters.table",sep="")
	graph_file <- paste(out_dir,"mst.gml",sep="")
	SPADE.FCSToTree(sampled_files, cells_file, graph_file, clust_file, 
					cols=cluster_cols, 
					arcsinh_cofactor=arcsinh_cofactor, 
					k=k, 
					desired_samples=clustering_samples)

	sampled_files <- c()
	for (f in density_files) {
		cat("Upsampling file:",f,"\n")
		f_sampled <- paste(f,".cluster.fcs",sep="")
		SPADE.addClusterToFCS(f, f_sampled, cells_file, cols=cluster_cols, arcsinh_cofactor=arcsinh_cofactor)
		sampled_files <- c(sampled_files, f_sampled)
	}

	graph  <- read.graph(graph_file, format="gml")
	layout_table <- layout(graph)
	write.table(layout_table,paste(out_dir,file="layout.table",sep=""),row.names = FALSE,col.names = FALSE)

	reference_medians <- NULL
	if (!is.null(reference_files)) {
		reference_files   <- sapply(as.vector(reference_files), function(rf) { paste(out_dir, rf, ".density.fcs.cluster.fcs",sep=""); })
		reference_medians <- SPADE.markerMedians(reference_files, cols=fold_cols, arcsinh_cofactor=arcsinh_cofactor, cluster_cols=cluster_cols)
	}

	for (f in sampled_files) {
		if (!is.null(reference_medians)) {	# If a reference file is specified		
			cat("Computing medians for file:",f,"\n")
			# Compute the median marker intensities in each node
			a <- SPADE.markerMedians(f, cols=median_cols, arcsinh_cofactor=arcsinh_cofactor, cluster_cols=cluster_cols)
			
			# Compute the overall cell frequency per node
			a[["percent"]] <- a$count / sum(a$count) * 100; colnames(a[["percent"]]) <- c("total");
			
			cat("Computing fold change for file:",f,"\n")
			b <- SPADE.markerMedians(f, cols=fold_cols, arcsinh_cofactor=arcsinh_cofactor, cluster_cols=cluster_cols)	

			# Compute the fold change compared to reference medians
			cc <- rownames(a$medians)[!is.na(match(rownames(b$medians),rownames(reference_medians$medians)))]  # Common clusters
			fold <- b$medians[cc,,drop=FALSE] - reference_medians$medians[cc,,drop=FALSE]

			dimnames(fold) <- list(cc, colnames(b$medians))
			b <- list(count=b$count, fold=fold)		

			# Merge the fold-change columns with the count, frequency, and median columns
			ab <- list(count = a$count, percent = a$percent, median = a$median, fold = b$fold)

			SPADE.write.graph(SPADE.annotateGraph(graph, layout=layout_table, anno=ab), paste(f,".medians.gml",sep=""), format="gml")
		} else {
			cat("Computing medians for file:",f,"\n")
			a <- SPADE.markerMedians(f, cols=median_cols, arcsinh_cofactor=arcsinh_cofactor, cluster_cols=cluster_cols)
			
			# Compute the overall cell frequency per node
			a[["percent"]] <- a$count / sum(a$count) * 100; colnames(a[["percent"]]) <- c("total");
			
			SPADE.write.graph(SPADE.annotateGraph(graph, layout=layout_table, anno=a), paste(f,".medians.gml",sep=""), format="gml")
		}
	}

	invisible(NULL)
}

# The following functions are copied from the TeachingDemos package, 
# authored by Greg Snow <greg.snow at imail.org>, and licensed under
# the Artistic-2.0 license.

cnvrt.coords <- function(x,y=NULL,input=c('usr','plt','fig','dev','tdev')) {

  input <- match.arg(input)
  xy <- xy.coords(x,y, recycle=TRUE)

  cusr <- par('usr')
  cplt <- par('plt')
  cfig <- par('fig')
  cdin <- par('din')
  comi <- par('omi')
  cdev <- c(comi[2]/cdin[1],(cdin[1]-comi[4])/cdin[1],
            comi[1]/cdin[2],(cdin[2]-comi[3])/cdin[2])

  if(input=='usr'){
    usr <- xy

    plt <- list()
    plt$x <- (xy$x-cusr[1])/(cusr[2]-cusr[1])
    plt$y <- (xy$y-cusr[3])/(cusr[4]-cusr[3])

    fig <- list()
    fig$x <- plt$x*(cplt[2]-cplt[1])+cplt[1]
    fig$y <- plt$y*(cplt[4]-cplt[3])+cplt[3]

    dev <- list()
    dev$x <- fig$x*(cfig[2]-cfig[1])+cfig[1]
    dev$y <- fig$y*(cfig[4]-cfig[3])+cfig[3]

    tdev <- list()
    tdev$x <- dev$x*(cdev[2]-cdev[1])+cdev[1]
    tdev$y <- dev$y*(cdev[4]-cdev[3])+cdev[3]

    return( list( usr=usr, plt=plt, fig=fig, dev=dev, tdev=tdev ) )
  }

  if(input=='plt') {

    plt <- xy

    usr <- list()
    usr$x <- plt$x*(cusr[2]-cusr[1])+cusr[1]
    usr$y <- plt$y*(cusr[4]-cusr[3])+cusr[3]

    fig <- list()
    fig$x <- plt$x*(cplt[2]-cplt[1])+cplt[1]
    fig$y <- plt$y*(cplt[4]-cplt[3])+cplt[3]

    dev <- list()
    dev$x <- fig$x*(cfig[2]-cfig[1])+cfig[1]
    dev$y <- fig$y*(cfig[4]-cfig[3])+cfig[3]

    tdev <- list()
    tdev$x <- dev$x*(cdev[2]-cdev[1])+cdev[1]
    tdev$y <- dev$y*(cdev[4]-cdev[3])+cdev[3]

    return( list( usr=usr, plt=plt, fig=fig, dev=dev, tdev=tdev ) )
  }

  if(input=='fig') {

    fig <- xy

    plt <- list()
    plt$x <- (fig$x-cplt[1])/(cplt[2]-cplt[1])
    plt$y <- (fig$y-cplt[3])/(cplt[4]-cplt[3])

    usr <- list()
    usr$x <- plt$x*(cusr[2]-cusr[1])+cusr[1]
    usr$y <- plt$y*(cusr[4]-cusr[3])+cusr[3]

    dev <- list()
    dev$x <- fig$x*(cfig[2]-cfig[1])+cfig[1]
    dev$y <- fig$y*(cfig[4]-cfig[3])+cfig[3]

    tdev <- list()
    tdev$x <- dev$x*(cdev[2]-cdev[1])+cdev[1]
    tdev$y <- dev$y*(cdev[4]-cdev[3])+cdev[3]

    return( list( usr=usr, plt=plt, fig=fig, dev=dev, tdev=tdev ) )
  }

  if(input=='dev'){
    dev <- xy

    fig <- list()
    fig$x <- (dev$x-cfig[1])/(cfig[2]-cfig[1])
    fig$y <- (dev$y-cfig[3])/(cfig[4]-cfig[3])

    plt <- list()
    plt$x <- (fig$x-cplt[1])/(cplt[2]-cplt[1])
    plt$y <- (fig$y-cplt[3])/(cplt[4]-cplt[3])

    usr <- list()
    usr$x <- plt$x*(cusr[2]-cusr[1])+cusr[1]
    usr$y <- plt$y*(cusr[4]-cusr[3])+cusr[3]

    tdev <- list()
    tdev$x <- dev$x*(cdev[2]-cdev[1])+cdev[1]
    tdev$y <- dev$y*(cdev[4]-cdev[3])+cdev[3]

    return( list( usr=usr, plt=plt, fig=fig, dev=dev, tdev=tdev ) )
  }

  if(input=='tdev'){
    tdev <- xy

    dev <- list()
    dev$x <- (tdev$x-cdev[1])/(cdev[2]-cdev[1])
    dev$y <- (tdev$y-cdev[3])/(cdev[4]-cdev[3])

    fig <- list()
    fig$x <- (dev$x-cfig[1])/(cfig[2]-cfig[1])
    fig$y <- (dev$y-cfig[3])/(cfig[4]-cfig[3])

    plt <- list()
    plt$x <- (fig$x-cplt[1])/(cplt[2]-cplt[1])
    plt$y <- (fig$y-cplt[3])/(cplt[4]-cplt[3])

    usr <- list()
    usr$x <- plt$x*(cusr[2]-cusr[1])+cusr[1]
    usr$y <- plt$y*(cusr[4]-cusr[3])+cusr[3]

    tdev <- list()
    tdev$x <- dev$x*(cdev[2]-cdev[1])+cdev[1]
    tdev$y <- dev$y*(cdev[4]-cdev[3])+cdev[3]

    return( list( usr=usr, plt=plt, fig=fig, dev=dev, tdev=tdev ) )
  }

}

subplot <- function(fun, x, y=NULL, size=c(1,1), vadj=0.5, hadj=0.5,
                    inset=c(0,0), type=c('plt','fig'), pars=NULL){

  old.par <- par(no.readonly=TRUE)
  on.exit(par(old.par))

  type <- match.arg(type)

  if(missing(x)) x <- locator(2)

  if(is.character(x)) {
      if(length(inset) == 1) inset <- rep(inset,2)
      x.char <- x
      tmp <- par('usr')
      x <- (tmp[1]+tmp[2])/2
      y <- (tmp[3]+tmp[4])/2

      if( length(grep('left',x.char, ignore.case=TRUE))) {
          x <- tmp[1] + inset[1]*(tmp[2]-tmp[1])
          if(missing(hadj)) hadj <- 0
      }
      if( length(grep('right',x.char, ignore.case=TRUE))) {
          x <- tmp[2] - inset[1]*(tmp[2]-tmp[1])
          if(missing(hadj)) hadj <- 1
      }
      if( length(grep('top',x.char, ignore.case=TRUE))) {
          y <- tmp[4] - inset[2]*(tmp[4]-tmp[3])
          if(missing(vadj)) vadj <- 1
      }
      if( length(grep('bottom',x.char, ignore.case=TRUE))) {
          y <- tmp[3] + inset[2]*(tmp[4]-tmp[3])
          if(missing(vadj)) vadj <- 0
      }
  }

  xy <- xy.coords(x,y)

  if(length(xy$x) != 2){
    pin <- par('pin')
    tmp <- cnvrt.coords(xy$x[1],xy$y[1],'usr')$plt

    x <- c( tmp$x - hadj*size[1]/pin[1],
            tmp$x + (1-hadj)*size[1]/pin[1] )
    y <- c( tmp$y - vadj*size[2]/pin[2],
            tmp$y + (1-vadj)*size[2]/pin[2] )

    xy <- cnvrt.coords(x,y,'plt')$fig
  } else {
    xy <- cnvrt.coords(xy,,'usr')$fig
  }

  par(pars)
  if(type=='fig'){
      par(fig=c(xy$x,xy$y), new=TRUE)
  } else {
      par(plt=c(xy$x,xy$y), new=TRUE)
  }
  fun
  tmp.par <- par(no.readonly=TRUE)

  return(invisible(tmp.par))
}

SPADE.normalize.trees <- function(files, file_pattern="*.gml", out_dir=".", layout=SPADE.layout.arch, attr_pattern="percent|median|fold", normalize="global") {
	if (length(files) == 1 && file.info(files)$isdir) {
		files <- dir(SPADE.strip.sep(files),full.names=TRUE,pattern=glob2rx(file_pattern))    
    }
	out_dir <- SPADE.normalize.out_dir(out_dir)


	clean_attr <- function(attr) {
		attr[attr == Inf | attr == -Inf] <- NA  # Clean up bogus values ...
		attr
	}

	attr_ranges <- function(file) {
		ar <- c()	
		graph <- read.graph(f, format="gml")
		attrs <- list.vertex.attributes(graph)
		for (i in grep(attr_pattern, attrs)) {
			attr <- clean_attr(get.vertex.attribute(graph,attrs[i]))			
			ar[[attrs[i]]] <- range(ar[[attrs[i]]], attr, na.rm=TRUE)
		}
		ar
	}

	sf_fn <- NULL
	if (normalize == "global") {
		gr <- c()  # Ranges of all encountered attributes
		for (f in files) {
			ar <- attr_ranges(f)
			for (i in seq_along(ar)) { n <- names(ar[i]); gr[[n]] <- range(gr[[n]], ar[i], na.rm=TRUE); }
		}			
		sf_fn <- function(a_name, a_vals) { 
			ifelse(is.null(gr[[a_name]]),max(abs(range(a_vals, na.rm=TRUE))), max(abs(gr[[a_name]])))
		}
	} else if (normalize == "local") {
		sf_fn <- function(a_name, a_vals) {
			max(abs(range(a_vals, na.rm=TRUE)))
		}
	} else
		stop("Unsupported kind of normalization")
			

	for (f in files) {
		graph <- read.graph(f, format="gml")
		attrs <- list.vertex.attributes(graph)
				
		if (is.function(layout))
			graph_l <- layout(graph)
		else
			graph_l <- layout	
		
		for (i in grep(attr_pattern, attrs)) {
			attr <- clean_attr(get.vertex.attribute(graph,attrs[i]))
			sf <- sf_fn(attrs[i], attr)
			if (sf != 0.0)
				attr <- scale(attr, center=0.0, scale=sf)
			graph <- set.vertex.attribute(graph, attrs[i], value=attr)
		}
		SPADE.write.graph(SPADE.annotateGraph(graph, layout=graph_l), paste(out_dir,basename(f),".",normalize,".norm.gml",sep=""), format="gml")
	}

}

SPADE.plot.trees <- function(files, file_pattern="*.gml", out_dir=".", layout=SPADE.layout.arch, attr_pattern="percent|median|fold", scale=NULL, pctile_color=c(0.02,0.98), normalize="global") {
    if (length(files) == 1 && file.info(files)$isdir) {
		files <- dir(SPADE.strip.sep(files),full.names=TRUE,pattern=glob2rx(file_pattern))    
    }
	out_dir <- SPADE.normalize.out_dir(out_dir)

	if (!is.null(scale) && (!is.vector(scale) || length(scale) !=2))
		stop("scale must be a two element vector")

	if (!is.vector(pctile_color) || length(pctile_color) != 2) {
		stop("pctile_color must be a two element vector with values in [0,1]")
	}

	clean_attr <- function(attr) {
		attr[attr == Inf | attr == -Inf] <- NA  # Clean up bogus values ...
		attr
	}

	boundaries <- NULL
	if (normalize == "global") {
		boundaries <- c()  # Calculate ranges of all encountered attributes with trimmed outliers
		all_attrs <- c()
		for (f in files) {
			graph <- read.graph(f, format="gml")
			attrs <- list.vertex.attributes(graph)
			for (i in grep(attr_pattern, attrs)) {
				all_attrs[[attrs[i]]] <- c(all_attrs[[attrs[i]]], clean_attr(get.vertex.attribute(graph,attrs[i])))
			}
		}
		for (i in seq_along(all_attrs)) {
			boundaries[[names(all_attrs)[i]]] <- quantile(all_attrs[[i]], probs=pctile_color, na.rm=TRUE)
		}
	}

    jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    colorscale <- jet.colors(100)

    for (f in files) {
		graph <- read.graph(f, format="gml")
		attrs <- list.vertex.attributes(graph)

		if (is.function(layout))
			graph_l <- layout(graph)
		else
			graph_l <- layout	

		vsize <- V(graph)$percenttotal
		vsize[vsize == Inf | vsize == -Inf] <- NA  # Clean up bogus values	
		vsize <- vsize/max(vsize,na.rm=TRUE) * 3 + 2
		vsize[is.na(vsize)] <- 1

		for (i in grep(attr_pattern, attrs)) {
			# Compute the color for each vertex using color gradient scaled from min to max
			attr <- get.vertex.attribute(graph,attrs[i])
			
			attr[attr == Inf | attr == -Inf] <- NA  # Clean up bogus values ...
			if (!is.null(scale))
				boundary <- scale
			else if (normalize == "global") { # Scale to global min/max
					boundary <- boundaries[[attrs[i]]] # Recall trimmed global boundary for this attribute
			} else # Scale to local min/max
					boundary <- quantile(attr, probs=pctile_color, na.rm=TRUE)  # Trim outliers for this attribtue
				
			if (boundary[1] == boundary[2]) {  boundary <- c(boundary[1]-1, boundary[2]+1); }  # Prevent "zero" width gradients
			if (length(grep("median|percent", attrs[i])))
				boundary <- c(min(boundary), max(boundary))  # Dont make range symmetric for median or percent values
			else
				boundary <- c(-max(abs(boundary)), max(abs(boundary)))  # Make range symmetric for fold-change values
			
			boundary <- round(boundary, 2) # Round boundary values to 2 digits of precision so scale looks nice
				
			grad <- seq(boundary[1], boundary[2], length.out=length(colorscale))
		
			color <- colorscale[findInterval(attr, grad,all.inside=TRUE)]
			color[is.na(attr)] <- "grey"
			V(graph)$color <- color
	    
			# Plot the tree, with legend showing the gradient
			pdf(paste(out_dir,basename(f),".",attrs[i],".pdf",sep=""))
	    
			plot(graph, layout=graph_l, vertex.shape="circle", edge.color="grey", vertex.size=vsize, vertex.frame.color=NA, vertex.label=NA) 
			
			# Substitute pretty attribute names
			if (length(grep("median", attrs[i])))
				attrs[i] <- sub("median", "Median of ", attrs[i])
			else if (length(grep("fold", attrs[i])))
				attrs[i] <- sub("fold", "Arcsinh diff. of ", attrs[i])
			else
				attrs[i] <- sub("percent", "Percent freq. of ", attrs[i])
						
			# Make parameters used for clustering obvious
			if (length(grep("_clust$", attrs[i])))
				attrs[i] <- sub("_clust", "\n(Used for tree-building)", attrs[i])
			
			title(main=paste(strsplit(basename(f),".fcs")[[1]][1], sub=attrs[i], sep="\n"))
			subplot(
				image(
					grad, c(1), matrix(1:length(colorscale),ncol=1), col=colorscale,
					xlab=ifelse(is.null(scale),paste("Range:",pctile_color[1],"to",pctile_color[2],"pctile"),""),
					ylab="", yaxt="n", xaxp=c(boundary,1)
				),
				x="right,bottom",size=c(1,.20)
			)

			dev.off()
		}
    }
}

SPADE.workflow.concat.FCS <- function(files, file_pattern="*.density.fcs.cluster.fcs", out_dir=".", cols=NULL, layout=SPADE.layout.arch, arcsinh_cofactor=5.0, in_graph_file="mst.gml", out_graph_file="concat.gml", cluster_cols=NULL)  {
	if (length(files) == 1 && file.info(files)$isdir) {
		in_graph_file <- paste(SPADE.strip.sep(files),in_graph_file,sep=.Platform$file)
		files <- dir(SPADE.strip.sep(files),full.names=TRUE,pattern=glob2rx(file_pattern))
    }
	out_dir <- SPADE.normalize.out_dir(out_dir)

	anno  <- SPADE.markerMedians(files,cols=cols,arcsinh_cofactor=arcsinh_cofactor,cluster_cols=cluster_cols)
	graph <- read.graph(in_graph_file, format="gml")

	# Compute the overall cell frequency per node
	anno[["percent"]] <- anno$count / sum(anno$count) * 100; colnames(anno[["percent"]]) <- c("total");

	out_graph_file <- paste(out_dir,out_graph_file,sep="")

	SPADE.write.graph(SPADE.annotateGraph(graph, layout=layout, anno=anno), out_graph_file, format="gml")
}


