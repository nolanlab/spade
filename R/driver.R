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

SPADE.flattenAnnotations <- function(annotations) {
	stopifnot(is.list(annotations))
	flat <- NULL
	for (i in seq_along(annotations)) {	
		df <- data.frame(annotations[[i]])
		
		to_paste <- names(annotations)[i] != colnames(df)
		if (any(to_paste))
			colnames(df) <- paste(names(annotations)[i],colnames(df),sep="")
	
		if (is.null(flat))
			flat <- df
		else
			flat <- cbind(flat, df)
	}
	flat
}

SPADE.driver <- function(files, file_pattern="*.fcs", out_dir=".", cluster_cols=NULL, panels=NULL, comp=TRUE, arcsinh_cofactor=5.0, downsampling_samples=20000, downsampling_exclude_pctile=0.01, downsampling_target_pctile=0.05, k=200, clustering_samples=50000, layout=layout.kamada.kawai, pctile_color=c(0.02,0.98)) {

	if (length(files) == 1 && file.info(files)$isdir) {
		files <- dir(SPADE.strip.sep(files),full.names=TRUE,pattern=glob2rx(file_pattern))
	}
	if (length(files) == 0) {
		stop("No input files found")
	} 
	out_dir <- SPADE.normalize.out_dir(out_dir)

	# Validate structure of panels
	if (!is.null(panels))  {
		if (!is.list(panels))
			stop("Invalid panels argument, see function documentation")
		lapply(panels, function(x) {
			if (!is.list(x) || !all(c("panel_files", "median_cols") %in% names(x)))
				stop("Invalid panel found, see function documentation for proper panel structure")
			if (!all(x$panel_files %in% basename(files)))
				stop("Panel files must be a subset of analysis files")
			if (!is.null(x$reference_files) && !all(x$reference_files %in% x$panel_files))
				stop("Panel reference files must be a subset of panel files")
		})
	}

	# Run downsampling/clustering/upsampling on all specified files 
	density_files <- c()
	sampled_files <- c()
	for (f in files) {
		cat("Downsampling file:",f,"\n")
		f_density <- paste(out_dir,basename(f),".density.fcs",sep="")
		f_sampled <- paste(out_dir,basename(f),".downsample.fcs",sep="")

		SPADE.addDensityToFCS(f, f_density, cols=cluster_cols, arcsinh_cofactor=arcsinh_cofactor, comp=comp)
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
					desired_samples=clustering_samples,
					comp=comp)

	sampled_files <- c()
	for (f in density_files) {
		cat("Upsampling file:",f,"\n")
		f_sampled <- paste(f,".cluster.fcs",sep="")
		SPADE.addClusterToFCS(f, f_sampled, cells_file, cols=cluster_cols, arcsinh_cofactor=arcsinh_cofactor,comp=comp)
		sampled_files <- c(sampled_files, f_sampled)
	}

	graph  <- read.graph(graph_file, format="gml")

	# Compute the layout once for the MST, write out for use by other functions
	layout_table <- layout(graph)
	if (deparse(substitute(layout)) != "SPADE.layout.arch")  # The igraph internal layouts are much more compact than arch.layout
		layout_table = layout_table * 50
	write.table(layout_table,paste(out_dir,file="layout.table",sep=""),row.names = FALSE,col.names = FALSE)
	
	# Track all attributes to compute global limits
	attr_values <- list()

	if (is.null(panels)) {  # Initialize panels if NULL
		panels <- list( list(panel_files=sampled_files, median_cols=NULL) )
	}
	
	for (p in panels) {
	
		reference_medians <- NULL
		if (!is.null(p$reference_files)) {
			reference_files   <- sapply(as.vector(p$reference_files), function(f) { sampled_files[grep(f, sampled_files)[1]] })
			reference_medians <- SPADE.markerMedians(reference_files, vcount(graph), cols=p$fold_cols, arcsinh_cofactor=arcsinh_cofactor, cluster_cols=cluster_cols, comp=comp)
		}

		for (f in as.vector(p$panel_files)) {
			f <- sampled_files[grep(f, sampled_files)[1]]

			# Compute the median marker intensities in each node, including the overall cell frequency per node	
			cat("Computing medians for file:",f,"\n")
			anno <- SPADE.markerMedians(f, vcount(graph), cols=p$median_cols, arcsinh_cofactor=arcsinh_cofactor, cluster_cols=cluster_cols, comp=comp)
					
			if (!is.null(reference_medians)) {	# If a reference file is specified								
				# Compute the fold change compared to reference medians
				cat("Computing fold change for file:",f,"\n")
				fold_anno <- SPADE.markerMedians(f, vcount(graph), cols=p$fold_cols, arcsinh_cofactor=arcsinh_cofactor, cluster_cols=cluster_cols, comp=comp)
				fold <- fold_anno$medians - reference_medians$medians
				
				ratio <- log10(fold_anno$percenttotal / reference_medians$percenttotal); 
				colnames(ratio) <- c("percenttotalratiolog")
				is.na(ratio) <- fold_anno$count == 0 | reference_medians$count == 0

				# Merge the fold-change columns with the count, frequency, and median columns
				anno <- c(anno, list(percenttotalratiolog = ratio, fold = fold))	
			}

			SPADE.write.graph(SPADE.annotateGraph(graph, layout=layout_table, anno=anno), paste(f,".medians.gml",sep=""), format="gml")

			# We save an R native version of the annotations to simpify plotting, and other downstream operations
			anno <- SPADE.flattenAnnotations(anno)
			for (c in colnames(anno)) { attr_values[[c]] <- c(attr_values[[c]], anno[,c]) }
			save(anno, file=paste(f,"anno.Rsave",sep="."))	
		}
	}
	
	# Compute the global limits (cleaning up attribute names to match those in GML files)
	attr_ranges <- t(sapply(attr_values, function(x) { quantile(x, probs=c(0.00, pctile_color, 1.00), na.rm=TRUE) }))
	rownames(attr_ranges) <- sapply(rownames(attr_ranges), function(x) { gsub("[^A-Za-z0-9_]","",x) })
	write.table(attr_ranges, paste(out_dir,"global_boundaries.table",sep=""), col.names=FALSE)

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


SPADE.plot.trees <- function(graph, files, file_pattern="*anno.Rsave", out_dir=".", layout=SPADE.layout.arch, attr_pattern="percent|median|fold", scale=NULL, pctile_color=c(0.02,0.98), normalize="global",size_scale_factor=1, edge.color="grey") {
    
	if (!is.igraph(graph)) {
		stop("Not a graph object")
    }	
	if (!is.null(scale) && (!is.vector(scale) || length(scale) !=2)) {
		stop("scale must be a two element vector")
	}
	if (!is.vector(pctile_color) || length(pctile_color) != 2) {
		stop("pctile_color must be a two element vector with values in [0,1]")
	}

	if (length(files) == 1 && file.info(files)$isdir) {
		files <- dir(SPADE.strip.sep(files),full.names=TRUE,pattern=glob2rx(file_pattern))    
  }
	out_dir <- SPADE.normalize.out_dir(out_dir)

	load_attr <- function(save_file) {
		l <- load(save_file)
		stopifnot(l == "anno")
		return(anno)
	}

	boundaries <- NULL
	if (normalize == "global") {
		boundaries <- c()  # Calculate ranges of all encountered attributes with trimmed outliers
		all_attrs  <- c()
		for (f in files) {
			attrs <- load_attr(f)
			for (i in grep(attr_pattern, colnames(attrs))) {
				n <- colnames(attrs)[i]
				all_attrs[[n]] <- c(all_attrs[[n]], attrs[,i])
			}
		}
		for (i in seq_along(all_attrs)) {
			boundaries[[names(all_attrs)[i]]] <- quantile(all_attrs[[i]], probs=pctile_color, na.rm=TRUE)
		}
	}

	if (is.function(layout))
		graph_l <- layout(graph)
	else
		graph_l <- layout	

    jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    colorscale <- jet.colors(100)

    for (f in files) {
		attrs <- load_attr(f)
	
		vsize <- attrs$percenttotal
		vsize <- vsize/(max(vsize,na.rm=TRUE)^(1/size_scale_factor)) * 3 + 2
		vsize[is.na(vsize)] <- 1

		for (i in grep(attr_pattern, colnames(attrs))) {
			name <- colnames(attrs)[i]
			
			# Compute the color for each vertex using color gradient scaled from min to max
			attr <- attrs[,i]
			if (!is.null(scale))
				boundary <- scale
			else if (normalize == "global") { # Scale to global min/max
				boundary <- boundaries[[name]] # Recall trimmed global boundary for this attribute
			} else # Scale to local min/max
				boundary <- quantile(attr, probs=pctile_color, na.rm=TRUE)  # Trim outliers for this attribtue
				
			if (boundary[1] == boundary[2]) {  boundary <- c(boundary[1]-1, boundary[2]+1); }  # Prevent "zero" width gradients
			if (length(grep("^median|percent", name)))
				boundary <- c(min(boundary), max(boundary))  # Dont make range symmetric for median or percent values
			else
				boundary <- c(-max(abs(boundary)), max(abs(boundary)))  # Make range symmetric for fold-change and ratio values
			
			boundary <- round(boundary, 2) # Round boundary values to 2 digits of precision so scale looks nice
				
			grad <- seq(boundary[1], boundary[2], length.out=length(colorscale))
		
			color <- colorscale[findInterval(attr, grad,all.inside=TRUE)]
			color[is.na(attr)] <- "grey"
			if (grepl("^logratio", name)) {
				# Color nodes with "infinite" ratios black
				color[is.na(attr) & attrs$count > 0] <- "black"
			}
				
			V(graph)$color <- color
	    
			# Plot the tree, with legend showing the gradient
			pdf(paste(out_dir,basename(f),".",name,".pdf",sep=""))
	    
			plot(graph, layout=graph_l, vertex.shape="circle", edge.color=edge.color, vertex.size=vsize, vertex.frame.color=NA, vertex.label=NA, edge.arrow.size=.25, edge.arrow.width=1) 
			
			# Substitute pretty attribute names
			if (length(grep("^median", name)))
				name <- sub("median", "Median of ", name)
			else if (length(grep("^fold", name)))
				name <- sub("fold", "Arcsinh diff. of ", name)
			else if (grepl("^percenttotal$", name))
				name <- sub("percent", "Percent freq. of ", name)
			else if (grepl("^percenttotalratiolog$", name))
				name <- "Log10 of Ratio of Percent Total of Cells in Each Cluster"

			# Make parameters used for clustering obvious
			if (grepl("_clust$", name))
				name <- sub("_clust", "\n(Used for tree-building)", name)
			
			title(main=paste(strsplit(basename(f),".fcs")[[1]][1], sub=name, sep="\n"))
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

SPADE.workflow.concat.FCS <- function(files, file_pattern="*.density.fcs.cluster.fcs", out_dir=".", cols=NULL, layout=SPADE.layout.arch, arcsinh_cofactor=5.0, in_graph_file="mst.gml", out_graph_file="concat.gml", cluster_cols=NULL, comp=TRUE)  {
	stop("Temporarily out of service. A candidate for deprecation.")
	if (length(files) == 1 && file.info(files)$isdir) {
		in_graph_file <- paste(SPADE.strip.sep(files),in_graph_file,sep=.Platform$file)
		files <- dir(SPADE.strip.sep(files),full.names=TRUE,pattern=glob2rx(file_pattern))
    }
	out_dir <- SPADE.normalize.out_dir(out_dir)

	anno  <- SPADE.markerMedians(files,cols=cols,arcsinh_cofactor=arcsinh_cofactor,cluster_cols=cluster_cols, comp=comp)
	graph <- read.graph(in_graph_file, format="gml")

	# Compute the overall cell frequency per node
	anno[["percent"]] <- anno$count / sum(anno$count) * 100; colnames(anno[["percent"]]) <- c("total");

	out_graph_file <- paste(out_dir,out_graph_file,sep="")

	SPADE.write.graph(SPADE.annotateGraph(graph, layout=layout, anno=anno), out_graph_file, format="gml")
}


