SPADE.strip.sep <- function(name) {
    ifelse(substr(name,nchar(name),nchar(name))==.Platform$file,substr(name,1,nchar(name)-1),name)
}

SPADE.driver <- function(files, file_pattern="*.fcs", out_dir=".", cluster_cols=NULL, arcsinh_cofactor=5.0, layout=SPADE.layout.arch, median_cols=NULL, reference_file = NULL, fold_cols=NULL, downsampling_samples=20000, k=200) {
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
	SPADE.downsampleFCS(f_density, f_sampled, desired_samples=downsampling_samples)
	
	density_files <- c(density_files, f_density)
	sampled_files <- c(sampled_files, f_sampled)	
    }
    
    cat("Clustering files...\n")
    cells_file <- paste(out_dir,"clusters.fcs",sep="")
    clust_file <- paste(out_dir,"clusters.table",sep="")
    graph_file <- paste(out_dir,"mst.gml",sep="")
    SPADE.FCSToTree(sampled_files, cells_file, graph_file, clust_file, cols=cluster_cols, arcsinh_cofactor=arcsinh_cofactor, k=k)

    sampled_files <- c()
    for (f in density_files) {
	cat("Upsampling file:",f,"\n")
	f_sampled <- paste(f,".cluster.fcs",sep="")
	SPADE.addClusterToFCS(f, f_sampled, cells_file, cols=cluster_cols, arcsinh_cofactor=arcsinh_cofactor)
	sampled_files <- c(sampled_files, f_sampled)
    }

    graph  <- read.graph(graph_file, format="gml")
    layout <- layout(graph)
    
    reference_medians = NULL
    if (!is.null(reference_file)) {
	reference_file <- paste(out_dir,reference_file,".density.fcs.cluster.fcs",sep="");
	reference_medians <- SPADE.markerMedians(reference_file, cols=fold_cols, arcsinh_cofactor=arcsinh_cofactor)
    }

    for (f in sampled_files) {

	cat("Computing medians for file:",f,"\n")
	a <- SPADE.markerMedians(f, cols=median_cols, arcsinh_cofactor=arcsinh_cofactor)
	g <- SPADE.annotateGraph(graph, layout=layout, a)
	SPADE.write.graph(g, paste(f,".medians.gml",sep=""), format="gml")
	
	if (!is.null(reference_medians) && f != reference_file) {
	    a <- SPADE.markerMedians(f, cols=fold_cols, arcsinh_cofactor=arcsinh_cofactor)
	    # Not all files might have all clusters represented
	    cc <- match(rownames(a$medians),rownames(reference_medians$medians))  # Common clusters
	    a <- list(count=a$count, fold=(a$medians-reference_medians$medians[cc,]))
	    g <- SPADE.annotateGraph(graph, layout=layout, a)
	    SPADE.write.graph(g, paste(f,".fold.gml",sep=""), format="gml")
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


SPADE.plot.trees <- function(files, file_pattern="*.gml", out_dir=".", layout=SPADE.layout.arch, attr_pattern="median|fold") {
    if (length(files) == 1 && file.info(files)$isdir) {
	files <- dir(SPADE.strip.sep(files),full.names=TRUE,pattern=glob2rx(file_pattern))    
    }

    out_dir_info <- file.info(out_dir)
    if (is.na(out_dir_info$isdir)) {
	dir.create(out_dir)
    }
    if (!file.info(out_dir)$isdir) {
	stop(paste("out_dir:",out_dir,"is not a directory"))
    }
    out_dir <- paste(SPADE.strip.sep(out_dir),.Platform$file,sep="")

    jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    colorscale <- jet.colors(100)

    for (f in files) {
	graph <- read.graph(f, format="gml")
	attrs <- list.vertex.attributes(graph)

	if (is.function(layout))
	    graph_l <- layout(graph)
	else
	    graph_l <- layout	

	vsize <- V(graph)$count
	vsize[vsize == Inf | vsize == -Inf] <- NA  # Clean up bogus values	
	vsize <- vsize/max(vsize,na.rm=TRUE) * 3 + 2
	vsize[is.na(vsize)] <- 1

	for (i in grep(attr_pattern, attrs)) {
	    # Compute the color for each vertex using color gradient
	    attr <- get.vertex.attribute(graph,attrs[i])
	    attr[attr == Inf | attr == -Inf] <- NA  # Clean up bogus values ...
	    
	    grad <- seq(min(attr,na.rm=TRUE),max(attr,na.rm=TRUE),length.out=length(colorscale))
	    color <- colorscale[findInterval(attr, grad)]
	    color[is.na(attr)] <- "grey"
	    V(graph)$color <- color
	    
	    # Plot the tree, with legend showing the gradient
	    pdf(paste(out_dir,basename(f),".",attrs[i],".pdf",sep=""))
	    
	    plot(graph, layout=graph_l, vertex.shape="csquare", edge.color="grey", vertex.size=vsize, vertex.frame.color=NA, vertex.label=NA)
	    
	    title(main=attrs[i])
	    subplot(image(grad,c(1),matrix(1:length(colorscale),ncol=1),col=colorscale,xlab="",ylab="",yaxt="n",xaxp=c(round(min(attr,na.rm=TRUE),2),round(max(attr,na.rm=TRUE),2),1)),x="right,bottom",size=c(1,.20))
	    
	    dev.off()
	}
    }

}
