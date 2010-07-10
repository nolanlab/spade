# Transpose table before call to put in row major order
FlowSPD.density <- function(tbl, kernel_mult=5.0, apprx_mult=1.5, med_samples=2000)
    .Call("FSPD_density",t(tbl),kernel_mult,apprx_mult,med_samples)

FlowSPD.addDensityToFCS <- function(infilename, outfilename, cols=NULL,
arcsinh_cofactor=5.0, kernel_mult=5.0, apprx_mult=1.5, med_samples=2000)
{
    # Load in FCS file
    in_fcs  <- read.FCS(infilename);
    in_data <- exprs(in_fcs);

    params <- parameters(in_fcs);
    pd     <- pData(params);

    # Select out the desired columns
    idxs <- c();
    if (!is.null(cols)) {
	for (cl in cols) {
	    idxs <- c(idxs,which(pd[,"desc"] == cl,arr.ind<-TRUE));
	}
    } else {
	idxs <- c(1:nrow(pd));
    }  

    # Compute the density
    density <- FlowSPD.density(asinh(in_data[,idxs]/arcsinh_cofactor),kernel_mult<-kernel_mult,apprx_mult<-apprx_mult,med_samples<-med_samples);

    # Reload FCS file without transformation, so it can be accurately rewritten...
    in_fcs <- read.FCS(infilename,transform<-FALSE)
    in_data <- exprs(in_fcs);

    # Add column named "density" to the FCS file
    channel_number <- ncol(in_fcs)+1;
    channel_id     <- paste("$P",channel_number,sep<-"");
    channel_name   <- "density";
    channel_range  <- max(density)+1;

    plist <- matrix(c(channel_name,channel_name,channel_range,0,channel_range-1));
    rownames(plist) <- c("name","desc","range","minRange","maxRange");
    colnames(plist) <- c(channel_id);

    pd <- rbind(pd,t(plist));
    pData(params) <- pd;

    out_data <- cbind(in_data,"density"<-density);
    out_frame <- flowFrame(out_data,params,description<-description(in_fcs));

    keyval<-list();
    keyval[[paste("$P",channel_number,"B",sep<-"")]] <- "32";
    keyval[[paste("$P",channel_number,"R",sep<-"")]] <- toString(channel_range);
    keyval[[paste("$P",channel_number,"E",sep<-"")]] <- "0,0";
    keyval[[paste("$P",channel_number,"N",sep<-"")]] <- channel_name;
    keyword(out_frame) <- keyval;

    write.FCS(out_frame,outfilename);
}
