percentCluster <- function(clusterId,parameter,data,medians,clusterIndex){
	percent = sum(data[data[,clusterIndex]==clusterId,parameter]>medians[parameter])/sum(data[,clusterIndex]==clusterId)
		return(percent)
}

markerPercents <- function(parameter,clusterIds,data,medians,clusterIndex){
	sapply(clusterIds,percentCluster,parameter=parameter,data=data,medians=medians,clusterIndex=clusterIndex)
}

indexOf <- function(x,data){
	which(x==data)
}

cleanStr <- function(x){
	x = sub(" ","_",x);
	x = sub(".txt","",x)
}


SPADE.evaluateCellTypeRule <- function(out_dir,fcsFileName,ruleCols,ruleDir,ruleFile) {
	fcs = SPADE.read.FCS(paste(fcsFileName,sep=""))

	data = exprs(fcs);
	data[,ruleCols] = asinh(data[,ruleCols]/5)
	clusterIndex = which(colnames(data)=="cluster")


	medians = apply(data,2,median)


	clusterIds = sort(unique(data[,clusterIndex]))

	percents = sapply(ruleCols,markerPercents,clusterIds=clusterIds,data=data,medians=medians,clusterIndex=clusterIndex)
	rownames(percents) = clusterIds
	message(paste("GML FILE:",paste(out_dir,"mst.gml",sep=""),"\n"))
	keyMST = read.graph(file=paste(out_dir,"mst.gml",sep=""),format="gml")
	layout_table = as.matrix(read.table(file=paste(out_dir,"layout.table",sep="")))
	message(paste(fcsFileName,".medians.gml",sep=""))
	fileMST = read.graph(file=paste(fcsFileName,".medians.gml",sep=""),format="gml")


	cat(paste("RULE:",ruleFile,"\n"));
	populationRules = read.table(paste(ruleDir,ruleFile,sep=""),sep="\t",header=F,col.names=c("parameter","hilo","desc"),comment.char="#")
	probability=1;
	for (ruleId in 1:nrow(populationRules)){
		parameter = as.vector(populationRules[ruleId,"parameter"])
		if (populationRules[ruleId,"hilo"]=="+"){
			print(paste(parameter,"Hi"))
				probability = probability*percents[,ruleId];
		} else {
			print(paste(parameter,"Low"))
			probability = probability*(1-percents[,ruleId]);
		}
	}
	cat("\n");

	attributeName=paste(cleanStr(ruleFile),"probability",sep="_");
	if (!exists("probabilityMatrix")){
		probabilityMatrix = matrix(probability,ncol=1,dimnames=list(clusterIds,attributeName));
	} else {
		probabilityMatrix = cbind(probabilityMatrix,probability);
		colnames(probabilityMatrix)[ncol(probabilityMatrix)]=attributeName
	}

	annotation=list();
	annotation[["autoassign_probability"]]=probabilityMatrix
	SPADE.write.graph(SPADE.annotateGraph(fileMST, layout=layout_table, anno=annotation),file=paste(fcsFileName,".medians.gml",sep=""),format="gml")
}
