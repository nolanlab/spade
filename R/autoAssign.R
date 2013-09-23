percentCluster = function(clusterId,parameter,data,medians,clusterIndex){
	percent = sum(data[data[,clusterIndex]==clusterId,parameter]>medians[parameter])/sum(data[,clusterIndex]==clusterId)
	return(percent)
}

markerPercents = function(parameter,clusterIds,data,medians,clusterIndex){
	sapply(clusterIds,percentCluster,parameter=parameter,data=data,medians=medians,clusterIndex=clusterIndex)
}

indexOf = function(x,data){
	which(x==data)
}

cleanStr = function(x){
	x = gsub("-","_",x)
	x = gsub(" ","_",x)
	x = gsub(".txt","",x)
}

SPADE.evaluateCellTypeRule = function(out_dir,fcsFileName,ruleCols,ruleDir,ruleFile) {
	#fcsFileName = "Marrow1_01_Unstim1_with_Manual_cell_type_ID.fcs.density.fcs.cluster.fcs"
	#ruleCols = c("Ir(190.960)-Dual","Cd(110,111,112,114)","Er(166.932)-Dual","Er(169.935)-Dual","Gd(157.924)-Dual","Gd(159.927)-Dual","In(114.903)-Dual","La(138.906)-Dual","Nd(141.907)-Dual","Nd(143.910)-Dual","Nd(144.912)-Dual","Nd(145.913)-Dual","Nd(147.916)-Dual","Sm(146.914)-Dual")


	#dataDirectory = "/Users/rbruggner/Desktop/clusterForGarry/spadeTest/output/"
	#key = read.table(paste(dataDirectory,"../../manualKey.txt",sep=""),sep="\t",header=1)
	fcs = SPADE.read.FCS(paste(fcsFileName,sep=""))


	data = exprs(fcs);
	data[,ruleCols] = asinh(data[,ruleCols]/5)
	clusterIndex = which(colnames(data)=="cluster")


	#dataRandomized = data
	#dataRandomized[,clusterIndex] = dataRandomized[sample(1:nrow(dataRandomized),nrow(dataRandomized)),clusterIndex]

	medians = apply(data,2,median)


	clusterIds = sort(unique(data[,clusterIndex]))


	#type=list()
	#for (i in clusterIds){
	#	type[[i]] = table(key[data[data[,clusterIndex]==i,"cell_type_index"],2])
	#}
	#tableNames=names(type[[1]])
	#type = do.call("rbind",type)
	#colnames(type)=tableNames

	percents = sapply(ruleCols,markerPercents,clusterIds=clusterIds,data=data,medians=medians,clusterIndex=clusterIndex)
	rownames(percents) = clusterIds
	#randomPercents = sapply(ruleCols,markerPercents,clusterIds=clusterIds,data=dataRandomized,medians=medians,clusterIndex=clusterIndex)
	#rownames(percents) = clusterIds
	message(paste("GML FILE:",paste(out_dir,"mst.gml",sep=""),"\n"))
	keyMST = read.graph(file=paste(out_dir,"mst.gml",sep=""),format="gml")
	layout_table = as.matrix(read.table(file=paste(out_dir,"layout.table",sep="")))
	message(paste(fcsFileName,".medians.gml",sep=""))
	fileMST = read.graph(file=paste(fcsFileName,".medians.gml",sep=""),format="gml")

	#ruleDir = "/Users/rbruggner/Desktop/clusterForGarry/populationRules/";

	#for (ruleFile in list.files(ruleDir,pattern=".txt")){
		cat(paste("RULE:",ruleFile,"\n"));
		populationRules = read.table(paste(ruleDir,ruleFile,sep=""),sep="\t",header=F,col.names=c("parameter","hilo","protein"),comment.char="#")
		probability=1;
		#randomProbability=1;
		for (ruleId in 1:nrow(populationRules)){
			parameter = as.vector(populationRules[ruleId,"protein"])
			if (populationRules[ruleId,"hilo"]=="+"){
				print(paste(parameter,"Hi"))
				probability = probability*percents[,ruleId];
				#randomProbability = randomProbability*randomPercents[,parameter];
			} else {
				print(paste(parameter,"Low"))
				probability = probability*(1-percents[,ruleId]);
				#randomProbability = randomProbability*(1-randomPercents[,parameter]);
		  	}
		}
		cat("\n");
	  #sort(probability)
	# 	highestScorers = order(probability,decreasing=T)
	  
	#  scores=type[highestScorers[1:10],]
	#  rownames(scores)=paste("Cluster",highestScorers[1:10])
	#  write.table(scores,file=paste(dataDirectory,"scores/scores_",ruleFile,".csv",sep=""),sep=",")
	 
		
		attributeName=paste(cleanStr(ruleFile),"probability",sep="_");
		if (!exists("probabilityMatrix")){
			probabilityMatrix = matrix(probability,ncol=1,dimnames=list(clusterIds,attributeName));
		} else {
			probabilityMatrix = cbind(probabilityMatrix,probability);
			colnames(probabilityMatrix)[ncol(probabilityMatrix)]=attributeName
		}
	#}

	annotation=list();
	annotation[["autoassign_probability"]]=probabilityMatrix
	gmlFilename = paste(fcsFileName,".medians.gml",sep="")
	SPADE.write.graph(SPADE.annotateGraph(fileMST, layout=layout_table, anno=annotation),file=gmlFilename, format="gml")

	cat("Create merge order file.\n")
	cellType=cleanStr(ruleFile)
	SPADE.createMergeOrderByCellProbability(out_dir=out_dir, gmlFilename=gmlFilename, cellType=cellType, outputFileName=paste("MergeOrder", gsub(" ", "_", ruleFile), sep="_"))
}

SPADE.createMergeOrderByCellProbability = function(out_dir, gmlFilename, cellType, outputFileName) {
		
	cat(paste("out_dir:", out_dir, "gmlFilename:", gmlFilename, "cellType:", cellType, "outputFileName:", outputFileName, "\n"))
	
	#Read gml file
	g=read.graph(paste(gmlFilename,sep=""),format="gml")
	
	#Create dataframe with node name and the cell type probability
	probabilityParameter = paste("autoassign_probability", cellType, "_probability", sep="")
	probabilities=get.vertex.attribute(g, probabilityParameter)
	df=data.frame(node=V(g)$id,probability=probabilities)

	#Sort by descending probability
	df=df[order(df$probability,decreasing=TRUE),]
	
	#Cluster the nodes
	cluster_count = 0
	output = list()
	cluster_count = 1
	src = as.numeric(as.character(df$node[[1]]))	# this is the source "id" and not "name"
	tgt = as.numeric(as.character(df$node[[2]]))	# this is the target "id" and not "name"
	output[[cluster_count]] = c(-src,-tgt)
	cluster_count = cluster_count + 1
	for (i in 3:nrow(df)) {
		src = -1 * as.numeric(as.character(df$node[[i]]))	# this is the source "id" and not "name"
		tgt = cluster_count - 1
		output[[cluster_count]] = c(src, tgt)
		cluster_count = cluster_count + 1
	}
	
	#Write out the file
	output_df=as.data.frame(do.call("rbind", output), stringsAsFactors = FALSE)
	write.table(output_df, paste(out_dir, outputFileName, sep=""))
}
