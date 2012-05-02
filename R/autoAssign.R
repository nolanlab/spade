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
	x = sub(" ","_",x);
  x = sub(".txt","",x)
}


evaluateCellTypeRule = function(fcsFileName,transformCols,dataDirectory,ruleDir,ruleFilename) {
#fcsFileName = "Marrow1_01_Unstim1_with_Manual_cell_type_ID.fcs.density.fcs.cluster.fcs"
#transformCols = c("Ir(190.960)-Dual","Cd(110,111,112,114)","Er(166.932)-Dual","Er(169.935)-Dual","Gd(157.924)-Dual","Gd(159.927)-Dual","In(114.903)-Dual","La(138.906)-Dual","Nd(141.907)-Dual","Nd(143.910)-Dual","Nd(144.912)-Dual","Nd(145.913)-Dual","Nd(147.916)-Dual","Sm(146.914)-Dual")


#dataDirectory = "/Users/rbruggner/Desktop/clusterForGarry/spadeTest/output/"
#key = read.table(paste(dataDirectory,"../../manualKey.txt",sep=""),sep="\t",header=1)
fcs = read.FCS(paste(dataDirectory,fcsFileName,sep=""))


data = exprs(fcs);
data[,transformCols] = asinh(data[,transformCols]/5)
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

percents = sapply(transformCols,markerPercents,clusterIds=clusterIds,data=data,medians=medians,clusterIndex=clusterIndex)
rownames(percents) = clusterIds
randomPercents = sapply(transformCols,markerPercents,clusterIds=clusterIds,data=dataRandomized,medians=medians,clusterIndex=clusterIndex)
rownames(percents) = clusterIds

keyMST = read.graph(file=paste(dataDirectory,"mst.gml",sep=""),format="gml")
layout_table = as.matrix(read.table(file=paste(dataDirectory,"layout.table",sep="")))

fileMST = read.graph(file=paste(dataDirectory,fcsFileName,".medians.gml",sep=""),format="gml")

#ruleDir = "/Users/rbruggner/Desktop/clusterForGarry/populationRules/";

#for (ruleFile in list.files(ruleDir,pattern=".txt")){
	cat(paste("RULE:",ruleFile,"\n"));
	populationRules = read.table(paste(ruleDir,ruleFile,sep=""),sep="\t",header=F,col.names=c("parameter","hilo","desc"),comment.char="#")
	probability=1;
	randomProbability=1;
	for (ruleId in 1:nrow(populationRules)){
		parameter = as.vector(populationRules[ruleId,"parameter"])
	  if (populationRules[ruleId,"hilo"]=="+"){
	  	print(paste(parameter,"Hi"))
	  	probability = probability*percents[,parameter];
		  randomProbability = randomProbability*randomPercents[,parameter];
	  } else {
	  	print(paste(parameter,"Low"))
	  	probability = probability*(1-percents[,parameter]);
		  randomProbability = randomProbability*(1-randomPercents[,parameter]);
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

SPADE.write.graph(SPADE.annotateGraph(fileMST, layout=layout_table, anno=annotation),file=paste(dataDirectory,fcsFileName,".medians.gml",sep=""),format="gml")
}
