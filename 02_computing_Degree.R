#############
# GabrieleT #
#   2020    #
#############

#load libraries
library("igraph")
library("gtools")

#set seed
set.seed(131)

#set path
patH <- "/home/gabriele/final/"

#load data
filteredData <- read.csv(paste0(patH, "node_table.csv"), sep="\t", header=T, stringsAsFactors=T)

#read network
stringInfo <- read.csv(paste0(patH, "string_over400_experimental.tsv"), sep="\t", header=T, stringsAsFactors=T)

#create edge list
edgeL <- as.matrix(cbind.data.frame(as.character(stringInfo$X.node1), as.character(stringInfo$node2)))
#and build the network
networkTmp <- graph_from_edgelist(edgeL)
#get biggest component
network <- decompose(networkTmp)[[1]]

#remove the F's from filtered data
#colnames(filteredData[, c(14:19, 32:37)])
#colnames(filteredData[, -c(14:19, 32:37)])
noFs <- filteredData[, -c(14:19, 32:37)]

#adding a 1, to those letters like N, U e H
colnames(noFs) <- c("proteinName", "N.1", "N.2", "N.3", "N.4", "N.5", "N.6", "N.7", "N.8", "N.9", "N.10", "N.11", "N.12", "U.1", "U.2", "U.3", "U.4", "U.5", "U.6", "H.1", "H.2", "H.3", "H.4", "H.5", "H.6")

degrees <- vector()
#computing degree for each node in the network, using experimental values
for (vrtx in V(network)$name) {
	#get the neighbours of vertex
	neighList <- neighbors(network, vrtx, mode="all")
	#initialise array for the Degree of the vertex, for each columns, i.e. each sample
	sumUp <- rep(0, length(noFs[, -1]))
	#for each neighbour of the selected vertex
	for (ngh in neighList$name) {
		#get the array of experimental values of that neighbour for all the samples
		#then sum to the previous existing values
		sumUp <- sumUp + noFs[which(noFs$proteinName==ngh), -1]
	}
	#attach the new values in the degree table: another vertex's degree is computed
	degrees <- rbind.data.frame(degrees, sumUp)
}
#the degree is a small, floating point value since all the original values were normalised by the molecular weight!
degreeS <- cbind.data.frame(V(network)$name, degrees)
colnames(degreeS) <- colnames(noFs)

#export degree table
write.table(degreeS, paste0(patH, "degreeS.csv"), sep="\t", col.names=T, row.names=F, quote=F)

####### TRY OUT WITH A VERTEX
#vrtx <- "ETFA"
#neighList <- neighbors(network, vrtx, mode="all")
#sumUp <- rep(0, length(noFs[, -1]))
#for (ngh in neighList$name) {
#	#get the row storing its exp values to compute the Degree, and sum to the previous existing values
#	sumUp <- sumUp + noFs[which(noFs$proteinName==ngh), -1]
#}
#degreeS[which(degreeS$proteinName==vrtx), -1] == sumUp
#############################

#TESTING HIGH FOLD-CHANGE PROTEINS

#find proteins belonging to the network: not all the proteins were included since not part of the biggest connected component
inRows <- which((noFs$proteinName %in% V(network)$name)==T)
#get the proteins from the table
proteins <- noFs$proteinName[inRows]

#compute the mean experimental value, by treatment
onlyUs <- rowMeans(noFs[inRows, grep("U\\.", colnames(noFs))])
onlyNs <- rowMeans(noFs[inRows, grep("N\\.", colnames(noFs))])
onlyHs <- rowMeans(noFs[inRows, grep("H\\.", colnames(noFs))])

#compute fold change over the row means, using the control as a reference
NUfold <- log2(onlyUs/onlyNs)
NHfold <- log2(onlyHs/onlyNs)

#set threshold equal to 1, to compute foldChange
fcThreshU <- fcThreshH <- read.csv(paste0(patH, "FC_threshold.csv"), sep="\t", header=T, stringsAsFactors=F)

#return all proteins fold change
write.table(cbind.data.frame(proteins, NUfold), paste0(patH, "NU_foldchange_values_LOG.csv"), sep="\t", col.names=T, row.names=F, quote=F)
write.table(cbind.data.frame(proteins, NHfold), paste0(patH, "NH_foldchange_values_LOG.csv"), sep="\t", col.names=T, row.names=F, quote=F)

#return interesting proteins
write.table(proteins[abs(NUfold)>=fcThreshU$th], paste0(patH, "NU_foldchange_LOG.csv"), sep="\t", col.names=T, row.names=F, quote=F)
write.table(proteins[abs(NHfold)>=fcThreshH$th], paste0(patH, "NH_foldchange_LOG.csv"), sep="\t", col.names=T, row.names=F, quote=F)
