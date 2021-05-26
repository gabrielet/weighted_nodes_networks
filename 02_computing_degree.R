#############
# GabrieleT #
#   2021    #
#############

# load libraries
library("igraph")
library("gtools")

# set seed
set.seed(131)

# set path
patH <- "~/WNNets/"

# load weights data set
filteredData <- read.csv(paste0(patH, "node_table.csv"), sep="\t", header=T, stringsAsFactors=T)

# load STRING network
stringInfo <- read.csv(paste0(patH, "string_interactions_over400_experimental.tsv"), sep="\t", header=T, stringsAsFactors=T)

# create edge list
edgeL <- as.matrix(cbind.data.frame(as.character(stringInfo$X.node1), as.character(stringInfo$node2)))
# build the network
networkTmp <- graph_from_edgelist(edgeL)
# and get biggest component
network <- decompose(networkTmp)[[1]]

# export the network, for Cytoscape
el <- as_edgelist(network)
toCyto <- data.frame(source=V(network)[el[, 1]]$name, target=V(network)[el[, 2]]$name)
write.table(toCyto, "cytoscape_network.csv", row.names=F, quote=F, sep="\t")


######################### THIS IS DONE TO OBTAIN THE BEFORE AND AFTER NORMALISATION TABLES ONLY
# export table before and after normalisation
beforeTmp <- read.csv(paste0(patH, "before_temporary.csv"), sep="\t", header=T, stringsAsFactors=T)
before <- beforeTmp[which(beforeTmp$proteins %in% unique(c(V(network)[el[,1]]$name, V(network)[el[,2]]$name))), ]
write.table(before, paste0(patH, "before_non_zero_and_normalisation.csv"), sep="\t", col.names=T, row.names=F, quote=F)

after <- filteredData[which(filteredData$proteinName %in% unique(c(V(network)[el[,1]]$name, V(network)[el[,2]]$name))), ]
write.table(after, paste0(patH, "after_non_zero_and_normalisation.csv"), sep="\t", col.names=T, row.names=F, quote=F)
###############################################################################################

# adding a 1, to those letters like N, U e H
colnames(filteredData) <- c("proteinName", "N.1", "N.2", "N.3", "N.4", "N.5", "N.6", "N.7", "N.8", "N.9", "N.10", "N.11", "N.12", "U.1", "U.2", "U.3", "U.4", "U.5", "U.6", "H.1", "H.2", "H.3", "H.4", "H.5", "H.6")

# initialise degree storage
degrees <- vector()
# computing degree for each node in the network, using experimental values
for (vrtx in V(network)$name) {
	# get the neighbours of vertex
	neighList <- neighbors(network, vrtx, mode="all")
	# initialise array for the Degree of the vertex, for each columns, i.e. each sample
	sumUp <- rep(0, length(filteredData[, -1]))
	# for each neighbour of the selected vertex
	for (ngh in neighList$name) {
		# get the array of experimental values of that neighbour for all the samples
		# then sum to the previous existing values
		sumUp <- sumUp + filteredData[which(filteredData$proteinName==ngh), -1]
	}
	# attach the new values in the degree table: another vertex degree is computed
	degrees <- rbind.data.frame(degrees, sumUp)
}
# the degree is a small, floating point value since all the original values were normalised by the molecular weight!
degreeS <- cbind.data.frame(V(network)$name, degrees)
colnames(degreeS) <- colnames(filteredData)

# export degree table
write.table(degreeS, paste0(patH, "degreeS.csv"), sep="\t", col.names=T, row.names=F, quote=F)

# TRY OUT WITH A VERTEX
# vrtx <- "ETFA"
# neighList <- neighbors(network, vrtx, mode="all")
# sumUp <- rep(0, length(filteredData[, -1]))
# for (ngh in neighList$name) {
#  	# get the row storing its exp values to compute the Degree, and sum to the previous existing values
#  	sumUp <- sumUp + filteredData[which(filteredData$proteinName==ngh), -1]
# }
# degreeS[which(degreeS$proteinName==vrtx), -1] == sumUp

# TESTING SIGNIFICANT FOLD-CHANGE PROTEINS, i.e. SFCPs

# find proteins that belong to the network: not all the proteins were included since not part of the biggest connected component
inRows <- which((filteredData$proteinName %in% V(network)$name)==T)

# compute the mean experimental value, by row

# for the U treatment
onlyUs <- rowMeans(filteredData[inRows, grep("U\\.", colnames(filteredData))])
# and corresponding healthy samples
# colnames(filteredData[inRows, c(2:7)])
onlyNUs <- rowMeans(filteredData[inRows, c(2:7)])

# for the H treatment
onlyHs <- rowMeans(filteredData[inRows, grep("H\\.", colnames(filteredData))])
# and corresponding healthy samples
# colnames(filteredData[inRows, c(8:13)])
onlyNHs <- rowMeans(filteredData[inRows, c(8:13)])

# compute fold change over the row means, using the control samples as a reference
NUfold <- log2(onlyUs/onlyNUs)
NHfold <- log2(onlyHs/onlyNHs)

# now compute t-test over the average experimental values, to test which proteins are significant with respect to fold change
# get the proteins names from filteredData
proteins <- filteredData$proteinName[inRows]

# get Us
pvalNU <- vector()
for (prt in seq(1, length(proteins), by=1)) {
	# get SpC for a protein, all the 6 treated samples
	subSetU <- filteredData[inRows, grep("U\\.", colnames(filteredData))]
	# same for healthy subs
	subSetNU <- filteredData[inRows, c(2:7)]
	# compute the statistics
	pvalNU <- append(pvalNU, t.test(subSetU[prt, ], subSetNU[prt, ])$p.value)
}

# get Hs
pvalNH <- vector()
for (prt in seq(1, length(proteins), by=1)) {
	# get SpC for a protein, all the 6 treated samples
	subSetH <- filteredData[inRows, grep("H\\.", colnames(filteredData))]
	# same for healthy subs
	subSetNH <- filteredData[inRows, c(8:13)]
	# compute the statistics
	pvalNH <- append(pvalNH, t.test(subSetH[prt, ], subSetNH[prt, ])$p.value)
}

# compute adjusted pvalues using false discovery rate
# load thresholds file
fcThresh <- read.csv("FC_FDR_threshold.csv", header=T, sep="\t")
# set threshold for FDR
tH <- fcThresh$FDR
# set FDR as method for correction
mth <- "fdr"
# perform correction using package stats
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust
padj_H <- p.adjust(pvalNH, method=mth)
padj_U <- p.adjust(pvalNU, method=mth)

# now put everything together
allSFCP <- cbind.data.frame(proteins, NUfold, padj_U, NHfold, padj_H)

# get SFCPs names
NH_names <- allSFCP$proteins[which(allSFCP[abs(allSFCP$NHfold) > fcThresh$foldC, 3] < 0.1)]
NU_names <- allSFCP$proteins[which(allSFCP[abs(allSFCP$NUfold) > fcThresh$foldC, 3] < 0.1)]

# and return their names
write.table(NU_names, paste0(patH, "NU_SFCPs.csv"), sep="\t", col.names=T, row.names=F, quote=F)
write.table(NH_names, paste0(patH, "NH_SFCPs.csv"), sep="\t", col.names=T, row.names=F, quote=F)
