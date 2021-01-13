#############
# GabrieleT #
#   2020    #
#############

#load libraries
library("igraph")
library("gtools")
library("ggplot2")
library("RColorBrewer")
library("extrafont")
loadfonts()

#set seed
set.seed(131)

#set paths
patH <- "/home/cbmc/wnnets/"
saveIn <- "/home/cbmc/wnnets/figures/"

#set palette with 8 colours for colouring plots
colorsPal <- brewer.pal(8, "Set2")

#load data
filteredData <- read.csv(paste0(patH, "node_table.csv"), sep="\t", header=T, stringsAsFactors=F)
#remove the F's from filtered data
noFs <- filteredData[, -c(14:19, 32:37)]

#read network
stringInfo <- read.csv(paste0(patH, "string_over400_experimental.tsv"), sep="\t", header=T, stringsAsFactors=T)

#create edge list
edgeL <- as.matrix(cbind.data.frame(as.character(stringInfo$X.node1), as.character(stringInfo$node2)))
#and build the network
networkTmp <- graph_from_edgelist(edgeL)
#get biggest component
network <- decompose(networkTmp)[[1]]

#initialise the storage
allDegsVect <- vector()
allDegsList <- list(list())

#vectors of pvalues
pvalNH <- vector()
pvalNU <- vector()
infoPNH <- vector()
infoPNU <- vector()

#generating 100 trials but using only 10 to plot the heatmap
trials <- 100
plotted <- 10

#test n=trials different random datasets
for (r in seq(1, trials, by=1)) {

	#compute random values for the random weighting matrix
	randomData <- cbind.data.frame(noFs$proteinName, matrix(runif(dim(noFs[, -1])[1] * dim(noFs[, -1])[2], min(noFs[, -1]), max(noFs[, -1])), ncol=dim(noFs[, -1])[2], nrow=dim(noFs[, -1])[1],), stringsAsFactors=F)
	colnames(randomData) <- colnames(noFs)

	print(paste0("iteration ", r))
	degrees <- vector()
	#computing degree for each node in the network, using experimental values
	for (vrtx in V(network)$name) {
		#get the neighbours of vertex
		neighList <- neighbors(network, vrtx, mode="all")
		#initialise array for the Degree of the vertex, for each columns, i.e. each sample
		sumUp <- rep(0, length(randomData[, -1]))
		#for each neighbour of the selected vertex
		for (ngh in neighList$name) {
			#get the array of experimental values of that neighbour for all the samples
			#then sum to the previous existing values
			sumUp <- sumUp + randomData[which(randomData$proteinName==ngh), -1]
		}
		#attach the new values in the degree table: another vertex's degree is computed
		degrees <- rbind.data.frame(degrees, sumUp)
	}
	#the degree is a small, floating point value since all the original values were normalised by the molecular weight!
	degreesGabri <- cbind.data.frame(V(network)$name, degrees)
	colnames(degreesGabri) <- colnames(randomData)
	#create net label
	nlabl <- rep(r, nrow(degrees))
	#unlist degrees and merge
	allDegsVect <- rbind.data.frame(allDegsVect, cbind.data.frame(protein=degreesGabri[, 1], rNet=nlabl, value=unlist(degreesGabri[, -1]), stringsAsFactors=F))
	allDegsList[[r]] <- degreesGabri

	#get only degree values
	degOnly <- degreesGabri[, -1]
	#and protein names
	proteins <- as.factor(degreesGabri$proteinName)

	# NOTE:
	# SUBSETTING CONSIDERING THAT THE COMPARISON MUST FOLLOW THESE RULES:

	#N7-N12 VERSUS H1-H6
	subSampleForH <- cbind.data.frame(proteins, degOnly[, colnames(degOnly[, grep("N", colnames(degOnly))])[c(7:12)]], degOnly[, grep("H", colnames(degOnly))])

	#N1-N6 VERSUS U1-U6
	subSampleForU <- cbind.data.frame(proteins, degOnly[, colnames(degOnly[, grep("N", colnames(degOnly))])[c(1:6)]], degOnly[, grep("U", colnames(degOnly))])

	#computing pvalues
	for (prt in proteins) {
		subsampleH <- as.matrix(subSampleForH[which(subSampleForH$protein==prt), -1])
		#get the N values from sumsampleForH
		valN <- subsampleH[, grep("N", colnames(subsampleH))]
		valH <- subsampleH[, grep("H", colnames(subsampleH))]
		#compute t.test
		pv <- t.test(valN, valH)$p.value
		infoPNH <- rbind.data.frame(infoPNH, cbind.data.frame(protein=prt, pval=pv, stringsAsFactors=F))
		
		subsampleU <- as.matrix(subSampleForU[which(subSampleForU$protein==prt), -1])
		#get the N values from sumsampleForU
		valN <- subsampleU[, grep("N", colnames(subsampleU))]
		valU <- subsampleU[, grep("U", colnames(subsampleU))]
		
		pv <- t.test(valN, valU)$p.value
		infoPNU <- rbind.data.frame(infoPNU, cbind.data.frame(protein=prt, pval=pv, stringsAsFactors=F))
	}
}

##now that it's truly random, it works!

#create the frequency table for the proteins, U
iPNU <- as.data.frame(table(infoPNU$protein[which(infoPNU$pval<=0.05)]))
colnames(iPNU) <- c("protein", "frequency")

#and H
iPNH <- as.data.frame(table(infoPNH$protein[which(infoPNH$pval<=0.05)]))
colnames(iPNH) <- c("protein", "frequency")

#and finally export the number of times a protein is found SSP
write.table(iPNU, paste0(patH, "summaryInfoPNU.csv"), quote=F, row.names=F, col.names=F)
write.table(iPNH, paste0(patH, "summaryInfoPNH.csv"), quote=F, row.names=F, col.names=F)

#set some textual parameters
sizeText <- 20

#lolliplotting for U
orderedIPNU <- iPNU[order(iPNU$frequency), ]
orderedIPNU$protein <- factor(orderedIPNU$protein, levels=unique(orderedIPNU$protein))
#compute the average value
avgPNU <- mean(table(infoPNU$protein[which(infoPNU$pval<=0.05)]))

#generate the plot
plotL <- ggplot(orderedIPNU, aes(x=protein, y=frequency)) +
	geom_segment(aes(x=protein, xend=protein, y=0, yend=frequency)) +
	ylim(c(0, 11)) +
	geom_point(size=3, color=colorsPal[2], fill=alpha(colorsPal[6], 0.3), alpha=0.9, shape=21, stroke=2) +
	coord_flip() +	
	annotate(geom="text", x=4, y=4, label=paste0("proteins count=", length(unique(infoPNU$protein[which(infoPNU$pval<=0.05)]))), hjust="left", size=8) +
	annotate(geom="text", x=2, y=4, label=paste0("average frequency=", round(avgPNU, digits=2)), hjust="left", size=8) +
	theme_bw() +
	theme(plot.title=element_text(size=sizeText, face="bold", hjust=0.5),
		axis.text.x=element_text(hjust=1, size=sizeText, face="bold"),
		axis.text.y=element_text(size=sizeText, face="bold"),
		axis.text=element_text(size=sizeText, face="bold"),
		axis.title=element_text(size=sizeText, face="bold"),
		legend.position="top",
		legend.title=element_text(size=sizeText),
		legend.text=element_text(size=sizeText))

#and export it
pdf(paste0(saveIn, "frequency_lolli_U.pdf"), width=9, height=12, font="Arial")
plot(plotL)
dev.off()

#lolliplotting for H
orderedIPNH <- iPNH[order(iPNH$frequency), ]
orderedIPNH$protein <- factor(orderedIPNH$protein, levels=unique(orderedIPNH$protein))
#compute the average value
avgPNH <- mean(table(infoPNH$protein[which(infoPNH$pval<=0.05)]))

#generate the plot
plotL <- ggplot(orderedIPNH, aes(x=protein, y=frequency)) +
	geom_segment(aes(x=protein, xend=protein, y=0, yend=frequency)) +
	ylim(c(0, 11)) +
	geom_point(size=3, color=colorsPal[3], fill=alpha(colorsPal[4], 0.3), alpha=0.9, shape=21, stroke=2) +
	coord_flip() +
	annotate(geom="text", x=4, y=4, label=paste0("proteins count=", length(unique(infoPNH$protein[which(infoPNH$pval<=0.05)]))), hjust="left", size=8) +
	annotate(geom="text", x=2, y=4, label=paste0("average frequency=", round(avgPNH, digits=2)), hjust="left", size=8) +
	theme_bw() +
	theme(plot.title=element_text(size=sizeText, face="bold", hjust=0.5),
		axis.text.x=element_text(hjust=1, size=sizeText, face="bold"),
		axis.text.y=element_text(size=sizeText, face="bold"),
		axis.text=element_text(size=sizeText, face="bold"),
		axis.title=element_text(size=sizeText, face="bold"),
		legend.position="top",
		legend.title=element_text(size=sizeText),
		legend.text=element_text(size=sizeText))

#and export it
pdf(paste0(saveIn, "frequency_lolli_H.pdf"), width=9, height=12, font="Arial")
plot(plotL)
dev.off()

##############################################################

#average degree comparison
originalRes <- read.table(paste0(patH, "degreeS.csv"), header=T, sep="\t", stringsAsFactors=T)

#test original degrees with n=plotted random degrees selected from the 100 generated!
randomChos <- sample(seq(1, 100, by=1), 10, replace=F)

#create label vector
subSet <- list(list())
for (rC in seq(1, plotted, by=1)) {
	subSet[[rC]] <- allDegsList[[randomChos[rC]]]
}

#and compute some statistics between 10 randomly chosen random dataset and the original dataset
pvals <- vector()
for (i in seq(1, plotted, by=1)) {
	pvals <- append(pvals, t.test(unlist(subSet[[i]][, -1]), unlist(originalRes[, -1]))$p.value)
}

#test 10 randomly chosen datasets degrees against the same 10 randomly datasets degrees. the result is a matrix whose diagonal will be equal to 1
pvalsRandom <- data.frame(matrix(nrow=plotted, ncol=plotted))
for (i in seq(1, plotted, by=1)) {
	for (j in seq(1, plotted, by=1)) {
		#if i==j then the result should be 1
		pvalsRandom[i ,j] <- t.test(unlist(subSet[[i]][, -1]), unlist(subSet[[j]][, -1]))$p.value
	}
}

#merge everything into a data structure to create a heatmap

#biological comparison
realPVs <- cbind.data.frame(rowS=rep("master", plotted), colS=paste0("random ", randomChos), Similarity=pvals, stringsAsFactors=F)
realPVs$rowS <- factor(realPVs$rowS, levels=unique(realPVs$rowS))
realPVs$colS <- factor(realPVs$colS, levels=unique(realPVs$colS))

#now create another label vector for the randomly chosen datasets
rowS <- vector()
#for each column of comparisonMat
for (r in seq(1, plotted, by=1)) {
	#create a vector of the first rowname, then the second, the third, etc etc which is long as the number of columns
	rowS <- append(rowS, rep(paste0("random ", randomChos[r]), plotted))
}

#random comparison
randomPVs <- cbind.data.frame(rowS=rowS, colS=rep(paste0("random ", randomChos), plotted), Similarity=unlist(pvalsRandom), stringsAsFactors=F)
randomPVs$rowS <- factor(randomPVs$rowS, levels=unique(randomPVs$rowS))
randomPVs$colS <- factor(randomPVs$colS, levels=unique(randomPVs$colS))

#put random and biological datasets together
chosenFrame <- rbind.data.frame(randomPVs, realPVs)

#and finally create the plot
heaT <- ggplot(chosenFrame, aes(rowS, colS, fill=Similarity)) + 
	geom_tile() +
	scale_fill_distiller(palette="YlGn", direction=2) +
	geom_text(aes(label=formatC(Similarity, format="e", digits=2)), size=5) +
	xlab("") +
	ylab("") +
	theme_bw() +
	theme(plot.title=element_text(size=sizeText, face="bold", hjust=0.5),
		axis.text.x=element_text(angle=45, hjust=1, size=sizeText, face="bold"),
		axis.text.y=element_text(size=sizeText, face="bold"),
		axis.text=element_text(size=sizeText, face="bold"),
		axis.title=element_text(size=sizeText, face="bold"),
		legend.position="right",
		legend.title=element_text(size=sizeText),
		legend.text=element_text(size=sizeText))
	
#and export it
pdf(paste0(saveIn, "heatmap_Similarity.pdf"), width=15, height=12)
plot(heaT)
dev.off()

