#############
# GabrieleT #
#   2021    #
#############

# load libraries
library("igraph")
library("gtools")
library("ggplot2")
library("ggridges")
library("RColorBrewer")
library("gridExtra")
library("extrafont")
loadfonts()

# set paths
patH <- "~/WNNets/"
saveIn <- "~/WNNets/figures/"

# set palette with 8 colours for colouring plots
colorsPal <- brewer.pal(8, "Set2")

# load data
filteredData <- read.csv(paste0(patH, "node_table.csv"), sep="\t", header=T, stringsAsFactors=F)

# read network
stringInfo <- read.csv(paste0(patH, "string_interactions_over400_experimental.tsv"), sep="\t", header=T, stringsAsFactors=T)

# create edge list
edgeL <- as.matrix(cbind.data.frame(as.character(stringInfo$X.node1), as.character(stringInfo$node2)))
# and build the network
networkTmp <- graph_from_edgelist(edgeL)
# get biggest component
network <- decompose(networkTmp)[[1]]

# initialise the storage
allDegsVect <- vector()
allDegsList <- list(list())

# vectors of pvalues
pvalNH <- vector()
pvalNU <- vector()
infoPNH <- vector()
infoPNU <- vector()

# generating 100 trials but using only 10 to plot the heatmap
trials <- 100
plotted <- 10

# set seed
set.seed(131)

# test n=trials different random datasets
for (r in seq(1, trials, by=1)) {

	# compute random values for the random weighting matrix
	randomData <- cbind.data.frame(filteredData$proteinName, matrix(runif(dim(filteredData[, -1])[1] * dim(filteredData[, -1])[2], min(filteredData[, -1]), max(filteredData[, -1])), ncol=dim(filteredData[, -1])[2], nrow=dim(filteredData[, -1])[1],), stringsAsFactors=F)
	colnames(randomData) <- colnames(filteredData)

	print(paste0("iteration ", r))
	degreesTmp <- vector()
	# computing degree for each node in the network, using experimental values
	for (vrtx in V(network)$name) {
		# get the neighbours of vertex
		neighList <- neighbors(network, vrtx, mode="all")
		# initialise array for the Degree of the vertex, for each columns, i.e. each sample
		sumUp <- rep(0, length(randomData[, -1]))
		# for each neighbour of the selected vertex compute degree
		for (ngh in neighList$name) {
			# get the array of experimental values of that neighbour for all the samples
			# then sum to the previous existing values
			sumUp <- sumUp + randomData[which(randomData$proteinName==ngh), -1]
		}
		# attach the new values in the degree table: another vertex's degree is computed
		degreesTmp <- rbind.data.frame(degreesTmp, sumUp)
	}
	# the degree is a small, floating point value since all the original values were normalised by the molecular weight!
	degreesAll <- cbind.data.frame(V(network)$name, degreesTmp)
	colnames(degreesAll) <- colnames(randomData)
	# create net label
	nlabl <- rep(r, nrow(degreesTmp))
	# unlist degrees and merge
	allDegsVect <- rbind.data.frame(allDegsVect, cbind.data.frame(protein=degreesAll[, 1], rNet=nlabl, value=unlist(degreesAll[, -1]), stringsAsFactors=F))
	allDegsList[[r]] <- degreesAll

	# get only degree values
	degOnly <- degreesAll[, -1]
	# and protein names
	proteins <- as.factor(degreesAll$proteinName)

	# NOTE:
	# SUBSETTING BY SAMPLE GROUPS

	# N7-N12 VERSUS H1-H6
	subSampleForH <- cbind.data.frame(proteins, degOnly[, colnames(degOnly[, grep("N", colnames(degOnly))])[c(7:12)]], degOnly[, grep("H", colnames(degOnly))])
	# colnames(subSampleForH)
	
	# N1-N6 VERSUS U1-U6
	subSampleForU <- cbind.data.frame(proteins, degOnly[, colnames(degOnly[, grep("N", colnames(degOnly))])[c(1:6)]], degOnly[, grep("U", colnames(degOnly))])
	# colnames(subSampleForU)

	# computing pvalues
	for (prt in proteins) {
		subsampleH <- as.matrix(subSampleForH[which(subSampleForH$protein==prt), -1])
		# get the N values from sumsampleForH
		valN <- subsampleH[, grep("N", colnames(subsampleH))]
		valH <- subsampleH[, grep("H", colnames(subsampleH))]
		# compute t.test
		pv <- t.test(valN, valH)$p.value
		infoPNH <- rbind.data.frame(infoPNH, cbind.data.frame(protein=prt, pval=pv, stringsAsFactors=F))
		
		subsampleU <- as.matrix(subSampleForU[which(subSampleForU$protein==prt), -1])
		# get the N values from sumsampleForU
		valN <- subsampleU[, grep("N", colnames(subsampleU))]
		valU <- subsampleU[, grep("U", colnames(subsampleU))]
		
		pv <- t.test(valN, valU)$p.value
		infoPNU <- rbind.data.frame(infoPNU, cbind.data.frame(protein=prt, pval=pv, stringsAsFactors=F))
	}
}

# RIDGE PLOTTING
# set FDR as method for correction
mth <- "fdr"
# perform correction using package stats
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust
infoPNH <- data.frame(infoPNH, padj=p.adjust(infoPNH$pval, method=mth))
infoPNU <- data.frame(infoPNU, padj=p.adjust(infoPNU$pval, method=mth))

# and plot p-adjusted
plotL <- grid.arrange(
	ggplot(infoPNH, aes(x = padj, y = protein, fill = protein)) + geom_density_ridges(from=0, to=0.5) + theme_ridges() + theme(legend.position = "none") + labs(x = "p-adjusted"),
	ggplot(infoPNH, aes(x = padj, y = protein, fill = protein)) + geom_density_ridges(from=0.5, to=1) + theme_ridges() + theme(legend.position = "none") + labs(x = "p-adjusted"),
	ncol=2)
		
# and export it
pdf(paste0(saveIn, "padj_ridge_H.pdf"), width=9, height=15, font="Arial")
plot(plotL)
dev.off()

# and plot p-adjusted
plotL <- grid.arrange(
	ggplot(infoPNU, aes(x = padj, y = protein, fill = protein)) + geom_density_ridges(from=0.5, to=0.82) + theme_ridges() + theme(legend.position = "none") + labs(x = "p-adjusted"),
	ggplot(infoPNU, aes(x = padj, y = protein, fill = protein)) + geom_density_ridges(from=0.95, to=1) + theme_ridges() + theme(legend.position = "none") + labs(x = "p-adjusted"),
	ncol=2)

# and export it
pdf(paste0(saveIn, "padj_ridge_U.pdf"), width=9, height=15, font="Arial")
plot(plotL)
dev.off()
