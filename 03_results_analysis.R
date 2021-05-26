#############
# GabrieleT #
#   2021    #
#############

# load libraries
library("ggplot2")
library("reshape2")
library("UpSetR")
library("grid")
library("ggrepel")
library("RColorBrewer")
library("igraph")
library("extrafont")
library("stats")
loadfonts()

# set seed
set.seed(131)

# set paths
patH <- "~/WNNets/"
saveIn <- "~/WNNets/figures/"

# read data from R
originalRes <- read.table(paste0(patH, "degreeS.csv"), header=T, sep="\t", stringsAsFactors=T)

# load thresholds file
fcThresh <- read.csv("FC_FDR_threshold.csv", header=T, sep="\t")

# set palette with 8 colours for colouring plots
coloursPal <- brewer.pal(8, "Set2")

# get only degree values
degOnly <- originalRes[, -1]
# and protein names
proteins <- as.factor(originalRes$proteinName)

# SUBSETTING BY SAMPLE GROUPS

# N1-N6 VERSUS U1-U6
subSampleForU <- cbind.data.frame(proteins, degOnly[, colnames(degOnly[, grep("N", colnames(degOnly))])[c(1:6)]], degOnly[, grep("U", colnames(degOnly))])
# get degree only
degreesU <- t(subSampleForU[, -1])
# and samples name
subjectsU <- as.factor(rownames(degreesU))
rownames(degreesU) <- subjectsU
colnames(degreesU) <- proteins

# N7-N12 VERSUS H1-H6
subSampleForH <- cbind.data.frame(proteins, degOnly[, colnames(degOnly[, grep("N", colnames(degOnly))])[c(7:12)]], degOnly[, grep("H", colnames(degOnly))])
# get degree only
degreesH <- t(subSampleForH[, -1])
# and samples name
subjectsH <- as.factor(rownames(degreesH))
rownames(degreesH) <- subjectsH
colnames(degreesH) <- proteins

allDegU <- melt(degreesU)
allDegH <- melt(degreesH)
colnames(allDegU) <- colnames(allDegH) <- c("subject", "proteins", "degree")

# these were values but now they become vectors of values
pvalNH <- vector()
pvalNU <- vector()
proteinNamesH <- vector()
proteinNamesU <- vector()

# computing statistics for H
for (prt in unique(allDegH$proteins)) {
	# append the protein name
	proteinNamesH <- append(proteinNamesH, prt)
	# subsample on the degree 
	subsample <- as.matrix(subSampleForH[which(subSampleForH$proteins==prt), -1])
	# get the H values from subsample
	valN <- subsample[, grep("N", colnames(subsample))]
	# and for N
	valH <- subsample[, grep("H", colnames(subsample))]
	# perform the t-test using the two samples
	pvalNH <- append(pvalNH, t.test(valN, valH)$p.value)
}

# computing statistics for U
for (prt in unique(allDegU$proteins)) {
	# append the protein name
	proteinNamesU <- append(proteinNamesU, prt)
	# subsample on the degree 
	subsample <- as.matrix(subSampleForU[which(subSampleForU$proteins==prt), -1])
	# get the H values from subsample
	valN <- subsample[, grep("N", colnames(subsample))]
	# and for N
	valU <- subsample[, grep("U", colnames(subsample))]
	# perform the t-test using the two samples
	pvalNU <- append(pvalNU, t.test(valN, valU)$p.value)
}

# verify that both tables contain the same set of proteins and generate a vector of names
ifelse(proteinNamesH==proteinNamesU, proteinNames<-proteinNamesH, print("ERRORE GRAVE ATTENZIONE!!!!!!!!"))

# compute adjusted pvalues using false discovery rate
# set threshold for FDR
tH <- fcThresh$FDR
# set FDR as method for corretcion
mth <- "fdr"
# perform correction using package stats
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust
padj_H <- p.adjust(pvalNH, method=mth)
padj_U <- p.adjust(pvalNU, method=mth)

print(paste0("p-adjusted < ", tH, " significant proteins, for NH-H: ", paste(proteins[padj_H < tH], collapse=", ")))
print(paste0("p-adjusted < ", tH, " significant proteins, for NU-U: ", paste(proteins[padj_U < tH], collapse=", ")))

# bubble plotting using adjusted p.values

# set some parameters
sizeText <- 20
bubblePal <- brewer.pal(4, "BuPu")

# merge everything!
tableAdj <- cbind.data.frame(proteinNames, padj_H, padj_U)
colnames(tableAdj) <- c("proteins", "padjNH", "padjNU")

# find significant proteins, i.e. SSPs
signifU <- as.character(tableAdj$proteins[tableAdj$padjNU < tH])
signifH <- as.character(tableAdj$proteins[tableAdj$padjNH < tH])

# export SSPs lists
write.table(signifU, paste0(patH, "U_SSPs.csv"), quote=F, row.names=F, col.names=F)
write.table(signifH, paste0(patH, "H_SSPs.csv"), quote=F, row.names=F, col.names=F)

# prepare for plotting
bubblePlot <- ggplot(tableAdj, aes(x=padjNU, y=padjNH)) +
	geom_point(alpha=0.8, shape=21, color="black", fill=coloursPal[3], size=5) +
	xlim(0, 1) +
	ylim(0, 1) +
	expand_limits(x=0) +
	labs(x="adjusted p-value NU-U", y="adjusted p-value NH-H") +
	ggtitle("") +
	geom_hline(yintercept=tH, color=coloursPal[2], lwd=0.8, linetype="dashed") +
	geom_vline(xintercept=tH, color=coloursPal[2], lwd=0.8, linetype="dashed") + 
	theme_bw() +
	theme(plot.title=element_text(size=sizeText, face="bold", hjust=0.5),
		axis.text.x=element_text(hjust=1, size=sizeText, face="bold"),
		axis.text.y=element_text(size=sizeText, face="bold"),
		axis.text=element_text(size=sizeText, face="bold"),
		axis.title=element_text(size=sizeText, face="bold"),
		legend.position="top",
		legend.title=element_text(size=sizeText),
		legend.text=element_text(size=sizeText)) +
	geom_text_repel(label=tableAdj$proteins, fontface="italic", size=4, colour="black", max.overlaps=50)

# and plot
print("bubblePlot adjusted")
pdf(paste(saveIn, "bubblePlot_adjusted.pdf", sep=""), width=20, height=15, font="Arial")
plot(bubblePlot)
dev.off()

# PLOTTING BEST PROTEIN LINEPLOT

# COMPUTING THE AVERAGE DIFFERENCE BETWEEN DEGREES IN THE CONDITIONS

# for H
avgDegH <- vector()
avgDegNH <- vector()

# loop through significant H, using the U proteins as a reference to compare the same proteins that were found for the NU-U comparison
for (prt in signifH) {
	# get the degree for a specific protein, using U proteins
	subsampleH <- allDegH[which(allDegH$proteins==prt), ]
	# compute the mean for H samples
	avgDegNH <- append(avgDegNH, mean(subsampleH[grep("N", subsampleH$subject), ]$degree))
	# compute the mean for NH samples
	avgDegH <- append(avgDegH, mean(subsampleH[grep("H", subsampleH$subject), ]$degree))
}

# prepare to plot
signifHPos <- vector()

# comparing proteins
for (prt in signifH) {
	pTmp <- grep(prt, allDegH$proteins)
	# store positions
	signifHPos <- append(signifHPos, pTmp)
}

# subset allDeg using the positions of the significant proteins
ssADegH <- allDegH[signifHPos, ]
# substitue N with NH, for plotting purposes
ssADegH$subject <- gsub("N", "NH", ssADegH$subject)
# removing unwanted characters and preparing new columns for the plot
ssADegH$group <- as.factor(gsub("[0-9]", "", ssADegH$subject))
# transforming the group in number
ssADegH$subgroup <- as.numeric(as.factor(gsub("[a-z]", "", ssADegH$subject)))
# finally, get only N or H
ssADegH <- ssADegH[grep("N|H", ssADegH$subject), ]

# computing the differences in terms of degree for the proteins that are considered significant
# first we compute the distance between each degree, then the average distance is rounded to 4 digits, to be plotted.
diffH <- round(mean(abs(avgDegNH - avgDegH)), digits=6)

# finally, plot
lPlot <- ggplot(data=ssADegH, aes(x=proteins, y=log(degree), group=subject, color=group, alpha=subgroup)) + 
	geom_line(lwd=0.8) +
	geom_point()+
	scale_alpha_continuous(range=c(1, 1)) +
	ylim(-15, -2) +
	scale_color_manual(values=c(coloursPal[1], coloursPal[2])) +
	guides(alpha=FALSE) +
	labs(x="", y="log(degree)") +
	ggtitle("") +
	theme_bw() +
	theme(plot.title=element_text(size=sizeText, face="bold", hjust=0.5),
	axis.text.x=element_text(angle=45, hjust=1, size=sizeText, face="bold"),
	axis.text.y=element_text(size=sizeText, face="bold"),
	axis.text=element_text(size=sizeText, face="bold"),
	axis.title=element_text(size=sizeText, face="bold"),
	legend.position="top",
	legend.title=element_text(size=sizeText),
	legend.text=element_text(size=sizeText)) +
	annotate(geom="text", x=1, y=-12, label=paste0("Average difference = ", diffH), hjust="left", size=6)

print("linePlotNH-ing")
pdf(paste0(saveIn, "linePlotNH_adjusted.pdf"), width=12, height=6, font="Arial")
plot(lPlot)
dev.off()

# REPEAT FOR U

# for U
avgDegU <- vector()
avgDegNU <- vector()

# loop through significant U
for (prt in signifU) {
	# get the degree for a specific protein
	subsampleU <- allDegU[which(allDegU$proteins==prt), ]
	# compute the mean for U samples
	avgDegU <- append(avgDegU, mean(subsampleU[grep("U", subsampleU$subject), ]$degree))
	# compute the mean for NU samples
	avgDegNU <- append(avgDegNU, mean(subsampleU[grep("N", subsampleU$subject), ]$degree))
}

# prepare to plot
signifUPos <- vector()

# comparing proteins
for (prt in signifU) {
	pTmp <- grep(prt, allDegU$proteins)
	signifUPos <- append(signifUPos, pTmp)
}

# subset allDeg
ssADegU <- allDegU[signifUPos, ]
# substitue N with NU, for plotting purposes
ssADegU$subject <- gsub("N", "NU", ssADegU$subject)
# removing unwanted characters and preparing new columns for the plot
ssADegU$group <- as.factor(gsub("[0-9]", "", ssADegU$subject))
# transforming the group in number
ssADegU$subgroup <- as.numeric(as.factor(gsub("[a-z]", "", ssADegU$subject)))
# finally, get only N or H
ssADegU <- ssADegU[grep("N|U", ssADegU$subject), ]

# computing the differences in terms of degree for the proteins that are considered significant
# first we compute the distance between each degree, then the average distance is rounded to 4 digits, to be plotted.
diffU <- round(mean(abs(avgDegNU - avgDegU)), digits=6)
# we expect that the average distance for the U samples is greater than the average distance for the H sample
# diffU > diffH

# finally, plot
lPlot <- ggplot(data=ssADegU, aes(x=proteins, y=log(degree), group=subject, color=group, alpha=subgroup)) + 
	geom_line(lwd=0.8) +
	geom_point()+
	scale_alpha_continuous(range=c(1, 1)) +
	ylim(-15, -2) +
	# invert colours since it plots by alphabetical order!
	scale_color_manual(values=c(coloursPal[2], coloursPal[1])) +
	guides(alpha=FALSE) +
	labs(x="", y="log(degree)") +
	ggtitle("") +
	theme_bw() +
	theme(plot.title=element_text(size=sizeText, face="bold", hjust=0.5),
	axis.text.x=element_text(angle=45, hjust=1, size=sizeText, face="bold"),
	axis.text.y=element_text(size=sizeText, face="bold"),
	axis.text=element_text(size=sizeText, face="bold"),
	axis.title=element_text(size=sizeText, face="bold"),
	legend.position="top",
	legend.title=element_text(size=sizeText),
	legend.text=element_text(size=sizeText)) +
	annotate(geom="text", x=20, y=-12, label=paste0("Average difference = ", diffU), hjust="left", size=6)

print("linePlotNU-ing")
pdf(paste0(saveIn, "linePlotNU_adjusted.pdf"), width=12, height=6, font="Arial")
plot(lPlot)
dev.off()

# COMPUTING BEST, AVERAGE DEGREE

# compute the average degree for N and U and put everything together
averageDegU <- cbind.data.frame(proteinName=subSampleForU$proteins, meanN=rowMeans(subSampleForU[, grep("N", colnames(subSampleForU))]), meanU=rowMeans(subSampleForU[, grep("U", colnames(subSampleForU))]))
# compute the average degree for N and H and put everything together
averageDegH <- cbind.data.frame(proteinName=subSampleForH$proteins, meanN=rowMeans(subSampleForH[, grep("N", colnames(subSampleForH))]), meanH=rowMeans(subSampleForH[, grep("H", colnames(subSampleForH))]))

# UPSET PLOTTING

# set palette with 8 colours for colouring plots
coloursPal <- brewer.pal(12, "Paired")

# get relevant fold change proteins for H
fcNH <- read.csv(paste0(patH, "NH_SFCPs.csv"), sep="\t", header=T, stringsAsFactors=T)
# get relevant fold change proteins for U
fcNU <- read.csv(paste0(patH, "NU_SFCPs.csv"), sep="\t", header=T, stringsAsFactors=T)

# UPSET for NH-H
listInputAll_H <-list()

# create the list with controls for H
listInputAll_H$ControlPs_NH <- as.character(averageDegH$proteinName[which(averageDegH$meanN > quantile(averageDegH$meanN)[4])])
# treated for H
listInputAll_H$TreatedPs_H <- as.character(averageDegH$proteinName[which(averageDegH$meanH > quantile(averageDegH$meanH)[4])])

# SSPs
listInputAll_H$SSPs_H <- signifH
# HFCPs
listInputAll_H$SFCPs_H <- as.character(fcNH$x)

# rename sublists to plot
names(listInputAll_H) <- c("NH-ControlPs", "H-TreatedPs", "H-SSPs", "H-SFCPs")

print("upsetALL NH-H")
pdf(paste0(saveIn, "upsetALL_NH_H.pdf"), width=8, height=6, font="Arial")
upset(fromList(listInputAll_H), order.by="freq", nsets=4, decreasing=T, sets=c("NH-ControlPs", "H-TreatedPs", "H-SSPs", "H-SFCPs"), keep.order=T, text.scale=1.5, sets.x.label="Set size", mainbar.y.label="Shared proteins", mainbar.y.max=25, set_size.scale_max=35,
	queries=list(
	list(query=intersects, params=list("NH-ControlPs", "H-TreatedPs"), color=coloursPal[4], active=T),
	list(query=intersects, params=list("NH-ControlPs", "H-TreatedPs", "H-SFCPs"), color=coloursPal[8], active=T),
	list(query=intersects, params=list("H-SSPs"), color=coloursPal[10], active=T)
	))
dev.off()

# UPSET for NU-U
listInputAll_U <-list()

# create the list with controls for U
listInputAll_U$ControlPs_NU <- as.character(averageDegU$proteinName[which(averageDegU$meanN > quantile(averageDegU$meanN)[4])])
# treated for U
listInputAll_U$TreatedPs_U <- as.character(averageDegU$proteinName[which(averageDegU$meanU > quantile(averageDegU$meanU)[4])])

# SSPs
listInputAll_U$SSPs_U <- signifU
# HFCPs
listInputAll_U$SFCPs_U <- as.character(fcNU$x)

# rename sublists to plot
names(listInputAll_U) <- c("NU-ControlPs", "U-TreatedPs", "U-SSPs", "U-SFCPs")

print("upsetALL NU-U")
pdf(paste0(saveIn, "upsetALL_NU_U.pdf"), width=8, height=6, font="Arial")
upset(fromList(listInputAll_U), order.by="freq", nsets=4, decreasing=T, sets=c("NU-ControlPs", "U-TreatedPs", "U-SSPs", "U-SFCPs"), keep.order=T, text.scale=1.5, sets.x.label="Set size", mainbar.y.label="Shared proteins", mainbar.y.max=25, set_size.scale_max=35,
	queries=list(
	list(query=intersects, params=list("NU-ControlPs", "U-TreatedPs"), color=coloursPal[4], active=T),
	list(query=intersects, params=list("U-SSPs", "U-SFCPs"), color=coloursPal[5], active=T),
	list(query=intersects, params=list("U-SSPs", "U-SFCPs", "NU-ControlPs"), color=coloursPal[1], active=T),	
	list(query=intersects, params=list("U-TreatedPs", "U-SSPs"), color=coloursPal[2], active=T),
	list(query=intersects, params=list("U-SSPs", "U-TreatedPs", "U-SFCPs"), color=coloursPal[12], active=T),
	list(query=intersects, params=list("NU-ControlPs", "U-SSPs"), color=coloursPal[6], active=T),
	list(query=intersects, params=list("NU-ControlPs", "U-TreatedPs", "U-SFCPs"), color=coloursPal[10], active=T),
	list(query=intersects, params=list("NU-ControlPs", "U-TreatedPs", "U-SSPs"), color=coloursPal[8], active=T)
	))
dev.off()

# UPSET WITH EVERYTHING
listInputAll <-list()

# create the list with controls for H
listInputAll$ControlPs_NH <- as.character(averageDegH$proteinName[which(averageDegH$meanN > quantile(averageDegH$meanN)[4])])
# treated for H
listInputAll$TreatedPs_H <- as.character(averageDegH$proteinName[which(averageDegH$meanH > quantile(averageDegH$meanH)[4])])
# controls for U
listInputAll$ControlPs_NU <- as.character(averageDegU$proteinName[which(averageDegU$meanN > quantile(averageDegU$meanN)[4])])
# treated for U
listInputAll$TreatedPs_U <- as.character(averageDegU$proteinName[which(averageDegU$meanU > quantile(averageDegU$meanU)[4])])

# SSPs
listInputAll$SSPs_U <- signifU
listInputAll$SSPs_H <- signifH
# HFCPs
listInputAll$SFCPs_H <- as.character(fcNH$x)
listInputAll$SFCPs_U <- as.character(fcNU$x)

# rename sublists to plot
names(listInputAll) <- c("NH-ControlPs", "H-TreatedPs", "NU-ControlPs", "U-TreatedPs", "U-SSPs", "H-SSPs", "H-SFCPs", "U-SFCPs")

print("upset full version")
pdf(paste0(saveIn, "upsetALL_fullversion.pdf"), width=8, height=6, font="Arial")
upset(fromList(listInputAll), order.by="freq", nsets=4, decreasing=T, sets=c("NH-ControlPs", "H-TreatedPs", "NU-ControlPs", "U-TreatedPs", "H-SSPs", "H-SFCPs", "U-SSPs", "U-SFCPs"), keep.order=T, text.scale=1.5, sets.x.label="Set size", mainbar.y.label="Shared proteins", mainbar.y.max=20, set_size.scale_max=35,
	queries=list(
	list(query=intersects, params=list("NH-ControlPs", "H-TreatedPs", "NU-ControlPs", "U-TreatedPs"), color=coloursPal[12], active=T),
	list(query=intersects, params=list("NH-ControlPs", "H-TreatedPs", "NU-ControlPs", "U-TreatedPs", "U-SSPs"), color=coloursPal[8], active=T),
	list(query=intersects, params=list("NH-ControlPs", "H-TreatedPs", "NU-ControlPs", "U-TreatedPs", "U-SFCPs"), color=coloursPal[5], active=T),
	list(query=intersects, params=list("NH-ControlPs", "H-TreatedPs", "NU-ControlPs", "U-TreatedPs", "U-SFCPs", "H-SFCPs"), color=coloursPal[4], active=T),
	list(query=intersects, params=list("U-SFCPs", "U-SSPs"), color=coloursPal[10], active=T),
	list(query=intersects, params=list("U-SFCPs"), color=coloursPal[2], active=T),
	list(query=intersects, params=list("U-SSPs"), color=coloursPal[6], active=T)
	))
dev.off()

# ############## PLOT NETWORK WITH MOST SIGNIFICANT PROTEINS AND THEIR FIRST NEIGHBOURS

# set palette with 8 colours for colouring plots
coloursPal <- brewer.pal(8, "Set2")

stringInfo <- read.csv(paste0(patH, "string_interactions_over400_experimental.tsv"), sep="\t", header=T, stringsAsFactors=T)

# create edge list
edgeL <- as.matrix(cbind.data.frame(as.character(stringInfo$X.node1), as.character(stringInfo$node2)))
# and build the network
networkTmp <- graph_from_edgelist(edgeL)
# get biggest component
network <- decompose(networkTmp)[[1]]

# for U samples

neighListU <- signifU
# get first neighbours of signifU proteins
for (vrtx in signifU) {
	# get the neighbours of vertex
	neighListU <- append(neighListU, neighbors(network, vrtx, mode="all")$name)
}

# export the neighbours of SSPs
neighsOfSignifU <- unique(neighListU)[!unique(neighListU) %in% signifU]
write.table(neighsOfSignifU, paste0(patH, "neighsOf_U_SFCPs.csv"), quote=F, row.names=F, col.names=F)

# create the network for SSPs and their neighbours
firstNeighsU <- induced.subgraph(network, unique(neighListU))

for (v in c(1:length(V(firstNeighsU)))) {
	# in pink both SSPs and SFCPs	
	if (V(firstNeighsU)$name[v] %in% signifU && V(firstNeighsU)$name[v] %in% fcNU$x) {
		V(firstNeighsU)$colour[v] <- coloursPal[4]
	# in green SSPs only
	} else if (V(firstNeighsU)$name[v] %in% signifU) {
		V(firstNeighsU)$colour[v] <- coloursPal[5]
	# in yellow SFCPs only
	} else if (V(firstNeighsU)$name[v] %in% fcNU$x) {
		V(firstNeighsU)$colour[v] <- coloursPal[6]
	} else if (V(firstNeighsU)$name[v] %in% c("CFL2", "MT-ATP6", "ACTN1")) {
		V(firstNeighsU)$colour[v] <- coloursPal[3]
	} else {
		V(firstNeighsU)$colour[v] <- "#ffffff"
	}
}

# U NET with vertexes labels
set.seed(108)
lytNicely <- layout_nicely(firstNeighsU)
print("signifU_orig")
pdf(paste0(saveIn, "signifU_orig.pdf"), font="Arial")
plot(firstNeighsU, vertex.label=V(firstNeighsU)$name, vertex.label.cex=0.6, vertex.size=16, edge.arrow.size=0, edge.color="dimgrey", edge.size=4, vertex.color=V(firstNeighsU)$colour, vertex.label.color="black", layout=lytNicely)
dev.off()

######### and for H samples

neighListH <- signifH
# get first neighbours of signifU proteins
for (vrtx in signifH) {
	# get the neighbours of vertex
	neighListH <- append(neighListH, neighbors(network, vrtx, mode="all")$name)
}

# export the neighbours of SSPs
neighsOfSignifH <- unique(neighListH)[!unique(neighListH) %in% signifH]
write.table(neighsOfSignifH, paste0(patH, "neighsOf_H_SFCPs.csv"), quote=F, row.names=F, col.names=F)

# create the network for SSPs and their neighbours
firstNeighsH <- induced.subgraph(network, unique(neighListH))

for (v in 1:length(V(firstNeighsH))) {
	# in pink both SSPs and SFCPs
	if (V(firstNeighsH)$name[v] %in% signifH && V(firstNeighsH)$name[v] %in% fcNH$x) {
		V(firstNeighsH)$colour[v] <- coloursPal[4]
	# in green SSPs only
	} else if (V(firstNeighsH)$name[v] %in% signifH) {
		V(firstNeighsH)$colour[v] <- coloursPal[5]
	# in yellow SFCPs only
	} else if (V(firstNeighsH)$name[v] %in% fcNH$x) {
		V(firstNeighsH)$colour[v] <- coloursPal[6]
	} else if (V(firstNeighsH)$name[v] %in% c("CLTC", "COX5A", "HSPA1A", "HP")) {
		V(firstNeighsH)$colour[v] <- coloursPal[3]
	} else {
		V(firstNeighsH)$colour[v] <- "#ffffff"
	}
}

# H NET with vertexes labels
set.seed(131)
lytNicely <- layout_nicely(firstNeighsH)
print("signifH_orig")
pdf(paste0(saveIn, "signifH_orig.pdf"), font="Arial")
plot(firstNeighsH, vertex.label=V(firstNeighsH)$name, vertex.label.cex=0.6, vertex.size=16, edge.arrow.size=0, edge.color="dimgrey", edge.size=4, vertex.color=V(firstNeighsH)$colour, vertex.label.color="black", layout=lytNicely)
dev.off()

