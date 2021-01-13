#############
# GabrieleT #
#   2020    #
#############

#load libraries
library("ggplot2")
library("reshape")
library("UpSetR")
library("grid")
library("ggrepel")
library("RColorBrewer")
library("igraph")
library("extrafont")
library("stats")
loadfonts()

#set seed
set.seed(131)

#set paths
patH <- "/home/cbmc/wnnets/"
saveIn <- "/home/cbmc/wnnets/figures/"

#read data from R
originalRes <- read.table(paste0(patH, "degreeS.csv"), header=T, sep="\t", stringsAsFactors=T)

#set palette with 8 colours for colouring plots
coloursPal <- brewer.pal(8, "Set2")

#get only degree values
degOnly <- originalRes[, -1]
#and protein names
proteins <- as.factor(originalRes$proteinName)

# NOTE:
# SUBSETTING CONSIDERING THAT THE COMPARISON MUST FOLLOW THESE RULES:

#N7-N12 VERSUS H1-H6
subSampleForH <- cbind.data.frame(proteins, degOnly[, colnames(degOnly[, grep("N", colnames(degOnly))])[c(7:12)]], degOnly[, grep("H", colnames(degOnly))])
#get degree only
degreesH <- t(subSampleForH[, -1])
#and samples name
subjectsH <- as.factor(rownames(degreesH))
rownames(degreesH) <- subjectsH
colnames(degreesH) <- proteins

#N1-N6 VERSUS U1-U6
subSampleForU <- cbind.data.frame(proteins, degOnly[, colnames(degOnly[, grep("N", colnames(degOnly))])[c(1:6)]], degOnly[, grep("U", colnames(degOnly))])
#get degree only
degreesU <- t(subSampleForU[, -1])
#and samples name
subjectsU <- as.factor(rownames(degreesU))
rownames(degreesU) <- subjectsU
colnames(degreesU) <- proteins

allDegU <- melt(degreesU)
allDegH <- melt(degreesH)
colnames(allDegU) <- colnames(allDegH) <- c("subject", "proteins", "degree")

#these were values but now they become vectors of values
pvalNH <- vector()
pvalNU <- vector()
proteinNamesH <- vector()
proteinNamesU <- vector()

#computing statistics for H
for (prt in unique(allDegH$proteins)) {
	proteinNamesH <- append(proteinNamesH, prt)
	subsample <- as.matrix(subSampleForH[which(subSampleForH$proteins==prt), -1])
	#get the N values from sumsampleForH
	valN <- subsample[, grep("N", colnames(subsample))]
	valH <- subsample[, grep("H", colnames(subsample))]
	
	pvalNH <- append(pvalNH, t.test(valN, valH)$p.value)
}

#computing statistics for U
for (prt in unique(allDegU$proteins)) {
	proteinNamesU <- append(proteinNamesU, prt)
	subsample <- as.matrix(subSampleForU[which(subSampleForU$proteins==prt), -1])
	#get the N values from sumsampleForU
	valN <- subsample[, grep("N", colnames(subsample))]
	valU <- subsample[, grep("U", colnames(subsample))]

	pvalNU <- append(pvalNU, t.test(valN, valU)$p.value)
}

#merge information for NH
tableNH <- cbind.data.frame(proteinNamesH, pvalNH)
#merge information for NU
tableNU <- cbind.data.frame(proteinNamesU, pvalNU)
colnames(tableNH) <- colnames(tableNU) <- c("proteins", "pvalue")

#verify that both tables contain the same set of proteins and generate a vector of names
ifelse(proteinNamesH==proteinNamesU, proteinNames<-proteinNamesH, print("ERRORE GRAVE ATTENZIONE!!!!!!!!"))

#merge everything!
tableBoth <- cbind.data.frame(proteinNames, pvalNH, pvalNU)
colnames(tableBoth) <- c("proteins", "pvalNH", "pvalNU")

#show significant pvalues
#tableBoth[tableBoth$pvalNU < 0.05, ]
#tableBoth[tableBoth$pvalNH < 0.05, ]

#get names of significant proteins, i.e. those proteins with pval lower than 0.05
signifU <- as.character(tableBoth$proteins[tableBoth$pvalNU < 0.05])
signifH <- as.character(tableBoth$proteins[tableBoth$pvalNH < 0.05])
#and the unique set
signifProt <- unique(c(signifU, signifH))

#and their pvalues
signifUVals <- tableBoth$pvalNU[tableBoth$pvalNU < 0.05]
signifHVals <- tableBoth$pvalNH[tableBoth$pvalNH < 0.05]

#print significant proteins names
write.table(signifU, paste0(patH, "signifU.csv"), quote=F, row.names=F, col.names=F)
write.table(signifH, paste0(patH, "signifH.csv"), quote=F, row.names=F, col.names=F)

#compute adjusted pvalues using false discovery rate
mth <- "fdr"
padj_H <- p.adjust(pvalNH, method=mth, n=length(pvalNH))
padj_U <- p.adjust(pvalNU, method=mth, n=length(pvalNU))

print(paste0("p-adjusted <= 0.1 significant proteins, for NH-H: ", paste(proteins[padj_H<=0.1], collapse=", ")))
print(paste0("p-adjusted <= 0.1 significant proteins, for NU-U: ", paste(proteins[padj_U<=0.1], collapse=", ")))

#bubble plotting

#set some parameters
sizeText <- 20
bubblePal <- brewer.pal(4, "BuPu")

bubblePlot <- ggplot(tableBoth, aes(x=pvalNU, y=pvalNH)) +
	geom_point(alpha=0.8, shape=21, color="black", fill=coloursPal[3], size=5) +
	xlim(0, 1) +
	ylim(0, 1) +
	expand_limits(x=0) +
	labs(x="P-value NU-U", y="P-value NH-H") +
	ggtitle("") +
	geom_hline(yintercept=0.05, color=coloursPal[2], lwd=0.8, linetype="dashed") +
	geom_vline(xintercept=0.05, color=coloursPal[2], lwd=0.8, linetype="dashed") + 
	theme_bw() +
	theme(plot.title=element_text(size=sizeText, face="bold", hjust=0.5),
		axis.text.x=element_text(hjust=1, size=sizeText, face="bold"),
		axis.text.y=element_text(size=sizeText, face="bold"),
		axis.text=element_text(size=sizeText, face="bold"),
		axis.title=element_text(size=sizeText, face="bold"),
		legend.position="top",
		legend.title=element_text(size=sizeText),
		legend.text=element_text(size=sizeText)) +
	geom_text_repel(label=tableBoth$proteins, fontface="italic", size=4, colour="black")

print("bubblePlot")
pdf(paste(saveIn, "bubblePlot.pdf", sep=""), width=12, height=9, font="Arial")
plot(bubblePlot)
dev.off()

############################################################# PLOTTING BEST PROTEIN LINEPLOT

######################## COMPUTING THE AVERAGE DIFFERENCE BETWEEN DEGREES IN THE CONDITIONS
#for U
avgDegU <- vector()
avgDegNU <- vector()
#loop through siginificant U
for (prt in signifU) {
  #get the degree for a specific protein
  subsampleU <- allDegU[which(allDegU$proteins==prt), ]
  avgDegU <- append(avgDegU, mean(subsampleU[grep("U", subsampleU$subject), ]$degree))
  avgDegNU <- append(avgDegNU, mean(subsampleU[grep("N", subsampleU$subject), ]$degree))
}

#for H
avgDegH <- vector()
avgDegNH <- vector()
#loop through siginificant H
for (prt in signifH) {
  subsampleH <- allDegH[which(allDegH$proteins==prt), ]
  avgDegNH <- append(avgDegNH, mean(subsampleH[grep("N", subsampleH$subject), ]$degree))
  avgDegH <- append(avgDegH, mean(subsampleH[grep("H", subsampleH$subject), ]$degree))
}

#computing the differences in terms of degree for the proteins that are considered significant
#first we compute the distance between each degree, then the average distance is rounded to 4 digits, to be plotted.
diffH <- round(mean(abs(avgDegNH - avgDegH)), digits=4)
diffU <- round(mean(abs(avgDegNU - avgDegU)), digits=4)
#we expect that the average distance for the U samples is greater than the average distance for the H sample
#diffU > diffH

######### Prepare to plot
signifHPos <- vector()
#find position of significant proteins!
for (prt in signifH) {
  pTmp <- grep(prt, allDegH$proteins)
  #store positions
  signifHPos <- append(signifHPos, pTmp)
}
#subset allDeg using the positions of the significant proteins
ssADegH <- allDegH[signifHPos, ]
#substitue N with NH, for plotting purposes
ssADegH$subject <- gsub("N", "NH", ssADegH$subject)

#removing unwanted characters and preparing new columns for the plot
ssADegH$group <- as.factor(gsub("[0-9]", "", ssADegH$subject))
#transforming the group in number
ssADegH$subgroup <- as.numeric(as.factor(gsub("[a-z]", "", ssADegH$subject)))
#finally, get only N or H
ssADegH <- ssADegH[grep("N|H", ssADegH$subject), ]

lPlot <- ggplot(data=ssADegH, aes(x=proteins, y=log(degree), group=subject, color=group, alpha=subgroup)) + 
  geom_line(lwd=0.8) +
  geom_point()+
  scale_alpha_continuous(range=c(1, 1)) +
  ylim(-11, -4) +
  scale_color_manual(values=c(colorsPal[1], colorsPal[2])) +
  guides(alpha=FALSE) +
  labs(x="", y="log(Degree)") +
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
  annotate(geom="text", x=3, y=-10, label=paste0("log(Average difference)=", diffH), hjust="left", size=6)

print("linePlotNH-ing")
pdf(paste0(saveIn, "linePlotNH.pdf"), width=6, height=4.5, font="Arial")
plot(lPlot)
dev.off()

########## REPEAT FOR U

signifUPos <- vector()
for (prt in signifU) {
  pTmp <- grep(prt, allDegU$proteins)
  signifUPos <- append(signifUPos, pTmp)
}
#subset allDeg
ssADegU <- allDegU[signifUPos, ]
#substitue N with NU, for plotting purposes
ssADegU$subject <- gsub("N", "NU", ssADegU$subject)

ssADegU$group <- as.factor(gsub("[0-9]", "", ssADegU$subject))
ssADegU$subgroup <- as.numeric(as.factor(gsub("[a-z]", "", ssADegU$subject)))
ssADegU <- ssADegU[grep("N|U", ssADegU$subject), ]

lPlot <- ggplot(data=ssADegU, aes(x=proteins, y=log(degree), group=subject, color=group, alpha=subgroup)) + 
  geom_line(lwd=0.8) +
  geom_point()+
  scale_alpha_continuous(range=c(1, 1)) +
  ylim(-11, -4) +
  #invert colours since it plots by alphabetical order!
  scale_color_manual(values=c(colorsPal[2], colorsPal[1])) +
  guides(alpha=FALSE) +
  labs(x="", y="log(Degree)") +
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
  annotate(geom="text", x=4, y=-10, label=paste0("log(Average difference)=", diffU), hjust="left", size=6)

print("linePlotNU-ing")
pdf(paste0(saveIn, "linePlotNU.pdf"), width=6, height=4.5, font="Arial")
plot(lPlot)
dev.off()

##################################### COMPUTING BEST, AVERAGE DEGREE

#compute the average degree for N and U and put everything together
averageDegU <- cbind.data.frame(proteinName=subSampleForU[, 1],  meanN=rowMeans(subSampleForU[, grep("N", colnames(subSampleForU))]), meanU=rowMeans(subSampleForU[, grep("U", colnames(subSampleForU))]) )
#order by meanN
bestDegU <- averageDegU[order(averageDegU$meanN, decreasing=T), ]

#compute the average degree for N and H and put everything together
averageDegH <- cbind.data.frame(proteinName=subSampleForH[, 1],  meanN=rowMeans(subSampleForH[, grep("N", colnames(subSampleForH))]), meanH=rowMeans(subSampleForH[, grep("H", colnames(subSampleForH))]) )
#order by meanN
bestDegH <- averageDegH[order(averageDegH$meanN, decreasing=T), ]

############################################### UPSET PLOT USING FOLD CHANGE

############################################### 5)
#upset fold change H
fcNH <- read.csv(paste0(patH, "NH_foldchange_LOG.csv"), sep="\t", header=T, stringsAsFactors=T)

#upset fold change U
fcNU <- read.csv(paste0(patH, "NU_foldchange_LOG.csv"), sep="\t", header=T, stringsAsFactors=T)

#determine the length of the set of proteins to compare using SFCPs
if (length(fcNU$x) >=length(fcNH$x)) {
	pMatch <- length(fcNU$x)
} else {
	pMatch <- length(fcNH$x)
}
#and SSPs
if (pMatch < length(signifH)) {
	pMatch <- length(signifH)
} else if (pMatch < length(signifU)){
	pMatch <- length(signifU)
}

listInputH <- list(list())
#they were ordered by N, so best healthy are the first n proteins
listInputH$ControlPs <- as.character(head(bestDegH$proteinName, pMatch))
#then order the table according to meanH and get the first n proteins
listInputH$TreatedPs <- as.character(head(averageDegH$proteinName[order(averageDegH$meanH, decreasing=T)], pMatch))
#and finally add the set of those proteins with relevant fold change
listInputH$SFCPs <- as.character(fcNH$x)

print("upsetFoldC_H")
pdf(paste0(saveIn, "upsetFoldC_H.pdf"), width=8, height=6, font="Arial")
upset(fromList(listInputH), order.by="freq", nsets=4, decreasing=T, sets=c("TreatedPs", "ControlPs", "SFCPs"), keep.order=T, text.scale=3,  sets.x.label="Set size", mainbar.y.label="Shared proteins", mainbar.y.max=15, set_size.scale_max=15, 
	queries=list(
		list(query=intersects, params=list("TreatedPs", "ControlPs"), color=coloursPal[5], active=T),
		list(query=intersects, params=list("TreatedPs", "ControlPs", "SFCPs"), color=coloursPal[2], active=T)
	))
grid.text("NH-H",x=0.65, y=0.95, gp=gpar(fontsize=20))
dev.off()

############################################### 6)

listInputU <- list(list())
#they were ordered by N, so best healthy are the first n proteins
listInputU$ControlPs <- as.character(head(bestDegU$proteinName, pMatch))
#then order the table according to meanH and get the first n proteins
listInputU$TreatedPs <- as.character(head(averageDegU$proteinName[order(averageDegU$meanU, decreasing=T)], pMatch))
#and finally add the set of those proteins with relevant fold change
listInputU$SFCPs <- as.character(fcNU$x)

print("upsetFoldC_U")
pdf(paste0(saveIn, "upsetFoldC_U.pdf"), width=8, height=6, font="Arial")
upset(fromList(listInputU), order.by="freq", nsets=4, decreasing=T, sets=c("TreatedPs", "ControlPs", "SFCPs"), keep.order=T, text.scale=3,  sets.x.label="Set size", mainbar.y.label="Shared proteins", mainbar.y.max=15, set_size.scale_max=15, 
	queries=list(
	list(query=intersects, params=list("TreatedPs", "ControlPs"), color=coloursPal[6], active=T),
	list(query=intersects, params=list("TreatedPs", "SFCPs"), color=coloursPal[5], active=T),
	list(query=intersects, params=list("ControlPs", "SFCPs"), color=coloursPal[2], active=T)
	))
grid.text("NU-U",x=0.65, y=0.95, gp=gpar(fontsize=20))
dev.off()

############################################### 7)

listInputH <- list(list())
#they were ordered by N, so best healthy are the first n proteins
listInputH$ControlPs <- as.character(head(bestDegH$proteinName, pMatch))
#then order the table according to meanH and get the first n proteins
listInputH$SSPs <- signifH
#and finally add the set of those proteins with relevant fold change
listInputH$SFCPs <- as.character(fcNH$x)

print("upsetFoldC_H_SIGNIF")
pdf(paste0(saveIn, "upsetFoldC_H_SIGNIF.pdf"), width=8, height=6, font="Arial")
upset(fromList(listInputH), order.by="freq", nsets=4, decreasing=T, sets=c("SSPs", "ControlPs", "SFCPs"), keep.order=T, text.scale=3,  sets.x.label="Set size", mainbar.y.label="Shared proteins", mainbar.y.max=15, set_size.scale_max=15,
	queries=list(
	list(query=intersects, params=list("SSPs", "SFCPs"), color=coloursPal[2], active=T),
	list(query=intersects, params=list("SSPs", "ControlPs"), color=coloursPal[5], active=T)))
grid.text("NH-H",x=0.65, y=0.95, gp=gpar(fontsize=20))
dev.off()

############################################### 8)

listInputU <- list(list())
#they were ordered by N, so best healthy are the first n proteins
listInputU$ControlPs <- as.character(head(bestDegU$proteinName, pMatch))
#then order the table according to meanH and get the first n proteins
listInputU$SSPs <- signifU
#and finally add the set of those proteins with relevant fold change
listInputU$SFCPs <- as.character(fcNU$x)

print("upsetFoldC_U_SIGNIF")
pdf(paste0(saveIn, "upsetFoldC_U_SIGNIF.pdf"), width=8, height=6, font="Arial")
upset(fromList(listInputU), order.by="freq", nsets=4, decreasing=T, sets=c("SSPs", "ControlPs", "SFCPs"), keep.order=T, text.scale=3,  sets.x.label="Set size", mainbar.y.label="Shared proteins", mainbar.y.max=15, set_size.scale_max=15, 
	queries=list(
	list(query=intersects, params=list("SSPs", "SFCPs"), color=coloursPal[2], active=T)
	))
grid.text("NU-U",x=0.65, y=0.95, gp=gpar(fontsize=20))
dev.off()

############################################# UPSET WITH EVERYTHING #############################################

listInputAll <-list(list())

listInputAll$ControlPs_H <- as.character(head(bestDegH$proteinName, pMatch))
listInputAll$TreatedPs_H <- as.character(head(averageDegH$proteinName[order(averageDegH$meanH, decreasing=T)], pMatch))
listInputAll$ControlPs_U <- as.character(head(bestDegU$proteinName, pMatch))
listInputAll$TreatedPs_U <- as.character(head(averageDegU$proteinName[order(averageDegU$meanU, decreasing=T)], pMatch))

listInputAll$SSPs_U <- signifU
listInputAll$SSPs_H <- signifH
listInputAll$SFCPs_H <- as.character(fcNH$x)
listInputAll$SFCPs_U <- as.character(fcNU$x)

print("upsetALL_fullversion")
pdf(paste0(saveIn, "upsetALL_fullversion.pdf"), width=8, height=6, font="Arial")
upset(fromList(listInputAll), order.by="freq", nsets=4, decreasing=T, sets=c("ControlPs_H", "TreatedPs_H", "ControlPs_U", "TreatedPs_U", "SSPs_H", "SFCPs_H", "SSPs_U", "SFCPs_U"), keep.order=T, text.scale=1.5, sets.x.label="Set size", mainbar.y.label="Shared proteins", mainbar.y.max=15, set_size.scale_max=15,
	queries=list(
	list(query=intersects, params=list("ControlPs_H", "TreatedPs_H", "ControlPs_U", "TreatedPs_U"), color=coloursPal[2], active=T),
	list(query=intersects, params=list("ControlPs_H", "TreatedPs_H", "ControlPs_U", "TreatedPs_U", "SSPs_H"), color=coloursPal[6], active=T),
	list(query=intersects, params=list("ControlPs_H", "TreatedPs_H", "ControlPs_U", "SFCPs_H"), color=coloursPal[5], active=T)
	))
dev.off()

############### PLOT NETWORK WITH MOST SIGNIFICANT PROTEINS AND FIRST NEIGHBOURS

stringInfo <- read.csv(paste0(patH, "string_over400_experimental.tsv"), sep="\t", header=T, stringsAsFactors=T)

#create edge list
edgeL <- as.matrix(cbind.data.frame(as.character(stringInfo$X.node1), as.character(stringInfo$node2)))
#and build the network
networkTmp <- graph_from_edgelist(edgeL)
#get biggest component
network <- decompose(networkTmp)[[1]]

######### for U samples
neighListU <- signifU
#get first neighbours of signifU proteins
for (vrtx in signifU) {
	#get the neighbours of vertex
	neighListU <- append(neighListU, neighbors(network, vrtx, mode="all")$name)
}

#export the neighbours of SSPs
neighsOfSignifU <- unique(neighListU)[!unique(neighListU) %in% signifU]
write.table(neighsOfSignifU, paste0(patH, "neighsOfSignifU.csv"), quote=F, row.names=F, col.names=F)

#create the network for SSPs and their neighbours
firstNeighsU <- induced.subgraph(network, unique(neighListU))

for (v in c(1:length(V(firstNeighsU)))) {
	#in pink both SSPs and SFCPs	
	if (V(firstNeighsU)$name[v] %in% signifU && V(firstNeighsU)$name[v] %in% fcNU$x) {
		V(firstNeighsU)$colour[v] <- coloursPal[4]
	#in green SSPs only
	} else if (V(firstNeighsU)$name[v] %in% signifU) {
		V(firstNeighsU)$colour[v] <- coloursPal[5]
	#in yellow SFCPs only
	} else if (V(firstNeighsU)$name[v] %in% fcNU$x) {
		V(firstNeighsU)$colour[v] <- coloursPal[6]
	} else if (V(firstNeighsU)$name[v] %in% c("ACTN2", "TPM1", "TPM2")) {
		V(firstNeighsU)$colour[v] <- coloursPal[2]
	} else {
		V(firstNeighsU)$colour[v] <- "#ffffff"
	}
}

# U NET
lytNicely <- layout_nicely(firstNeighsU)
print("signifU_orig")
pdf(paste0(saveIn, "signifU_orig.pdf"), font="Arial")
plot(firstNeighsU, vertex.label="", vertex.size=15, edge.arrow.size=0, edge.color="dimgrey", edge.size=4, vertex.color=V(firstNeighsU)$colour, layout=lytNicely)
dev.off()

######### and for H samples
neighListH <- signifH
#get first neighbours of signifU proteins
for (vrtx in signifH) {
	#get the neighbours of vertex
	neighListH <- append(neighListH, neighbors(network, vrtx, mode="all")$name)
}

#export the neighbours of SSPs
neighsOfSignifH <- unique(neighListH)[!unique(neighListH) %in% signifH]
write.table(neighsOfSignifH, paste0(patH, "neighsOfSignifH.csv"), quote=F, row.names=F, col.names=F)

#create the network for SSPs and their neighbours
firstNeighsH <- induced.subgraph(network, unique(neighListH))

for (v in 1:length(V(firstNeighsH))) {
	#in pink both SSPs and SFCPs
	if (V(firstNeighsH)$name[v] %in% signifH && V(firstNeighsH)$name[v] %in% fcNH$x) {
		V(firstNeighsH)$colour[v] <- coloursPal[4]
	#in green SSPs only
	} else if (V(firstNeighsH)$name[v] %in% signifH) {
		V(firstNeighsH)$colour[v] <- coloursPal[5]
	#in yellow SFCPs only
	} else if (V(firstNeighsH)$name[v] %in% fcNH$x) {
		V(firstNeighsH)$colour[v] <- coloursPal[6]
	} else if (V(firstNeighsH)$name[v] %in% c("ACTN2", "VIM", "TPM1", "ACTN1")) {
		V(firstNeighsH)$colour[v] <- coloursPal[2]
	} else {
		V(firstNeighsH)$colour[v] <- "#ffffff"
	}
}

# H NET
lytLgl <- layout_with_lgl(firstNeighsH)
print("signifH_orig")
pdf(paste0(saveIn, "signifH_orig.pdf"), font="Arial")
plot(firstNeighsH, vertex.label="", vertex.size=15, edge.arrow.size=0, edge.color="dimgrey", edge.size=4, vertex.color=V(firstNeighsH)$colour, layout=lytLgl)
dev.off()

############################################################# volcano plotting as shown here:
#https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html

#set palette with 8 colours for colouring plots
coloursPal <- brewer.pal(8, "Paired")

#now import fold change values
fcNUVals <- read.csv(paste0(patH, "NU_foldchange_values_LOG.csv"), sep="\t", header=T, stringsAsFactors=T)
fcNHVals <- read.csv(paste0(patH, "NH_foldchange_values_LOG.csv"), sep="\t", header=T, stringsAsFactors=T)

tmpV <- merge(tableBoth, fcNUVals, by="proteins")

volcanoData_H <- merge(tmpV, fcNHVals, by="proteins")

#set threshold equal to 1, to compute foldChange
fcThreshU <- fcThreshH <- read.csv(paste0(patH, "FC_threshold.csv"), sep="\t", header=T, stringsAsFactors=F)

############################################################# NH

volcanoData_H$expr <- "NO"
# if NHfold > mean(NHfold) and pvalNH < 0.05, set as "UP" 
volcanoData_H$expr[volcanoData_H$NHfold > fcThreshH$th & volcanoData_H$pvalNH < 0.05] <- "UP"
# if NHfold < mean(NHfold) and pvalNH < 0.05, set as "DOWN"
volcanoData_H$expr[volcanoData_H$NHfold < -fcThreshH$th & volcanoData_H$pvalNH < 0.05] <- "DOWN"

volcanoData_H$label <- ""
volcanoData_H$label[volcanoData_H$expr !="NO"] <- volcanoData_H$proteins[volcanoData_H$expr !="NO"]

# plot adding up all layers we have seen so far
p <- ggplot(data=volcanoData_H, aes(x=NHfold, y=-log10(pvalNH), col=expr, label=label)) +
        geom_point(size=2) +
	      geom_text_repel(label=volcanoData_H$label, size=5, fontface="bold", force=1, point.padding=unit(1,'lines'), direction='y', nudge_x=0.2, segment.size=0.2) +
        # scale_color_manual(values=c(coloursPal[8], coloursPal[2]), breaks=c("DOWN", "NO"), guide=F) +
        scale_color_manual(values=c(coloursPal[3], coloursPal[6]), guide=F) +
        geom_vline(xintercept=c(-fcThreshH$th, fcThreshH$th), col=coloursPal[5], lwd=0.4) +
        geom_hline(yintercept=-log10(0.05), col=coloursPal[5], lwd=0.4) +
        xlim(c(-3,3)) +
        ylim(c(0, 3.5)) +
       	labs(x="NH-H log2(fold Change)", y="-log10(P-value)", color="Fold Change") +
        theme_bw() +
	theme(plot.title=element_text(size=sizeText, face="bold", hjust=0.5),
		axis.text.x=element_text(hjust=1, size=sizeText, face="bold"),
		axis.text.y=element_text(size=sizeText, face="bold"),
		axis.text=element_text(size=sizeText, face="bold"),
		axis.title=element_text(size=sizeText, face="bold"),
		legend.position="top",
		legend.title=element_text(size=sizeText),
		legend.text=element_text(size=sizeText))
#		+ guides(colour=guide_legend(override.aes=list(size=c(3, 3))))

		
pdf(paste(saveIn, "volcanoPlot_NH.pdf", sep=""), width=6, height=9, font="Arial")
plot(p)
dev.off()

############################################################# NU
#now import fold change values
fcNUVals <- read.csv(paste0(patH, "NU_foldchange_values_LOG.csv"), sep="\t", header=T, stringsAsFactors=T)
fcNHVals <- read.csv(paste0(patH, "NH_foldchange_values_LOG.csv"), sep="\t", header=T, stringsAsFactors=T)

tmpV <- merge(tableBoth, fcNUVals, by="proteins")

volcanoData_U <- merge(tmpV, fcNHVals, by="proteins")

volcanoData_U$expr <- "NO"
# if NUfold > mean(NUfold) and pvalNU < 0.05, set as "UP" 
volcanoData_U$expr[volcanoData_U$NUfold > fcThreshU$th & volcanoData_U$pvalNU < 0.05] <- "UP"
# if NUfold < mean(NUfold) and pvalNU < 0.05, set as "DOWN"
volcanoData_U$expr[volcanoData_U$NUfold < -fcThreshU$th & volcanoData_U$pvalNU < 0.05] <- "DOWN"

volcanoData_U$label <- ""
volcanoData_U$label[volcanoData_U$expr !="NO"] <- volcanoData_U$proteins[volcanoData_U$expr !="NO"]

# plot adding up all layers we have seen so far
p <- ggplot(data=volcanoData_U, aes(x=NUfold, y=-log10(pvalNU), col=expr, label=label)) +
        geom_point(size=2) +
	geom_text_repel(label=volcanoData_U$label, size=5, fontface="bold", force=1, point.padding=unit(1,'lines'), direction='y', nudge_x=0.1, segment.size=0.3) +
#        scale_color_manual(values=c(coloursPal[1], coloursPal[8], coloursPal[2]), breaks=c("UP", "DOWN"), guide=F) +
        scale_color_manual(values=c(coloursPal[2], coloursPal[3], coloursPal[6]), guide=F) +
        geom_vline(xintercept=c(-fcThreshU$th, fcThreshU$th), col=coloursPal[5], lwd=0.4) +
        geom_hline(yintercept=-log10(0.05), col=coloursPal[5], lwd=0.4) +
        xlim(c(-3,3)) +
        ylim(c(0, 3.5)) +
       	labs(x="NU-U log2(fold Change)", y="-log10(P-value)", color="Fold Change") +
        theme_bw() +
	theme(plot.title=element_text(size=sizeText, face="bold", hjust=0.5),
		axis.text.x=element_text(hjust=1, size=sizeText, face="bold"),
		axis.text.y=element_text(size=sizeText, face="bold"),
		axis.text=element_text(size=sizeText, face="bold"),
		axis.title=element_text(size=sizeText, face="bold"),
		legend.position="top",
		legend.title=element_text(size=sizeText),
		legend.text=element_text(size=sizeText))
#		+ guides(colour=guide_legend(override.aes=list(size=c(3, 3))))
	
pdf(paste(saveIn, "volcanoPlot_NU.pdf", sep=""), width=6, height=9, font="Arial")
plot(p)
dev.off()

