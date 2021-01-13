#############
# GabrieleT #
#   2020    #
#############

#load libraries
library("ggplot2")
library("reshape")
library("ggrepel")
library("RColorBrewer")
library("igraph")
library("extrafont")
loadfonts()

#set seed
set.seed(131)

#set paths
patH <- "/home/cbmc/wnnets/"
saveIn <- "/home/cbmc/wnnets/figures/"

############################################################# PLOTTING DEGREE EXAMPLE BARPLOT

#set palette with 8 colours for colouring plots
coloursPal <- brewer.pal(8, "Set2")

data <- read.csv(paste0(patH, "tableExampleDeg.csv"), header=T, sep=",")
filtered <- data[, c(1, 3, 5, 7)]
#rename for fitting latex document
colnames(filtered) <- c("node", "(a)", "(b)", "(c)")
#set some parameters
sizeText <- 20

plotIt <- cbind.data.frame(node=rep(c("1", "2", "3", "4", "5"), 3), melt(filtered[, -1]))

plotIt$node <- factor(plotIt$node, levels=unique(plotIt$node))

p <- ggplot(plotIt, aes(fill=node, y=value, x=variable)) + 
	geom_bar(position=position_dodge(width=0.6), stat="identity", color="black", lwd=0.6, width=0.6) +
	scale_fill_manual(values=coloursPal[c(1:5)]) +
	xlab("") +
	ylab("degree") +
	ggtitle("") +
	theme_bw() +
	theme(plot.title=element_text(size=sizeText, face="bold", hjust=0.5),
		axis.text.x=element_text(hjust=1, size=sizeText, face="bold"),
		axis.text.y=element_text(size=sizeText, face="bold"),
		axis.text=element_text(size=sizeText, face="bold"),
		axis.title=element_text(size=sizeText, face="bold"),
		legend.position="top",
		legend.title=element_text(size=sizeText),
		legend.text=element_text(size=sizeText))

pdf(paste0(saveIn, "columnsExample.pdf"), width=8, height=16, font="Arial")
plot(p)
dev.off()

########################################## plot networks with coloured nodes ORIGINAL

df <- data.frame("from"=c("1", "2", "3", "3", "4"), 
	"to"=c("2", "3", "4", "5", "5"))

meta <- data.frame("name"=c("1", "2", "3", "4", "5"), 
	"x"=c(5, 8, 12, 9, 7), 
	"y"=c(7, 5, 8, 12, 9))

g <- graph.data.frame(df, directed=FALSE, vertices=meta)
lo <- layout.norm(as.matrix(meta[,2:3]))

for (v in 1:length(V(g))) {
	V(g)$colour[v] <- coloursPal[v]
}

E(g)$colour <- c("dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey")

E(g)$width <- c(2, 2, 2, 2, 2)

pdf(paste0(saveIn, "original.pdf"), font="Arial")
plot.igraph(g,
	layout=lo,
	vertex.color=V(g)$colour,
	vertex.size=50,
	vertex.label.cex=3.5,
	vertex.label.font=2,
	vertex.label.color="black",
	edges.color=E(g)$colour,
	edges.size=E(g)$width
)
dev.off()

########################################## plot networks with coloured nodes MULTI FIVE PARTIAL

df <- data.frame("from"=c("1", "2", "3", "3", "4", "5", "5", "5.1"), 
	"to"=c("2", "3", "4", "5", "5", "5.1", "5.2", "5.2"))

meta <- data.frame("name"=c("1", "2", "3", "4", "5", "5.1", "5.2"), 
	"x"=c(5, 8, 12, 9, 7, 5, 6),
	"y"=c(7, 5, 8, 12, 9, 11, 13))

g <- graph.data.frame(df, directed=FALSE, vertices=meta)
lo <- layout.norm(as.matrix(meta[,2:3]))

V(g)$colour[1] <- coloursPal[1]
V(g)$colour[2] <- coloursPal[2]
V(g)$colour[3] <- coloursPal[3]
V(g)$colour[4] <- coloursPal[4]
V(g)$colour[5] <- coloursPal[5]
V(g)$colour[6] <- coloursPal[5]
V(g)$colour[7] <- coloursPal[5]

E(g)$width <- c(2, 2, 2, 2, 2, 6, 6, 6)

E(g)$colour <- c("dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "black", "black", "black")

pdf(paste0(saveIn, "multiFive_no_neighs.pdf"), font="Arial")
plot.igraph(g,
	layout=lo,
	vertex.color=V(g)$colour,
	vertex.size=50,
	edges.color=E(g)$colour,
	edges.size=E(g)$width,
	vertex.label.cex=3.5,
	vertex.label.font=2,
	vertex.label.color="black"
)
dev.off()

########################################## plot networks with coloured nodes MULTI FIVE COMPLETE

df <- data.frame("from"=c("1", "2", "3", "3", "4", "5", "5", "5.1", "5.1", "5.1", "5.2", "5.2"), 
	"to"=c("2", "3", "4", "5", "5", "5.1", "5.2", "5.2", "3", "4", "3", "4"))

meta <- data.frame("name"=c("1", "2", "3", "4", "5", "5.1", "5.2"), 
	"x"=c(5, 8, 12, 9, 7, 5, 6),
	"y"=c(7, 5, 8, 12, 9, 11, 13))

g <- graph.data.frame(df, directed=FALSE, vertices=meta)
lo <- layout.norm(as.matrix(meta[,2:3]))

V(g)$colour[1] <- coloursPal[1]
V(g)$colour[2] <- coloursPal[2]
V(g)$colour[3] <- coloursPal[3]
V(g)$colour[4] <- coloursPal[4]
V(g)$colour[5] <- coloursPal[5]
V(g)$colour[6] <- coloursPal[5]
V(g)$colour[7] <- coloursPal[5]

E(g)$width <- c(2, 2, 2, 2, 2, 6, 6, 6, 6, 6, 6, 6)

E(g)$colour <- c("dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "black", "black", "black", "black", "black", "black", "black")

pdf(paste0(saveIn, "multiFive.pdf"), font="Arial")
plot.igraph(g,
	layout=lo,
	vertex.color=V(g)$colour,
	vertex.size=50,
	edges.color=E(g)$colour,
	edges.size=E(g)$width,
	vertex.label.cex=3.5,
	vertex.label.font=2,
	vertex.label.color="black"
)
dev.off()

########################################## plot networks with coloured nodes MULTI THREE COMPLETE

df <- data.frame("from"=c("1", "2", "3", "3", "4", "3", "3", "3.1", "3.1", "3.1", "3.2", "3.2", "3.1", "3.2"), 
	"to"=c("2", "3", "4", "5", "5", "3.2", "3.1", "3.2", "5", "4", "5", "4", "2", "2"))

meta <- data.frame("name"=c("1", "2", "3", "4", "5", "3.1", "3.2"),
	"x"=c(5, 8, 12, 9, 7, 8, 11),
	"y"=c(7, 5, 8, 12, 9, 7, 6))

g <- graph.data.frame(df, directed=FALSE, vertices=meta)
lo <- layout.norm(as.matrix(meta[,2:3]))

V(g)$colour[1] <- coloursPal[1]
V(g)$colour[2] <- coloursPal[2]
V(g)$colour[3] <- coloursPal[3]
V(g)$colour[4] <- coloursPal[4]
V(g)$colour[5] <- coloursPal[5]
V(g)$colour[6] <- coloursPal[3]
V(g)$colour[7] <- coloursPal[3]

E(g)$width <- c(2, 2, 2, 2, 2, 2, 6, 6, 6, 6, 6, 6, 6, 6)

E(g)$colour <- c("dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "black", "black", "black", "black", "black", "black", "black", "black")

pdf(paste0(saveIn, "multiThree.pdf"), font="Arial")
plot.igraph(g,
	layout=lo,
	vertex.color=V(g)$colour,
	vertex.size=50,
	edges.color=E(g)$colour,
	edges.size=E(g)$width,
	vertex.label.cex=3.5,
	vertex.label.font=2,
	vertex.label.color="black"
)

dev.off()

###################################################################################################################
########################################## plot networks with coloured nodes MULTI FIVE COMPLETE FOR DEGREE EXAMPLE

df <- data.frame("from"=c("1", "2", "3", "3", "4", "5", "5", "5.1", "5.1", "5.1", "5.2", "5.2"), 
	"to"=c("2", "3", "4", "5", "5", "5.1", "5.2", "5.2", "3", "4", "3", "4"))

meta <- data.frame("name"=c("1", "2", "3", "4", "5", "5.1", "5.2"), 
	"x"=c(5, 8, 12, 9, 7, 5, 6),
	"y"=c(7, 5, 8, 12, 9, 11, 13))

g <- graph.data.frame(df, directed=FALSE, vertices=meta)
lo <- layout.norm(as.matrix(meta[,2:3]))

V(g)$colour[1] <- coloursPal[1]
V(g)$colour[2] <- coloursPal[2]
V(g)$colour[3] <- coloursPal[3]
V(g)$colour[4] <- coloursPal[4]
V(g)$colour[5] <- coloursPal[5]
V(g)$colour[6] <- coloursPal[5]
V(g)$colour[7] <- coloursPal[5]

E(g)$width <- c(2, 6, 6, 6, 2, 2, 2, 2, 6, 2, 6, 2)

E(g)$colour <- c("dimgrey", "black", "black", "black", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "black", "dimgrey", "black", "dimgrey")

pdf(paste0(saveIn, "multiFive_degreeExample.pdf"), font="Arial")
plot.igraph(g,
	layout=lo,
	vertex.color=V(g)$colour,
	vertex.size=50,
	edges.color=E(g)$colour,
	edges.size=E(g)$width,
	vertex.label.cex=3.5,
	vertex.label.font=2,
	vertex.label.color="black"
)
dev.off()

########################################## plot networks with coloured nodes MULTI THREE COMPLETE FOR DEGREE EXAMPLE

df <- data.frame("from"=c("1", "4", "3", "3", "3.1", "3.1", "3.1", "3.2", "3.2", "3.1", "3.2", "2", "3", "3"), 
	           "to"=c("2", "5", "3.1", "3.2", "3.2", "5", "4", "5", "4", "2", "2", "3", "4", "5"))

meta <- data.frame("name"=c("1", "2", "3", "4", "5", "3.1", "3.2"),
	"x"=c(5, 8, 12, 9, 7, 8, 11),
	"y"=c(7, 5, 8, 12, 9, 7, 6))

g <- graph.data.frame(df, directed=FALSE, vertices=meta)
lo <- layout.norm(as.matrix(meta[,2:3]))

V(g)$colour[1] <- coloursPal[1]
V(g)$colour[2] <- coloursPal[2]
V(g)$colour[3] <- coloursPal[3]
V(g)$colour[4] <- coloursPal[4]
V(g)$colour[5] <- coloursPal[5]
V(g)$colour[6] <- coloursPal[3]
V(g)$colour[7] <- coloursPal[3]

E(g)$width <- c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 6, 6, 6)

E(g)$colour <- c("dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "black", "black", "black")

pdf(paste0(saveIn, "multiThree_degreeExample.pdf"), font="Arial")
plot.igraph(g,
	layout=lo,
	vertex.color=V(g)$colour,
	vertex.size=50,
	edges.color=E(g)$colour,
	edges.size=E(g)$width,
	vertex.label.cex=3.5,
	vertex.label.font=2,
	vertex.label.color="black"
)

dev.off()

########################################## plot networks with coloured nodes THE EDGE ORIGINAL

#set a new palette
coloursPal <- brewer.pal(5, "Paired")

df <- data.frame("from"=c("1", "2"), 
	"to"=c("2", "3"))

meta <- data.frame("name"=c("1", "2", "3"), 
	"x"=c(5, 8, 9), 
	"y"=c(7, 5, 6.5))

g <- graph.data.frame(df, directed=FALSE, vertices=meta)
lo <- layout.norm(as.matrix(meta[,2:3]))

V(g)$colour[1] <- coloursPal[3]
V(g)$colour[2] <- coloursPal[5]
V(g)$colour[3] <- coloursPal[1]

E(g)$colour <- c("dimgrey", "dimgrey")

E(g)$width <- c(2, 2)

pdf(paste0(saveIn, "theEdgePlain.pdf"), font="Arial")
plot.igraph(g,
	layout=lo,
	vertex.color=V(g)$colour,
	vertex.size=50,
	edges.color=E(g)$colour,
	edges.size=E(g)$width,
	vertex.label.cex=3.5,
	vertex.label.font=2,
	vertex.label.color="black"
)
dev.off()

########################################## plot networks with coloured nodes WRONGH PATH

df <- data.frame("from"=c("1", "1.1", "2"), 
	"to"=c("2", "2", "3"))

meta <- data.frame("name"=c("1", "2", "1.1", "3"), 
	"x"=c(5, 8, 5, 9),
	"y"=c(7, 5, 4, 6.5))

g <- graph.data.frame(df, directed=FALSE, vertices=meta)
lo <- layout.norm(as.matrix(meta[,2:3]))

V(g)$colour[1] <- coloursPal[3]
V(g)$colour[2] <- coloursPal[5]
V(g)$colour[3] <- coloursPal[3]
V(g)$colour[4] <- coloursPal[1]

E(g)$width <- c(6, 6, 2)

E(g)$colour <- c("black", "black", "dimgrey")

pdf(paste0(saveIn, "wrongPath.pdf"), font="Arial")
plot.igraph(g,
	layout=lo,
	vertex.color=V(g)$colour,
	vertex.size=50,
	edges.color=E(g)$colour,
	edges.size=E(g)$width,
	vertex.label.cex=3.5,
	vertex.label.font=2,
	vertex.label.color="black"
)
dev.off()

########################################## plot networks with coloured nodes ADDING THAT EDGE

df <- data.frame("from"=c("1", "1.1", "1.1", "2"), 
	"to"=c("2", "2", "1", "3"))

meta <- data.frame("name"=c("1", "2", "1.1", 3), 
	"x"=c(5, 8, 5, 9),
	"y"=c(7, 5, 4, 6.5))

g <- graph.data.frame(df, directed=FALSE, vertices=meta)
lo <- layout.norm(as.matrix(meta[,2:3]))

V(g)$colour[1] <- coloursPal[3]
V(g)$colour[2] <- coloursPal[5]
V(g)$colour[3] <- coloursPal[3]
V(g)$colour[4] <- coloursPal[1]

E(g)$width <- c(2, 2, 6, 2)

E(g)$colour <- c("dimgrey", "dimgrey", "black", "dimgrey")

pdf(paste0(saveIn, "addThatEdge.pdf"), font="Arial")
plot.igraph(g,
	layout=lo,
	vertex.color=V(g)$colour,
	vertex.size=50,
	edges.color=E(g)$colour,
	edges.size=2,
	vertex.label.cex=3.5,
	vertex.label.font=2,
	vertex.label.color="black"
)
dev.off()

########################################################################## PLOTTING FIGURES APPENDIX B ##########################################################################

df <- data.frame("from"=c("1", "2"), 
	"to"=c("2", "3"))
	
meta <- data.frame("name"=c("1", "2", "3"),
	"x"=c(4, 8, 12),
	"y"=c(30, 50, 20))
	
g <- graph.data.frame(df, directed=FALSE, vertices=meta)
lo <- layout.norm(as.matrix(meta[,2:3]))

V(g)$colour[1] <- coloursPal[3]
V(g)$colour[2] <- coloursPal[5]
V(g)$colour[3] <- coloursPal[1]
	
V(g)$label <- c("1", "1", "1")

E(g)$width <- c(2, 2)

E(g)$colour <- c("dimgrey", "dimgrey")

pdf(paste0(saveIn, "real_int_plain.pdf"), font="Arial")
plot.igraph(g,
	layout=lo,
	vertex.color=V(g)$colour,
	vertex.size=50,
	vertex.label.cex=3,
	edges.color=E(g)$colour,
	edges.size=E(g)$width,
	vertex.label.font=2,
	vertex.label.color="black"
)
dev.off()

##########################################################################

dfWeight <- data.frame("from"=c("1", "2"), 
	"to"=c("2", "3"))
	
metaWeight <- data.frame("name"=c("1", "2", "3"),
	"x"=c(4, 8, 12),
	"y"=c(30, 50, 20))
	
g <- graph.data.frame(df, directed=FALSE, vertices=metaWeight)
lo <- layout.norm(as.matrix(metaWeight[,2:3]))

V(g)$colour[1] <- coloursPal[3]
V(g)$colour[2] <- coloursPal[5]
V(g)$colour[3] <- coloursPal[1]
	
V(g)$label <- c("4.5", "1.5", "3.5")

E(g)$width <- c(2, 2)

E(g)$colour <- c("dimgrey", "dimgrey")

pdf(paste0(saveIn, "real_int_weight.pdf"), font="Arial")
plot.igraph(g,
	layout=lo,
	vertex.color=V(g)$colour,
	vertex.size=50,
	vertex.label.cex=3,
	edges.color=E(g)$colour,
	edges.size=E(g)$width,
	vertex.label.font=2,
	vertex.label.color="black"
)
dev.off()
	
##########################################################################

dfMulti <- data.frame("from"=c("1", "2", "1.1", "1.2", "1.3", "1.1", "1.2", "1.3", "1.1", "1.1", "1.2", "3.1", "3.2", "3.1", "3.2", "3.1"), 
	"to"=c("2", "3", "1", "1", "1", "2", "2", "2", "1.2", "1.3", "1.3", "3", "3", "2", "2", "3.2"))
	
metaMulti <- data.frame("name"=c("1", "2", "3", "1.1", "1.2", "1.3", "3.1", "3.2"),
	"x"=c(4, 24, 50, -5, -20, -10, 55, 32),
	"y"=c(28, 32, 30, 18, 25, 35, 24, 26))
	

g <- graph.data.frame(dfMulti, directed=FALSE, vertices=metaMulti)
lo <- layout.norm(as.matrix(metaMulti[,2:3]))

V(g)$colour[1] <- coloursPal[3]
V(g)$colour[2] <- coloursPal[5]
V(g)$colour[3] <- coloursPal[1]
V(g)$colour[4] <- coloursPal[3]
V(g)$colour[5] <- coloursPal[3]
V(g)$colour[6] <- coloursPal[3]
V(g)$colour[7] <- coloursPal[1]
V(g)$colour[8] <- coloursPal[1]

V(g)$label <- c("1", "1", "1", "1", "1", "1", "1", "1")

E(g)$colour <- c("dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey")

E(g)$width <- c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)

pdf(paste0(saveIn, "real_int_multi.pdf"), font="Arial")
plot.igraph(g,
	layout=lo,
	vertex.color=V(g)$colour,
	vertex.size=50,
	vertex.label.cex=3,
	edges.color=E(g)$colour,
	edges.size=E(g)$width,
	vertex.label.font=2,
	vertex.label.color="black"
)
dev.off()

data <- read.csv(paste0(patH, "tableExampleDeg_APPENDIX_B.csv"), header=T, sep=",", stringsAsFactors=F)
filtered <- data[, c(1, 3, 5, 7)]
#rename for fitting latex document
colnames(filtered) <- c("node", "(a)", "(b)", "(c)")
#set some parameters
sizeText <- 18
sizeText <- 20
sizeText <- 18

plotIt <- cbind.data.frame(node=rep(c("1", "2", "3"), 3), melt(filtered[, -1]))

plotIt$node <- factor(plotIt$node, levels=unique(plotIt$node))

p <- ggplot(plotIt, aes(fill=node, y=value, x=variable)) + 
	geom_bar(stat="identity", position=position_dodge(width=0.6), color="black", lwd=0.6, width=0.4) +
	scale_fill_manual(values=coloursPal[c(3, 5, 1)]) +
	xlab("") +
	ylab("degree") +
	ggtitle("") +
	scale_y_continuous(limits=c(0, 8), breaks=seq(0, 8, by=1)) +	
	theme_bw() +
	theme(plot.title=element_text(size=sizeText, face="bold", hjust=0.5),
		axis.text.x=element_text(hjust=1, size=sizeText, face="bold"),
		axis.text.y=element_text(size=sizeText, face="bold"),
		axis.text=element_text(size=sizeText, face="bold"),
		axis.title=element_text(size=sizeText, face="bold"),
		legend.position="top",
		legend.title=element_text(size=sizeText),
		legend.text=element_text(size=sizeText)) +
	guides(colour=guide_legend(reverse=F))


pdf(paste0(saveIn, "real_int_multi_BARPLOT.pdf"), width=8, height=16, font="Arial")
plot(p)
dev.off()

########################################## plot networks FIGURE 1 overview

coloursPal <- brewer.pal(8, "Paired")

df <- data.frame("from"=c("1", "2", "3", "3", "4"), 
	"to"=c("2", "3", "4", "5", "5"))

meta <- data.frame("name"=c("1", "2", "3", "4", "5"), 
	"x"=c(5, 8, 12, 9, 7), 
	"y"=c(7, 5, 8, 12, 9))

g <- graph.data.frame(df, directed=FALSE, vertices=meta)
lo <- layout.norm(as.matrix(meta[,2:3]))

for (v in 1:length(V(g))) {
	V(g)$colour[v] <- coloursPal[v]
}

E(g)$colour <- c("dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey")

E(g)$width <- c(2, 2, 2, 2, 2)

pdf(paste0(saveIn, "original_figureOne.pdf"), font="Arial")
plot.igraph(g,
	layout=lo,
	vertex.color=V(g)$colour,
	vertex.size=50,
	vertex.label.cex=3.5,
	vertex.label.font=2,
	edges.color=E(g)$colour,
	edges.size=E(g)$width,
	vertex.label.color="black"
)
dev.off()

###################################################################################################

df <- data.frame("from"=c("1", "2", "3", "3", "4", "5", "5", "5.1", "5.1", "5.1", "5.2", "5.2"), 
	"to"=c("2", "3", "4", "5", "5", "5.1", "5.2", "5.2", "3", "4", "3", "4"))

meta <- data.frame("name"=c("1", "2", "3", "4", "5", "5.1", "5.2"), 
	"x"=c(5, 8, 12, 9, 7, 5, 6),
	"y"=c(7, 5, 8, 12, 9, 11, 13))

g <- graph.data.frame(df, directed=FALSE, vertices=meta)
lo <- layout.norm(as.matrix(meta[,2:3]))

V(g)$colour[1] <- coloursPal[1]
V(g)$colour[2] <- coloursPal[2]
V(g)$colour[3] <- coloursPal[3]
V(g)$colour[4] <- coloursPal[4]
V(g)$colour[5] <- coloursPal[5]
V(g)$colour[6] <- coloursPal[5]
V(g)$colour[7] <- coloursPal[5]

E(g)$colour <- c("dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey")

E(g)$width <- c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)

pdf(paste0(saveIn, "multiFive_figureOne.pdf"), font="Arial")
plot.igraph(g,
	layout=lo,
	vertex.color=V(g)$colour,
	vertex.size=50,
	vertex.label.cex=3.5,
	edges.color=E(g)$colour,
	edges.size=E(g)$width,
	vertex.label.font=2,
	vertex.label.color="black"
)
dev.off()

###################################################################################################

df <- data.frame("from"=c("1", "2", "3", "3", "4", "3", "3", "3.1", "3.1", "3.1", "3.2", "3.2", "3.1", "3.2"), 
	           "to"=c("2", "3", "4", "5", "5", "3.1", "3.2", "3.2", "5", "4", "5", "4", "2", "2"))

meta <- data.frame("name"=c("1", "2", "3", "4", "5", "3.1", "3.2"),
	"x"=c(5, 8, 12, 9, 7, 8, 11),
	"y"=c(7, 5, 8, 12, 9, 7, 6))

g <- graph.data.frame(df, directed=FALSE, vertices=meta)
lo <- layout.norm(as.matrix(meta[,2:3]))

V(g)$colour[1] <- coloursPal[1]
V(g)$colour[2] <- coloursPal[2]
V(g)$colour[3] <- coloursPal[3]
V(g)$colour[4] <- coloursPal[4]
V(g)$colour[5] <- coloursPal[5]
V(g)$colour[6] <- coloursPal[3]
V(g)$colour[7] <- coloursPal[3]

E(g)$colour <- c("dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey", "dimgrey")

E(g)$width <- c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)

pdf(paste0(saveIn, "multiThree_degreeExample_figureOne.pdf"), font="Arial")
plot.igraph(g,
	layout=lo,
	vertex.color=V(g)$colour,
	vertex.size=50,
	vertex.label.cex=3.5,
	edges.color=E(g)$colour,
	edges.size=E(g)$width,
	vertex.label.font=2,
	vertex.label.color="black"
)

dev.off()
