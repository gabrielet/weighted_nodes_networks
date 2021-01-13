#############
# GabrieleT #
#   2020    #
#############

#set seed
set.seed(131)

#set paths
patH <- "/home/cbmc/wnnets/"

#read data
originalNOT <- read.table(paste0(patH, "original_data_no_normalisation.csv"), header=T, sep="\t", stringsAsFactors=T, dec=",")

#get gene names
pNames <- originalNOT[, 1]

#normalise according to molecular weight
original <- originalNOT[, -c(1,2)]/originalNOT$MW
original <- cbind.data.frame(pNames, original)

##################### CONSIDER ONLY THOSE ROWS THAT CONTAINS NO ZEROES

#find those rows that contains all not-zero values 
rowSub <- apply(original[, -1], 1, function(row) all(row!=0))

#TEST
##find those rows that contains at least a zero
#invRowSub <- apply(original[, -1], 1, function(row) any(row==0))
##then check if it worked: this should be true
#all(!rowSub == invRowSub)

#no rounding no multiplying, keep the values as they are.
filteredData <- cbind.data.frame(toupper(original[rowSub, 1]), original[rowSub, -1])

#colnames the dataset with the normalised values
colnames(filteredData) <- c("proteinName", colnames(original[, -1]))
#assign the proper type to the column "proteinName"
filteredData$proteinName <- levels(factor(unique(filteredData$proteinName)))

#substitute the wrong STRING names  

#PROBLEM: PKM AND PKM2 are the same gene! we must merge the two lines and conside both as PKM
#compute the averaged values
avgVals <- colMeans(filteredData[grep("PKM", filteredData$proteinName), -1])
#remove PKM
filteredNoPKM <- filteredData[-grep("PKM", filteredData$proteinName), ]

#create character vector of protein names
pNames <- append(filteredNoPKM$proteinName, "PKM")

#add the new averaged line
filteredData <- cbind.data.frame(pNames, rbind.data.frame(filteredNoPKM[, -1], avgVals), stringsAsFactors=F)
colnames(filteredData) <- c("proteinName", colnames(original[, -1]))
 
#and then modify two names that are found in STRING
filteredData$proteinName[grep("TUBB5", filteredData$proteinName)] <- "TUBB4A"
filteredData$proteinName[grep("SERPINA3", filteredData$proteinName)] <- factor("AACT", levels="AACT")

#export the new dataset
write.table(filteredData, paste0(patH, "node_table.csv"), sep="\t", col.names=T, row.names=F, quote=F)
write.table(filteredData$proteinName, paste0(patH, "gene_list.txt"), sep="\t", col.names=F, row.names=F, quote=F)

###################### NOW WHAT ######################
#submit lista_geni.txt to STRING to get the network
