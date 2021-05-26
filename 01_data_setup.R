#############
# GabrieleT #
#   2021    #
#############

# set seed
set.seed(131)

# set paths
patH <- "~/WNNets/"

# read data
originalNOT <- read.table(paste0(patH, "SpC_data_no_normalisation.csv"), header=T, sep="\t", stringsAsFactors=T, dec=",")

# remove F samples
# colnames(originalNOT[, c(15:20, 33:38)])
originalNOFs <- originalNOT[, -c(15:20, 33:38)]

###################### THIS IS DONE TO OBTAIN THE BEFORE-NORMALISATION DATA TABLE ONLY
# create the table before finding non-zero values and before the normalisation
beforeTmp <- data.frame(as.character(originalNOFs$gene.sus), originalNOFs[, -1])
colnames(beforeTmp) <- c("proteins", colnames(originalNOFs[, -1]))
# remove PKM2 because of issues with STRING
beforeTmp <- beforeTmp[which(beforeTmp$proteins != "PKM2"), ]
# remove TUBB4A because it is recognised as Tubb5 from STRING
beforeTmp <- beforeTmp[which(beforeTmp$proteins != "TUBB4A"), ]
# remove TUBA1B which is a duplicated (with SpC always at zero) for Tuba1B
beforeTmp <- beforeTmp[which(beforeTmp$proteins != "TUBA1B"), ]
# now, modify those names that are found by STRING under other names
# to be consistent with the network that will be downloaded
beforeTmp$proteins[grep("Tubb5", beforeTmp$proteins)] <- "TUBB4A"
beforeTmp$proteins[grep("SERPINA3", beforeTmp$proteins)] <- "AACT"
beforeTmp$proteins[grep("ENSG00000125166", beforeTmp$proteins)] <- "GOT2"
beforeTmp$proteins[grep("COX2", beforeTmp$proteins)] <- "MT-CO2"
beforeTmp$proteins[grep("ATPase 6", beforeTmp$proteins)] <- "MT-ATP6"
# be careful with IGVH, because there are many. search for IGVH$ to be sure
beforeTmp$proteins[grep("IGHV$", beforeTmp$proteins)] <- "IGHV4-38-2"
# plus those proteins that are lowercase
beforeTmp$proteins[grep("Tuba1b", beforeTmp$proteins)] <- "TUBA1B"

# and export it
write.table(beforeTmp, "before_temporary.csv", row.names=F, quote=F, sep="\t")

######################################################################################

# get gene names
pNames <- originalNOFs[, 1]

# normalise according to molecular weight without considering protein names and MW itself, i.e. columns 1 and 2
#  colnames(originalNOFs[, c(1,2)])
original <- originalNOFs[, -c(1,2)]/originalNOFs$MW
# re-attach proteins names to the normalised data set
original <- cbind.data.frame(pNames, original)

# find those rows that contains non-zero values only. don't consider the first column (protein names)
# colnames(original[, -1])
rowSub <- apply(original[, -1], 1, function(row) all(row!=0))

# TEST
# #find those rows that contains at least a zero
# invRowSub <- apply(original[, -1], 1, function(row) any(row==0))
# #then check if it worked: this should be true
# all(!rowSub == invRowSub)

# get only non-zero rows from the original, normalised matrix. these values will be used to weigh the master nodes
nonZeros <- cbind.data.frame(toupper(original[rowSub, 1]), original[rowSub, -1])
# colnames the dataset with the normalised values
colnames(nonZeros) <- c("proteinName", colnames(original[, -1]))
# assign the proper type to the column "proteinName"
nonZeros$proteinName <- as.character(nonZeros$proteinName)

# PROBLEM: PKM AND PKM2 are the same gene, i.e. the two names are synonyms
# however PKM2 is synonym and UniProt has info for PKM only

# remove PKM2 which is not found by STRING
cleanNames <- nonZeros[which(nonZeros$proteinName != "PKM2"), ]

# now, modify those names that are found by STRING under a synonym

# to be consistent with the network that will be downloaded
cleanNames$proteinName[grep("TUBB5", cleanNames$proteinName)] <- "TUBB4A"
cleanNames$proteinName[grep("SERPINA3", cleanNames$proteinName)] <- "AACT"
cleanNames$proteinName[grep("ENSG00000125166", cleanNames$proteinName)] <- "GOT2"
cleanNames$proteinName[grep("COX2", cleanNames$proteinName)] <- "MT-CO2"
cleanNames$proteinName[grep("ATPASE 6", cleanNames$proteinName)] <- "MT-ATP6"
# be careful with IGVH, because there are many. search for IGVH$ to be sure
cleanNames$proteinName[grep("IGHV$", cleanNames$proteinName)] <- "IGHV4-38-2"

# export the new dataset
write.table(cleanNames, paste0(patH, "node_table.csv"), sep="\t", col.names=T, row.names=F, quote=F)
write.table(cleanNames$proteinName, paste0(patH, "shared_proteins.txt"), sep="\t", col.names=F, row.names=F, quote=F)

# at this point
# submit shared_proteins.txt to STRING to get the network. Use these options:
# experimental interactions only
# threshold = 0.400 (STRING default threshold)
