library(TCGAbiolinks)
library(maftools)
library(NMF)
library(BSgenome.Hsapiens.UCSC.hg38)



# Downloading somatic mutation data for all patients 

skcm.muse.maf <- GDCquery_Maf("SKCM", pipelines = "muse")

coad.muse.maf <- GDCquery_Maf("COAD", pipelines = "muse")

luad.muse.maf <- GDCquery_Maf("LUAD", pipelines = "muse")

lusc.muse.maf <- GDCquery_Maf("LUSC", pipelines = "muse")




#reading and storing data in variable using maftools

skcm= read.maf(maf = skcm.muse.maf)

coad= read.maf(maf = coad.muse.maf)

luad= read.maf(maf = luad.muse.maf)

lusc= read.maf(maf = lusc.muse.maf)


#Getting overall mutation summarry and Pathway enrichments in the cancer datasets

plotmafSummary(maf = skcm, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

plotmafSummary(maf = coad, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

plotmafSummary(maf = luad, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

plotmafSummary(maf = lusc, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)





#Plotting/seeing enrichment of pathways
#Available pathways: "Cell_Cycle", "Hippo", "MYC", "NOTCH", "NRF2", "PI3K", "RTK-RAS" , "TGF-Beta", "TP53", "WNT"

PlotOncogenicPathways(skcm, pathways = "MYC", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(skcm, pathways = "TP53", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(coad, pathways = "MYC", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(coad, pathways = "TP53", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(lusc, pathways = "MYC", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(lusc, pathways = "TP53", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(luad, pathways = "MYC", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(luad, pathways = "TP53", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)


# Getting trinucleotide matrix by using a reference genome, the downloaded file were based on hg38
#so we will use the BSgenome.Hsapiens.UCSC.hg38 



skcmTNM <- trinucleotideMatrix(maf =skcm, ref_genome = "BSgenome.Hsapiens.UCSC.hg38", prefix = NULL,
                               add = TRUE, ignoreChr = NULL, useSyn = FALSE, fn = NULL)

coadTNM <- trinucleotideMatrix(maf =coad, ref_genome = "BSgenome.Hsapiens.UCSC.hg38", prefix = NULL,
                               add = TRUE, ignoreChr = NULL, useSyn = FALSE, fn = NULL)

luadTNM <- trinucleotideMatrix(maf =luad, ref_genome = "BSgenome.Hsapiens.UCSC.hg38", prefix = NULL,
                               add = TRUE, ignoreChr = NULL, useSyn = FALSE, fn = NULL)

luscTNM <- trinucleotideMatrix(maf =lusc, ref_genome = "BSgenome.Hsapiens.UCSC.hg38", prefix = NULL,
                               add = TRUE, ignoreChr = NULL, useSyn = FALSE, fn = NULL)



#Getting Overall Signatures in the cancers

SigsSKCM = extractSignatures(mat= skcmTNM , n = NULL, nTry = 10, plotBestFitRes = TRUE, parallel = NULL, pConstant = NULL)


plotSignatures(nmfRes = SigsSKCM, contributions = FALSE, color = NULL, patient_order = NULL, font_size = 1.2, show_title = TRUE, axis_lwd = 2, title_size = 0.9, show_barcodes = FALSE, yaxisLim = 0.3)


SigsCOAD = extractSignatures(mat= coadTNM , n = NULL, nTry = 10, plotBestFitRes = TRUE, parallel = NULL, pConstant = NULL)


plotSignatures(nmfRes = SigsCOAD, contributions = FALSE, color = NULL, patient_order = NULL, font_size = 1.2, show_title = TRUE, axis_lwd = 2, title_size = 0.9, show_barcodes = FALSE, yaxisLim = 0.3)


SigsLUAD = extractSignatures(mat= luadTNM , n = NULL, nTry = 10, plotBestFitRes = TRUE, parallel = NULL, pConstant = NULL)


plotSignatures(nmfRes = SigsLUAD, contributions = FALSE, color = NULL, patient_order = NULL, font_size = 1.2, show_title = TRUE, axis_lwd = 2, title_size = 0.9, show_barcodes = FALSE, yaxisLim = 0.3)


SigsLUSC = extractSignatures(mat= luadTNM , n = NULL, nTry = 10, plotBestFitRes = TRUE, parallel = NULL, pConstant = NULL)


plotSignatures(nmfRes = SigsLUSC, contributions = FALSE, color = NULL, patient_order = NULL, font_size = 1.2, show_title = TRUE, axis_lwd = 2, title_size = 0.9, show_barcodes = FALSE, yaxisLim = 0.3)


#Saving trinucleotide matrix with mutation counts to ananlyze using MATLAB code (clustering, seperating by mutation burden)
#see MATLAB  code for determining muutation burden or clustering

skcmtnm <- skcmTNM[["nmf_matrix"]]
write.csv(skcmtnm, file = "SKCMTNM.csv")

coadtnm <- coadTNM[["nmf_matrix"]]
write.csv(coadtnm, file = "coadTNM.csv")

luadtnm <- luadTNM[["nmf_matrix"]]
write.csv(skcmtnm, file = "LUADTNM.csv")

lusctnm <- luscTNM[["nmf_matrix"]]
write.csv(skcmtnm, file = "LUSCTNM.csv")




#Saving contribution scores of each sample
#determining what is dominant signature for each patient was done in MATLAB


skcm_s <- sigsSKCM[["contributions"]]
write.csv(skcm_s, file = "SKCMscores.csv")

coad_s <- sigsCOAD[["contributions"]]
write.csv(coad_s, file = "Coadscores.csv")

luad_s <- sigsLUAD[["contributions"]]
write.csv(luad_s, file = "LUADscores.csv")

lusc_s <- luscTNM[["contributions"]]
write.csv(lusc_s, file = "LUSCscores.csv")


# Creating a Subset MAF file based on a given list patient IDS (see MATLAB  code for determining muutation burden or clustering)
#IDS in the CSV were seperated into 4 columns (one per cluster)
#4 MAFs were created (one per cluster)



dat1 <- read.csv("SKCM_cluster.csv",header=F)$V1
dat1 <- as.character(dat1)
skcmsubsetmaf1 <- subsetMaf(skcm, tsb = dat1, genes = NULL, fields = NULL, query = NULL, mafObj = TRUE, includeSyn = TRUE, isTCGA = FALSE, dropLevels = TRUE, restrictTo = "all")
dat2 <- read.csv("SKCM_cluster.csv",header=F)$V2
dat2 <- as.character(dat2)
skcmsubsetmaf2 <- subsetMaf(skcm, tsb = dat2, genes = NULL, fields = NULL, query = NULL, mafObj = TRUE, includeSyn = TRUE, isTCGA = FALSE, dropLevels = TRUE, restrictTo = "all")
dat3 <- read.csv("SKCM_cluster.csv",header=F)$V3
dat3 <- as.character(dat3)
skcmsubsetmaf3 <- subsetMaf(skcm, tsb = dat3, genes = NULL, fields = NULL, query = NULL, mafObj = TRUE, includeSyn = TRUE, isTCGA = FALSE, dropLevels = TRUE, restrictTo = "all")
dat4 <- read.csv("SKCM_cluster.csv",header=F)$V4
dat4 <- as.character(dat4)
skcmsubsetmaf4 <- subsetMaf(skcm, tsb = dat4, genes = NULL, fields = NULL, query = NULL, mafObj = TRUE, includeSyn = TRUE, isTCGA = FALSE, dropLevels = TRUE, restrictTo = "all")

dat1 <- read.csv("COAD_cluster.csv",header=F)$V1
dat1 <- as.character(dat1)
coadsubsetmaf1 <- subsetMaf(coad, tsb = dat1, genes = NULL, fields = NULL, query = NULL, mafObj = TRUE, includeSyn = TRUE, isTCGA = FALSE, dropLevels = TRUE, restrictTo = "all")
dat2 <- read.csv("COAD_cluster.csv",header=F)$V2
dat2 <- as.character(dat2)
coadsubsetmaf2 <- subsetMaf(coad, tsb = dat2, genes = NULL, fields = NULL, query = NULL, mafObj = TRUE, includeSyn = TRUE, isTCGA = FALSE, dropLevels = TRUE, restrictTo = "all")
dat3 <- read.csv("COAD_cluster.csv",header=F)$V3
dat3 <- as.character(dat3)
coadsubsetmaf3 <- subsetMaf(coad, tsb = dat3, genes = NULL, fields = NULL, query = NULL, mafObj = TRUE, includeSyn = TRUE, isTCGA = FALSE, dropLevels = TRUE, restrictTo = "all")
dat4 <- read.csv("COAD_cluster.csv",header=F)$V4
dat4 <- as.character(dat4)
coadsubsetmaf4 <- subsetMaf(coad, tsb = dat4, genes = NULL, fields = NULL, query = NULL, mafObj = TRUE, includeSyn = TRUE, isTCGA = FALSE, dropLevels = TRUE, restrictTo = "all")


dat1 <- read.csv("LUAD_cluster.csv",header=F)$V1
dat1 <- as.character(dat1)
luadsubsetmaf1 <- subsetMaf(luad, tsb = dat1, genes = NULL, fields = NULL, query = NULL, mafObj = TRUE, includeSyn = TRUE, isTCGA = FALSE, dropLevels = TRUE, restrictTo = "all")
dat2 <- read.csv("LUAD_cluster.csv",header=F)$V2
dat2 <- as.character(dat2)
luadsubsetmaf2 <- subsetMaf(luad, tsb = dat2, genes = NULL, fields = NULL, query = NULL, mafObj = TRUE, includeSyn = TRUE, isTCGA = FALSE, dropLevels = TRUE, restrictTo = "all")
dat3 <- read.csv("LUAD_cluster.csv",header=F)$V3
dat3 <- as.character(dat3)
Luadsubsetmaf3 <- subsetMaf(luad, tsb = dat3, genes = NULL, fields = NULL, query = NULL, mafObj = TRUE, includeSyn = TRUE, isTCGA = FALSE, dropLevels = TRUE, restrictTo = "all")
dat4 <- read.csv("LUAD_cluster.csv",header=F)$V4
dat4 <- as.character(dat4)
luadsubsetmaf4 <- subsetMaf(luad, tsb = dat4, genes = NULL, fields = NULL, query = NULL, mafObj = TRUE, includeSyn = TRUE, isTCGA = FALSE, dropLevels = TRUE, restrictTo = "all")



dat1 <- read.csv("LUSC_cluster.csv",header=F)$V1
dat1 <- as.character(dat1)
luscsubsetmaf1 <- subsetMaf(lusc, tsb = dat1, genes = NULL, fields = NULL, query = NULL, mafObj = TRUE, includeSyn = TRUE, isTCGA = FALSE, dropLevels = TRUE, restrictTo = "all")
dat2 <- read.csv("LUSC_cluster.csv",header=F)$V2
dat2 <- as.character(dat2)
luscsubsetmaf2 <- subsetMaf(lusc, tsb = dat2, genes = NULL, fields = NULL, query = NULL, mafObj = TRUE, includeSyn = TRUE, isTCGA = FALSE, dropLevels = TRUE, restrictTo = "all")
dat3 <- read.csv("LUSC_cluster.csv",header=F)$V3
dat3 <- as.character(dat3)
Luscsubsetmaf3 <- subsetMaf(lusc, tsb = dat3, genes = NULL, fields = NULL, query = NULL, mafObj = TRUE, includeSyn = TRUE, isTCGA = FALSE, dropLevels = TRUE, restrictTo = "all")
dat4 <- read.csv("LUSC_cluster.csv",header=F)$V4
dat4 <- as.character(dat4)
Luscsubsetmaf4 <- subsetMaf(lusc, tsb = dat4, genes = NULL, fields = NULL, query = NULL, mafObj = TRUE, includeSyn = TRUE, isTCGA = FALSE, dropLevels = TRUE, restrictTo = "all")

#Getting TNM for subset MAF files

skcmTNM1 <- trinucleotideMatrix(maf =skcmsubsetmaf1, ref_genome = "BSgenome.Hsapiens.UCSC.hg38", prefix = NULL,
                               add = TRUE, ignoreChr = NULL, useSyn = FALSE, fn = NULL)
skcmTNM2 <- trinucleotideMatrix(maf =skcmsubsetmaf2, ref_genome = "BSgenome.Hsapiens.UCSC.hg38", prefix = NULL,
                                add = TRUE, ignoreChr = NULL, useSyn = FALSE, fn = NULL)
skcmTNM3 <- trinucleotideMatrix(maf =skcmsubsetmaf3, ref_genome = "BSgenome.Hsapiens.UCSC.hg38", prefix = NULL,
                                add = TRUE, ignoreChr = NULL, useSyn = FALSE, fn = NULL)
skcmTNM4 <- trinucleotideMatrix(maf =skcmsubsetmaf4, ref_genome = "BSgenome.Hsapiens.UCSC.hg38", prefix = NULL,
                                add = TRUE, ignoreChr = NULL, useSyn = FALSE, fn = NULL)

coadTNM1 <- trinucleotideMatrix(maf =coadubsetmaf1, ref_genome = "BSgenome.Hsapiens.UCSC.hg38", prefix = NULL,
                                add = TRUE, ignoreChr = NULL, useSyn = FALSE, fn = NULL)
coadTNM2 <- trinucleotideMatrix(maf =coadubsetmaf2, ref_genome = "BSgenome.Hsapiens.UCSC.hg38", prefix = NULL,
                                add = TRUE, ignoreChr = NULL, useSyn = FALSE, fn = NULL)
coadTNM3 <- trinucleotideMatrix(maf =coadubsetmaf3, ref_genome = "BSgenome.Hsapiens.UCSC.hg38", prefix = NULL,
                                add = TRUE, ignoreChr = NULL, useSyn = FALSE, fn = NULL)
coadTNM4 <- trinucleotideMatrix(maf =coadubsetmaf4, ref_genome = "BSgenome.Hsapiens.UCSC.hg38", prefix = NULL,
                                add = TRUE, ignoreChr = NULL, useSyn = FALSE, fn = NULL)


luadTNM1 <- trinucleotideMatrix(maf =luadsetmaf1, ref_genome = "BSgenome.Hsapiens.UCSC.hg38", prefix = NULL,
                                add = TRUE, ignoreChr = NULL, useSyn = FALSE, fn = NULL)

luadTNM2 <- trinucleotideMatrix(maf =luadsetmaf2, ref_genome = "BSgenome.Hsapiens.UCSC.hg38", prefix = NULL,
                                add = TRUE, ignoreChr = NULL, useSyn = FALSE, fn = NULL)

luadTNM3 <- trinucleotideMatrix(maf =luadsetmaf3, ref_genome = "BSgenome.Hsapiens.UCSC.hg38", prefix = NULL,
                                add = TRUE, ignoreChr = NULL, useSyn = FALSE, fn = NULL)

luadTNM4 <- trinucleotideMatrix(maf =luadsetmaf4, ref_genome = "BSgenome.Hsapiens.UCSC.hg38", prefix = NULL,
                                add = TRUE, ignoreChr = NULL, useSyn = FALSE, fn = NULL)


#Comparing mutation load in subsets 

skcmsubsetmaf1.mutload = tcgaCompare(maf = skcmsubsetmaf1, cohortName = 'subset1')
skcmsubsetmaf2.mutload = tcgaCompare(maf = skcmsubsetmaf2, cohortName = 'subset2')
skcmsubsetmaf3.mutload = tcgaCompare(maf = skcmsubsetmaf3, cohortName = 'subset3')
skcmsubsetmaf1.mutload = tcgaCompare(maf = skcmsubsetmaf4, cohortName = 'subset4')

coadsubsetmaf1.mutload = tcgaCompare(maf = coadsubsetmaf1, cohortName = 'subset1')
coadsubsetmaf2.mutload = tcgaCompare(maf = coadsubsetmaf2, cohortName = 'subset2')
coadsubsetmaf3.mutload = tcgaCompare(maf = coadsubsetmaf3, cohortName = 'subset3')
coadsubsetmaf1.mutload = tcgaCompare(maf = coadsubsetmaf4, cohortName = 'subset4')

luscsubsetmaf1.mutload = tcgaCompare(maf = luscsubsetmaf1, cohortName = 'subset1')
luscsubsetmaf2.mutload = tcgaCompare(maf = luscsubsetmaf2, cohortName = 'subset2')
luscsubsetmaf3.mutload = tcgaCompare(maf = luscsubsetmaf3, cohortName = 'subset3')
luscsubsetmaf1.mutload = tcgaCompare(maf = luscsubsetmaf4, cohortName = 'subset4')

luadsubsetmaf1.mutload = tcgaCompare(maf = luadsubsetmaf1, cohortName = 'subset1')
luadsubsetmaf2.mutload = tcgaCompare(maf = luadsubsetmaf2, cohortName = 'subset2')
luadsubsetmaf3.mutload = tcgaCompare(maf = luadsubsetmaf3, cohortName = 'subset3')
luadsubsetmaf1.mutload = tcgaCompare(maf = luadsubsetmaf4, cohortName = 'subset4')

#To find differences in TP53 and MYC pathway 


#SKCM 

PlotOncogenicPathways(skcmsubsetmaf1, pathways = "MYC", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(skcmsubsetmaf1, pathways = "TP53", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)


PlotOncogenicPathways(skcmsubsetmaf2, pathways = "MYC", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(skcmsubsetmaf2, pathways = "TP53", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(skcmsubsetmaf3, pathways = "MYC", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(skcmsubsetmaf3, pathways = "TP53", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(skcmsubsetmaf4, pathways = "MYC", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(skcmsubsetmaf4, pathways = "TP53", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

#COAD

PlotOncogenicPathways(coadsubsetmaf1, pathways = "MYC", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(coadsubsetmaf1, pathways = "TP53", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(coadsubsetmaf2, pathways = "MYC", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(coadsubsetmaf2, pathways = "TP53", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(coadsubsetmaf3, pathways = "MYC", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(coadsubsetmaf4, pathways = "TP53", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(coadsubsetmaf4, pathways = "MYC", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(coadsubsetmaf3, pathways = "TP53", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

#LUAD

PlotOncogenicPathways(luadsubsetmaf1, pathways = "MYC", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(luadsubsetmaf1, pathways = "TP53", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(luadsubsetmaf2, pathways = "MYC", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(luadsubsetmaf2, pathways = "TP53", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(luadsubsetmaf3, pathways = "MYC", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(luadsubsetmaf4, pathways = "TP53", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(luadsubsetmaf4, pathways = "MYC", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

#LUSC

PlotOncogenicPathways(luscsubsetmaf1, pathways = "MYC", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(luscubsetmaf1, pathways = "TP53", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(luscsubsetmaf2, pathways = "MYC", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(luscsubsetmaf2, pathways = "TP53", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(luscsubsetmaf3, pathways = "MYC", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(luscsubsetmaf4, pathways = "TP53", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)

PlotOncogenicPathways(luscsubsetmaf4, pathways = "MYC", fullPathway = FALSE,
                      removeNonMutated = TRUE, tsgCol = "red", ogCol = "royalblue",
                      fontSize = 0.6, showTumorSampleBarcodes = FALSE,
                      SampleNamefontSize = 0.6)





#Getting Mutation Signatures from subsetMAF files
#this step  requies obtaining Tri nucleotide matrices using the code above


#SKCM

SigSKCM1 = extractSignatures(mat= skcmTNM1, n = NULL, nTry = 10, plotBestFitRes = TRUE, parallel = 'P4', pConstant = NULL)
plotSignatures(nmfRes = SigSKCM, contributions = FALSE, color = NULL, patient_order = NULL, font_size = 1.2, show_title = TRUE, axis_lwd = 2, title_size = 0.9, show_barcodes = FALSE, yaxisLim = 0.3)

SigSKCM2 = extractSignatures(mat= skcmTNM2, n = NULL, nTry = 10, plotBestFitRes = TRUE, parallel = 'P4', pConstant = NULL)
plotSignatures(nmfRes = SigSKCM2, contributions = FALSE, color = NULL, patient_order = NULL, font_size = 1.2, show_title = TRUE, axis_lwd = 2, title_size = 0.9, show_barcodes = FALSE, yaxisLim = 0.3)

SigSKCM3 = extractSignatures(mat= skcmTNM3, n = NULL, nTry = 10, plotBestFitRes = TRUE, parallel = 'P4', pConstant = NULL)
plotSignatures(nmfRes = SigSKCM3, contributions = FALSE, color = NULL, patient_order = NULL, font_size = 1.2, show_title = TRUE, axis_lwd = 2, title_size = 0.9, show_barcodes = FALSE, yaxisLim = 0.3)

SigSKCM4 = extractSignatures(mat= skcmTNM4, n = NULL, nTry = 10, plotBestFitRes = TRUE, parallel = 'P4', pConstant = NULL)
plotSignatures(nmfRes = SigSKCM4, contributions = FALSE, color = NULL, patient_order = NULL, font_size = 1.2, show_title = TRUE, axis_lwd = 2, title_size = 0.9, show_barcodes = FALSE, yaxisLim = 0.3)

#COAD
SigCOAD1 = extractSignatures(mat= coadTNM1, n = NULL, nTry = 10, plotBestFitRes = TRUE, parallel = 'P4', pConstant = NULL)
plotSignatures(nmfRes = SigCOAD1, contributions = FALSE, color = NULL, patient_order = NULL, font_size = 1.2, show_title = TRUE, axis_lwd = 2, title_size = 0.9, show_barcodes = FALSE, yaxisLim = 0.3)

SigCOAD2 = extractSignatures(mat= coadTNM2, n = NULL, nTry = 10, plotBestFitRes = TRUE, parallel = 'P4', pConstant = NULL)
plotSignatures(nmfRes = SigCOAD2, contributions = FALSE, color = NULL, patient_order = NULL, font_size = 1.2, show_title = TRUE, axis_lwd = 2, title_size = 0.9, show_barcodes = FALSE, yaxisLim = 0.3)

SigCOAD3 = extractSignatures(mat= coadTNM3, n = NULL, nTry = 10, plotBestFitRes = TRUE, parallel = 'P4', pConstant = NULL)
plotSignatures(nmfRes = SigCOAD3, contributions = FALSE, color = NULL, patient_order = NULL, font_size = 1.2, show_title = TRUE, axis_lwd = 2, title_size = 0.9, show_barcodes = FALSE, yaxisLim = 0.3)

SigCOAD4 = extractSignatures(mat= coadTNM4, n = NULL, nTry = 10, plotBestFitRes = TRUE, parallel = 'P4', pConstant = NULL)
plotSignatures(nmfRes = SigCOAD4, contributions = FALSE, color = NULL, patient_order = NULL, font_size = 1.2, show_title = TRUE, axis_lwd = 2, title_size = 0.9, show_barcodes = FALSE, yaxisLim = 0.3)

#LUSC

SigLUSC1 = extractSignatures(mat= luscTNM1, n = NULL, nTry = 10, plotBestFitRes = TRUE, parallel = 'P4', pConstant = NULL)
plotSignatures(nmfRes = SigLUSC1, contributions = FALSE, color = NULL, patient_order = NULL, font_size = 1.2, show_title = TRUE, axis_lwd = 2, title_size = 0.9, show_barcodes = FALSE, yaxisLim = 0.3)

SigLUSC2 = extractSignatures(mat= luscTNM2, n = NULL, nTry = 10, plotBestFitRes = TRUE, parallel = 'P4', pConstant = NULL)
plotSignatures(nmfRes = SigLUSC2, contributions = FALSE, color = NULL, patient_order = NULL, font_size = 1.2, show_title = TRUE, axis_lwd = 2, title_size = 0.9, show_barcodes = FALSE, yaxisLim = 0.3)

SigLUSC3 = extractSignatures(mat= luscTNM3, n = NULL, nTry = 10, plotBestFitRes = TRUE, parallel = 'P4', pConstant = NULL)
plotSignatures(nmfRes = SigLUSC3, contributions = FALSE, color = NULL, patient_order = NULL, font_size = 1.2, show_title = TRUE, axis_lwd = 2, title_size = 0.9, show_barcodes = FALSE, yaxisLim = 0.3)

SigLUSC4 = extractSignatures(mat= luscTNM4, n = NULL, nTry = 10, plotBestFitRes = TRUE, parallel = 'P4', pConstant = NULL)
plotSignatures(nmfRes = SigLUSC4, contributions = FALSE, color = NULL, patient_order = NULL, font_size = 1.2, show_title = TRUE, axis_lwd = 2, title_size = 0.9, show_barcodes = FALSE, yaxisLim = 0.3)

#LUAD
SigLUAD1 = extractSignatures(mat= luadTNM1, n = NULL, nTry = 10, plotBestFitRes = TRUE, parallel = 'P4', pConstant = NULL)
plotSignatures(nmfRes = SigLUAD1, contributions = FALSE, color = NULL, patient_order = NULL, font_size = 1.2, show_title = TRUE, axis_lwd = 2, title_size = 0.9, show_barcodes = FALSE, yaxisLim = 0.3)

SigLUAD2 = extractSignatures(mat= luadTNM2, n = NULL, nTry = 10, plotBestFitRes = TRUE, parallel = 'P4', pConstant = NULL)
plotSignatures(nmfRes = SigLUAD2, contributions = FALSE, color = NULL, patient_order = NULL, font_size = 1.2, show_title = TRUE, axis_lwd = 2, title_size = 0.9, show_barcodes = FALSE, yaxisLim = 0.3)

SigLUAD3 = extractSignatures(mat= luadTNM3, n = NULL, nTry = 10, plotBestFitRes = TRUE, parallel = 'P4', pConstant = NULL)
plotSignatures(nmfRes = SigLUAD3, contributions = FALSE, color = NULL, patient_order = NULL, font_size = 1.2, show_title = TRUE, axis_lwd = 2, title_size = 0.9, show_barcodes = FALSE, yaxisLim = 0.3)

SigLUAD4 = extractSignatures(mat= luadTNM4, n = NULL, nTry = 10, plotBestFitRes = TRUE, parallel = 'P4', pConstant = NULL)
plotSignatures(nmfRes = SigLUAD4, contributions = FALSE, color = NULL, patient_order = NULL, font_size = 1.2, show_title = TRUE, axis_lwd = 2, title_size = 0.9, show_barcodes = FALSE, yaxisLim = 0.3)








