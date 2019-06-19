# WGCNA pCR Basal-like

#######################################################################################
# SETTING UP R SESSION
#######################################################################################

workingDir = '.'
setwd(workingDir)

set.seed(12345)
options(stringsAsFactors = FALSE)

# Load required libraries
source("RCode/0_loadLibraries.R")
loadpkg("NanoStringQCPro")
loadpkg("dplyr")
loadpkg("WGCNA")
loadpkg("biomaRt")
loadpkg("clusterProfiler")
loadpkg("linkcomm")
loadpkg('ggplot2')
loadpkg("igraph")

allowWGCNAThreads() # multi-threading within WGCNA. RStudio. 
# For r, Rscript use enableWGCNAThreads()
# enableWGCNAThreads()

dir.create(path = 'Data/RData')
dir.create(path = 'Results')
dir.create(path = 'Results/WGCNAplots')
dir.create(path = 'Results/linkcomm')
dir.create(path = 'Results/modules_enrichments')
dir.create(path = 'Results/expanded_modules')
dir.create(path = 'Results/expanded_modules_enrichments')
dir.create(path = 'Results/expanded_modules_networks')

#######################################################################################
# DATA INPUT
#######################################################################################

### Load RCC files

DataDir <- 'Data'
rccDir <- file.path(DataDir, "RCCData")
## Process file names: <sample>_<key1>_<key2>_<subtype>.RCC
## We have 204 files. There are 6 files without <subtype>
## Key 1 contains pre/post info. 0 = PRE 
## Key 2 contains pCR info 0 = no pCR; 1 = pCR
### We want basal-like pre NAC samples 
pre_samples = data.frame(fname=list.files(rccDir, pattern="*.RCC")) %>%
  tidyr::separate(fname,c("sample","key1","key2","subtype"),sep="_",remove=F) %>%
  tidyr::separate(subtype, c("subtype", "ext"), sep="\\.") %>%
  dplyr::filter(subtype=="BL") %>%
  dplyr::select (fname,sample, key1,key2) %>%
  dplyr::filter(key1 == 0)


## Files without subtype has key2 = 'key2.RCC' so we have to remove '.RCC' to avoid further errors
for (i in grep(pattern = '\\.' ,pre_samples$key2)) {
  pre_samples$key2[i] <- unlist(strsplit(pre_samples$key2[i], '\\.'))[1]
}

rccSet <- newRccSet(
  rccFiles               = file.path(rccDir, pre_samples$fname)
  #,rccCollectorToolExport = file.path(exampleDataDir, "nSolver", "RCC_collector_tool_export.csv")
  ,rlf                    = file.path(DataDir,"RLFData/NS_CancerPath_C2535.rlf") #rlf file in RLFData subdir
  #,cdrDesignData          = file.path(exampleDataDir, "CDR", "CDR-DesignData.csv")
  #,extraPdata             = file.path(exampleDataDir, "extraPdata", "SampleType.txt")
  #,blankLabel             = "blank"
  ,experimentData.name    = "PRE_pCR_BL"
  ,experimentData.lab     = "LBMC"
  ,experimentData.contact = "amoyag@uma.es"
  ,experimentData.title   = "WGCNA correlation with pCR in basal-like samples"
  ,experimentData.abstract= "Gene co-expression network to find genes associated with pCR in basal-like pre-NAC samples."
)

### Expression data preprocessing and cleaning

## Preprocessing and normalization
norm_rccSet <- preprocRccSet(rccSet = rccSet, doPosCtrlNorm = T, doContentNorm =T, normMethod = "housekeeping", bgReference = "negative")

## QC report
qc_rccSet <- makeQCReport(norm_rccSet, "QC_report"
                          , outputDir = paste(getwd(),'Results',sep='/'))

## Expression matrix with normalised data
ncounts <- assayData(qc_rccSet)$normData

## Apply QC filter

keepFiles <- dplyr::filter(qc_rccSet@phenoData@data, TechnicalFlags == FALSE)$FileName
ncounts <- ncounts[, keepFiles]

## Remove control genes (positive, negative, housekeeping)
rm_control <- (qc_rccSet@featureData@data %>% 
                 tibble::rownames_to_column(var = "FullName") %>%
                 dplyr::select(FullName, CodeClass) %>%
                 dplyr::filter(!grepl("Pos",CodeClass)) %>%
                 dplyr::filter(!grepl("Neg", CodeClass)) %>%
                 dplyr::filter(!grepl("Spike", CodeClass)) %>%
                 dplyr::filter(!grepl("Ligati", CodeClass)) %>%
                 dplyr::filter(!grepl("Housekee", CodeClass)))$FullName

ncounts <- ncounts[rm_control, ]

## Change columns from file name to sample name
pre_samples <- filter(pre_samples, fname %in% colnames(ncounts))
colnames(ncounts) <- pre_samples$sample

## Transpose the expression data
datExpr0 = as.data.frame(t(ncounts))

## Check samples for missing expression values
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

## Clustering to detect outliers
sampleTree = hclust(dist(datExpr0), method = "average");

## Plot the sample tree

pdf(file = "Results/WGCNAplots/sampleClustering.pdf");
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 65, col = "red");
dev.off()
#################################

## Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 65, minSize = 10)
table(clust)

keepSamples = (clust==1)
datExpr = datExpr0[keepSamples,]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

### Load clinical data

pCRDat <- pre_samples %>% dplyr::select(sample,key2)

names(pCRDat) <- c("sample", "pCR")
pCRDat = pCRDat[keepSamples,]

pCR.to.chemores <- function(pcr){
  if(pcr == 0){return(1)}
  if(pcr == 1){return(0)}
}

chemores <- pCRDat
chemores$ChemoRes<- sapply(chemores$pCR, pCR.to.chemores)
chemores <- dplyr::select(chemores, sample, ChemoRes)
chemores$ChemoRes <- as.numeric(chemores$ChemoRes)
datTrait <- dplyr::select(chemores, ChemoRes)
rownames(datTrait) <- chemores$sample

### Plot a cluster of the samples and the pCR values as color

# White means low (pCR =0) and red means pCR =1
sampleTree2 = hclust(dist(datExpr), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTrait, signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
pdf("Results/WGCNAplots/dendogram.pdf")
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTrait),
                    main = "Sample dendrogram chemoresistance heatmap")
dev.off()
save(datExpr, datTrait, file = "Data/RData/chemoresNAC_wgcna-01-dataInput.RData")


#######################################################################################
# NETWORK CONSTRUCTION AND MODULE DETECTION
#######################################################################################

### Soft-thresholding power selection

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2));
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5);

# Plot the results:
#sizeGrWindow(9, 5);
par(mfrow = c(1,2));
pdf("Results/WGCNAplots/softpower_threshold.pdf")
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# This line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=0,col='red')
dev.off()


### Build a Topological Overlaping Matrix from adjacency matrix  

softPower = 4
adjacency = adjacency(datExpr, power = softPower);

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM


### Gene clustering on TOM-based dissimilarity  

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");

# Plot the resulting clustering tree (dendrogram)

plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

### Module identification  

# Since I get few modules, I'll try with smaller modules:
minModuleSize = 20;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


### Gene clustering modules  

# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath

pdf("Results/WGCNAplots/gene_clustering_modules.pdf")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
  

### Modules clustering  

# The Dynamic Tree Cut may identify modules whose expression profiles are very similar. 
# It may be prudent to merge such modules since their genes are highly co-expressed. 
# To quantify co-expression similarity of entire modules, we calculate their eigengenes and cluster them on their correlation:  

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

# Plot the result
pdf("Results/WGCNAplots/modules-merge.pdf")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.25 # Corresponds with a correlation of 0.75
# Plot the cut line into the dendrogram
abline(h = MEDissThres, col = "red")
dev.off()

### Merge modules with highly co-expressed genes

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors = merge$colors;

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


pdf(file = "Results/WGCNAplots/gene_clustering_modules-merge.pdf")
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

### Rename, label and save  

# Rename to moduleColors
moduleColors = mergedColors

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "Data/RData/chemoresNAC_wgcna-02-networkConstruction.RData")


#######################################################################################
# RELATE MODULES TO CLINICAL CHEMORESISTANCE
#######################################################################################

# Load data (datExpr and datTrait) if already saved:
#load("RData/chemoresNAC_wgcna-01-dataInput.RData")
#load("RData/chemoresNAC_wgcna-02-networkConstruction.RData")

### Quantifing module-trait associations  

# In this analysis we would like to identify modules that are significantly associated with the measured clinical traits. 
# Since we already have a summary profile (eigengene) for each module, we simply correlate eigengenes with external traits and look for the most significant associations:

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
#moduleTraitCor = cor(MEs, datTrait, use = "p");
moduleTraitCor = bicor(MEs, datTrait, use="p",robustY=FALSE, maxPOutliers=1)
# moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
moduleTraitPvalue = corPvalueFisher(moduleTraitCor, nSamples, twoSided = F)
#bicorAndPvalue(x = MEs, y = datTrait, use = "p",alternative = "g")


# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

par(mar = c(6.5, 8, 3, 3));
# Display the correlation values within a heatmap plot
pdf(file = "Results/WGCNAplots/module-trait.pdf", wi = 8, he = 10)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTrait),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.0,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()


# Module eigengene yellow shows significant anti-correlation with chemoresistance (-0.28 pval=0.0495201).



### Gene relationship to trait and important modules:Gene significance and Module Membership  

# We quantify associations of individual genes with our trait of interest (pCR) by defining Gene Significance GS as (the absolute value of) the correlation between the gene and the trait. 
# For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile. 
# This allows us to quantify the similarity of all genes on the array to every module.  

# Define variable weight containing the weight column of datTrait
ChemoRes = as.data.frame(datTrait$ChemoRes);
names(ChemoRes) = "ChemoRes"

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, ChemoRes, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(ChemoRes), sep="");
names(GSPvalue) = paste("p.GS.", names(ChemoRes), sep="");


### Intramodular analysis: identfying genes with GS and MM  

# Using the GS and MM measures, we can identify genes that have a high significance for pCR as well as high module membership in interesting modules.  

MM_vs_GS <- function(module, modNames, moduleColors){
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  pdf(file = 'Results/WGCNAplots/',module,'_MMvsGS.pdf')
  par(mfrow = c(1,1));
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for Chemoresistance",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module, pch=19)
  abline(h=0.2, col = "red")
  abline(v=0.7, col = "red")
  dev.off()
  return(NULL)
}

MM_vs_GS(module = "yellow", modNames, moduleColors)


### Summary output of network analysis result  

# We have found modules with high association with our trait of interest, and have identified their central players by the Module Membership measure. 
# We now merge this statistical information with gene annotation and perform an enrichment analysis. 
# After that we write out several files that summarizes the most important results.  

# Load gene annotation info from NS_CancerPath_C2535.rlf file
annot <- readRlf("Data/RLFData/NS_CancerPath_C2535.rlf")

get.probe <- function(complexid){
  parts = strsplit(complexid, split = "_")[[1]]
  paste("NM", parts[4], sep="_")
}

probes = unname(sapply(names(datExpr), get.probe))
probes2annot <-  match(probes, annot$Accession)
sum(is.na(probes2annot)) # number or probes without annotation. Should be 0

annot.geneid <- bitr(annot$GeneName, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
names(annot.geneid) <- c("GeneName", "ENTREZID")
annot <- merge(annot, annot.geneid, by="GeneName")
annot <- annot %>% dplyr::select(Accession, GeneName, ENTREZID, ProbeID)

## Define the universe backgroud for the enrichments
universe = unique(annot$ENTREZID) # Note that we are not using the universe in the enrichments

# Create the starting data frame
geneInfo0 = data.frame(Accession = probes,
                       geneSymbol = annot$GeneName[probes2annot],
                       entrezID = annot$ENTREZID[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, ChemoRes, use = "p")));

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)){
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(-abs(geneInfo0$GS.ChemoRes),geneInfo0$moduleColor);
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "Results/geneInfo.csv", quote = F,row.names = F)

## Write interesting modules gene info in a file
intModules <- c("yellow");
module_gene_info <- function(module, geneInfo){
  
  if(module=="yellow"){
    ChemoRes.genes <- geneInfo[which(geneInfo$moduleColor==module), c(1,2,3,4,5,6,7,8)]
    ChemoRes.genes <- ChemoRes.genes %>% filter(abs(MM.yellow) >= 0.70 & abs(GS.ChemoRes) >= 0.1)
  }
  write.csv(ChemoRes.genes, file = paste("Results/",module,"-genesInfo.csv", sep=""), quote = F,row.names = F)
  return(NULL)
}
lapply(intModules, module_gene_info, geneInfo)

#######################################################################################
# WGCNA MODULES ENRICHMENT ANALYSIS
#######################################################################################

### Loading data   

# Load genes info
geneInfo <- read.csv(file = "Results/geneInfo.csv", header = T)

# Separate modules
yellow_module <- geneInfo %>% 
  dplyr::filter(abs(MM.yellow) >=0.7 & abs(GS.ChemoRes) >= 0.1 ) %>%
  dplyr::select(Accession, geneSymbol, entrezID, GS.ChemoRes, p.GS.ChemoRes, MM.yellow, p.MM.yellow)

### Modules enrichment analysis

source("RCode/enrichment.R")
module_enrichment_analysis <- function(geneInfo, module, GS, MM){
  
  MM.module <- paste("MM.",module,sep="")
  
  res <- geneInfo %>% 
    dplyr::filter(moduleColor == module & !is.na(entrezID) & abs(GS.ChemoRes) >= GS & abs(geneInfo[,MM.module]) >= MM) %>%
    dplyr::select(Accession, geneSymbol, entrezID, matches("GS"), matches("MM")) %>% 
    dplyr::mutate(entrezID = as.character(entrezID))
  
  genes <- res$entrezID
  
  enrich_rs <- enrich_cp(genes, comparison = "ChemoRes")
  enrich_summary <- enrich_rs$summary %>% arrange(qvalue)
  enrich_summary <- convert_enriched_ids(enrich_summary, entrezsymbol = entrezsymbol) %>% filter(qvalue <=0.05) %>% arrange(-Count)
  
  write.csv(res, file = paste("Results/modules_enrichments/",module,"_module_genes.csv", sep="") )
  write.csv(enrich_summary, file = paste("Results/modules_enrichments/",module,"_module_enrichment.csv", sep=""))
}

module_enrichment_analysis(geneInfo, "yellow", GS=0.1, MM=0.7)

### Expanding modules with DIAMOnD and building networks

# Before using DIAMOnD, we have to download several high confidence human interactomes in order to run the algorithm. With the interactome and the seed genes we can expand our initial modules and construct bigger networks from them. 
# We will use:  
  
#  * The DIAMOnD paper interactome provided by Ghiassian et al.(2015)[1].  
#  * STRINGdb human interactome with a score threshold > 0.7.  
#  * iRef human interactome  
#  * A metabolism centered interactome  

# The first one doesn't need any processing before run DIAMOnD. Nevertheless, we have to process the others.  
# To get the STRINGdb interactome we download the whole interactome of Homo sapiens from https://string-db.org/cgi/download.pl?sessionId=ED2PAnsg3O9r&species_text=Homo+sapiens. After that we apply a filter to preserve only those interactions with a score higher than 0.7 (high confidence interactome) and we map the protein ids to its appropriate gene id and we keep only the gene ids and the score.  
# To get the iRef interactome we have tu use the iRefR package to load a human interactome from iRef using iRefR::get_irefindex(tax_id = "9606", iref_version = "13.0"). The loaded table has different charasteristics. We keep iROGid´s of the interactors and the confidence of the interaction.  
# On the one hand, iRefIndex guarantees the non-redundancy of protein information by assigning a different protein identifier (called Redundant Object Group, ROG) to every different protein sequence.  
# On the other hand, the column 'confidence' refers to three confidence scores: lpr, hpr and np. We are interested only in np which is the total number of unique PMIDs(PubMeds identifiers) used to support the interaction described in the row. (More information is available on http://irefindex.org/wiki/index.php?title=README_MITAB2.6_for_iRefIndex#Column_number:_15_.28confidence.29)  
# We map the protein iROGids to gene ids using the function iRefR::create_id_conversion_table() and we keep only those interactions with an np confidence score higher than ??????.
# Finally, we download the metabolism centered interactome from https://bioinfo.uth.edu/ccmGDB/download/2072genes_PathwayCommons.txt?csrt=5379671582968607161. This interactome has Pathway Commons format (see http://www.pathwaycommons.org/pc/sif_interaction_rules.do). We have to process it to keep only the INTERACTS_WITH relationships and to map the symbol ids to gene ids.  
# Once we have the interactomes ready we can run DIAMOnD.

source('RCode/expansion_functions.R')
# We have to set the path to our python source
# To see your current python source type:
# system('python --version')
pybin <- "/anaconda3/bin/python"

# We set the path to DIAMOnD
diamond_path <- "DIAMOnD.py"

expand_module(pybin, diamond_path, "yellow", network_path = "Interactomes/Ghiassian_interactome.tsv")
expand_module(pybin, diamond_path, "yellow", network_path = "Interactomes/Metabolism-centered_interactome.tsv")
expand_module(pybin, diamond_path, "yellow", network_path = "Interactomes/STRINGdb_interactome.tsv")

lapply(list.files(path = "Results/expanded_modules", full.names = T), network_building)

### Expanded modules enrichments
lapply(list.files(path = "Results/expanded_modules", full.names = T), expanded_module_enrichment)

#######################################################################################
# LINK COMMUNITIES ANALYSIS SUBNETWORKS PROCCESSING
#######################################################################################

source('RCode/process_subnetworks.R');
subnetworks = list.files(path = 'Results/expanded_modules_networks', full.names = F);
lapply(subnetworks, get_gml)

## References  
# 1. Ghiassian SD, Menche J, Barabási A-L (2015) A DIseAse MOdule Detection (DIAMOnD) Algorithm Derived from a Systematic Analysis of Connectivity Patterns of Disease Proteins in the Human Interactome. PLoS Comput Biol 11(4): e1004120. https://doi.org/10.1371/journal.pcbi.1004120

