######### --Libraries--#########
setwd("/vol08/ngs/P51/Aging/Aging02/Cornelius_analysis/de_work")
library(stringr)
library(limma)
library(edgeR)
library(PCAtools)
library(plyr)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(WebGestaltR)
source("/vol08/ngs/P51/Aging/Aging02/Cornelius_analysis/de_work/heatmap3LW_function.R")

# Helper file
rhesus2human <- read.csv(
  file = "/vol08/ngs/P51/Aging/Aging02/Cornelius_analysis/de_work/Macaca_mulatta_Mmul_10.109/rhesus2human109v2.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)


# --Read in target files
message("STATUS: Load tables")
cm <- read.table("./count_matrix.txt", header = TRUE, sep = "\t", row.names = 1, as.is = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
target <- read.csv("./Aging02_Target_File_Slim.csv", sep = ",", row.names = 1, as.is = TRUE, check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)


# Rename samples
newsampleIDs <- c()
for (i in colnames(cm)) {
  i <- str_remove(i, "_RNA\\d+_Lib\\d+\\S*$")
  i <- str_replace_all(i, "-", "_")
  newsampleIDs <- c(newsampleIDs, i)
}
colnames(cm) <- newsampleIDs

# --generate figure of all counts
png("de_intensities_raw_counts.png", res = 100)
par(xpd = TRUE)
if (length(rownames(target)) > 10) {
  plotDensities(log2(cm + 0.1), legend = FALSE)
} else {
  plotDensities(log2(cm + 0.1),
                legend = "topright",
                inset = c(-0.2, 0), levels(rownames(target))
  )
}
dev.off()


## Generate model matrix ##

target$bioreps <- paste(target$Time_Point, target$Aged, sep = "_")
target$bioreps <- str_replace_all(target$bioreps, "-", "_")
bioreps <- factor(target$bioreps)
mm <- model.matrix(~ 0 + bioreps)

# id <- factor(target[, "Animal_ID"])
# tp <- factor(target[, "Time_Point"])
# age <- factor(target[, "Aged"])
# tr <- factor(target[, "Treatement"])
# mm <- model.matrix(~ 0 + id:tp:age)

rownames(mm) <- rownames(target)
colnames(mm) <- make.names(colnames(mm))
mm <- mm[, colnames(mm)[order(tolower(colnames(mm[, ])))]]
mm <- mm[, colSums(mm) > 0]

excludeAll <- nonEstimable(mm)
if (length(excludeAll) > 0) {
  message("WARNING: These samples are nonEstimatable, design matrix ", excludeAll)
}

if ("ti" %in% excludeAll) {
  return("interactions term non estimable")
}
mm <- mm[, !colnames(mm) %in% excludeAll]
if (!is.fullrank(mm)) {
  return("not full rank")
}

# Save model matrix for later analysis
#saveRDS(mm, file = "model.matrix.RDS")

# -- Normalize data
# order target and count matrix so they are the same (THIS IS IMPORTANT)
cm <- cm[, rownames(target)]

# CHECK IF ORDER IS THE SAME
if (all.equal(colnames(cm), rownames(target)) != TRUE) {
  print("MASSIVE WARNING: RESULTS WILL BE WRONG IF THIS IS NOT EQUAL!!!!!!!!")
  print(rownames(target))
  print(colnames(cm))
}

# normalize
cm2 <- DGEList(counts = cm)
cm2 <- calcNormFactors(cm2, method = "TMM") # TMM normalization
corfit <- duplicateCorrelation(cm2$counts, block = factor(target$Animal_ID))
png("mean_variance_norm.png")
Pi.CPM <- voom(counts = cm2, design = mm, correlation = corfit, normalize.method = "none", plot = T, span = 0.1)
dev.off()
write.csv(Pi.CPM$E, "Pi.CPM$E_all.csv")


# box plot of unfiltered data
png("boxplot_vnorm_all.png", width = 10, height = 8, units = "in", res = 100)
# par(mar=c(1,1,1,1))
minvalue <- min(Pi.CPM$E)
maxvalue <- max(Pi.CPM$E)
boxplot(Pi.CPM$E,
        labels = target$GaleID, ylim = c(minvalue - 1, maxvalue + 1),
        ylab = "voom expression", main = "Count matrix", cex.axis = .6, las = 2,
        frame = FALSE
)
dev.off()


#-- filter out genes from each group that are below mean count of 3 across samples 
#Iteratively adjusted thresholds and decided that 3 was the best cutoff to get rid of the 
# mean-variance hook shown in the initial mean-variance plot (iterative results not saved just ran in R)
A <- rowMeans(cm)
isexpr <- A >= 3
cmfl_counts <- cm[isexpr, ]
write.csv(cmfl_counts, "count_matrix_renamed_fl.csv")

#Normalize again with filtering
cm2 <- DGEList(counts = cmfl_counts)
cm2 <- calcNormFactors(cm2, method = "TMM") # TMM normalization
corfit <- duplicateCorrelation(cm2$counts, block = factor(target$Animal_ID))
png("mean_variance_norm_fl.png")
Pi.CPM <- voom(counts = cm2, design = mm, correlation = corfit, normalize.method = "none", plot = T, span = 0.1)
dev.off()
write.csv(Pi.CPM$E, "Pi.CPM$E_all_fl.csv")

# Save voom normalized object to resume analysis
saveRDS(Pi.CPM, file = "Pi.CPM.rds")
Pi.CPM <- readRDS("Pi.CPM.rds")

# Get sig genes with gene names
sig_HGNC <- merge(rhesus2human, Pi.CPM$E,
                  by.x = "Gene.stable.ID",
                  by.y = "row.names",
                  all.X = T, all.Y = T
)

sig_HGNC <- sig_HGNC[, !(names(sig_HGNC) %in% c("Gene.stable.ID"))]
sig_HGNC <- avereps(sig_HGNC,
                    ID = sig_HGNC$HGNC.symbol
)
rownames(sig_HGNC) <- sig_HGNC[, "HGNC.symbol"]
sig_HGNC <- sig_HGNC[, !(colnames(sig_HGNC) %in% c("HGNC.symbol"))]
sig_HGNC <- as.matrix(data.frame(sig_HGNC))
write.csv(sig_HGNC, "norm_matrix_HGNC_Aging_02.csv", quote = FALSE)


# box plot of filtered data
png("boxplot_vnorm_all_fl.png", width = 10, height = 8, units = "in", res = 100)
# par(mar=c(1,1,1,1))
minvalue <- min(Pi.CPM$E)
maxvalue <- max(Pi.CPM$E)
boxplot(Pi.CPM$E,
        labels = target$GaleID, ylim = c(minvalue - 1, maxvalue + 1),
        ylab = "voom normalized expression", main = "Normalized count matrix", cex.axis = .6, las = 2,
        frame = FALSE
)
dev.off()

#Feature reduction
p <- PCAtools::pca(Pi.CPM$E, 
                   metadata = target)

PCAtools::biplot(p, x='PC1', y='PC2', 
                 lab=NULL, 
                 shape='Aged',
                 colby='Time_Point',
                 pointSize = 3,
                 legendPosition = 'right',
                 #hline=0,
                 #vline=0,
                 axisLabSize = 10,
                 legendLabSize = 10,
                 legendTitleSize = 10,
                 showLoadings = FALSE,
                 legendIconSize = 3,
                 gridlines.major = FALSE,
                 gridlines.minor = FALSE)
dev.off()

screeplot(p, xlab = "Principal component", title = "SCREE plot", ylim = c(0,100), components = getComponents(p, 1:10),
          colBar = "dodgerblue")

PCAtools::biplot(p, x='PC1', y='PC2', 
                 lab=NULL, 
                 shape='Aged',
                 colby='Animal_ID',
                 pointSize = 3,
                 legendPosition = 'right',
                 #hline=0,
                 #vline=0,
                 axisLabSize = 10,
                 legendLabSize = 10,
                 legendTitleSize = 10,
                 showLoadings = FALSE,
                 legendIconSize = 3,
                 gridlines.major = FALSE,
                 gridlines.minor = FALSE)









## Generate lmfit object ##

Pi.lmfit <- lmFit(Pi.CPM, design = mm, block = target$Animal_ID, correlation = corfit$consensus)
# saveRDS(Pi.lmfit, file = "lmfitobject.RDS")
#Pi.lmfit <- readRDS("lmfitobject.RDS")

contrastsmatrix <- c(
  "bioreps1_8_A-bioreps1_8_Y",
  "bioreps2_26_A-bioreps2_26_Y",
  "bioreps3_22_A-bioreps3_22_Y",
  "bioreps4_19_A-bioreps4_19_Y",
  "bioreps4_26_A-bioreps4_26_Y",
  "bioreps5_3_A-bioreps5_3_Y"
)

contr <- makeContrasts(contrasts = contrastsmatrix, levels = mm)

fit <- contrasts.fit(Pi.lmfit, contr) # This results in a 3-dimensional array, which is not entirely analogous to fit2 from above. The next two lines extract the coefficients and p-value matrices from the array, which is all you need to use the decideTests() function.
fit2 <- eBayes(fit, robust = TRUE, trend = TRUE)
#fit2 <- treat(fit, lfc = 0.58, robust = TRUE, trend = TRUE)
results <- decideTests(fit2, lfc = (.58), method = "separate", adjust.method = "BH", p.value = 0.05)
sum_results  <- summary(results)
write.csv(sum_results,file = "sum_results.csv")

n <- 0
for (i in contrastsmatrix) {
  n <- n+1
  topTable_data <- topTable(fit2, coef=n, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
  
  topTable_data$id <- 1:nrow(topTable_data)
  
  topTable_data_genes <- merge(rhesus2human, topTable_data,
                               by.x = "Gene.stable.ID",
                               by.y = "row.names",
                               all.X = T, all.Y = T
  )
  
  topTable_data_genes <- topTable_data_genes[order(topTable_data_genes$id),]
  write.csv(topTable_data_genes, file = paste0(i, "_gene_table.csv"))
}



dataMatrixde <- fit2$coefficients
sigMask <- dataMatrixde * (results**2) # 1 if significant, 0 otherwise

ExpressMatrixde <- subset(dataMatrixde, rowSums(sigMask) != 0)

ExpressMatrixde <- merge(rhesus2human, ExpressMatrixde,
                         by.x = "Gene.stable.ID",
                         by.y = "row.names",
                         all.X = T, all.Y = T
)

ExpressMatrixde$HGNC.symbol[ExpressMatrixde$HGNC.symbol==''] = ExpressMatrixde$Gene.stable.ID[ExpressMatrixde$HGNC.symbol ==""]
ExpressMatrixde$HGNC.symbol <- make.unique(ExpressMatrixde$HGNC.symbol)
rownames(ExpressMatrixde) <- ExpressMatrixde[, "HGNC.symbol"]

ExpressMatrixde <- ExpressMatrixde[,c(3,4,5,6,7,8)]

ExpressMatrixde <- as.matrix(data.frame(ExpressMatrixde, check.names = FALSE))
ExpressMatrixde

class(ExpressMatrixde) <- "numeric"
collabels <- c(
  "Aged - Young 1/8/21",
  "Aged - Young 2/26/21",
  "Aged - Young 3/22/21",
  "Aged - Young 4/19/21",
  "Aged - Young 4/26/21",
  "Aged - Young 5/3/21"
)

colnames(ExpressMatrixde) <- collabels

hmap_big <- heatmap.L.4(ExpressMatrixde,
            figmargins = c(20, 10), cutoffmethod = "number",
            cutoff = 4, distmethod = "euclidean", cexcol = 2,
            clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9, labRow = FALSE
)

###WebGestaltR analysis####
for (cluster in unique(hmap_big$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(hmap_big$modulesrows == cluster)
  #print(genes)
  WebGestaltR(enrichMethod="ORA",
              organism="hsapiens",
              interestGene = names(genes),
              interestGeneType="genesymbol",
              referenceGene = rhesus2human[,2],
              referenceGeneType = "genesymbol",
              projectName = paste0(cluster,"_ORA"),
              minNum=5,
              enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                 "geneontology_Molecular_Function_noRedundant",
                                 "geneontology_Cellular_Component_noRedundant",
                                 "pathway_KEGG"
              ),
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)
}





## Generate tables for Jenny ##
topTable_data <- topTable(fit2, coef=1, resort.by=NULL, number=length(fit2$coefficients))
dataMatrixde <- fit2$coefficients
sigMask <- dataMatrixde * (results**2) # 1 if significant, 0 otherwise

topTable_data_sort <- topTable_data[match(rownames(dataMatrixde), rownames(topTable_data)),]
topTable_data_sig <- subset(topTable_data_sort, rowSums(sigMask) != 0)
topTable_data_sig_genes <- merge(rhesus2human, topTable_data_sig,
                                 by.x = "Gene.stable.ID",
                                 by.y = "row.names",
                                 all.X = T, all.Y = T
)

topTable_data_sig_genes[duplicated(topTable_data_sig_genes$Gene.stable.ID),]


n <- 0
top_big <- data.frame(matrix(NA, nrow = 169, ncol = 1))
for (i in contrastsmatrix) {
  n <- n+1
  topTable_data <- topTable(fit2, coef=n, resort.by=NULL, number=length(fit2$coefficients))
  
  topTable_data_sort <- topTable_data[match(rownames(dataMatrixde), rownames(topTable_data)),]
  topTable_data_sig <- subset(topTable_data_sort, rowSums(sigMask) != 0)
  topTable_data_sig_genes <- merge(rhesus2human, topTable_data_sig,
                                   by.x = "Gene.stable.ID",
                                   by.y = "row.names",
                                   all.X = T, all.Y = T
  )
  top_big <- cbind(top_big, topTable_data_sig_genes)
}

write.csv(top_big, file = "aged-young_DEG_all.csv")




#### Look at loadings ####
loadings <- p$loadings

loadings <- merge(rhesus2human, loadings,
                         by.x = "Gene.stable.ID",
                         by.y = "row.names",
                         all.X = T, all.Y = T
)

write.csv(loadings, file = "loadings.csv")


# Generate PCA with gene names on the loadings
E <- Pi.CPM$E

E <- merge(rhesus2human, E,
           by.x = "Gene.stable.ID",
           by.y = "row.names",
           all.X = T, all.Y = T
)

E$HGNC.symbol[E$HGNC.symbol==''] = E$Gene.stable.ID[E$HGNC.symbol ==""]
E$HGNC.symbol <- make.unique(E$HGNC.symbol)
rownames(E) <- E[, "HGNC.symbol"]
E <- subset(E, select = -c(Gene.stable.ID, HGNC.symbol))


q <- PCAtools::pca(E, 
                   metadata = target)

PCAtools::biplot(q, x='PC1', y='PC2', 
                 lab=NULL, 
                 shape='Aged',
                 colby='Time_Point',
                 pointSize = 3,
                 legendPosition = 'right',
                 #hline=0,
                 #vline=0,
                 axisLabSize = 10,
                 legendLabSize = 10,
                 legendTitleSize = 10,
                 showLoadings = TRUE,
                 legendIconSize = 3,
                 gridlines.major = FALSE,
                 gridlines.minor = FALSE)


PCAtools::biplot(q, x='PC1', y='PC2', 
                 lab=NULL, 
                 shape='Aged',
                 colby='Animal_ID',
                 pointSize = 3,
                 legendPosition = 'right',
                 #hline=0,
                 #vline=0,
                 axisLabSize = 10,
                 legendLabSize = 10,
                 legendTitleSize = 10,
                 showLoadings = TRUE,
                 legendIconSize = 3,
                 gridlines.major = FALSE,
                 gridlines.minor = FALSE)


##### Young monkeys only #####


## Generate lmfit object ##

Pi.lmfit <- lmFit(Pi.CPM, design = mm, block = target$Animal_ID, correlation = corfit$consensus)
# saveRDS(Pi.lmfit, file = "lmfitobject.RDS")
#Pi.lmfit <- readRDS("lmfitobject.RDS")

contrastsmatrix <- c(
  "bioreps3_22_Y-(bioreps1_8_Y+bioreps2_26_Y)",
  "bioreps4_19_Y-(bioreps1_8_Y+bioreps2_26_Y)",
  "bioreps4_26_Y-(bioreps1_8_Y+bioreps2_26_Y)",
  "bioreps5_3_Y-(bioreps1_8_Y+bioreps2_26_Y)"
)

contr <- makeContrasts(contrasts = contrastsmatrix, levels = mm)

fit <- contrasts.fit(Pi.lmfit, contr) # This results in a 3-dimensional array, which is not entirely analogous to fit2 from above. The next two lines extract the coefficients and p-value matrices from the array, which is all you need to use the decideTests() function.
fit2 <- eBayes(fit, robust = TRUE, trend = TRUE)
#fit2 <- treat(fit, lfc = 0.58, robust = TRUE, trend = TRUE)
results <- decideTests(fit2, lfc = (.58), method = "separate", adjust.method = "BH", p.value = 0.05)
sum_results  <- summary(results)
write.csv(sum_results,file = "young_sum_results.csv")

n <- 0
for (i in contrastsmatrix) {
  n <- n+1
  topTable_data <- topTable(fit2, coef=n, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
  
  topTable_data$id <- 1:nrow(topTable_data)
  
  topTable_data_genes <- merge(rhesus2human, topTable_data,
                               by.x = "Gene.stable.ID",
                               by.y = "row.names",
                               all.X = T, all.Y = T
  )
  
  topTable_data_genes <- topTable_data_genes[order(topTable_data_genes$id),]
  write.csv(topTable_data_genes, file = paste0(i, "_gene_table.csv"))
}



dataMatrixde <- fit2$coefficients
sigMask <- dataMatrixde * (results**2) # 1 if significant, 0 otherwise

ExpressMatrixde <- subset(dataMatrixde, rowSums(sigMask) != 0)

ExpressMatrixde <- merge(rhesus2human, ExpressMatrixde,
                         by.x = "Gene.stable.ID",
                         by.y = "row.names",
                         all.X = T, all.Y = T
)

ExpressMatrixde$HGNC.symbol[ExpressMatrixde$HGNC.symbol==''] = ExpressMatrixde$Gene.stable.ID[ExpressMatrixde$HGNC.symbol ==""]
ExpressMatrixde$HGNC.symbol <- make.unique(ExpressMatrixde$HGNC.symbol)
rownames(ExpressMatrixde) <- ExpressMatrixde[, "HGNC.symbol"]
ExpressMatrixde <- subset(ExpressMatrixde, select = -c(Gene.stable.ID, HGNC.symbol))

ExpressMatrixde <- as.matrix(data.frame(ExpressMatrixde, check.names = FALSE))
ExpressMatrixde

class(ExpressMatrixde) <- "numeric"
collabels <- c(
  "Young 3/22/21 - Baseline",
  "Young 4/19/21 - Baseline",
  "Young 4/26/21 - Baseline",
  "Young 5/3/21 - Baseline"
)

colnames(ExpressMatrixde) <- collabels

heatmap.L.4(ExpressMatrixde,
            figmargins = c(20, 10),
            cutoff = 1, distmethod = "euclidean", cexcol = 2,
            clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9, labRow = FALSE
)

nrow(ExpressMatrixde)


##### Aged monkeys only #####


## Generate lmfit object ##

Pi.lmfit <- lmFit(Pi.CPM, design = mm, block = target$Animal_ID, correlation = corfit$consensus)
# saveRDS(Pi.lmfit, file = "lmfitobject.RDS")

contrastsmatrix <- c(
  "bioreps3_22_A-(bioreps1_8_A+bioreps2_26_A)",
  "bioreps4_19_A-(bioreps1_8_A+bioreps2_26_A)",
  "bioreps4_26_A-(bioreps1_8_A+bioreps2_26_A)",
  "bioreps5_3_A-(bioreps1_8_A+bioreps2_26_A)"
)

contr <- makeContrasts(contrasts = contrastsmatrix, levels = mm)

fit <- contrasts.fit(Pi.lmfit, contr) # This results in a 3-dimensional array, which is not entirely analogous to fit2 from above. The next two lines extract the coefficients and p-value matrices from the array, which is all you need to use the decideTests() function.
fit2 <- eBayes(fit, robust = TRUE, trend = TRUE)
#fit2 <- treat(fit, lfc = 0.58, robust = TRUE, trend = TRUE)
results <- decideTests(fit2, lfc = (.58), method = "separate", adjust.method = "BH", p.value = 0.05)
sum_results  <- summary(results)
write.csv(sum_results,file = "aged_sum_results.csv")

n <- 0
for (i in contrastsmatrix) {
  n <- n+1
  topTable_data <- topTable(fit2, coef=n, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
  
  topTable_data$id <- 1:nrow(topTable_data)
  
  topTable_data_genes <- merge(rhesus2human, topTable_data,
                               by.x = "Gene.stable.ID",
                               by.y = "row.names",
                               all.X = T, all.Y = T
  )
  
  topTable_data_genes <- topTable_data_genes[order(topTable_data_genes$id),]
  write.csv(topTable_data_genes, file = paste0(i, "_gene_table.csv"))
}



dataMatrixde <- fit2$coefficients
sigMask <- dataMatrixde * (results**2) # 1 if significant, 0 otherwise

ExpressMatrixde <- subset(dataMatrixde, rowSums(sigMask) != 0)

ExpressMatrixde <- merge(rhesus2human, ExpressMatrixde,
                         by.x = "Gene.stable.ID",
                         by.y = "row.names",
                         all.X = T, all.Y = T
)

ExpressMatrixde$HGNC.symbol[ExpressMatrixde$HGNC.symbol==''] = ExpressMatrixde$Gene.stable.ID[ExpressMatrixde$HGNC.symbol ==""]
ExpressMatrixde$HGNC.symbol <- make.unique(ExpressMatrixde$HGNC.symbol)
rownames(ExpressMatrixde) <- ExpressMatrixde[, "HGNC.symbol"]
ExpressMatrixde <- subset(ExpressMatrixde, select = -c(Gene.stable.ID, HGNC.symbol))

ExpressMatrixde <- as.matrix(data.frame(ExpressMatrixde, check.names = FALSE))
ExpressMatrixde

class(ExpressMatrixde) <- "numeric"
collabels <- c(
  "Aged 3/22/21 - Baseline",
  "Aged 4/19/21 - Baseline",
  "Aged 4/26/21 - Baseline",
  "Aged 5/3/21 - Baseline"
)

colnames(ExpressMatrixde) <- collabels

hmap <- heatmap.L.4(ExpressMatrixde,
            figmargins = c(20, 10),
            cutoff = 1, distmethod = "euclidean", cexcol = 2,
            clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9, labRow = FALSE
)

hmap

nrow(ExpressMatrixde)


#### Look at just first two timepoints ####

Pi.lmfit <- lmFit(Pi.CPM, design = mm, block = target$Animal_ID, correlation = corfit$consensus)
# saveRDS(Pi.lmfit, file = "lmfitobject.RDS")

contrastsmatrix <- c(
  "bioreps1_8_A-bioreps1_8_Y",
  "bioreps2_26_A-bioreps2_26_Y"
)

contr <- makeContrasts(contrasts = contrastsmatrix, levels = mm)

fit <- contrasts.fit(Pi.lmfit, contr) # This results in a 3-dimensional array, which is not entirely analogous to fit2 from above. The next two lines extract the coefficients and p-value matrices from the array, which is all you need to use the decideTests() function.
fit2 <- eBayes(fit, robust = TRUE, trend = TRUE)
#fit2 <- treat(fit, lfc = 0.58, robust = TRUE, trend = TRUE)
results <- decideTests(fit2, lfc = (.58), method = "separate", adjust.method = "BH", p.value = 0.05)
sum_results  <- summary(results)
write.csv(sum_results,file = "sum_results_first_two.csv")

n <- 0
for (i in contrastsmatrix) {
  n <- n+1
  topTable_data <- topTable(fit2, coef=n, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
  
  topTable_data$id <- 1:nrow(topTable_data)
  
  topTable_data_genes <- merge(rhesus2human, topTable_data,
                               by.x = "Gene.stable.ID",
                               by.y = "row.names",
                               all.X = T, all.Y = T
  )
  
  topTable_data_genes <- topTable_data_genes[order(topTable_data_genes$id),]
  write.csv(topTable_data_genes, file = paste0(i, "_gene_table_first_two.csv"))
}



dataMatrixde <- fit2$coefficients
sigMask <- dataMatrixde * (results**2) # 1 if significant, 0 otherwise

ExpressMatrixde <- subset(dataMatrixde, rowSums(sigMask) != 0)

ExpressMatrixde <- merge(rhesus2human, ExpressMatrixde,
                         by.x = "Gene.stable.ID",
                         by.y = "row.names",
                         all.X = T, all.Y = T
)

ExpressMatrixde$HGNC.symbol[ExpressMatrixde$HGNC.symbol==''] = ExpressMatrixde$Gene.stable.ID[ExpressMatrixde$HGNC.symbol ==""]
rownames(ExpressMatrixde) <- ExpressMatrixde[, "HGNC.symbol"]
ExpressMatrixde <- ExpressMatrixde[,c(3,4)]

ExpressMatrixde <- as.matrix(data.frame(ExpressMatrixde, check.names = FALSE))
ExpressMatrixde

class(ExpressMatrixde) <- "numeric"
collabels <- c(
  "Aged - Young 1/8/21",
  "Aged - Young 2/26/21"
)

colnames(ExpressMatrixde) <- collabels

hmap <- heatmap.L.4(ExpressMatrixde,
            figmargins = c(20, 10),
            cutoff = 0, distmethod = "euclidean", cexcol = 2,
            clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9, labRow = FALSE
)

dev.off()


## Generate list of de genes per cluster and LFCs ##

for (cluster in unique(hmap$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(hmap$modulesrows == cluster)
  # write.csv(names(genes), paste0(cluster,"_genes_first_two.csv"))
  
  temp <- ExpressMatrixde[row.names(ExpressMatrixde)%in%names(genes),]
  write.csv(temp, paste0(cluster,"_gene_lfcs_first_two.csv"))
}


## Generate tables for Jenny ##
topTable_data <- topTable(fit2, coef=1, resort.by=NULL, number=length(fit2$coefficients))
dataMatrixde <- fit2$coefficients
sigMask <- dataMatrixde * (results**2) # 1 if significant, 0 otherwise

topTable_data_sort <- topTable_data[match(rownames(dataMatrixde), rownames(topTable_data)),]
topTable_data_sig <- subset(topTable_data_sort, rowSums(sigMask) != 0)
topTable_data_sig_genes <- merge(rhesus2human, topTable_data_sig,
                                 by.x = "Gene.stable.ID",
                                 by.y = "row.names",
                                 all.X = T, all.Y = T
)

topTable_data_sig_genes[duplicated(topTable_data_sig_genes$Gene.stable.ID),]




topTable_data_2 <- topTable(fit2, coef=2, resort.by=NULL, number=length(fit2$coefficients))
dataMatrixde <- fit2$coefficients
sigMask <- dataMatrixde * (results**2) # 1 if significant, 0 otherwise

topTable_data_sort_2 <- topTable_data_2[match(rownames(dataMatrixde), rownames(topTable_data_2)),]
topTable_data_sig_2 <- subset(topTable_data_sort_2, rowSums(sigMask) != 0)
topTable_data_sig_genes_2 <- merge(rhesus2human, topTable_data_sig_2,
                                 by.x = "Gene.stable.ID",
                                 by.y = "row.names",
                                 all.X = T, all.Y = T
)

topTable_data_sig_genes_2[duplicated(topTable_data_sig_genes_2$Gene.stable.ID),]

topTable_data_sig_genes_all <- cbind(topTable_data_sig_genes, topTable_data_sig_genes_2)

write.csv(topTable_data_sig_genes_all, file = "aged-young_DEG.csv")





#### Aged Post-rap comparison to 2/26 baseline ####

Pi.lmfit <- lmFit(Pi.CPM, design = mm, block = target$Animal_ID, correlation = corfit$consensus)
# saveRDS(Pi.lmfit, file = "lmfitobject.RDS")

contrastsmatrix <- c(
  "bioreps4_26_A-bioreps2_26_A",
  "bioreps5_3_A-bioreps2_26_A"
)

contr <- makeContrasts(contrasts = contrastsmatrix, levels = mm)

fit <- contrasts.fit(Pi.lmfit, contr) # This results in a 3-dimensional array, which is not entirely analogous to fit2 from above. The next two lines extract the coefficients and p-value matrices from the array, which is all you need to use the decideTests() function.
fit2 <- eBayes(fit, robust = TRUE, trend = TRUE)
#fit2 <- treat(fit, lfc = 0.58, robust = TRUE, trend = TRUE)
results <- decideTests(fit2, lfc = (.58), method = "separate", adjust.method = "BH", p.value = 0.05)
sum_results  <- summary(results)
write.csv(sum_results,file = "sum_results_aged_post_rap.csv")

n <- 0
for (i in contrastsmatrix) {
  n <- n+1
  topTable_data <- topTable(fit2, coef=n, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
  
  topTable_data$id <- 1:nrow(topTable_data)
  
  topTable_data_genes <- merge(rhesus2human, topTable_data,
                               by.x = "Gene.stable.ID",
                               by.y = "row.names",
                               all.X = T, all.Y = T
  )
  
  topTable_data_genes <- topTable_data_genes[order(topTable_data_genes$id),]
  write.csv(topTable_data_genes, file = paste0(i, "_gene_table_aged_post_rap.csv"))
}



dataMatrixde <- fit2$coefficients
sigMask <- dataMatrixde * (results**2) # 1 if significant, 0 otherwise

ExpressMatrixde <- subset(dataMatrixde, rowSums(sigMask) != 0)

ExpressMatrixde <- merge(rhesus2human, ExpressMatrixde,
                         by.x = "Gene.stable.ID",
                         by.y = "row.names",
                         all.X = T, all.Y = T
)

ExpressMatrixde$HGNC.symbol[ExpressMatrixde$HGNC.symbol==''] = ExpressMatrixde$Gene.stable.ID[ExpressMatrixde$HGNC.symbol ==""]
ExpressMatrixde$HGNC.symbol <- make.unique(ExpressMatrixde$HGNC.symbol)
rownames(ExpressMatrixde) <- ExpressMatrixde[, "HGNC.symbol"]
ExpressMatrixde <- ExpressMatrixde[,c(3,4)]

ExpressMatrixde <- as.matrix(data.frame(ExpressMatrixde, check.names = FALSE))
ExpressMatrixde

class(ExpressMatrixde) <- "numeric"
collabels <- c(
  "Aged 4/26 - Baseline",
  "Aged 5/3 - Baseline"
)

colnames(ExpressMatrixde) <- collabels

hmap <- heatmap.L.4(ExpressMatrixde,
                    figmargins = c(20, 10),
                    cutoff = 0, distmethod = "euclidean", cexcol = 2,
                    clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9, labRow = FALSE
)

dev.off()


## Generate list of de genes per cluster and LFCs ##

for (cluster in unique(hmap2$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(hmap2$modulesrows == cluster)
  # write.csv(names(genes), paste0(cluster,"_genes_first_two.csv"))
  
  temp <- ExpressMatrixde[row.names(ExpressMatrixde)%in%names(genes),]
  write.csv(temp, paste0(cluster,"_gene_lfcs_aged_post_rap.csv"))
}













#Heatmap(ExpressMatrixde, clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D2", cluster_columns = FALSE, cluster_rows = TRUE, row_labels = rownames(ExpressMatrixde), row_names_gp = gpar(fontsize = 4))




###WebGestaltR analysis####
for (cluster in unique(hmap$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(hmap$modulesrows == cluster)
  #print(genes)
  WebGestaltR(enrichMethod="ORA",
              organism="hsapiens",
              interestGene = names(genes),
              interestGeneType="genesymbol",
              referenceGene = rhesus2human[,2],
              referenceGeneType = "genesymbol",
              projectName = paste0(cluster,"_ORA_small"),
              minNum=5,
              enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                 "geneontology_Molecular_Function_noRedundant",
                                 "geneontology_Cellular_Component_noRedundant",
                                 "pathway_KEGG"
              ),
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)
}

## Generate list of de genes per cluster and LFCs ##

for (cluster in unique(global_modulesalllung$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(global_modulesalllung$modulesrows == cluster)
  write.csv(names(genes), paste0(cluster,"_genes_10wk.csv"))
  
  temp <- ExpressMatrixde_genes[row.names(ExpressMatrixde_genes)%in%names(genes),]
  write.csv(temp, paste0(cluster,"_gene_lfcs_all_lung.csv"))
}

## Cibersort
cibersort <- Pi.CPM$E
cibersort <- merge(rhesus2human, cibersort,
                   by.x = "Gene.stable.ID",
                   by.y = "row.names",
                   all.X = T, all.Y = T
)
cibersort <- cibersort[, !(names(cibersort) %in% c("Gene.stable.ID"))]
cibersort <- avereps(cibersort,
                     ID = cibersort$HGNC.symbol
)
rownames(cibersort) <- cibersort[, "HGNC.symbol"]
cibersort <- cibersort[, !(colnames(cibersort) %in% c("HGNC.symbol"))]
cibersort <- as.matrix(data.frame(cibersort))
# cibersort <- make.unique(rownames(cibersort))
class(cibersort) <- "numeric"

cibersortnotlog <- 2^cibersort
cibersortnotlog <- cibersortnotlog[!(rownames(cibersortnotlog) %in% c("")), ]
write.table(cibersortnotlog, file.path("HGNC_cibersort_allnotlog.txt"), sep="\t", quote=F)
write.table(cibersort, file.path("HGNC_cibersort_all.txt"), sep="\t", quote=F)




