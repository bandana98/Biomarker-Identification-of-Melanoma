#1. DEGs analysis for TCGA data
sampleinfo <- read.csv("data/input/sampleinfoT.csv")
sampleinfo <- data.frame(sampleinfo, row.names = T)

# load expression matrix
load("data/input/tcga_gtex_exprsT.RData")

# check if no.of samples and names are the same in both files
all(colnames(tcga_gtex_exprsT) == rownames(sampleinfo))
all(colnames(tcga_gtex_exprsT) %in% rownames(sampleinfo))

# group samples
group = factor(sampleinfo$cancer_status)
group = relevel(group, ref="Normal")

# Create DGE object
dge <- DGEList(counts = tcga_gtex_exprsT, group = group)

# Filter out lowly expressed genes
dge <- dge[rowSums(cpm(dge) > 1) >= 3, , keep.lib.sizes = FALSE]

# Normalization
dge <- calcNormFactors(estimateDisp(dge), method = "TMM")

# Testing for Differential Expression Analysis
test <- exactTest(dge)
tcga_degs_TvN <- topTags(test, n = Inf)$table

# save degs file
save(tcga_degs_TvN, file="data/output/tcga_degs_TvsN.RData")


#2. DEGs analysis for GEO data
phenoData <- read.csv("data/input/GSE15605_pheno.csv", row.names = 1)

# load expression matrix
load("data/input/GSE15605_exprs.RData")

# # check if no.of samples and names are the same in both files
all(rownames(phenoData) %in% colnames(GSE15605_exprs))
all(rownames(phenoData) == colnames(GSE15605_exprs))

# Design and fit model
design <- model.matrix(~0+Tumor_status, phenoData)
colnames(design) <- c("Normal", "Tumor")
fit <- lmFit(GSE15605_exprs, design)

# Make contrast and fit 
fit <- eBayes(contrasts.fit(lmFit(GSE15605_exprs, design), makeContrasts(Tumor - Normal, levels = design)), 0.01)

# Get top table
gse15605_degs_TvN <- topTable(fit, adjust.method = "BH", sort.by="B", number=Inf)

# save degs file
save(gse15605_degs_TvN, file="data/output/gse15605_degs_TvN.RData")
