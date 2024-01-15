################################################################################
# Pathway prioritization using MAGMA based on BLR model summary statistics
################################################################################
#
#library(devtools)
#devtools::install_github("psoerensen/gact")
#devtools::install_github("psoerensen/qgg")
#
################################################################################


library(qgg)
library(gact)
library(msigdbr)


# Load GAlist with information on gact database
GAlist <- readRDS(file="/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds")

# Select study
studyID <- "GWAS2"

# Load BLR result
fit <- readRDS(file=file.path(GAlist$dirs["gbayes"],paste0("fit_blr_pruned_",studyID,".rds")))
head(fit$post)
head(fit$conv)
head(fit$stat)

# Extract gene-marker sets
msets <- getMarkerSetsDB(GAlist = GAlist, feature = "Genesplus")
msets <- mapSets(msets,fit$stat$rsids,index=TRUE)

# Compute gene-level statistics from BLR fit object
bm <- sapply(msets,function(x){sum(abs(fit$stat$bm[x]))})
dm <- sapply(msets,function(x){sum(abs(fit$stat$dm[x]))})
vm <- sapply(msets,function(x){sum(abs(fit$stat$vm[x]))})

# Get gene sets for KEGG pathways from mSigDB
msigdb <- msigdbr(species = "human", category = "C2", subcategory = "CP:KEGG")
psets <- split(msigdb$ensembl_gene, f=msigdb$gs_name)
psets <- mapSets(psets,names(bm), index=FALSE)

# Get gene sets for GWAS catalog from gact database
psets <- getSetsDB(GAlist = GAlist, feature = "GWAScatalogPlus")
psets <- mapSets(psets,names(bm), index=FALSE)
psets <- psets[sapply(psets,length)>20]


# Fit linear model (i.e. MAGMA procedure) using BLR gene-level statistics
resBm <- qgg:::magma(stat=bm, sets=psets)
resDm <- qgg:::magma(stat=dm, sets=psets)
resVm <- qgg:::magma(stat=vm, sets=psets)


# Compute BRL derived gene set statistics (e.g. sum of gene-level statistics)
BmSet <- sapply(psets,function(x){sum(bm[x])})
DmSet <- sapply(psets,function(x){sum(dm[x])})
VmSet <- sapply(psets,function(x){sum(vm[x])})


pairs(cbind(BmSet,DmSet,VmSet))
cor(cbind(BmSet,DmSet,VmSet))


# QQ-plot of MAGMA results for different BLR derived statistics
qplot <- function(p=NULL, main = "") {
  mlogObs <- -log10(p)
  m <- length(mlogObs)
  o <- order(mlogObs, decreasing = TRUE)
  mlogExp <- -log10((1:m) / m)
  plot( y = mlogObs[o], x = mlogExp, col = 2, pch = "+",
        frame.plot = FALSE, main = main, xlab = "Expected -log10(p)", ylab = "Observed -log10(p)")
  abline(a = 0, b = 1)
}

layout(matrix(1:3, ncol=3))
qplot(p=resBm$p, main="Bm")
qplot(p=resDm$p, main="Dm")
qplot(p=resVm$p, main="Vm")


