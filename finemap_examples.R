library(qgg)
library(data.table)
library(gact)

# Glist for 1000G genetic data used in the GSEA
Glist <- readRDS(file="C:/Users/au223366/Dropbox/Projects/1000G/Glist_MAF1_LD1K_1000G.rds")

# GAlist with information about GDT database
GAlist <- readRDS(file="C:/Users/au223366/Dropbox/Projects/balder/gdtdb/GAlist_t2dm-gact-0.0.2.rds")


stat <- getMarkerStat(GAlist=GAlist, studies="GWAS1")

fit <- gmap(Glist=Glist, stat=stat, 
            trait=1, msize=5000,   
            method="bayesC", nit=500, pi=0.001,
            updateB=TRUE, updatePi=TRUE, updateG=TRUE, updateE=TRUE)

ves <- sapply(fit$ves, function(x){ mean(x) })
vbs <- sapply(fit$vbs, function(x){ mean(x) })
vgs <- sapply(fit$vgs, function(x){ mean(x) })
sum(vgs)

# create some marker sets
msize <- 5000
chr <- 20
rsidsLD <- Glist$rsidsLD[[chr]]
rsidsLD <- rsidsLD[rsidsLD%in%rownames(stat$b)]
sets20 <- split(rsidsLD, ceiling(seq_along(rsidsLD) / msize))
chr <- 22
rsidsLD <- Glist$rsidsLD[[chr]]
rsidsLD <- rsidsLD[rsidsLD%in%rownames(stat$b)]
sets22 <- split(rsidsLD, ceiling(seq_along(rsidsLD) / msize))

sets <- list(sets20[[1]],sets22[[1]],sets20[[2]])

fit <- gmap(Glist=Glist, stat=stat, 
            sets=sets, formatLD="dense", shrinkLD=TRUE,   
            method="bayesR", nit=1000, pi=0.001,
            updateB=TRUE, updatePi=TRUE, updateG=TRUE, updateE=TRUE)

ves <- sapply(fit$ves, function(x){ mean(x) })
vbs <- sapply(fit$vbs, function(x){ mean(x) })
vgs <- sapply(fit$vgs, function(x){ mean(x) })
sum(vgs)


fit <- gmap(Glist=Glist, stat=stat, 
                  sets=sets, formatLD="sparse", shrinkLD=TRUE,   
                  method="bayesC", nit=1000, pi=0.001,
                  updateB=TRUE, updatePi=TRUE, updateG=TRUE, updateE=TRUE)

