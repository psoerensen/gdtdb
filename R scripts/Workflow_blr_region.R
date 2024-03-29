################################################################################
# Fine-mapping using BLR models for LD regions
################################################################################
#
#library(devtools)
#devtools::install_github("psoerensen/gact")
#devtools::install_github("psoerensen/qgg")
#
################################################################################

library(qgg)
library(gact)


# Load GAlist with information on gact database
GAlist <- readRDS(file="/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds")

# Check studies in gact database
GAlist$studies

# List BLR results files allready in gact database
list.files(GAlist$dirs["gbayes"])

# Load Glist with information on 1000G
Glist <- readRDS(file=file.path(GAlist$dirs["marker"],"Glist_1000G.rds"))

# Select study
studyID <- "GWAS1"

# Extract GWAS summary statistics for GWAS1
stat <- getMarkerStatDB(GAlist=GAlist, studyID=studyID)

# Check and align summary statistics based on marker information in Glist
stat <- checkStat(Glist=Glist, stat=stat,
                  excludeMAF=0.05,
                  excludeMAFDIFF=0.05,
                  excludeINFO=0.8,
                  excludeCGAT=TRUE,
                  excludeINDEL=TRUE,
                  excludeDUPS=TRUE,
                  excludeMHC=FALSE,
                  excludeMISS=0.05,
                  excludeHWE=1e-12)

# Fit BLR model across all chromosomes
for ( chr in 22:1) {

  # Create LD sets  
  sets <- qgg:::createLDsets(ldscores=Glist$ldscores[[chr]], 
                             maxsize=2000, msize=50, verbose=TRUE)
  
  # Map LD sets to rsids in stat
  sets <- qgg:::mapSets(sets, rsids=stat$rsids, index=FALSE)

  
  # Fine-mapping using BLR model on selected regions
  fit <- gmap(Glist=Glist, stat=stat, sets=sets,
              method="bayesC", pi=0.001, h2=0.1,
              vb=0.01, ssb_prior=0.01, nub=5,
              checkConvergence=TRUE, ntrial=3, 
              pruneLD=TRUE, r2=0.95, checkLD=FALSE,
              nit=1000, nburn=500, verbose=TRUE,
              updateB=TRUE, updatePi=FALSE, updateG=TRUE, updateE=TRUE)
  
  # Posterior estimates of hyper-parameters (ve,vg,vb,pi,pip) for every fine-mapped region
  head(fit$post)
  
  # Convergence statistics for every fine-mapped region
  head(fit$conv)
  
  # Posterior estimates of marker effects for every fine-mapped region
  head(fit$stat)
  
  # Save fit object for all sets on chromosome in directory of choice (here gbayes folder in database)
  #saveRDS(fit, file=file.path(GAlist$dirs["gbayes"],paste0("fit_blr_pruned_chr",chr,"_",studyID,".rds")))
  
}



