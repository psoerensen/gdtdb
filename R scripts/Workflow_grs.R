################################################################################
# Genomic risk scoring based on BLR summary statistics
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

# Load Glist with information on 1000G (i.e. genotypes for the target/test population)
Glist <- readRDS(file=file.path(GAlist$dirs["marker"],"Glist_1000G.rds"))


# Check studies in gact database
GAlist$studies

# List BLR results files allready in gact database
list.files(GAlist$dirs["gbayes"])

# Select study
studyID <- "GWAS1"

# Extract GWAS summary statistics for GWAS1
stat <- fread(file.path(GAlist$dirs["gbayes"],paste0(studyID,"_stat_BayesC.txt.gz")), data.table=FALSE)

# Check and align summary statistics based on marker information in Glist
stat <- checkStat(Glist=Glist, stat=stat)

# Compute genomic risk score for GWAS1
grs <- gscore(Glist=Glist,stat=stat)

# Save grs object for all genes in directory of choice (here gbayes folder in database)
saveRDS(fit, file=file.path(GAlist$dirs["grisk"],"grs_blr_genes_GWAS1.rds"))

