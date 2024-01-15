################################################################################
# Install new version of gact and qgg
################################################################################
#
#library(devtools)
#devtools::install_github("psoerensen/gact")
#devtools::install_github("psoerensen/qgg")
#
################################################################################

library(devtools)
devtools::install_github("psoerensen/gact")
devtools::install_github("psoerensen/qgg")

################################################################################
# Download and install gact database
################################################################################

# Load libraries
library(qgg)
library(gact)

# Define working directory for storing the data base
dbdir <- "/faststorage/project/ukbiobank/projects/gact"


# Create infrastructure and download database (can take several minutes)
GAlist <- gact(version="hsa.0.0.1", dbdir=dbdir, task="download")

# Save GAlist for use later
saveRDS(GAlist, file="/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds")

# Update GAlist for use later
GAlist <- readRDS(file="/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds")


################################################################################
# Overview of database content
################################################################################

# Overview of GWAS summary statistics in database
GAlist$studies

# Accessing the 'dirs' slot in the GAlist object
GAlist$dirs

# Listing all directories recursively under each path in GAlist$dirs
list.dirs(GAlist$dirs, recursive = TRUE)

# Listing files under the 'marker' directory specified in GAlist$dirs
list.files(GAlist$dirs["marker"])

# Listing files under the 'gstat' directory specified in GAlist$dirs
list.files(GAlist$dirs["gstat"])

# Listing files under the 'gsea' directory specified in GAlist$dirs
list.files(GAlist$dirs["gsea"])



