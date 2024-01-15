################################################################################
# Prepare Glist for 1000G data
################################################################################
#
#library(devtools)
#devtools::install_github("psoerensen/gact")
#devtools::install_github("psoerensen/qgg")
#
################################################################################


# Load libraries
library(qgg)
library(gact)
library(data.table)

# Load GAlist
GAlist <- readRDS(file="/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds")

# Download 1000G data (if not allready downloaded)
GAlist <- downloadDB(GAlist=GAlist, what="1000G")

# Define location of bed/bim/fam files
path <- file.path(GAlist$dirs["marker"],"1000G_EUR_Phase3_plink/1000G.EUR.QC.")
bedfiles <- paste(path,1:22,".bed",sep="")
bimfiles <- paste(path,1:22,".bim",sep="")
famfiles <- paste(path,1:22,".fam",sep="")

# Prepare summary (i.e. Glist) for 1000G genotype data
Glist <- gprep(study="1000G",
               bedfiles=bedfiles,
               bimfiles=bimfiles,
               famfiles=famfiles)

# Save Glist for use later
saveRDS(Glist, file=file.path(GAlist$dirs["marker"],"Glist_1000G.rds"))


################################################################################
# Prepare Glist for 1000G data for different ancestries
################################################################################

# Load GAlist
GAlist <- readRDS(file="/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds")

# Marker IDs in database
rsids <- GAlist$rsids


#####################################
# Process EUR 1000G data
#####################################

# Define the file paths for the original bed/bim/fam files to read
bedRead <- file.path(GAlist$dirs["marker"], "g1000_eur.bed")
bimRead <- file.path(GAlist$dirs["marker"], "g1000_eur.bim")
famRead <- file.path(GAlist$dirs["marker"], "g1000_eur.fam")

# Define the file paths for the filtered bed/bim/fam files to write
bedWrite <- file.path(GAlist$dirs["marker"], "g1000_eur_filtered.bed")
bimWrite <- file.path(GAlist$dirs["marker"], "g1000_eur_filtered.bim")
famWrite <- file.path(GAlist$dirs["marker"], "g1000_eur_filtered.fam")

# Call the writeBED function to filter and write the data
writeBED(
 bedRead = bedRead,
 bimRead = bimRead,
 famRead = famRead,
 bedWrite = bedWrite,
 bimWrite = bimWrite,
 famWrite = famWrite,
 rsids = rsids
)

# Define location of bed/bim/fam files
bedfiles <- file.path(GAlist$dirs["marker"],"g1000_eur_filtered.bed")
bimfiles <- file.path(GAlist$dirs["marker"],"g1000_eur_filtered.bim")
famfiles <- file.path(GAlist$dirs["marker"],"g1000_eur_filtered.fam")

# Prepare summary (i.e. Glist) for 1000G genotype data
Glist <- gprep(study="1000G EUR",
               bedfiles=bedfiles,
               bimfiles=bimfiles,
               famfiles=famfiles)

# Save Glist for use later
saveRDS(Glist, file=file.path(GAlist$dirs["marker"],"Glist_1000G_eur_filtered.rds"))



#####################################
# Process EAS
#####################################

# Define the file paths for the original bed/bim/fam files to read
bedRead <- file.path(GAlist$dirs["marker"],"g1000_eas.bed")
bimRead <- file.path(GAlist$dirs["marker"],"g1000_eas.bim")
famRead <- file.path(GAlist$dirs["marker"],"g1000_eas.fam")

# Define the file paths for the filtered bed/bim/fam files to write
bedWrite <- file.path(GAlist$dirs["marker"],"g1000_eas_filtered.bed")
bimWrite <- file.path(GAlist$dirs["marker"],"g1000_eas_filtered.bim")
famWrite <- file.path(GAlist$dirs["marker"],"g1000_eas_filtered.fam")

# Call the writeBED function to filter and write the data
writeBED(bedRead=bedRead, bimRead=bimRead, famRead=famRead,
         bedWrite=bedWrite, bimWrite=bimWrite, famWrite=famWrite,
         rsids=rsids)

# Define location of bed/bim/fam files
bedfiles <- file.path(GAlist$dirs["marker"],"g1000_eas_filtered.bed")
bimfiles <- file.path(GAlist$dirs["marker"],"g1000_eas_filtered.bim")
famfiles <- file.path(GAlist$dirs["marker"],"g1000_eas_filtered.fam")

# Prepare summary (i.e. Glist) for 1000G genotype data
Glist <- gprep(study="1000G EAS",
               bedfiles=bedfiles,
               bimfiles=bimfiles,
               famfiles=famfiles)

# Save Glist for use later
saveRDS(Glist, file=file.path(GAlist$dirs["marker"],"Glist_1000G_eas_filtered.rds"))

#####################################
# Process SAS
#####################################

# Define the file paths for the original bed/bim/fam files to read
bedRead <- file.path(GAlist$dirs["marker"],"g1000_sas.bed")
bimRead <- file.path(GAlist$dirs["marker"],"g1000_sas.bim")
famRead <- file.path(GAlist$dirs["marker"],"g1000_sas.fam")

# Define the file paths for the filtered bed/bim/fam files to write
bedWrite <- file.path(GAlist$dirs["marker"],"g1000_sas_filtered.bed")
bimWrite <- file.path(GAlist$dirs["marker"],"g1000_sas_filtered.bim")
famWrite <- file.path(GAlist$dirs["marker"],"g1000_sas_filtered.fam")

# Call the writeBED function to filter and write the data
writeBED(bedRead=bedRead, bimRead=bimRead, famRead=famRead,
         bedWrite=bedWrite, bimWrite=bimWrite, famWrite=famWrite,
         rsids=rsids)

# Define location of bed/bim/fam files
bedfiles <- file.path(GAlist$dirs["marker"],"g1000_sas_filtered.bed")
bimfiles <- file.path(GAlist$dirs["marker"],"g1000_sas_filtered.bim")
famfiles <- file.path(GAlist$dirs["marker"],"g1000_sas_filtered.fam")

# Prepare summary (i.e. Glist) for 1000G genotype data
Glist <- gprep(study="1000G SAS",
               bedfiles=bedfiles,
               bimfiles=bimfiles,
               famfiles=famfiles)

# Save Glist for use later
saveRDS(Glist, file=file.path(GAlist$dirs["marker"],"Glist_1000G_sas_filtered.rds"))







