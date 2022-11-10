########################################################################################################################################
# Script to create and insert data into GDT database
########################################################################################################################################

library(qgg)
library(data.table)
library(gact)

# The following steps are needed to create the GDT database:

################################################################################
# Create database infrastructure and GAlist
################################################################################
dbdir <- "C:/Users/au223366/Dropbox/Projects/balder/gdtdb"
GAlist <- gact(version="t2dm-gact-0.0.1", dbdir=dbdir, task="download")

# This step is only done once and the GAlist can be saved and used later
saveRDS(GAlist, file="C:/Users/au223366/Dropbox/Projects/balder/gdtdb/GAlist_t2dm-gact-0.0.1.rds")


################################################################################
# Process external GWAS summary statistic data and add to database
################################################################################

# Mahajan.NatGenet2018b.T2D-noUKBB.European.txt
fname_stat <- "C:/Users/au223366/Dropbox/Projects/balder/data/Mahajan.NatGenet2018b.T2D-noUKBB.European.txt"
stat <- fread(fname_stat, data.table=FALSE)
stat <- stat[, c("SNP","Chr","Pos","EA","NEA","EAF","Beta","SE","Pvalue")]
  colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p")

GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="Mahajan.NatGenet2018b.T2D-noUKBB.European.txt",
                       trait="T2DM",
                       type = "binary",
                       gender = "both",
                       reference = "PMID:30297969",
                       n = 898130,
                       ncase = 74124,
                       ncontrol = 824006)

# CARDIoGRAMplusC4D.txt.gz
fname_stat <- "C:/Users/au223366/Dropbox/Projects/balder/data/CARDIoGRAMplusC4D.txt.gz"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)
stat <- stat[, c(1:6,9:11)]
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p")

GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="CARDIoGRAMplusC4D.txt.gz",
                       trait="CAD",
                       type = "binary",
                       gender = "both",
                       reference = "PMID:26343387",
                       n = 184305,
                       ncase = 60801,
                       ncontrol = 123504)

# DIAMANTE-EUR.sumstat.txt.gz
fname_stat <- "C:/Users/au223366/Dropbox/Projects/balder/data/DIAMANTE-EUR.sumstat.txt.gz"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)
stat <- stat[, c(4,1:2,5:10)]
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p")

GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="DIAMANTE-EUR.sumstat.txt.gz",
                       trait="T2DM",
                       type = "binary",
                       gender = "both",
                       reference = "https://doi.org/10.1038/s41588-022-01058-3",
                       n = 180834,
                       ncase = 1159055,
                       ncontrol = 1339889)

# Mahajan.NatGenet2018b.T2D.MALE.European.txt
fname_stat <- "C:/Users/au223366/Dropbox/Projects/balder/data/Mahajan.NatGenet2018b.T2D.MALE.European.txt"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)
stat <- stat[, c("SNP","Chr","Pos","EA","NEA","EAF","Beta","SE","Pvalue")]
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p")

GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="Mahajan.NatGenet2018b.T2D.MALE.European.txt",
                       trait="T2DM",
                       type = "binary",
                       gender = "males",
                       reference = "PMID:30297969",
                       n = 898130,
                       ncase = 74124,
                       ncontrol = 824006)


# Mahajan.NatGenet2018b.T2D.FEMALE.European.txt

fname_stat <- "C:/Users/au223366/Dropbox/Projects/balder/data/Mahajan.NatGenet2018b.T2D.FEMALE.European.txt"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)
stat <- stat[, c("SNP","Chr","Pos","EA","NEA","EAF","Beta","SE","Pvalue")]
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p")

GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="Mahajan.NatGenet2018b.T2D.FEMALE.European.txt",
                       trait="T2DM",
                       type = "binary",
                       gender = "females",
                       reference = "PMID:30297969",
                       n = 898130,
                       ncase = 74124,
                       ncontrol = 824006)


# Mahajan.NatGenet2018b.T2D.European.txt

fname_stat <- "C:/Users/au223366/Dropbox/Projects/balder/data/Mahajan.NatGenet2018b.T2D.European.txt"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)
stat <- stat[, c("SNP","Chr","Pos","EA","NEA","EAF","Beta","SE","Pvalue")]
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p")

GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="Mahajan.NatGenet2018b.T2D.European.txt",
                       trait="T2DM",
                       type = "binary",
                       gender = "both",
                       reference = "PMID:30297969",
                       n = 898130,
                       ncase = 74124,
                       ncontrol = 824006)

# Save updated GAlist
saveRDS(GAlist, file="C:/Users/au223366/Dropbox/Projects/balder/gdtdb/GAlist_t2dm-gact-0.0.1.rds")

