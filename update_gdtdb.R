########################################################################################################################################
# Script to create and insert data into GDT database
########################################################################################################################################

library(qgg)
library(data.table)
library(gact)

# The following steps are needed to create the GDT database:




################################################################################
# Continue processing external GWAS summary statistic data and add to database
################################################################################

# Load  GAlist
GAlist <- readRDS(file="C:/Users/au223366/Dropbox/Projects/balder/gdtdb/GAlist_t2dm-gact-0.0.1.rds")



# DIAMANTE-EAS.sumstat.txt.gz
fname_stat <- "C:/Users/au223366/Dropbox/Projects/balder/data/DIAMANTE-EAS.sumstat.txt.gz"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)
stat <- stat[, c(4,1:2,5:10)]
head(stat)
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p")

GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="DIAMANTE-EAS.sumstat.txt.gz",
                       trait="T2DM",
                       type = "binary",
                       gender = "both",
                       reference = "PMID:35551307",
                       n = 283423,
                       ncase = 56268,
                       ncontrol = 227155,
                       writeStatDB=TRUE)


# DIAMANTE-SAS.sumstat.txt.gz
fname_stat <- "C:/Users/au223366/Dropbox/Projects/balder/data/DIAMANTE-SAS.sumstat.txt.gz"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)
stat <- stat[, c(4,1:2,5:10)]
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p")

GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="DIAMANTE-SAS.sumstat.txt.gz",
                       trait="T2DM",
                       type = "binary",
                       gender = "both",
                       reference = "PMID:35551307",
                       n = 49492,
                       ncase = 16540,
                       ncontrol = 32952,
                       writeStatDB=TRUE)

# Save updated GAlist
#saveRDS(GAlist, file="C:/Users/au223366/Dropbox/Projects/balder/gdtdb/GAlist_t2dm-gact-0.0.1.rds")


################################################################################
# Continue processing external GWAS summary statistic data and add to database
################################################################################

#GAlist <- readRDS(file="C:/Users/au223366/Dropbox/Projects/balder/gdtdb/GAlist_t2dm-gact-0.0.1.rds")



# Meta-analysis of BMI in UK Biobank and GIANT data.
# Combined set of samples, max N = 697,734.

fname_stat <- "C:/Users/au223366/Dropbox/Projects/balder/data/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz"
#fname_stat <- "https://portals.broadinstitute.org/collaboration/giant/images/1/14/Bmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt.gz"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)
dim(stat)
stat <- na.omit(stat)
dim(stat)
stat <- stat[, c("SNP","CHR","POS","Tested_Allele","Other_Allele",
                 "Freq_Tested_Allele", "BETA", "SE", "P")]
head(stat)
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p")

if(any(!sapply(stat, class)[c("eaf", "b", "seb", "p")]%in%"numeric")) warning("Column classes wrong")
stat$p <- as.numeric(stat$p)
stat$p[stat$p==0] <- .Machine$double.xmin

GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz",
                       trait="BMI",
                       type = "quantitative",
                       gender = "both",
                       reference = "PMID:30239722",
                       n = 505454,
                       ncase = 0,
                       ncontrol = 505454,
                       writeStatDB=TRUE)



# Meta-analysis of whr in UK Biobank and GIANT data.

fname_stat <- "C:/Users/au223366/Dropbox/Projects/balder/data/whr.giant-ukbb.meta-analysis.combined.23May2018.txt.gz"
#fname_stat <- "https://portals.broadinstitute.org/collaboration/giant/images/6/6e/Whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt.gz"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)
dim(stat)
stat <- na.omit(stat)
dim(stat)
stat <- stat[, c("SNP","CHR","POS","Tested_Allele","Other_Allele",
                 "Freq_Tested_Allele", "BETA", "SE", "P")]
head(stat)
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p")

if(any(!sapply(stat, class)[c("eaf", "b", "seb", "p")]%in%"numeric")) warning("Column classes wrong")
stat$p <- as.numeric(stat$p)
stat$p[stat$p==0] <- .Machine$double.xmin

GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="whr.giant-ukbb.meta-analysis.combined.23May2018.txt.gz",
                       trait="WHR",
                       type = "quantitative",
                       gender = "both",
                       reference = "PMID:30239722",
                       n = 498318,
                       ncase = 0,
                       ncontrol = 498318,
                       writeStatDB=TRUE)


# Save updated GAlist
saveRDS(GAlist, file="C:/Users/au223366/Dropbox/Projects/balder/gdtdb/GAlist_t2dm-gact-0.0.1.rds")




################################################################################
# Continue processing external GWAS summary statistic data and add to database
################################################################################

#GAlist <- readRDS(file="C:/Users/au223366/Dropbox/Projects/balder/gdtdb/GAlist_t2dm-gact-0.0.1.rds")



# CKD

fname_stat <- "C:/Users/au223366/Dropbox/Projects/balder/data/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)
dim(stat)
stat <- na.omit(stat)
dim(stat)
stat <- stat[, c("RSID","Chr","Pos_b37","Allele1","Allele2", "Freq1",
                 "Effect", "StdErr", "P-value")]

colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf","b", "seb", "p")

if(any(!sapply(stat, class)[c("eaf", "b", "seb", "p")]%in%"numeric")) warning("Column classes wrong")
stat$p <- as.numeric(stat$p)
stat$p[stat$p==0] <- .Machine$double.xmin

GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz",
                       trait="CKD",
                       type = "quantitative",
                       gender = "both",
                       reference = "PMID:31152163",
                       n = 1046070,
                       ncase = 0,
                       ncontrol = 1046070,
                       writeStatDB=TRUE)



# Atrial Fibrillation

fname_stat <- "C:/Users/au223366/Dropbox/Projects/balder/data/AF_HRC_GWAS_ALLv11.txt"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)
dim(stat)
stat <- na.omit(stat)
dim(stat)
stat <- stat[, c("MarkerName","chr","pos","Allele1","Allele2",
                 "Effect", "StdErr", "P-value")]
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "b", "seb", "p")

if(any(!sapply(stat, class)[c("b", "seb", "p")]%in%"numeric")) warning("Column classes wrong")
stat$p <- as.numeric(stat$p)
stat$p[stat$p==0] <- .Machine$double.xmin

GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="AF_HRC_GWAS_ALLv11.txt",
                       trait="AF",
                       type = "binary",
                       gender = "both",
                       reference = "PMID:29892015",
                       n = 588190,
                       ncase = 65446,
                       ncontrol = 522744,
                       writeStatDB=TRUE)

# Save updated GAlist
saveRDS(GAlist, file="C:/Users/au223366/Dropbox/Projects/balder/gdtdb/GAlist_t2dm-gact-0.0.1.rds")








# Meta-analysis of whr in UK Biobank and GIANT data.
# Female samples only, max N = 381,152.

fname_stat <- "https://portals.broadinstitute.org/collaboration/giant/images/6/62/Bmi.giant-ukbb.meta-analysis.females.23May2018.HapMap2_only.txt.gz"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)
stat <- stat[, c("SNP","CHR","POS","Tested_Allele","Other_Allele",
                 "Freq_Tested_Allele", "BETA", "SE", "P")]
head(stat)
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p")

if(any(!sapply(stat, class)[c("eaf", "b", "seb", "p")]%in%"numeric")) warning("Column classes wrong")
stat$p <- as.numeric(stat$p)
stat$p[stat$p==0] <- .Machine$double.xmin

GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="Bmi.giant-ukbb.meta-analysis.females.23May2018.HapMap2_only.txt.gz",
                       trait="BMI",
                       type = "quantitative",
                       gender = "females",
                       reference = "PMID:30239722",
                       n = 381152,
                       ncase = 0,
                       ncontrol = 381152,
                       writeStatDB=TRUE)


# Meta-analysis of whr in UK Biobank and GIANT data.
# Male samples only, max N = 316,772.

fname_stat <- "https://portals.broadinstitute.org/collaboration/giant/images/a/ab/Bmi.giant-ukbb.meta-analysis.males.23May2018.HapMap2_only.txt.gz"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)
stat <- stat[, c("SNP","CHR","POS","Tested_Allele","Other_Allele",
                 "Freq_Tested_Allele", "BETA", "SE", "P")]
head(stat)
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p")

GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="Bmi.giant-ukbb.meta-analysis.males.23May2018.HapMap2_only.txt.gz",
                       trait="BMI",
                       type = "quantitative",
                       gender = "males",
                       reference = "PMID:30239722",
                       n = 316772,
                       ncase = 0,
                       ncontrol = 316772,
                       writeStatDB=TRUE)


# Meta-analysis of whr in UK Biobank and GIANT data.

fname_stat <- "https://portals.broadinstitute.org/collaboration/giant/images/6/6e/Whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt.gz"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)
stat <- stat[, c("SNP","CHR","POS","Tested_Allele","Other_Allele",
                 "Freq_Tested_Allele", "BETA", "SE", "P")]
head(stat)
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p")

if(any(!sapply(stat, class)[c("eaf", "b", "seb", "p")]%in%"numeric")) warning("Column classes wrong")
stat$p <- as.numeric(stat$p)
stat$p[stat$p==0] <- .Machine$double.xmin

GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="Whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt.gz",
                       trait="WHR",
                       type = "quantitative",
                       gender = "both",
                       reference = "PMID:30239722",
                       n = 697734,
                       ncase = 0,
                       ncontrol = 697734,
                       writeStatDB=TRUE)


# Meta-analysis of whr in UK Biobank and GIANT data.
# Female samples only, max N = 381,152.

fname_stat <- "https://portals.broadinstitute.org/collaboration/giant/images/9/97/Whradjbmi.giant-ukbb.meta-analysis.females.23May2018.HapMap2_only.txt.gz"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)
stat <- stat[, c("SNP","CHR","POS","Tested_Allele","Other_Allele",
                 "Freq_Tested_Allele", "BETA", "SE", "P")]
head(stat)
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p")

if(any(!sapply(stat, class)[c("eaf", "b", "seb", "p")]%in%"numeric")) warning("Column classes wrong")
stat$p <- as.numeric(stat$p)
stat$p[stat$p==0] <- .Machine$double.xmin

GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="Whradjbmi.giant-ukbb.meta-analysis.females.23May2018.HapMap2_only.txt.gz",
                       trait="WHR",
                       type = "quantitative",
                       gender = "females",
                       reference = "PMID:30239722",
                       n = 381152,
                       ncase = 0,
                       ncontrol = 381152,
                       writeStatDB=TRUE)


# Meta-analysis of whr in UK Biobank and GIANT data.
# Male samples only, max N = 316,772.

fname_stat <- "https://portals.broadinstitute.org/collaboration/giant/images/7/71/Whradjbmi.giant-ukbb.meta-analysis.males.23May2018.HapMap2_only.txt.gz"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)
stat <- stat[, c("SNP","CHR","POS","Tested_Allele","Other_Allele",
                 "Freq_Tested_Allele", "BETA", "SE", "P")]
head(stat)
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p")

if(any(!sapply(stat, class)[c("eaf", "b", "seb", "p")]%in%"numeric")) warning("Column classes wrong")
stat$p <- as.numeric(stat$p)
stat$p[stat$p==0] <- .Machine$double.xmin

GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="Whradjbmi.giant-ukbb.meta-analysis.males.23May2018.HapMap2_only.txt.gz",
                       trait="WHR",
                       type = "quantitative",
                       gender = "males",
                       reference = "PMID:30239722",
                       n = 316772,
                       ncase = 0,
                       ncontrol = 316772,
                       writeStatDB=TRUE)


# Save updated GAlist
saveRDS(GAlist, file="C:/Users/au223366/Dropbox/Projects/balder/gdtdb/GAlist_t2dm-gact-0.0.1.rds")

