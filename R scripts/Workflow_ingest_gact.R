################################################################################
# Ingest new data into gact database
################################################################################
#
#library(devtools)
#devtools::install_github("psoerensen/gact")
#devtools::install_github("psoerensen/qgg")


# Load libraries
library(qgg)
library(gact)
library(data.table)


# Load GAlist
GAlist <- readRDS(file="/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds")



################################################################################
# UK Biobank 2021 HbA1c GWAS: European ancestry
################################################################################
# link:

# Load GWAS data
fname_stat <- "/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/download/Glycated_haemoglobin_HbA1c.imp.gz"
stat <- fread(fname_stat, data.table = FALSE)
head(stat)

# Modify columns according to required format


# Subset and rename columns according to required format
stat <- stat[, c("MarkerName", "#CHROM", "POS", "ALT", "REF", "Effect", "StdErr", "P-value")]
colnames(stat) <- c("marker", "chr", "pos", "ea", "nea", "b", "seb", "p")


# Update database
GAlist <- updateStatDB(GAlist = GAlist,
                       stat = stat,
                       source = "Glycated_haemoglobin_HbA1c.imp",
                       trait = "HbA1C",
                       type = "quantitative",
                       gender = "both",
                       ancestry = "EUR",
                       build = "GRCh37",
                       reference = "PMID:33462484",
                       n = 318779,
                       ncase = 0,
                       ncontrol = 0,
                       comments = "Just UK biobank",
                       writeStatDB = TRUE)

# Save updated database
saveRDS(GAlist, file = "/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds", compress = FALSE)


################################################################################
# Random glucose 2023 GWAS: European ancestry
################################################################################
# link: https://t2d.hugeamp.org/dinspector.html?dataset=Qiao2023_BloodGlucose_EU


# Load GWAS data
fname_stat <- "/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/download/random_glucose.txt"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)

# Modify columns according to required format


# Subset and rename columns according to required format
stat <- stat[, c("SNP","CHR","BP","A1","A2", "freq", "b", "se", "P", "N")]
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "n")

# Update database
GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="random_glucose.txt",
                       trait="Random glucose",
                       type = "quantitative",
                       gender = "both",
                       ancestry = "EUR",
                       build = "GRCh37",
                       reference = "PMID: 36707517",
                       n = 	518615,
                       ncase = 0,
                       ncontrol = 518615,
                       comments ="Just UK biobank",
                       writeStatDB=TRUE)


saveRDS(GAlist, file="/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds", compress = FALSE)


################################################################################
# EHR-defined NAFLD 2021 GWAS: European ancestry
################################################################################
# link: https://t2d.hugeamp.org/dinspector.html?dataset=Fairfield2022_NAFLD_EU

# Load GWAS data
fname_stat <- "/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/download/GCST90054782_buildGRCh37.tsv"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)

# Modify columns according to required format
beta = log(stat$odds_ratio)
stat$beta <- beta
se = sqrt(((beta)^2)/qchisq(stat$p_value,1,lower.tail=F))
stat$se <- se

# Subset and rename columns according to required format
stat <- stat[, c("variant_id","chromosome","base_pair_location","effect_allele","other_allele","effect_allele_frequency",
                 "beta", "se", "p_value", "info")]
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "info")


# Update database
GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="GCST90054782_buildGRCh37.tsv",
                       trait="NAFLD",
                       type = "binary",
                       gender = "both",
                       ancestry = "EUR",
                       build = "GRCh37",
                       reference = "PMID: 34535985",
                       n = 377998,
                       ncase = 4761,
                       ncontrol = 373227,
                       comments ="Include UK biobank",
                       writeStatDB=TRUE)


saveRDS(GAlist, file="/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds", compress = FALSE)

################################################################################
# CKDGen 2019 GWAS
################################################################################
# link: https://t2d.hugeamp.org/dinspector.html?dataset=GWAS_CKDGenConsortium


# Load GWAS data
fname_stat <- "/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/download/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)


# Modify columns according to required format

# Subset and rename columns according to required format
stat <- stat[, c("RSID","Chr","Pos_b37","Allele1","Allele2","Freq1","Effect", "StdErr", "P-value", "n_total_sum")]
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "n")


# Update database
GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="CKD_overall_EA_JW_20180223_nstud23.dbgap.txt",
                       trait="CKD",
                       type = "binary",
                       gender = "both",
                       ancestry = "EUR",
                       build = "GRCh37",
                       reference = "PMID: 31152163",
                       n = 480698,
                       ncase = 41395,
                       ncontrol = 439303,
                       comments ="Include UK biobank",
                       writeStatDB=TRUE)

saveRDS(GAlist, file="/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds", compress = FALSE)

################################################################################
# Heart rate and hypertension GWAS
################################################################################
# link: https://t2d.hugeamp.org/dinspector.html?dataset=Zhu2019_COPD_CVD_eu

# Load GWAS data
fname_stat <- "/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/download/ZhuZ_30940143_ukbb.bolt_460K_selfRepWhite.doctor_highbloodpressure.assoc.gz"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)


# Modify columns according to required format

# Subset and rename columns according to required format
stat <- stat[, c("SNP","CHR","BP","A1","A0", "BETA", "SE", "P", "INFO")]
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "b", "seb", "p", "info")

# Update database
GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="ZhuZ_30940143_ukbb.bolt_460K_selfRepWhite.doctor_highbloodpressure.assoc",
                       trait="HTN",
                       type = "binary",
                       gender = "both",
                       ancestry = "EUR",
                       build = "GRCh37",
                       reference = "PMID:30940143",
                       n = 458554,
                       ncase = 144793,
                       ncontrol = 313761,
                       comments ="Include UK biobank",
                       writeStatDB=TRUE)

saveRDS(GAlist, file="/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds", compress = FALSE)

################################################################################
#  Meta-analysis of body mass index (bmi) in UK Biobank and GIANT data.
################################################################################
# link: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6298238/

# Load GWAS data
fname_stat <- "/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/download/Bmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt.gz"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)

# Modify columns according to required format

# Subset and rename columns according to required format
stat <- stat[, c("SNP","CHR","POS","Tested_Allele","Other_Allele",
                 "Freq_Tested_Allele","BETA", "SE", "P", "N","INFO")]
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "n", "info")

# Update database
GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="Bmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt",
                       trait="BMI",
                       type = "quantitative",
                       gender = "both",
                       ancestry = "EUR",
                       build = "GRCh37",
                       reference = "PMID:30239722",
                       n = 806834,
                       ncase = 0,
                       ncontrol = 806834,
                       comments ="Include UK biobank",
                       writeStatDB=TRUE)
saveRDS(GAlist, file="/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds", compress = FALSE)


################################################################################
# CARDIoGRAMplusC4D-UK Biobank CAD 2022 GWAS Meta-analysis
################################################################################
# link: https://t2d.hugeamp.org/dinspector.html?dataset=GWAS_CADMETA_eu

# Load GWAS data
fname_stat <- "/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/download/GCST90132314_buildGRCh37.tsv"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)

# Modify columns according to required format

# Subset and rename columns according to required format
stat <- stat[, c("markername","chromosome","base_pair_location","effect_allele","other_allele","effect_allele_frequency","beta", "standard_error", "p_value")]
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p")

# Update database
GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="CAD_META.gz",
                       trait="CAD",
                       type = "binary",
                       gender = "both",
                       ancestry = "EUR",
                       build = "GRCh37",
                       reference = "PMID:29212778",
                       n = 547261,
                       ncase = 122733,
                       ncontrol = 424528,
                       comments ="Include UK biobank",
                       writeStatDB=TRUE)

saveRDS(GAlist, file="/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds", compress = FALSE)


################################################################################
# Meta-analysis of whr in UK Biobank and GIANT data.
################################################################################
# link: https://t2d.hugeamp.org/dinspector.html?dataset=GWAS_UKBiobankGIANT_eu

# Load GWAS data
fname_stat <- "/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/download/Whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt.gz"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)

# Modify columns according to required format

# Subset and rename columns according to required format
stat <- stat[, c("SNP","CHR","POS","Tested_Allele","Other_Allele",
                 "Freq_Tested_Allele","BETA", "SE", "P", "N","INFO")]
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "n", "info")

# Update database
GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="Whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt",
                       trait="WHR",
                       type = "quantitative",
                       gender = "both",
                       ancestry = "EUR",
                       build = "GRCh37",
                       reference = "PMID:30239722",
                       n = 694649,
                       ncase = 0,
                       ncontrol = 694649,
                       comments ="Include UK biobank",
                       writeStatDB=TRUE)

saveRDS(GAlist, file="/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds", compress = FALSE)

################################################################################
# Standing Height
################################################################################
# link:

# Load GWAS data
fname_stat <- "/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/download/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)

# Modify columns according to required format

# Subset and rename columns according to required format
stat <- stat[, c("RSID","CHR","POS","EFFECT_ALLELE","OTHER_ALLELE",
                 "EFFECT_ALLELE_FREQ","BETA", "SE", "P", "N")]
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "n")

# Update database
GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="Whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt",
                       trait="Height",
                       type = "quantitative",
                       gender = "both",
                       ancestry = "EUR",
                       build = "GRCh37",
                       reference = "PMID:36224396",
                       n = 4080687,
                       ncase = 0,
                       ncontrol = 4080687,
                       comments ="Include UK biobank",
                       writeStatDB=TRUE)

saveRDS(GAlist, file="/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds", compress = FALSE)

################################################################################
# MEGASTROKE GWAS: Europeans
################################################################################
# link: https://t2d.hugeamp.org/dinspector.html?dataset=GWAS_MEGASTROKE_eu

# Load GWAS data
fname_stat <- "/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/download/metastroke.all.chr.bp"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)

# Modify columns according to required format
stat$Allele1 <- toupper(stat$Allele1)
stat$Allele2 <- toupper(stat$Allele2)

# Subset and rename columns according to required format
stat <- stat[, c("MarkerName","chromosome","position","Allele1","Allele2",
                 "Freq1","Effect", "StdErr", "P.value")]
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p")

# Update database
GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="metastroke.all.chr.bp",
                       trait="Stroke",
                       type = "binary",
                       gender = "both",
                       ancestry = "EUR",
                       build = "GRCh37",
                       reference = "https://doi.org/10.1038/s41588-018-0058-3",
                       n = 446696,
                       ncase = 40585,
                       ncontrol = 406111,
                       comments ="Include UK biobank",
                       writeStatDB=TRUE)

saveRDS(GAlist, file="/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds", compress = FALSE)


################################################################################
# Dyslipidemia
################################################################################
# link: https://t2d.hugeamp.org/dinspector.html?dataset=GWAS_GERA_eu

# Load GWAS data
fname_stat <- "/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/download/Results_DYSLIPID_GERA_chr_1_to_22.txt.gz"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)

# Modify columns according to required format

# Subset and rename columns according to required format
stat <- stat[, c("rs_id_all","chr","position","alleleB","alleleA",
                 "freq_effect_allele","frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue", "all_total")]
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "n")

# Update database
GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="Results_DYSLIPID_GERA_chr_1_to_22.txt.gz",
                       trait="DLP",
                       type = "binary",
                       gender = "both",
                       ancestry = "EUR",
                       build = "GRCh37",
                       reference = "PMID:33893285",
                       n = 56637,
                       ncase = 30244,
                       ncontrol = 26393,
                       comments ="Include UK biobank",
                       writeStatDB=TRUE)

saveRDS(GAlist, file="/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds", compress = FALSE)


################################################################################
# CVD
################################################################################
# link: https://t2d.hugeamp.org/dinspector.html?dataset=GWAS_GERA_eu

# Load GWAS data
fname_stat <- "/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/download/Results_CARD_GERA_chr_1_to_22.txt.gz"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)

# Modify columns according to required format

# Subset and rename columns according to required format
stat <- stat[, c("rs_id_all","chr","position","alleleB","alleleA",
                 "freq_effect_allele","frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue", "all_total")]
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "n")

# Update database
GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="Results_CARD_GERA_chr_1_to_22.txt.txt.gz",
                       trait="CVD",
                       type = "binary",
                       gender = "both",
                       ancestry = "EUR",
                       build = "GRCh37",
                       reference = "PMID:33893285",
                       n = 56637,
                       ncase = 15009,
                       ncontrol = 41628,
                       comments ="Include UK biobank",
                       writeStatDB=TRUE)

saveRDS(GAlist, file="/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds", compress = FALSE)


#######################################################################################
# GLGC 2021 Lipids GWAS
################################################################################
# link: https://t2d.hugeamp.org/dinspector.html?dataset=Graham_2021_lipids_Mixed

# Load GWAS data
fname_stat <- "/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/download/HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)

# Modify columns according to required format

# Subset and rename columns according to required format
stat <- stat[, c("rsID","CHROM","POS_b37","ALT","REF", "POOLED_ALT_AF", "EFFECT_SIZE", "SE", "pvalue", "N")]
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "n")

# Update database
GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results",
                       trait="HDL",
                       type = "quantitative",
                       gender = "both",
                       ancestry = "EUR",
                       build = "GRCh37",
                       reference = "PMID:34887591",
                       n = 1320016,
                       ncase = 0,
                       ncontrol = 1320016,
                       comments ="Include UK biobank",
                       writeStatDB=TRUE)
saveRDS(GAlist, file="/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds", compress = FALSE)


################################################################################
# LDL
################################################################################
# link:

# Load database
fname_stat <- "/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/download/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)

# Modify columns according to required format

# Subset and rename columns according to required format
stat <- stat[, c("rsID","CHROM","POS_b37","ALT","REF", "POOLED_ALT_AF", "EFFECT_SIZE", "SE", "pvalue", "N")]
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "n")

# Update database
GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results",
                       trait="LDL",
                       type = "quantitative",
                       gender = "both",
                       ancestry = "EUR",
                       build = "GRCh37",
                       reference = "PMID:34887591",
                       n = 1320016,
                       ncase = 0,
                       ncontrol = 1320016,
                       comments ="Include UK biobank",
                       writeStatDB=TRUE)
saveRDS(GAlist, file="/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds", compress = FALSE)


################################################################################
# TC
################################################################################
# link:

# Load database
fname_stat <- "/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/download/TC_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)

# Modify columns according to required format

# Subset and rename columns according to required format
stat <- stat[, c("rsID","CHROM","POS_b37","ALT","REF", "POOLED_ALT_AF", "EFFECT_SIZE", "SE", "pvalue", "N")]
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "n")

# Update database
GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="TC_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results",
                       trait="TC",
                       type = "quantitative",
                       gender = "both",
                       ancestry = "EUR",
                       build = "GRCh37",
                       reference = "PMID:34887591",
                       n = 1320016,
                       ncase = 0,
                       ncontrol = 1320016,
                       comments ="Include UK biobank",
                       writeStatDB=TRUE)

saveRDS(GAlist, file="/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds", compress = FALSE)


################################################################################
# TG
################################################################################
# link:

# Load GWAS data
fname_stat <- "/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/download/logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)


# Modify columns according to required format

# Subset and rename columns according to required format
stat <- stat[, c("rsID","CHROM","POS_b37","ALT","REF", "POOLED_ALT_AF", "EFFECT_SIZE", "SE", "pvalue", "N")]
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "n")

# Update database
GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results",
                       trait="TG",
                       type = "quantitative",
                       gender = "both",
                       ancestry = "EUR",
                       build = "GRCh37",
                       reference = "PMID:34887591",
                       n = 1320016,
                       ncase = 0,
                       ncontrol = 1320016,
                       comments ="Include UK biobank",
                       writeStatDB=TRUE)

saveRDS(GAlist, file="/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds", compress = FALSE)


################################################################################
# Systolic blood presure
################################################################################
# link: https://t2d.hugeamp.org/dinspector.html?dataset=Evangelou2018_bp_eu

# Load GWAS data
fname_stat <- "/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/download/Evangelou_30224653_SBP.txt.gz"
stat <- fread(fname_stat, data.table=FALSE)
head(stat)

# Modify columns according to required format
cp <- strsplit(stat$MarkerName, ":")
stat$chr <- sapply(cp, function(x){x[1]})
stat$pos <- as.numeric(sapply(cp, function(x){x[2]}))
stat$Allele1 <- toupper(stat$Allele1)
stat$Allele2 <- toupper(stat$Allele2)

# Subset and rename columns according to required format
stat <- stat[, c("MarkerName", "chr","pos", "Allele1","Allele2","Freq1","Effect", "StdErr", "P", "N_effective")]
colnames(stat) <- c("marker","chr", "pos", "ea", "nea", "eaf", "b", "seb", "p", "n")

# Update database
GAlist <- updateStatDB(GAlist=GAlist,
                       stat=stat,
                       source="Evangelou_30224653_SBP.txt",
                       trait="SBP",
                       type = "quantitative",
                       gender = "both",
                       ancestry = "EUR",
                       build = "GRCh37",
                       reference = "PMID:302246538",
                       n = 757601,
                       ncase = 0,
                       ncontrol = 757601,
                       comments ="Include UK biobank",
                       writeStatDB=TRUE)

saveRDS(GAlist, file="/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds", compress = FALSE)


