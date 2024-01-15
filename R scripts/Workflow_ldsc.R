################################################################################
# LD score regression
################################################################################
#
#library(devtools)
#devtools::install_github("psoerensen/gact")
#devtools::install_github("psoerensen/qgg")
#
################################################################################


library(qgg)
library(gact)
library(ggcorrplot)


# Load GAlist with information on gact database
GAlist <- readRDS(file="/faststorage/project/ukbiobank/projects/gact/hsa.0.0.1/GAlist_hsa.0.0.1.rds")

# Check studies in gact database
GAlist$studies

# Get ldscores
ldscores <- getLDscoresDB(GAlist=GAlist, ancestry="EUR")

# Select study
studyIDs <- paste0("GWAS",c(1:2, 4:6))

# Get GWAS summary statistics from gact database
z <- getMarkerStatDB(GAlist=GAlist, studyID=studyIDs, what="z")
dim(z)

# Get number of samples for each GWAS study
n <- GAlist$studies[studyIDs,"neff"]
names(n) <- studyIDs

# Estimate narrow sense heritability using ldsc
fit <- ldsc(z=z, ldscores=ldscores, n=n, what="h2")
fit

# Estimate genetic correlation using using ldsc
fit <- ldsc(z=z, ldscores=ldscores, n=n, what="rg")
fit$h2
fit$rg

ggcorrplot::ggcorrplot(fit$rg,
                       title = "Genetic correlations",
                       hc.order = TRUE,
                       ggtheme = theme_bw(),
                       colors = c("#D73027", "#FFFFFF", "#1A9850"),
                       tl.cex = 10,
                       tl.srt = 90) +
  theme(axis.text.x = element_text(vjust = 0.5))

