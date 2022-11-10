####################################################################
# Prepare Glist for 1000G
####################################################################

#download.file(url="https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz",dest="./1000G_Phase3_plinkfiles.tgz")
#cmd <- "tar -xvzf 1000G_Phase3_plinkfiles.tgz"
#system(cmd)


library(qgg)
library(data.table)

bedfiles <- paste("C:\\Users\\au223366\\Dropbox\\Projects\\1000G\\1000G_EUR_Phase3_plink\\1000G.EUR.QC.",1:22,".bed",sep="")
bimfiles <- paste("C:\\Users\\au223366\\Dropbox\\Projects\\1000G\\1000G_EUR_Phase3_plink\\1000G.EUR.QC.",1:22,".bim",sep="")
famfiles <- paste("C:\\Users\\au223366\\Dropbox\\Projects\\1000G\\1000G_EUR_Phase3_plink\\1000G.EUR.QC.",1:22,".fam",sep="")

# Process bed/bim/fam files and gather information about samples and genotypes
Glist <- gprep(study="1000G", bedfiles=bedfiles, bimfiles=bimfiles, famfiles=famfiles)


# Compute sparse LD matrices and LD scores
# only use QC'ed markers for computation of sparse LD matrices

rsidsLD <-  gfilter(Glist = Glist,
                    excludeMAF=0.01,
                    excludeMISS=0.05,
                    excludeHWE=1e-12,
                    excludeMHC=TRUE,
                    excludeCGAT=TRUE,
                    excludeINDEL=TRUE,
                    excludeDUPS=TRUE,
                    assembly="GRCh37")

ldfiles <- paste("C:\\Users\\au223366\\Dropbox\\Projects\\1000G\\1000G_EUR_Phase3_plink\\LD_1000G.EUR.QC.",1:22,".ld",sep="")
Glist <- gprep(Glist, task="sparseld",msize=1000, rsids=rsidsLD, ldfiles=ldfiles, overwrite=TRUE)

saveRDS(Glist, file="C:\\Users\\au223366\\Dropbox\\Projects\\1000G\\Glist_MAF1_LD1K_1000G.rds")


