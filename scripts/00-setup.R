# Set up project

#--- Install renv for version tracking
install.packages("renv")

#--- Initialize renv
renv::init(bioconductor = TRUE)

#--- Install packages
renv::install(c("tidymodels",
                "vip",
                "tidyverse",
                "TCGAbiolinks",
                "DESeq2",
                "xgboost"))

#--- Snapshot
renv::snapshot()
