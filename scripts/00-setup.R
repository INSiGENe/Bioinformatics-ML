# Set up project

#--- Install renv for version tracking
install.packages("renv")

#--- Initialize renv
renv::init(bioconductor = TRUE)

#--- Install packages
renv::install(c("tidymodels",
                "rpart.plot",
                "vip",
                "GEOquery",
                "tidyverse",
                "viridis",
                "gt"))

#--- Snapshot
renv::snapshot()
