
library(HandyPack)

rm(list=ls())

source('../CradleWare.R')
## ####################################################
## ####################################################

assays = c('RNA','integrated')
resolutions = c(1,2)

for(assay in assays)
{
    for(res in resolutions)
    {
        fileIn = paste0(assay,
                        '_resolution_',
                        res,
                        '/geneClusters.txt')

        fileOut = paste0(assay,
                         '_resolution_',
                         res,
                         '/geneClusters_wide.txt')

        widenGeneClusterLists(fileIn,fileOut)
    }
}

    
