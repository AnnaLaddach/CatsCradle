
library(Seurat)
library(ggplot2)
library(HandyPack)

rm(list=ls())
graphics.off()

source('../CradleWare.R')
source('../PCAWare.R')
source('AccessObjects.R')

## ####################################################
## ####################################################

assays = c('integrated','RNA')
res = 1
cats = c(cells='f',genes='fPrime')

for(assay in assays)
{
    pair = getObjectPair(assay)

    for(n in names(cats))
    {
        dim = 10
        outDir = c(paste0(assay,'_resolution_',res),
                   paste0(n,'CentroidPCSeparation'))
        outDir = nameAndMakeDir(outDir)

        plotList = makePCSeparationPlotEachVsEach(pair[[cats[n]]],dim,outDir)
    }
}
