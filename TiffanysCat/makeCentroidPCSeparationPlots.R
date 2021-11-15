
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

cases = c('adult','combined')
cats = c(cells='f',genes='fPrime')

for(case in cases)
{
    pair = getObjectPair(case)
    for(n in names(cats))
    {
        dim = 10
        outDir = c('figures',
                   paste0(case,
                          '_',
                          n,
                          'CentroidPCSeparation'))
        outDir = nameAndMakeDir(outDir)
        
        plotList = makePCSeparationPlotEachVsEach(pair[[cats[n]]],dim,outDir)
    }
}

