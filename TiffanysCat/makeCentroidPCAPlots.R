
library(Seurat)
library(ggplot2)
library(plotly)
library(cowplot)

rm(list=ls())
graphics.off()

source('AccessObjects.R')
source('../CradleWare.R')
source('../PCAWare.R')
## ####################################################
## ####################################################

cases = c('adult','combined')

for(case in cases)
{
    pair = getObjectPair(case)

    objects = c('fPrime','f')
    tag = c()
    tag['fPrime'] = 'genes'
    tag['f'] = 'cells'
    dim = 5

    for(object in objects)
    {
        fThis = pair[[object]]
        outDir = c('figures',
                   paste0('PCA_',tag[object],'_',case,'_centroid',dim))
        outDir = nameAndMakeDir(outDir)

        ## Get the PCDF:
        if(object == 'fPrime')
        {
            PCADF = makeGenePCDF(fThis,dim)
        } else if(object == 'f') {
            PCADF = makePCDF(fThis,dim)
        }

        ## Get the centroid DF
        centroidDF = makeCentroidDF(PCADF)

        ## Do the plotting:
        dataCols = paste0('PC_',1:dim)
        coloring = 'seurat_clusters'
        labeling = 'seurat_clusters'
        outDir = outDir

        g = plotPCs(centroidDF,dataCols,coloring,labeling,outDir)
        dev.new()
        print(g)
    }
}
