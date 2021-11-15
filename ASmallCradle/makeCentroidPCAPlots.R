
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

assays = c('integrated','RNA')
res = 1

objects = c('fPrime','f')
tag = c()
tag['fPrime'] = 'genes'
tag['f'] = 'cells'
dim = 5

for(assay in assays)
{
    pair = getObjectPair(assay)

    for(object in objects)
    {
        fThis = pair[[object]]
        clusters = unique(fThis$seurat_clusters)
        clusters = clusters[order(as.numeric(clusters))]
        fThis$seurat_clusters = factor(fThis$seurat_clusters,
                                       levels=clusters)
                                       

        outDir = c(paste0(assay,'_resolution_',res),
                   paste0('PCA_',tag[object],'_centroid',dim))
        outDir = nameAndMakeDir(outDir)

        ## Get PCDF:
        PCADF = makePCDF(fThis,dim)

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
