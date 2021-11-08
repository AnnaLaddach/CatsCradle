
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
    fPrime = pair$fPrime

    outDir = c('figures',paste0('PCA_genes_',case))
    outDir = nameAndMakeDir(outDir)

    df = makeGenePCDF(fPrime,4)
    dims = paste0('PC_',1:4)

    g = plotPCs(df,
                dataCols=dims,
                coloring='seurat_clusters',
                labeling=c('gene','seurat_clusters'),
                outDir=outDir)
    
    dev.new()
    print(g)
}
