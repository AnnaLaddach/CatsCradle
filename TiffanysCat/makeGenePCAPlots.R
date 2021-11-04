
library(Seurat)
library(ggplot2)
library(plotly)
library(cowplot)

rm(list=ls())
graphics.off()

source('AccessObjects.R')
source('../CradleWare.R')
## ####################################################
## ####################################################

cases = c('adult','combined')

dims = paste0('PC_',1:4)

for(case in cases)
{
    pair = getObjectPair(case)
    fPrime = pair$fPrime

    outDir = c('figures',paste0('PCA_',case))
    outDir = nameAndMakeDir(outDir)

    df = FetchData(fPrime,dims)
    df$gene = colnames(fPrime)
    df$cluster = fPrime$seurat_clusters

    g = plotPCs(df,
                dataCols=dims,
                coloring='cluster',
                labeling=c('gene','cluster'),
                outDir=outDir)

    dev.new()
    print(g)
}
