
library(HandyPack)
library(ggplot2)
library(pheatmap)
library(stringr)

rm(list=ls())
graphics.off()

source('CradleWare.R')


assays = c('integrated','RNA')
geneSets = readLines('~/geneSets/mouse/h.all.v7.2.symbols.gmt')
geneSets = str_split(geneSets,'\t')
which = c('log','density')


for(assay in assays)
{
    for(whichFunction in which)
    {
        clusterDF = Read.Table(paste0(assay,'_resolution_2/geneClusters.txt'))
        background = 21281
        M = getGeneSetsVsClustersMatrix(geneSets,
                                        clusterDF,
                                        whichFunction,
                                        background)

        figName = paste0(assay,'_resolution_2/hallmarkHeatMap_',
                         whichFunction,'.jpg')
        jpeg(figName,
             height=12,width=12,units='in',res=200)
        pheatmap(M,
                 treeheight_row=0,
                 treeheight_col=0)
        dev.off()
    }

}
    
