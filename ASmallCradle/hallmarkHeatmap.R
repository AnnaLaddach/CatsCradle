
library(HandyPack)
library(tictoc)
library(ggplot2)
library(pheatmap)
library(stringr)

rm(list=ls())
graphics.off()

source('../CradleWare.R')
source('AccessObjects.R')


assays = c('integrated','RNA')
geneSets = readLines('~/geneSets/mouse/h.all.v7.2.symbols.gmt')
geneSets = str_split(geneSets,'\t')
which = c('log','density')

for(res in 1:2)
{
    for(assay in assays)
    {
        pair = getObjectPair(assay=assay,res=res)
        f = pair$f
        
        for(whichFunction in which)
        {
            Tic(paste('hallmark',assay,res))
                
            clusterDF = Read.Table(paste0(assay,'_resolution_',res,'/geneClusters.txt'))
            names(clusterDF)[2] = 'geneCluster'
            names(clusterDF)[1] = 'gene'
            
            M = getGeneSetsVsClustersMatrix(geneSets,
                                            clusterDF,
                                            whichFunction,
                                            backgroundGenes=rownames(f))
            
            if(whichFunction == 'log')
                M = log10(M+1)

            figName = paste0(assay,'_resolution_',res,'/hallmarkHeatMap_',
                             whichFunction,'.jpg')
            jpeg(figName,
                 height=12,width=12,units='in',res=200)
            pheatmap(M,
                     treeheight_row=0,
                     treeheight_col=0)
            dev.off()
            toc()
        }
    }
}
    
