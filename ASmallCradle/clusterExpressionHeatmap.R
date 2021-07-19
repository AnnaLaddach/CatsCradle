
library(Seurat)
library(pheatmap)
library(HandyPack)
library(tictoc)

rm(list=ls())
source('CradleWare.R')

for(res in 1:2)
{
    for(assay in c('integrated','RNA'))
    {
        Tic(paste(assay,res))
        objectPair = getObjectPair(assay,res)
        f = objectPair$f
        fPrime = objectPair$fPrime

        M = getExpressionTotalsMatrix(f,fPrime)

        fileName = paste0(assay,
                          '_resolution_',
                          res,
                          '/clusterExpressionHeatmap.jpg')

        jpeg(fileName,
             height=8,width=8,units='in',res=100)
        pheatmap(M,
                 treeheight_row=0,
                 treeheight_col=0)
        dev.off()

        
        toc()
    }
}
