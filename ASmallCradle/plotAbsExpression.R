
library(Seurat)
library(ggplot2)

rm(list=ls())
graphics.off()

source('~/FranzeSingleCell07/SeuratWare02.R')
source('CradleWare.R')

## ###################################################
## ###################################################

assays = c('integrated','RNA')

for(assay in assays)
{
    objectPair = getObjectPair(assay)
    f = objectPair$f
    fPrime = objectPair$fPrime

    df = getAbsoluteExpressionTotals(f,fPrime)
    
    g = ggplot(df,aes(x=geneCluster,y=shortName,fill=expression)) +
        geom_tile() +
        scale_fill_viridis_c() +
        ggtitle(paste('Absolute expression totals',assay))

    dev.new()
    print(g)
}


