
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

    df = getExpressionTotals(f,fPrime)
    
    g = ggplot(df,aes(x=geneCluster,y=shortName,fill=expression)) +
        geom_tile() +
        scale_fill_viridis_c() +
        ggtitle(paste('Expression totals',assay))

    dev.new()
    print(g)
    fileName = paste0(assay,
                      '_resolution_2/expressionTotals_',
                      assay,
                      '.jpg')
    ggsave(plot=g,
           filename=fileName,
           height=8,width=12,units='in')
}


