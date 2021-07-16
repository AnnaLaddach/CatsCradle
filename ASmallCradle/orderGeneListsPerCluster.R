
library(HandyPack)
library(stringr)

rm(list=ls())

source('CradleWare.R')

## ###################################################
## ###################################################


assays = c('integrated','RNA')
whichFunction = 'log'

geneSets = c(hallmark='~/geneSets/mouse/h.all.v7.2.symbols.gmt',
             c2.all='~/geneSets/mouseCheap/c2.all.v7.4.symbols.gmt')

background = 21281

for(assay in assays)
{
    clusterDF = Read.Table(paste0(assay,'_resolution_2/geneClusters.txt'))
    for(n in names(geneSets))
    {
        geneSet = geneSets[n]
        theseGeneSets = readLines(geneSet)
        theseGeneSets = str_split(theseGeneSets,'\t')

        M = getGeneSetsVsClustersMatrix(theseGeneSets,
                                        clusterDF,
                                        whichFunction,
                                        background)

        df = orderByEachColumn(M,extended=TRUE)
        dirOut = paste0(assay,'_resolution_2/')
        fileOut = paste0(dirOut,
                         n,
                         '_',
                         assay,
                         '.txt')
        Write.Table(df,
                    fileOut)
    }
}

