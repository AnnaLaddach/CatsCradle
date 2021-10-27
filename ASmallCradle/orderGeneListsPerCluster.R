
library(HandyPack)
library(tictoc)
library(stringr)

rm(list=ls())

source('../CradleWare.R')
source('AccessObjects.R')

## ###################################################
## ###################################################


assays = c('integrated','RNA')
whichFunction = 'log'

geneSets = c(hallmark='~/geneSets/mouse/h.all.v7.2.symbols.gmt',
             c2.all='~/geneSets/mouse/c2.all.v7.4.symbols.gmt',
             c5.go.bp='~/geneSets/mouse/c5.go.bp.v7.4.symbols.gmt')

res = 1

for(assay in assays)
{
    pair = getObjectPair(assay=assay,res=res)
    f = pair$f

    clusterDF = Read.Table(paste0(assay,'_resolution_',res,'/geneClusters.txt'))
    names(clusterDF)[2] = 'geneCluster'
    names(clusterDF)[1] = 'gene'
    
    for(n in names(geneSets))
    {
        Tic(paste(assay,res,n))
        geneSet = geneSets[n]
        theseGeneSets = readLines(geneSet)
        theseGeneSets = str_split(theseGeneSets,'\t')

        M = getGeneSetsVsClustersMatrix(theseGeneSets,
                                        clusterDF,
                                        whichFunction,
                                        backgroundGenes=rownames(f))

        df = orderByEachColumn(M,extended=TRUE)
        dirOut = paste0(assay,'_resolution_',res,'/')
        fileOut = paste0(dirOut,
                         n,
                         '_',
                         assay,
                         '.txt')
        Write.Table(df,
                    fileOut)
        toc()
    }
}

