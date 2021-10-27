
library(HandyPack)
library(stringr)

rm(list=ls())

source('AccessObjects.R')
source('../CradleWare.R')

## ###################################################
## ###################################################


cases = c('adult','combined')
whichFunction = 'log'

geneSets = c(hallmark='~/geneSets/fish/h.all.v7.2.symbols.gmt',
             c2.all='~/geneSets/fish/c2.all.v7.4.symbols.gmt',
             c5.go.bp='~/geneSets/fish/c5.go.bp.v7.4.symbols.gmt')



for(case in cases)
{
    pair = getObjectPair(case)
    f = pair$f
    
    backgroundGenes = rownames(f)
    
    clusterDF = Read.Table(paste0('tables/geneClusters_',
                                  case,
                                  '.txt'))
    names(clusterDF)[2] = 'geneCluster'
                           
    for(n in names(geneSets))
    {
        geneSet = geneSets[n]
        theseGeneSets = readLines(geneSet)
        theseGeneSets = str_split(theseGeneSets,'\t')

        M = getGeneSetsVsClustersMatrix(theseGeneSets,
                                        clusterDF,
                                        whichFunction,
                                        backgroundGenes=rownames(f))

        df = orderByEachColumn(M,extended=TRUE)
 
        fileOut = paste0('tables/',
                         n,
                         '_',
                         case,
                         '.txt')
        Write.Table(df,
                    fileOut)
    }
}

