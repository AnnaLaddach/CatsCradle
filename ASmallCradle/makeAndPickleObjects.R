
library(Seurat)

rm(list=ls())

source('CradleWare.R')

objectDir = 'SeuratObject'
if(! dir.exists(objectDir))
    dir.create(objectDir)

assays = c('integrated','RNA')

for(assay in assays)
{
    ## Read the original object:
    f = readRDS('~/FranzeSingleCell07/SeuratObject/SeuratObject.rds')
    f@active.assay = assay

    ## Remove rows with no variation:
    M = getM(f)
    delta = c()
    for(i in 1:nrow(M))
        delta[i] = max(M[i,]) - min(M[i,])

    idx = delta > 0

    print(assay)
    print(table(idx))
    f = f[idx,]

    ## Rename the cells:
    f = reviseF(f)

    ## Transpose:
    res = 1
    fPrime = makeFPrime(f,assay,res=res)

    ## Save the objects:
    fOut = paste0(objectDir,
                 '/f_',
                 assay,
                 '.rds')
    ## saveRDS(f,fOut)

    fPrimeOut = paste0(objectDir,
                      '/fPrime_',
                      assay,
                      '_',
                      res,
                      '.rds')
    saveRDS(fPrime,fPrimeOut)
}


    
    
    
