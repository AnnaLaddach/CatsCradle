
library(Seurat)
library(HandyPack)

rm(list=ls())

source('AccessObjects.R')
source('../CradleWare.R')

## ####################################################
## ####################################################
cases = c('adult','combined')
cellClusterName = c(adult='shortName',
                    combined='clusterName')
for(case in cases[2])
{
    pair = getObjectPair(case)
    f = pair$f
    fPrime = pair$fPrime

    ## A hack:
    if(case == 'combined')
        f = RenameCells(object=f,new.names=rownames(fPrime))
    
    DECells = findDECells(f,fPrime,
                          cellClusterName=cellClusterName[case])
    outDir = c('tables',
               paste0('/DECells_',case))
    outDir = nameAndMakeDir(outDir)

    for(tag in names(DECells))
    {
        fileOut = paste0(outDir,
                         '/DECells_',
                         tag,
                         '.txt')
        Write.Table(DECells[[tag]],
                    fileOut)
    }
}

                      
