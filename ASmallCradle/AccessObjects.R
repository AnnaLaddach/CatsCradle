
library(Seurat)

## ###################################################
getObjectPair = function(assay,res=1)
{
    objectPair = list()

    objectDir = '~/CatsCradle/ASmallCradle/SeuratObject'
    fName = paste0(objectDir,
                 '/f_',
                 assay,
                 '.rds')
    objectPair[['f']] = readRDS(fName)

    fPrimeName = paste0(objectDir,
                      '/fPrime_',
                      assay,
                      '_',
                      res,
                      '.rds')
    objectPair[['fPrime']] = readRDS(fPrimeName)

    return(objectPair)
}

## ###################################################
getSeuratObject = function()
{
    fileName = '~/FranzeSingleCell07/SeuratObject/SeuratObject.rds'
    f = readRDS(fileName)

    return(f)
}

