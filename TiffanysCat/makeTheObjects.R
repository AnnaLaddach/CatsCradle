
library(Seurat)
library(HandyPack)
library(stringr)
library(tictoc)

rm(list=ls())

source('../ASmallCradle/CradleWare.R')

## Note: the following overwrites getSeuratObject()
source('~/TiffanyZebraFish/ZebraWare.R')

## ####################################################
## ####################################################

cases = c('larval','adult','combined')
for(case in cases)
{
    Tic(case)
    ## Zebra fish object:
    f = getSeuratObject(case)

    if(case == 'combined')
    {
        active.assay = 'integrated'
        f@active.assay = active.assay
    }
    if(case == 'adult')
    {
        ## Subset to var features:
        f = f[f@assays$RNA@var.features,]
    }
    if(case == 'larval')
    {
        f = FindVariableFeatures(f)
        f = f[f@assays$RNA@var.features,]
    }
    
    ## Its transpose:
    fPrime = makeFPrime(f,active.assay=f@active.assay,res=1)
    
    objectDir = nameAndMakeDir('SeuratObject')

    fName = paste0('SeuratObject/f_',
                   case,
                   '.rds')
    fPrimeName = str_replace(fName,'f_','fPrime_')
    
    saveRDS(f,
            fName)
    saveRDS(fPrime,
            fPrimeName)
    toc()
}




