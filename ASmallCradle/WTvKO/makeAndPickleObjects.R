
library(Seurat)

rm(list=ls())

source('../CradleWare.R')

objectDir = 'SeuratObject'
if(! dir.exists(objectDir))
    dir.create(objectDir)

assays = c('integrated','RNA')
genotypes = c('WT','KO')

res = 1
for(assay in assays)
{
    ## Read the original object:
    f = readRDS('~/FranzeSingleCell07/SeuratObject/SeuratObject.rds')
    f@active.assay = assay

    ## Rename the cells:
    newNames = paste0('C',1:ncol(f))
    f = RenameCells(f,new.names=newNames)

    objectQuad = list()
    Delta = list()
    for(genotype in genotypes)
    {
        tag = paste0('f',genotype)
        objectQuad[[tag]] = f[,f@meta.data$genotype == genotype]

        ## We'll want genes that vary in both:
        M = getM(objectQuad[[tag]])

        delta = c()
        for(i in 1:nrow(M))
            delta[i] = max(M[i,]) - min(M[i,])
        Delta[[tag]] = delta > 0
    }

    idx = Delta[[1]] & Delta[[2]]
    for(genotype in genotypes)
    {
        tag = paste0('f',genotype)
        objectQuad[[tag]] = objectQuad[[tag]][idx,]
    }
    
    for(genotype in genotypes)
    {
        tag = paste0('f',genotype)
        tagPrime = paste0('fPrime',genotype)

        ## Transpose:
        objectQuad[[tagPrime]] = makeFPrime(objectQuad[[tag]],assay,res=res)
    }

    fOut = paste0(objectDir,
                  '/',
                  assay,
                  '_resolution_',
                  res,
                  '_objectQuad.RData')

    save(list='objectQuad',
         file=fOut)
}


    
    
    
