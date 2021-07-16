
library(Seurat)

rm(list=ls())

## ###################################################
reviseF = function(f)
{
    new.names = paste0('C',1:ncol(f))
    F = RenameCells(f,new.names=new.names)
    
    return(F)
}


## ###################################################
runTest = function(f)
{
    df = FindMarkers(f,
                     assay=f@active.assay,
                     only.pos=TRUE,
                     group.by=f@meta.data$seurat_clusters,
                     ident.1=5,
                     test='MAST')
    return(df)
}

## ###################################################
fff = readRDS('~/FranzeSingleCell07/SeuratObject/SeuratObject.rds')

df = runTest(fff)

newNames = paste0('C',1:ncol(fff))

dictionary = Read.Table('CellDictionary.txt')
newNames = dictionary$newNames[1:ncol(fff)]

## FFF = RenameCells(fff,new.names=newNames)
FFF = reviseF(fff)

DF = runTest(FFF)

