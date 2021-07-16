
library(Seurat)
library(HandyPack)

source('~/FranzeSingleCell07/SeuratWare02.R')

f = getSeuratObject()
N = ncol(f)

oldNames = colnames(f)
newNames = sprintf('C%.5d',1:N)

df = data.frame(newNames,oldNames)
Write.Table(df,
            'CellDictionary.txt')


