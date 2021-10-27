
library(Seurat)
library(ggplot2)
library(plotly)
library(HandyPack)

rm(list=ls())
graphics.off()

source('../CradleWare.R')
source('AccessObjects.R')

## ####################################################
assay = 'integrated'
res = 1
objectPair = getObjectPair(assay=assay,res=res)

fPrime = objectPair$fPrime

plots = makeUMAPPlot(fPrime,'integrated res=1 gene clusters','genes')
p = plots$p

fileName = paste0(assay,
                  '_resolution_',
                  res,
                  '/genesUMAP.html')

saveBrowseable(p,fileName)

