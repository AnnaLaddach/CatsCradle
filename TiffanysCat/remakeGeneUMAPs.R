
library(Seurat)
library(HandyPack)
library(tictoc)

rm(list=ls())

source('AccessObjects.R')
source('../CradleWare.R')

## ####################################################
## ####################################################
figDir = nameAndMakeDir('figures')
tableDir = nameAndMakeDir('tables')


cases = c('larval','adult','combined')
for(case in cases)
{
    Tic(case)

    pair = getObjectPair(case)
    fPrime = pair$fPrime

    a = makeUMAPPlot(fPrime,
                     title=paste(case,'Genes'),
                     which='genes',
                     size=1.5)
    pPrime = a$p
    print(pPrime)
    
    fileName = paste0('figures/umapOfGenes_',
                      case,
                      '.html')
    saveBrowseable(as_widget(pPrime),fileName)
}
