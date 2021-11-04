
library(HandyPack)

rm(list=ls())

source('../CradleWare.R')
## ####################################################
## ####################################################

cases = c('larval','adult','combined')

for(case in cases)
{
    fileIn = paste0('tables/geneClusters_',
                    case,
                    '.txt')
    fileOut = paste0('tables/geneClusters_',
                     case,
                     '_wide.txt')

    widenGeneClusterLists(fileIn,fileOut)
}
