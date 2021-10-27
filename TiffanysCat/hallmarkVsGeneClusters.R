
library(pheatmap)
library(HandyPack)
library(stringr)

rm(list=ls())
graphics.off()

stop('This is kind of a hack regarding mouse v fish')

source('../CradleWare.R')

## ####################################################
## ####################################################
getHallmark = function()
{
    fileName = '~/geneSets/mouse/h.all.v7.2.symbols.gmt'
    a = readLines(fileName)
    a = tolower(a)
    a = str_split(a,'\t')
    A = list()
    for(i in 1:length(a))
    {
        tag = a[[i]][1]
        A[[tag]] = a[[i]][3:length(a[[i]])]
    }
    return(A)
}
    
## ####################################################
getGeneClusters = function(case)
{
    fileName = paste0('tables/geneClusters_',
                     case,
                     '.txt')
    df = Read.Table(fileName)

    clusters = unique(df$cluster)
    A = list()
    for(cluster in clusters)
    {
        tag = paste0('Cluster_',cluster)
        idx = df$cluster == cluster
        A[[tag]] = df$gene[idx]
    }
    return(A)
}
    
## ####################################################
## ####################################################

hallmark = getHallmark()

case = 'adult'
geneClusters = getGeneClusters(case)

M = makeOverlapMatrix(hallmark,geneClusters,background=2000,
                      takeLog=FALSE)


