
library(Seurat)
library(ggplot2)

rm(list=ls())
graphics.off()

source('CradleWare.R')
source('~/FranzeSingleCell07/SeuratWare02.R')

## ###################################################
## ###################################################
orderMByCellType = function(f)
{
    shortNames = getShortNames()
    cluster = c()
    for(i in 1:length(shortNames))
    {
        shortName = shortNames[i]
        A = getM(f[,f@meta.data$shortName == shortName])
        cluster = c(cluster,rep(i,ncol(A)))
        if(i == 1)
        {
            M = A
        } else {
            M = cbind(M,A)
        }
    }
    
    return(list(M,cluster))
}
    
## ###################################################
    
    res = 1
assay = 'integrated'
objectPair = getObjectPair(assay,res)
f = objectPair$f
fPrime = objectPair$fPrime

stuff = orderMByCellType(f)
M = stuff[[1]]
cellCluster = stuff[[2]]
numCells = 100
whichOnes = sample(ncol(M),numCells)
whichOnes = whichOnes[order(whichOnes)]
M = M[,whichOnes]
cellCluster = cellCluster[whichOnes]

## get some genes:
startingFrom = c(0,4,8)
for(start in startingFrom)
{
    clusters = start:(start+3)
    howManyPerCluster = 6
    
    df = data.frame(cluster=numeric(0),
                    gene=character(0),
                    expression=numeric(0),
                    x=numeric(0),
                    tag=numeric(0),
                    cellCluster=numeric(0))
    
    for(cluster in clusters)
    {
        genes = colnames(fPrime[,fPrime@meta.data$seurat_clusters == cluster])
        genes = sample(genes,howManyPerCluster)
        
        for(i in 1:length(genes))
        {
            gene = genes[i]
            a = data.frame(cluster,
                           gene,
                           expression=M[gene,],
                           x=1:numCells,
                           tag=i,
                           cellCluster)
            df = rbind(df,a)
        }
    }
    
    df$tag = factor(df$tag)
    df$cluster = paste('Gene cluster',df$cluster)
    df$cellCluster = factor(df$cellCluster)

    title = paste('Gene clusters',start,'-',start+3)
    g = ggplot(df,aes(x=x,y=expression,color=tag)) +
        geom_line(size=.2) +
        geom_point(aes(color=cellCluster),size=.5) +
        facet_wrap(~cluster,ncol=1) +
        ggtitle(title) +
        theme(legend.position='none')

    dev.new()
    print(g)

}
