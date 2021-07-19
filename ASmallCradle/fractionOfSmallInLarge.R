
library(HandyPack)
library(BioHandy)
library(ggplot2)
library(pheatmap)

rm(list=ls())
graphics.off()

## ###################################################
matrixToDF = function(M)
{
    clusterSmall = c()
    clusterLarge = c()
    value = c()
    finger = 1
    for(i in 1:nrow(M))
        for(j in 1:ncol(M))
        {
            clusterLarge[finger] = rownames(M)[i]
            clusterSmall[finger] = colnames(M)[j]
            value[finger] = M[i,j]
            finger = finger + 1
        }
    df = data.frame(clusterLarge,clusterSmall,value)
    df$clusterLarge = factor(df$clusterLarge,levels=rownames(M))
    df$clusterSmall = factor(df$clusterSmall,levels=colnames(M))

    return(df)
}

## ###################################################
## ###################################################
res = 1
smallDF = Read.Table(paste0('integrated_resolution_',res,'/geneClusters.txt'))
largeDF = Read.Table(paste0('RNA_resolution_',res,'/geneClusters.txt'))

clustersSmall = unique(smallDF$geneCluster)
clustersLarge = unique(largeDF$geneCluster)

clustersSmall = clustersSmall[order(clustersSmall)]
clustersLarge = clustersLarge[order(clustersLarge)]

NSmall = length(clustersSmall)
NLarge = length(clustersLarge)

MFraction = matrix(0,nrow=NLarge,ncol=NSmall)
MPValue = matrix(0,nrow=NLarge,ncol=NSmall)

rownames(MFraction) = as.character(clustersLarge)
colnames(MFraction) = as.character(clustersSmall)

rownames(MPValue) = as.character(clustersLarge)
colnames(MPValue) = as.character(clustersSmall)

background = 21281

for(i in 1:NLarge)
{
    largeIdx = largeDF$geneCluster == clustersLarge[i]
    genesLarge = largeDF$genes[largeIdx]
    
    for(j in 1:NSmall)
    {
        smallIdx = smallDF$geneCluster == clustersSmall[j]
        genesSmall = smallDF$genes[smallIdx]

        A = length(genesLarge)
        B = length(genesSmall)
        C = length(intersect(genesLarge,genesSmall))

        MFraction[i,j] = C / (A * B)
        MPValue[i,j] = -log10(geneListPValue(A,B,C,background))
    }
}

fractionDF = matrixToDF(MFraction)
pValueDF = matrixToDF(MPValue)

jpeg(paste0('figures/fractionOfIntegratedInRNAClusters_res_',res,'.jpg'),
     height=8,width=8,units='in',res=100)
pheatmap(MFraction,
         treeheight_row=0,
         treeheight_col=0)
dev.off()


jpeg(paste0('figures/minusLogPForOverlapOfIntegratedAndRNAClusters_res_',res,'.jpg'),
     height=8,width=8,units='in',res=100)
pheatmap(MPValue,
         treeheight_row=0,
         treeheight_col=0)
dev.off()
