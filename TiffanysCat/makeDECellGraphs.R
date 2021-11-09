
library(ggplot2)
library(HandyPack)
library(stringr)

rm(list=ls())
graphics.off()

source('../CradleWare.R')
source('AccessObjects.R')

## ####################################################
## ####################################################

cases = c('adult','combined')
cellClusterName = c(adult='shortName',
                    combined='clusterName')

for(case in cases[2])
{
    pair = getObjectPair(case)
    f = pair$f
    fPrime = pair$fPrime

    df = data.frame(geneCluster=numeric(0),
                    cellCluster=character(0),
                    cellClusterCount=numeric(0),
                    cellClusterFraction=numeric(0))                    
    
    geneClusters = unique(fPrime$seurat_clusters)
    geneClusters = geneClusters[order(as.numeric(geneClusters))]

    cellClusters = unique(f@meta.data[,cellClusterName[case]])

    if(case == 'adult')
        cellClusters = cellClusters[order(cellClusters)]
    else
    {
        theirOrder = str_replace(cellClusters,'C','')
        cellClusters = cellClusters[order(as.numeric(theirOrder))]
    }

    tableDir = c('tables',
                 paste0('DECells_',case))
    tableDir = nameAndMakeDir(tableDir)

    for(geneCluster in geneClusters)
    {
        tag = paste0('seurat_cluster_',geneCluster)
        fileIn = paste0(tableDir,
                        '/DECells_',
                        tag,
                        '.txt')
        if(! file.exists(fileIn))
        {
            writeLines(paste('skipping',fileIn))
            next
        }
        dfIn = Read.Table(fileIn)

        count = rep(0,length(cellClusters))
        names(count) = cellClusters

        a = data.frame(table(dfIn$cellType))
        for(i in 1:nrow(a))
            count[a[i,1]] = a[i,2]

        b = data.frame(geneCluster=geneCluster,
                       cellCluster=names(count),
                       cellClusterCount=count)
        b$cellClusterFraction = b$cellClusterCount / sum(b$cellClusterCount)

        df = rbind(df,b)
    }

    if(case == 'combined')
        df$cellCluster = factor(df$cellCluster,levels=cellClusters)

    df$geneCluster = factor(df$geneCluster,levels=geneClusters)

    title = paste(case,'DE Cells per gene cluster')
    g = ggplot(df,aes(x=cellCluster,y=cellClusterFraction,fill=cellCluster)) +
        geom_col() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        facet_wrap(~geneCluster) +
        ggtitle(title)
    
    dev.new()
    print(g)

    fileName = paste0('figures/',
                      case,
                      '_DECellsPerGeneCluster.jpg')
    ggsave(plot=g,
           filename=fileName,
           width=15,height=10,unit='in')
}
