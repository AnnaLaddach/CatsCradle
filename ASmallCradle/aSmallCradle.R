
library(Seurat)
library(ggplot2)
library(stringr)
library(HandyPack)
library(cowplot)

graphics.off()
rm(list=ls())

source('~/FranzeSingleCell07/SeuratWare02.R')
source('CradleWare.R')

## ###################################################
## ###################################################

for(res in 1)
{
    for(active.assay in c('integrated','RNA'))
    {
        objectPair = getObjectPair(active.assay,res=res)
        f = objectPair$f
        fPrime = objectPair$fPrime

        ## f = initialAnalysis(f)
        ## cells = plotUMAP(f,'Cells')
        ## print(cells)
        
        genes = plotUMAP(fPrime,'Genes')
        dev.new()
        print(genes)
        
        df = data.frame(genes=colnames(fPrime),
                        geneCluster=fPrime@meta.data$seurat_clusters)
        df = df[order(df$geneCluster),]
        
        saveDir = paste0(active.assay,
                         '_resolution_',
                         res)
        
        if(!dir.exists(saveDir))
            dir.create(saveDir)
        
        ggsave(plot=genes,
               filename=paste0(saveDir,'/genes.pdf'),
               height=8,width=12,units='in')
        
        Write.Table(df,
                    paste0(saveDir,'/geneClusters.txt'))
    }
}
