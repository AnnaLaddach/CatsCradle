
library(Seurat)
library(HandyPack)
library(ggplot2)
library(stringr)
library(tictoc)

rm(list=ls())
graphics.off()

source('../ASmallCradle/CradleWare.R')
## ####################################################
## ####################################################
makeUMAPPlot = function(f,title='',size=1)
{
    a = FetchData(f,c('UMAP_1','UMAP_2'))
    df = data.frame(UMAP_1=as.numeric(a[,1]),
                    UMAP_2=as.numeric(a[,2]),
                    cluster=f@meta.data$seurat_clusters)
    clusters = unique(df$cluster)
    clusters = clusters[order(clusters)]
    df$cluster = factor(df$cluster,levels=clusters)
    df = df[order(df$cluster),]
    N = length(clusters)

    legendDotSize = 4
    g = ggplot(df,aes(x=UMAP_1,y=UMAP_2,color=cluster)) + 
        geom_point(size=size) +
        ggtitle(title) +
        scale_color_manual(breaks=unique(df$cluster),
                           values=getThePalette(N)) +
        guides(color = guide_legend(override.aes = list(size=legendDotSize)))
    return(g)
}

## ####################################################
figDir = nameAndMakeDir('figures')
tableDir = nameAndMakeDir('tables')

cases = c('larval','adult','combined')
for(case in cases)
{
    Tic(case)
    fName = paste0('SeuratObject/f_',
                   case,
                   '.rds')
    fPrimeName = str_replace(fName,'f_','fPrime_')
    
    f = readRDS(fName)
    if(case == 'combined')
        f@active.assay = 'integrated'
    fPrime = readRDS(fPrimeName)

    if(case != 'larval')
    {
        g = makeUMAPPlot(f,paste(case,'Cells'))
        print(g)
        ggsave(plot=g,
               filename=paste0('figures/umapOfCells_',
                               case,
                               '.jpg'))
    }

    gPrime = makeUMAPPlot(fPrime,paste(case,'Genes'))
    dev.new()
    print(gPrime)
    
 
    ggsave(plot=gPrime,
           filename=paste0('figures/umapOfGenes_',
                           case,
                           '.jpg'))
    
    geneDF = data.frame(gene=colnames(fPrime),
                        cluster=fPrime@meta.data$seurat_clusters)
    geneDF = geneDF[order(geneDF$cluster),]
    Write.Table(geneDF,
                paste0('tables/geneClusters_',
                       case,
                       '.txt'))
    
    ## ####################################################
    ## Expression totals:
    M = getExpressionTotalsMatrix(f,fPrime)
    disambiguation = c('Cells_','Genes_')
    sankeyPair = sankeyPairFromMatrix(M,disambiguation)
    print(sankeyPair$up)
    
    saveSankeyGraph(sankeyPair$up,
                    paste0('figures/CellsGenesExpressionUpSankeyGraph_',
                           case,
                           '.html'))
    toc()
}





