
library(Seurat)
library(HandyPack)
library(ggplot2)
library(stringr)
library(tictoc)
library(plotly)

rm(list=ls())
graphics.off()

source('../ASmallCradle/CradleWare.R')
## ####################################################
## ####################################################
makeUMAPPlot = function(f,title='',which,size=1)
{
    a = FetchData(f,c('UMAP_1','UMAP_2'))
    df = data.frame(UMAP_1=as.numeric(a[,1]),
                    UMAP_2=as.numeric(a[,2]),
                    cluster=f@meta.data$seurat_clusters)
    if(which == 'genes')
    {
        df$label = colnames(f)
    } 
    
    clusters = unique(df$cluster)
    clusters = clusters[order(clusters)]
    df$cluster = factor(df$cluster,levels=clusters)
    df = df[order(df$cluster),]
    N = length(clusters)

    legendDotSize = 4
    if(which == 'genes')
    {
        g = ggplot(df,aes(x=UMAP_1,y=UMAP_2,color=cluster,label=label))
    } else {
        g = ggplot(df,aes(x=UMAP_1,y=UMAP_2,color=cluster))
    }
    g = g  + 
        geom_point(size=size) +
        ggtitle(title) +
        scale_color_manual(breaks=unique(df$cluster),
                           values=getThePalette(N)) +
        guides(color = guide_legend(override.aes = list(size=legendDotSize)))
    return(g)
}

## ####################################################
labelRows = function(M,f,nameCol,cutoff=.3)
{
    ## Assumption here is that rownames(M) is seurat_clusters
    for(i in 1:nrow(M))
    {
        idx = f$seurat_clusters == rownames(M)[i]
        theNames = f@meta.data[idx,nameCol]
        df = data.frame(table(theNames))
        df = df[order(-df$Freq),]
        df$frac = df$Freq / length(theNames)
        df = df[df$frac >= cutoff,]
        rownames(M)[i] = paste(df$theNames,collapse='_')
    }
    return(M)
}

## ####################################################
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
        g = makeUMAPPlot(f,title=paste(case,'Cells'),which='cells')
        print(g)
        ggsave(plot=g,
               filename=paste0('figures/umapOfCells_',
                               case,
                               '.jpg'))
    }

    gPrime = makeUMAPPlot(fPrime,title=paste(case,'Genes'),which='genes',size=1.5)
    dev.new()
    print(gPrime)
 
    ggsave(plot=gPrime,
           filename=paste0('figures/umapOfGenes_',
                           case,
                           '.jpg'))

    pPrime = ggplotly(gPrime,tooltip='label')
    print(pPrime)

    fileName = paste0('figures/umapOfGenes_',
                      case,
                      '.html')
    saveBrowseable(as_widget(pPrime),fileName)

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

    nameCols = c(larval='FinalCat',adult='shortName',combined='shorterName')


    M = labelRows(M,f,nameCol=nameCols[case])

    disambiguation = c('Cells_','Genes_')
    sankeyPair = sankeyPairFromMatrix(M,disambiguation)
    print(sankeyPair$up)
    
    saveSankeyGraph(sankeyPair$up,
                    paste0('figures/CellsGenesExpressionUpSankeyGraph_',
                           case,
                           '.html'))
    toc()
}





