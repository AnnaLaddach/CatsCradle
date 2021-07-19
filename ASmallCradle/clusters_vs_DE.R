
library(HandyPack)
library(stringr)
library(networkD3)
library(plotly)
library(pheatmap)

library(htmlwidgets)
library(htmltools)

rm(list=ls())

source('CradleWare.R')
source('~/FranzeSingleCell07/SeuratWare02.R')

## ###################################################
## ###################################################
rename = function(M,toRename,clusterTag,DETag)
{
    stopifnot(toRename %in% c('rows','cols'))
    
    if(toRename == 'rows')
    {
        theNames = rownames(M)
    } else {
        theNames = colnames(M)
    }

    for(i in 1:length(theNames))
        theNames[i] = getShortName(as.numeric(theNames[i]))

    if(toRename == 'rows')
    {
        rownames(M) = theNames
        colnames(M) = paste0(DETag,colnames(M))
    } else {
        colnames(M) = theNames
        rownames(M) = paste0(clusterTag,rownames(M))
    }

    return(M)
}

## ###################################################
assays = c('integrated','RNA')
background = c(integrated=2000,
               RNA=11134)

whichComparison = c('cells','genes')

res = 1
for(assay in assays)
{
    objectPair = getObjectPair(assay,res)
    f = objectPair$f
    fPrime = objectPair$fPrime

    for(which in whichComparison)
    {
        if(which == 'cells')
        {
            obj = f
            DEdir = paste0(assay,'_resolution_',res,'/DE/fPrime')
            toRename = 'rows'
        } else {
            obj = fPrime
            DEdir = paste0(assay,'_resolution_',res,'/DE/f')
            toRename = 'cols'
        }
        clusters = unique(obj$seurat_clusters)
        clusters = clusters[order(clusters)]
        clusters = as.character(clusters)
        obj@meta.data$seurat_clusters = as.character(obj@meta.data$seurat_clusters)

        DEfiles = Sys.glob(paste0(DEdir,'/*txt'))
        DEitems = str_replace(DEfiles,'.*group_','')
        DEitems = str_replace(DEitems,'\\.txt','')

        M = matrix(0,nrow=length(clusters),ncol=length(DEitems))
        rownames(M) = as.character(clusters)
        colnames(M) = as.character(DEitems)

        for(i in 1:length(clusters))
        {
            a = colnames(obj)[obj@meta.data$seurat_clusters == clusters[i]]
            for(j in 1:length(DEitems))
            {
                DEDF = Read.Table(DEfiles[j])
                b = DEDF$id
                A = length(a)
                B = length(b)
                C = length(intersect(a,b))
                M[i,j] = -log10(geneListPValue(A,B,C,background[assay]))
            }

        }

        ## This is a hack:
        Max = 300
        M[M>Max] = Max

        
        clusterTag = paste0(which,'Cluster_')
        DETag = paste0('DE_',which,'_')

        M = rename(M,toRename,clusterTag,DETag)
        
        p = sankeyFromMatrix(M,disambiguation=c(clusterTag,DETag))

##        p =htmlwidgets::prependContent(p,htmltools::tags$h3("Title")) 
        
        print(p)
  
        fileName = paste0(assay,'_resolution_',res,'/',
                          assay,'_',which,
                          'Cluster_vs_DE',which,
                          '.html')
        htmlwidgets::saveWidget(as_widget(p),fileName)

        fileName = str_replace(fileName,'\\.html','Heatmap.jpg')
        jpeg(fileName,
             height=8,width=8,units='in',res=100)
        pheatmap(M,
                 treeheight_row=0,
                 treeheight_col=0)
        dev.off()

        
    }  ## whichComparison                     
} ## assay

                            
