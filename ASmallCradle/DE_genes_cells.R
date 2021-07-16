
library(Seurat)
library(HandyPack)
library(tictoc)

rm(list=ls())

source('CradleWare.R')

## ###################################################
## ###################################################
prepend = function(df)
{
    df = cbind(data.frame(id=rownames(df),
                          stringsAsFactors=FALSE),
               df)

    return(df)
}

## ###################################################
assays = c('integrated','RNA')

for(assay in assays)
{
    dirs = c(DEdir=paste0(assay,'_resolution_2/DE'),
             f=paste0(assay,'_resolution_2/DE/f'),
             fPrime=paste0(assay,'_resolution_2/DE/fPrime'))
    
    for(d in dirs)
        if(! dir.exists(d))
            dir.create(d)
    
    objectPair = getObjectPair(assay)
    f = objectPair$f
    fPrime = objectPair$fPrime
    
    ## ###################################################    
    ## Go after the DE stuff:
    
    ## ###################################################    
    ## For f:
    groups = unique(f@meta.data$seurat_clusters)
    groups = groups[order(groups)]

    for(g in groups)
    {
        writeLines('XXXXXXXXXXXXXXXXXXXXXXXXXXX')
        Tic(paste(assay,'f, group',g))
        
        groupMarkerDF = FindMarkers(f,
                                    only.pos=TRUE,
                                    group.by='seurat_clusters',
                                    ident.1=g)

        groupMarkerDF = prepend(groupMarkerDF)
        
        ## ###################################################
        ## Subset and save:
        cutoff = 0.05
        idx = groupMarkerDF$p_val_adj <= cutoff
        groupMarkerDF = groupMarkerDF[idx,]
        
        if(nrow(groupMarkerDF) == 0)
            next
        
        fileName = paste0(dirs[2],'/group_',g,'.txt')
        Write.Table(groupMarkerDF,
                    fileName)
        
        toc()
    }
    
    ## ###################################################    
    ## For fPrime:
    fPrime@meta.data$seurat_clusters = as.character(fPrime@meta.data$seurat_clusters)
    groups = unique(fPrime@meta.data$seurat_clusters)
    groups = groups[order(groups)]
    
    for(g in groups)
    {
        writeLines('XXXXXXXXXXXXXXXXXXXXXXXXXXX')
        Tic(paste(assay,'fPrime, group',g))
        
        groupMarkerDF = FindMarkers(fPrime,
                                    only.pos=TRUE,
                                    group.by='seurat_clusters',
                                    ident.1=g)

        groupMarkerDF = prepend(groupMarkerDF)        
        
        ## ###################################################
        ## Subset and save:
        cutoff = 0.05
        idx = groupMarkerDF$p_val_adj <= cutoff
        groupMarkerDF = groupMarkerDF[idx,]
        
        if(nrow(groupMarkerDF) == 0)
            next
        
        fileName = paste0(dirs[3],'/group_',g,'.txt')
        Write.Table(groupMarkerDF,
                    fileName)
        
        toc()        
    }
}
