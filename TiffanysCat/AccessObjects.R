
library(Seurat)

## ####################################################
getObjectPair = function(case)
{
    fName = paste0('SeuratObject/f_',
                   case,
                   '.rds')
    fPrimeName = str_replace(fName,'f_','fPrime_')

    f = readRDS(fName)
    if(case == 'combined')
        f@active.assay = 'integrated'
    fPrime = readRDS(fPrimeName)

    answer = list()
    answer$f = f
    answer$fPrime = fPrime

    return(answer)
}

## ####################################################
getSeuratObject = function(which,clustersAsCharacter=TRUE)
{
    if(which == 'combined')
    {
        fileName =  paste0('~/TiffanyZebraFish/',
                           'SC18139_and_GSE_152906/data/provisional.rds')
        f = readRDS(fileName)
        
        if(clustersAsCharacter)
            f@meta.data$seurat_clusters = as.character(f@meta.data$seurat_clusters)

        ## Eliminate the Unresolved:
        idx = str_detect(f$shorterName,'Unresolved')
        f = f[,!idx]

        return(f)
    }

    if(which == 'adult')
    {
        fileName = paste0('~/TiffanyZebraFish/adultAlone/',
                          '/data/adultSeuratObjectSubclustered.rds')
        f = readRDS(fileName)
        if(clustersAsCharacter)
            f@meta.data$seurat_clusters = as.character(f@meta.data$seurat_clusters)

        return(f)
    }
    if(which == 'larval')
    {
        fileName = paste0('~/TiffanyZebraFish/larvalZebraFish/data/',
                          'SeuratObject.rds')
        f = readRDS(fileName)
        if(clustersAsCharacter)
            f@meta.data$seurat_clusters = as.character(f@meta.data$seurat_clusters)
        return(f)
    }

    stop(paste('No Seurat object for',which))        
}
