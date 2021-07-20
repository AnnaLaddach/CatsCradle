
library(Seurat)

getObjectQuad = function(assay,res=1)
{
    objectDir = 'SeuratObject'

    fileName = paste0(objectDir,
                      '/',
                      assay,
                      '_resolution_',
                      res,
                      '_objectQuad.RData')

    load(fileName)

    return(objectQuad)
}
