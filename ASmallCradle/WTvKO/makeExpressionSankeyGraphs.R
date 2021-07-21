
library(Seurat)
library(networkD3)
library(plotly)
library(pheatmap)

source('../CradleWare.R')

res = 1
assays = c('integrated','RNA')

for(assay in assays)
{
    outDir = paste0(assay,'_resolution_',res)
    if(! dir.exists(outDir))
        dir.create(outDir)
    
    objectQuad = getObjectQuad(assay,res)
    fWT = objectQuad$fWT
    fKO = objectQuad$fKO
    fPrimeWT = objectQuad$fPrimeWT
    fPrimeKO = objectQuad$fPrimeKO

    MWT = getExpressionTotalsMatrix(fWT,fPrimeWT)
    MKO = getExpressionTotalsMatrix(fKO,fPrimeKO)

    pairWT = sankeyPairFromMatrix(MWT)
    pairKO = sankeyPairFromMatrix(MKO)

    print(pairWT[[1]])
    print(pairKO[[1]])

    print(pairWT[[2]])
    print(pairKO[[2]])    
}
    
