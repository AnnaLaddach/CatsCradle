
library(Seurat)
library(pheatmap)
library(networkD3)
library(plotly)
library(stringr)

rm(list=ls())
graphics.off()

source('../CradleWare.R')

## ###################################################
## ###################################################
ensureVariablityInBothCases = function(f)
{
    genotypes = c('WT','KO')
    Delta = list()
    for(genotype in genotypes)
    {
        tag = paste0('f',genotype)

        ## We'll want genes that vary in both:
        M = getM(f[,f@meta.data$genotype == genotype])

        delta = c()
        for(i in 1:nrow(M))
            delta[i] = max(M[i,]) - min(M[i,])
        Delta[[tag]] = delta > 0
    }

    idx = Delta[[1]] & Delta[[2]]

    return(idx)
}
        


res = 1
for(assay in c('integrated','RNA')[2])
{
    objectPair = getObjectPair(assay,res)
    f = objectPair$f
    fPrime = objectPair$fPrime

    ## Get rid of genes that don't have variation in both subsets:
    idx = ensureVariablityInBothCases(f)
    f = f[idx,]
    fPrime = fPrime[,idx]

    
    for(genotype in c('WT','KO'))
    {
        fGenotype = f[,f@meta.data$genotype == genotype]
        M = getExpressionTotalsMatrix(fGenotype,fPrime)

        sankeyPair = sankeyPairFromMatrix(M)
        p = sankeyPair$up

        print(p)
        fileName = paste0(assay,'_resolution_',res,
                          '/sankeyUp_',genotype,'.html')
        htmlwidgets::saveWidget(as_widget(p), fileName)


        fileName = str_replace(fileName,'sankeyUp','heatmap')
        fileName = str_replace(fileName,'html','jpg')

        jpeg(fileName,
             height=8,width=8,units='in',res=100)
        pheatmap(M,
                 treeheight_row=0,
                 treeheight_col=0)
        dev.off()
    }
}
