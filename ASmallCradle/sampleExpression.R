
library(Seurat)
library(ggplot2)

graphics.off()
rm(list=ls())

source('CradleWare.R')

## ###################################################
makeSampleHistograms = function(f,nrow=3,ncol=5,bins=100)
{
    M = getM(f)
    M = zScores(M)
    
    N = nrow * ncol

    select = sample(nrow(M),N)
    df = data.frame(gene=character(0),
                    expression=numeric(0))

    for(i in select)
        df = rbind(df,
                   data.frame(gene=rownames(M)[i],
                              expression=M[i,]))

    g = ggplot(df,aes(x=expression,color=gene)) +
        geom_histogram(bins=bins) +
        facet_wrap(~gene,nrow=nrow) +
        theme(legend.position='none')
    
    return(g)
}

## ###################################################
makeMHistogram = function(f,whichRM,bins=100)
{
    M = getM(f)

    if(whichRM == 'mean')
    {
        m = rowMeans(M)
    }
    else if(whichRM == 'median')
    {
        m = Biobase::rowMedians(M)
    }
    
    df = data.frame(m)

    g = ggplot(df,aes(x=m)) +
        geom_histogram(bins=bins)

    return(g)
}

## ###################################################


## ###################################################

for(active.assay in c('integrated','RNA'))
{
    f = getSeuratObject()
    
    f@active.assay = active.assay
    f = reviseF(f,active.assay)
    
    fPrime = makeFPrime(f,active.assay)

    title = paste('Gene expression',
                  active.assay,
                  'assay')
    
    g = makeSampleHistograms(f) +
        ggtitle(title)

    dev.new()
    print(g)
    
    fileName = paste0(active.assay,
                      '_resolution_2/sampledGeneExpressionHistogram.jpg')
    ggsave(plot=g,
           filename=fileName,
           height=8,width=12,units='in')
    
    title = paste('Cell expression',
                  active.assay,
                  'assay')
    
    h = makeSampleHistograms(fPrime) +
        ggtitle(title)
    dev.new()
    print(h)
    
    fileName = paste0(active.assay,
                      '_resolution_2/sampledCellExpressionHistogram.jpg')
    ggsave(plot=h,
           filename=fileName,
           height=8,width=12,units='in')
}
