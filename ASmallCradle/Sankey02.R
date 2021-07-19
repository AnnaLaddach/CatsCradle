
library(Seurat)
library(networkD3)
library(dplyr)
library(plotly)

source('CradleWare.R')

res = 1
assays = c('integrated','RNA')
for(assay in assays)
{
    objectPair = getObjectPair(assay,res)
    
    f = objectPair$f
    fPrime = objectPair$fPrime
    
    fClusters = unique(f@meta.data$shortName)
    fClusters = fClusters[order(fClusters)]
    fClusters = as.character(fClusters)
    
    fPrimeClusters = unique(fPrime@meta.data$seurat_clusters)
    fPrimeClusters = fPrimeClusters[order(fPrimeClusters)]
    fPrimeClusters = as.character(fPrimeClusters)
    
    source = rep(fClusters,each=length(fPrimeClusters))
    
    target = paste('genes',rep(fPrimeClusters,length(fClusters)))
    
    expressionDF = getExpressionTotals(f,fPrime)
    value = expressionDF$expression
    
    links = data.frame(source,target,value,
                       stringsAsFactors=FALSE)

    up = links$value >= 0
    linksUp = links[up,]
    linksDown = links[!up,]
    linksDown$value = - linksDown$value

    
    nodesUp = unique(c(linksUp$source,linksUp$target))
    nodesUp = data.frame(name=nodesUp,
                       stringsAsFactors=FALSE)
    
    linksUp$IDsource = match(linksUp$source,nodesUp$name) - 1
    linksUp$IDtarget = match(linksUp$target,nodesUp$name) - 1
    
    pUp = sankeyNetwork(Links = linksUp, Nodes = nodesUp,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name", 
                      sinksRight=FALSE,
                      fontSize=10)
    
    print(pUp)
    fileName = paste0(assay,'_resolution_',res,'/sankeyGraphExpressionUp.html')
    htmlwidgets::saveWidget(as_widget(pUp), fileName)

    ## Now down:
    
    nodesDown = unique(c(linksDown$source,linksDown$target))
    nodesDown = data.frame(name=nodesDown,
                       stringsAsFactors=FALSE)
    
    linksDown$IDsource = match(linksDown$source,nodesDown$name) - 1
    linksDown$IDtarget = match(linksDown$target,nodesDown$name) - 1
    
    pDown = sankeyNetwork(Links = linksDown, Nodes = nodesDown,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name", 
                      sinksRight=FALSE,
                      fontSize=10)
    
    print(pDown)
    fileName = paste0(assay,'_resolution_',res,'/sankeyGraphExpressionDown.html')
    htmlwidgets::saveWidget(as_widget(pDown), fileName)    
}
