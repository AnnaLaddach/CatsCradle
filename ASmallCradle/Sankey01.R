
library(Seurat)
library(networkD3)
library(dplyr)

source('CradleWare.R')

assays = c('integrated','RNA')
for(assay in assays)
{
    objectPair = getObjectPair(assay)
    
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
    
    expressionDF = getAbsoluteExpressionTotals(f,fPrime)
    value = expressionDF$expression
    
    links = data.frame(source,target,value,
                       stringsAsFactors=FALSE)
    
    nodes = unique(c(links$source,links$target))
    nodes = data.frame(name=nodes,
                       stringsAsFactors=FALSE)
    
    links$IDsource = match(links$source,nodes$name) - 1
    links$IDtarget = match(links$target,nodes$name) - 1
    
    p = sankeyNetwork(Links = links, Nodes = nodes,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name", 
                      sinksRight=FALSE,
                      fontSize=10)
    
    print(p)
    fileName = paste0(assay,'_resolution_2/sankeyGraphAbsoluteValue.html')
    htmlwidgets::saveWidget(as_widget(p), fileName)
}
