
library(BioHandy)
library(dplyr)
library(networkD3)
library(plotly)


source('~/FranzeSingleCell07/SeuratWare02.R')


## ###################################################
getM = function(f)
{
    M = f@assays[[f@active.assay]]@data
    M = as.matrix(M)
    rownames(M) = rownames(f)
    
    return(M)
}


## ###################################################
zScores = function(M)
{
    mu = rowMeans(M)
    M = M - mu

    for(i in 1:nrow(M))
        M[i,] = M[i,] / sd(M[i,])
    
    return(M)
}

## ###################################################
getSeuratObject = function()
{
    versionDir = '~/FranzeSingleCell07'
    fileName = paste0(versionDir,'/SeuratObject/SeuratObject.rds')
    SeuratObject = readRDS(fileName)
    
    return(SeuratObject)        
}

## ###################################################
getTranslation = function()
{
    df = Read.Table('CellDictionary.txt')
    dictionary = df$newNames
    names(dictionary) = df$oldNames
    
    return(dictionary)
}

## ###################################################
reviseF = function(f)
{
    new.names = paste0('C',1:ncol(f))
    F = RenameCells(f,new.names=new.names)
    
    return(F)
}

## ###################################################
makeFPrime = function(f,active.assay,npcs=30,dims=1:20,res=2)
{
    fPrime = transposeSeuratObject(f,active.assay)
    fPrime = initialAnalysis(fPrime,npcs,dim)
    fPrime = FindClusters(fPrime,res=res)

    return(fPrime)
}

## ###################################################
transposeSeuratObject = function(f,active.assay)
{
    f@active.assay = active.assay
    M = f@assays[[active.assay]]@data
    M = as.matrix(M)
    MPrime = t(M)
    rownames(MPrime) = colnames(f)
    colnames(MPrime) = rownames(f)

    fPrime = CreateSeuratObject(MPrime,assay=active.assay)
    return(fPrime)
}

## ###################################################
initialAnalysis = function(f,npcs=30,dims=1:20)
{
    f = FindVariableFeatures(f)
    f = ScaleData(f)
    f = RunPCA(f,npcs=30)
    f = RunUMAP(f,reduction='pca',dims=1:20)
    f = FindNeighbors(f,reduction='pca',dims=1:20)

    return(f)
}

## ###################################################
plotUMAP = function(f,title,legend=TRUE)
{
    df = f@reductions$umap@cell.embeddings
    df = data.frame(df)
    df$cluster = f@meta.data$seurat_clusters

    
    breaks = unique(df$cluster)
    breaks = breaks[order(breaks)]
    breaks = factor(breaks)
    df$cluster = factor(df$cluster,levels=breaks)
    
    values = getThePalette(ncolor=length(break))

    legendDotSize=4
    g = ggplot(df,aes(x=UMAP_1,y=UMAP_2,color=cluster)) +
        geom_point(size=.5) +
        scale_color_manual(values=values,
                           breaks=breaks) +
        guides(color = guide_legend(override.aes = list(size=legendDotSize))) +
        ggtitle(title)

    if(!legend)
        g = g +
            theme(legend.position='none')
        
    return(g)
}

## ###################################################
getAbsoluteExpressionTotals = function(f,fPrime)
{
    return(getExpressionTotals(f,fPrime,absolute=TRUE))
}


## ###################################################
getExpressionTotals = function(f,fPrime,absolute=FALSE)
{
    M = getM(f)
    M = zScores(M)

    if(absolute)
        M = abs(M)

    cellClusters = unique(f@meta.data$seurat_clusters)
    cellClusters = cellClusters[order(cellClusters)]

    geneClusters = unique(fPrime@meta.data$seurat_clusters)
    geneClusters = geneClusters[order(geneClusters)]

    N = length(cellClusters) * length(geneClusters)

    df = data.frame(cellCluster=numeric(N),
                    geneCluster=numeric(N),
                    expression=numeric(N))

    finger = 1
    for(i in cellClusters)
    {
        for(j in geneClusters)
        {
            cellIdx = f@meta.data$seurat_cluster == i
            geneIdx = fPrime@meta.data$seurat_cluster == j
            lambda = sum(cellIdx) * sum(geneIdx)
            df$cellCluster[finger] = i
            df$geneCluster[finger] = j
            if(lambda > 0)
            {
                df$expression[finger] = sum(M[geneIdx,cellIdx]) / lambda
            } else {
                df$expression[finger] = 0
            }

            if(is.na(df$expression[finger]))
                stop('NA NA')

            ## Too specific to previous case:
            ## df$shortName[finger] = getShortName(i)
                

            finger = finger + 1
        }
    }

    
    df$cellCluster = factor(df$cellCluster,levels=cellClusters)
    df$geneCluster = factor(df$geneCluster,levels=geneClusters)
    ## df$shortName = factor(df$shortName,levels=getShortNames())
    
    
    return(df)
}

## ###################################################
getExpressionTotalsMatrix = function(f,fPrime,absolute=FALSE)
{
    df = getExpressionTotals(f,fPrime,absolute=FALSE)

    cellCluster = as.character(unique(df$cellCluster))
    geneCluster = as.character(unique(df$geneCluster))

    M = matrix(0,nrow=length(cellCluster),ncol=length(geneCluster))
    rownames(M) = cellCluster
    colnames(M) = geneCluster

    for(i in 1:nrow(df))
        M[df$cellCluster[i],df$geneCluster[i]] = df$expression[i]

    return(M)
}
    

## ###################################################
getGeneSetsVsClustersMatrix = function(geneSets,
                                       clusterDF,
                                       whichFunction,
                                       background=NA,
                                       pathways=TRUE)
{
    stopifnot(whichFunction %in% c('log','density'))
    stopifnot(pathways)
    
    clusters = unique(clusterDF$geneCluster)
    clusters = clusters[order(clusters)]
    NClusters = length(clusters)
    NGeneSets = length(geneSets)

    names(geneSets) = unlist(lapply(geneSets,function(x) return(x[1])))
    for(i in 1:length(geneSets))
        geneSets[[i]] = geneSets[[i]][3:length(geneSets[[i]])]

    M = matrix(0,nrow=NGeneSets,ncol=NClusters)
    rownames(M) = names(geneSets)
    colnames(M) = as.character(clusters)

    
    for(i in 1:NGeneSets)
    {
        for(j in 1:NClusters)
        {
            cluster = clusters[j]
            idx = clusterDF$geneCluster == cluster
            clusterGenes = clusterDF$genes[idx]

            A = length(geneSets[[i]])
            B = length(clusterGenes)
            C = length(intersect(geneSets[[i]],clusterGenes))

            if(whichFunction == 'log')
                M[i,j] = -log10(geneListPValue(A,B,C,background))

            if(whichFunction == 'density')
                M[i,j] = C / (A * B)
        }
    }

    return(M)
}

## ###################################################
orderByEachColumn = function(M,extended=FALSE)
{
    for(i in 1:ncol(M))
    {
        MPrime = M[order(-M[,i]),]
        if(! extended)
        {
            a = data.frame(rownames(MPrime))
            names(a)[1] = colnames(M)[i]
        } else {
            a = data.frame(rownames(MPrime),MPrime[,i])
            col = colnames(M)[i]
            names(a) = c(col,paste0(col,'_value'))
        }
            
        if(i == 1)
            df = a
        else
            df = cbind(df,a)
    }
    return(df)
}

## ###################################################
getObjectPair = function(assay,res=2)
{
    objectPair = list()

    objectDir = '~/CatsCradle/ASmallCradle/SeuratObject'
    fName = paste0(objectDir,
                 '/f_',
                 assay,
                 '.rds')
    objectPair[['f']] = readRDS(fName)

    fPrimeName = paste0(objectDir,
                      '/fPrime_',
                      assay,
                      '_',
                      res,
                      '.rds')
    objectPair[['fPrime']] = readRDS(fPrimeName)

    return(objectPair)
}

## ###################################################
sankeyFromMatrix = function(M,disambiguation=c('R_','C_'),fontSize=20)
{
    stopifnot(min(M) >= 0)

    from = rownames(M)
    to = colnames(M)

    if(length(intersect(from,to)) > 0)
    {
        from = paste0(disambiguation[1],from)
        to = paste0(disambiguation[2],to)
    }

    source = rep(from,each=length(to))
    target = rep(to,length(from))

    value = as.numeric(t(M))

    links = data.frame(source,target,value,
                       stringsAsFactors=FALSE)

    nodes = unique(c(links$source,links$target))
    nodes = data.frame(name=nodes,
                       stringsAsFactors=FALSE)

    links$IDsource = match(links$source,nodes$name) - 1
    links$IDtarget = match(links$target,nodes$name) - 1

    p = sankeyNetwork(Links=links,Nodes=nodes,
                      Source='IDsource',Target='IDtarget',
                      Value='value',NodeID='name',
                      sinksRight=FALSE,
                      fontSize=fontSize)
    
    return(p)    
}

## ###################################################
sankeyPairFromMatrix = function(M,disambiguation=c('R_','C_'),fontSize=20)
{
    sankeyPair = list()

    M_up = M
    M_up[M_up < 0] = 0

    sankeyPair[['up']] = sankeyFromMatrix(M_up,disambiguation,fontSize)

    M_down = -M
    M_down[M_down < 0] = 0
    
    sankeyPair[['down']] = sankeyFromMatrix(M_down,disambiguation,fontSize)

    return(sankeyPair)
}

## ####################################################
saveSankeyGraph = function(p,fileName)
{
     htmlwidgets::saveWidget(as_widget(p), fileName)
}
    
