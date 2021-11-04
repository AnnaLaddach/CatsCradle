
library(Seurat)
library(BioHandy)
library(dplyr)
library(networkD3)
library(ggplot2)
library(plotly)
library(cowplot)
library(HandyPack)
library(stringr)


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

            finger = finger + 1
        }
    }

    
    df$cellCluster = factor(df$cellCluster,levels=cellClusters)
    df$geneCluster = factor(df$geneCluster,levels=geneClusters)
    
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
                                       backgroundGenes,
                                       pathways=TRUE)
{
    stopifnot(whichFunction %in% c('log','density'))
    stopifnot(pathways)

    background = length(backgroundGenes)
    
    clusters = unique(clusterDF$geneCluster)
    clusters = clusters[order(clusters)]
    NClusters = length(clusters)
    NGeneSets = length(geneSets)

    names(geneSets) = unlist(lapply(geneSets,function(x) return(x[1])))
    for(i in 1:length(geneSets))
    {
        geneSets[[i]] = geneSets[[i]][3:length(geneSets[[i]])]
        geneSets[[i]] = intersect(geneSets[[i]],backgroundGenes)
    }

    M = matrix(0,nrow=NGeneSets,ncol=NClusters)
    rownames(M) = names(geneSets)
    colnames(M) = as.character(clusters)

    
    for(i in 1:NGeneSets)
    {
        for(j in 1:NClusters)
        {
            cluster = clusters[j]
            idx = clusterDF$geneCluster == cluster
            clusterGenes = clusterDF$gene[idx]

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
    
## ####################################################
saveBrowseable = function(p,fileName)
{
     htmlwidgets::saveWidget(as_widget(p), fileName)
}

## ###################################################
makeOverlapMatrix = function(A,B,backgroundGenes,takeLog=TRUE)
{
    m = length(A)
    n = length(B)
    background = length(backgroundGenes)

    M = matrix(0,nrow=m,ncol=n)
    if(! is.null(names(A)))
        rownames(M) = names(A)
    if(! is.null(names(B)))
        colnames(M) = names(B)

    for(i in 1:m)
        A[[i]] = intersect(A[[i]],backgroundGenes)

    for(j in 1:n)
        B[[j]] = intersect(B[[j]],backgroundGenes)

    for(i in 1:m)
        for(j in 1:n)
        {
            aHat = length(A[[i]])
            bHat = length(B[[j]])
            cHat = length(intersect(A[[i]],B[[j]]))
            M[i,j] = geneListPValue(aHat,bHat,cHat,
                                    background)
        }

    if(takeLog)
        M = -log10(M)

    return(M)
}

## ####################################################
makeUMAPPlot = function(f,title='',which,size=1)
{
    a = FetchData(f,c('UMAP_1','UMAP_2'))
    df = data.frame(UMAP_1=as.numeric(a[,1]),
                    UMAP_2=as.numeric(a[,2]),
                    cluster=f@meta.data$seurat_clusters)
    if(which == 'genes')
    {
        df$label = colnames(f)
    } 
    
    clusters = unique(df$cluster)
    clusters = clusters[order(clusters)]
    df$cluster = factor(df$cluster,levels=clusters)
    df = df[order(df$cluster),]
    N = length(clusters)

    legendDotSize = 4
    if(which == 'genes')
    {
        g = ggplot(df,aes(x=UMAP_1,y=UMAP_2,color=cluster,label=label))
    } else {
        g = ggplot(df,aes(x=UMAP_1,y=UMAP_2,color=cluster))
    }
    g = g  + 
        geom_point(size=size) +
        ggtitle(title) +
        scale_color_manual(breaks=unique(df$cluster),
                           values=getThePalette(N)) +
        guides(color = guide_legend(override.aes = list(size=legendDotSize)))

    answer = list()
    answer$g = g

    if(which == 'genes')
    {
        answer$p = ggplotly(g,tooltip=c('label','cluster'))
    }

    return(answer)
}

## ####################################################
plotPCs = function(df,dataCols,coloring,labeling,outDir)
{
    dim = length(dataCols)

    nameAndMakeDir(outDir)

    plotList = list()
    finger = 1

    for(i in 1:dim)
    {
        for(j in 2:dim)
        {
            g = ggplot(df,aes_string(x=dataCols[j],
                                     y=dataCols[i],
                                     color=coloring,
                                     gene='gene')) +
                geom_point()

            if(i < j)
            {
                p = ggplotly(g,tooltip=labeling)
                fileName = paste0(outDir,
                               '/',
                               dataCols[i],
                               '_vs_',
                               dataCols[j],
                               '.html')
                
                htmlwidgets::saveWidget(as_widget(p),fileName)
            }

            if(max(i,j) > 4)
                next

            if(j <= i)
                plotList[[finger]] = NULL
            else
                plotList[[finger]] = g
            finger = finger + 1
        }
    }

    if(! dim %in% c(3,4))
    {
        writeLines('Currently supporting plots for 3 or 4 dimensional data')
        return(NULL)
    }
    
    
    if(dim == 3)
        h = plot_grid(plotList[[1]],plotList[[2]],plotList[[3]],plotList[[4]],
                      ncol=dim-1)
    
    if(dim == 4)
        h = plot_grid(plotList[[1]],plotList[[2]],plotList[[3]],plotList[[4]],
                      plotList[[5]],plotList[[6]],plotList[[7]],plotList[[8]],
                      plotList[[9]],
                      ncol=dim-1)

                 fileName = paste0(outDir,
                                   '/dimensions1-',
                                   length(dataCols),
                                   '.jpg')
    ggsave(plot=h,
           filename=fileName,
           width=15,height=15,units='in')
  
                
    return(h)
}

## ####################################################
widenGeneClusterLists = function(fileIn,fileOut)
{
    df = Read.Table(fileIn)
    headings = names(df)

    df[,2] = as.numeric(df[,2])
    clusters = unique(df[,2])
    clusters = clusters[order(clusters)]

    geneLists = list()
    for(i in 1:length(clusters))
    {
        idx = df[,2] == clusters[i]
        tag = paste0('geneCluster_',clusters[i])
        geneLists[[tag]] = df[idx,1]
    }
    longest = max(unlist(lapply(geneLists,length)))

    pad = function(genes,N)
    {
        n = max(0,N - length(genes))
        return(c(genes,rep('',n)))
    }
    for(tag in names(geneLists))
        geneLists[[tag]] = pad(geneLists[[tag]],longest)

    for(i in 1:length(geneLists))
    {
        tag = names(geneLists)[i]
        a = data.frame(genes=geneLists[[tag]])
        names(a) = tag
        if(i == 1)
            wide = a
        else
            wide = cbind(wide,a)
    }
    Write.Table(wide,
                fileOut)
}

## ####################################################
trimGeneListTable = function(fileIn,fileOut,cutoff=0.05)
{
    df = Read.Table(fileIn)
    names(df) = str_replace_all(names(df),'X','')

    ## Blank out the useless stuff:
    numCol = ncol(df)
    for(I in seq(from=1,to=numCol,by=2))
    {
        J = I + 1

        for(i in 1:nrow(df))
        {
            logPValue = df[i,J]
            if(logPValue < -log10(cutoff))
            {
                df[i,I] = ''
                df[i,J] = -1
            }
        }
    }

    ## Find the longest:
    longest = c()
    for(I in seq(from=1,to=numCol,by=2))
    {
        idx = df[,I] != ''
        longest = c(longest,sum(idx))
    }
    longest = max(longest)

    ## Trim:
    df = df[1:longest,]

    ## Get rid of the -1s:
    for(I in seq(from=1,to=numCol,by=2))
    {
        J = I + 1
        idx = df[,J] == -1
        df[,J] = as.character(df[,J])
        df[idx,J] = ''
    }
    Write.Table(df,
                fileOut)
}

    
