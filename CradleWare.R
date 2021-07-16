

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
reviseF = function(f,active.assay)
{
    M = f@assays[[active.assay]]@data
    M = as.matrix(M)
    
    dictionary = getTranslation()
    colnames(M) = dictionary[colnames(M)]
    
    F = CreateSeuratObject(M,assay=active.assay)
    F@meta.data = f@meta.data
    
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
    M = getM(f)
    M = zScores(M)
    M = abs(M)

    cellClusters = unique(f@meta.data$seurat_clusters)
    cellClusters = cellClusters[order(cellClusters)]

    geneClusters = unique(fPrime@meta.data$seurat_clusters)
    geneClusters = geneClusters[order(geneClusters)]

    N = length(cellClusters) * length(geneClusters)

    df = data.frame(cellCluster=numeric(N),
                    geneCluster=numeric(N),
                    expression=numeric(N),
                    shortName=character(N))

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
                df$expression[finger] = sum(M[geneIdx,cellIdx]) / lambda
            else
                df$expression[finger] = 0
            df$shortName[finger] = getShortName(i)
                

            finger = finger + 1
        }
    }

    
df$cellCluster = factor(df$cellCluster,levels=cellClusters)
df$geneCluster = factor(df$geneCluster,levels=geneClusters)
df$shortName = factor(df$shortName,levels=getShortNames())

  
    return(df)
}
