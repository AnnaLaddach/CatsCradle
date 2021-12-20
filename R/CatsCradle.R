



## ###################################################
#' Create the transpose of a Seurat object
#'
#' This takes a Seurat object f and creates a new Seurat object whose
#' expression matrix is the transpose of that of f.
#'
#' @param f - a Seurat object
#' @param active.assay - the assay to use.  Defaults to the
#'     active.assay
#' @param npcs - number of principal components, defaults to 30
#' @param dims - dimensions to use for umap and nearest neighbors,
#'     defaults to 1:20
#' @param res - the clustering resolution, defaults to 1
#' @return A Seurat object
#' @export
#' @import Seurat
transposeSeuratObject = function(f,active.assay=f@active.assay,
                                 npcs=30,dims=1:20,res=1)
{
    f@active.assay = active.assay
    M = f@assays[[active.assay]]@data
    M = as.matrix(M)
    MPrime = t(M)
    rownames(MPrime) = colnames(f)
    colnames(MPrime) = rownames(f)

    fPrime = CreateSeuratObject(MPrime,assay=active.assay)

    fPrime = FindVariableFeatures(fPrime)
    fPrime = ScaleData(fPrime)
    fPrime = RunPCA(fPrime,npcs=npcs)
    fPrime = RunUMAP(fPrime,reduction='pca',dims=dims)
    fPrime = FindNeighbors(fPrime)
    fPrime = FindClusters(fPrime,resolution=res)
    
    return(fPrime)
}

## ###################################################
## This retrieves the expression matrix:
##
## f - a Seurat object
getExpression = function(f)
{
    M = f@assays[[f@active.assay]]@data
    M = as.matrix(M)
    rownames(M) = rownames(f)
    
    return(M)
}

## ###################################################
## This computes z-scores:
##
## M - a matrix
zScores = function(M)
{
    mu = rowMeans(M)
    M = M - mu

    for(i in 1:nrow(M))
        M[i,] = M[i,] / sd(M[i,])
    
    return(M)
}


## ###################################################
#' This computes average expression of each gene cluster in
#' each cell cell cluster and returns the result as a matrix
#'
#' @param f - The Seurat object of cells
#' @param fPrime - The Seurat object of genes
#' @return A matrix of the average expression
#' @export
getAverageExpressionMatrix = function(f,fPrime)
{
    f$seurat_clusters = as.character(f$seurat_clusters)
    cellCluster = unique(f$seurat_clusters)
    cellCluster = cellCluster[order(as.numeric(cellCluster))]

    fPrime$seurat_clusters = as.character(fPrime$seurat_clusters)
    geneCluster = unique(fPrime$seurat_clusters)
    geneCluster = geneCluster[order(as.numeric(geneCluster))]
    
    ## Get z-scored expression:
    X = getExpression(f)
    X = zScores(X)
    
    M = matrix(0,nrow=length(cellCluster),ncol=length(geneCluster))
    rownames(M) = cellCluster
    colnames(M) = geneCluster

    for(i in cellCluster)
    {
        for(j in geneCluster)
        {
            idxI = f$seurat_clusters == i
            idxJ = fPrime$seurat_clusters == j

            M[i,j] = sum(X[idxJ,idxI]) / (sum(idxI) * sum(idxJ))
        }
    }
    return(M)
}

## ####################################################
#' This gussies up the rownames and colnames of M
#'
#' @param M - a matrix, typically the average expression matrix
#' @param ccTag - a prefix for the row (cell cluster) names
#' @param gcTag - a prefix for the column (gene cluster names
#' @export
tagRowAndColNames = function(M,ccTag='CC_',gcTag='GC_')
{
    rownames(M) = paste0(rownames(M),ccTag)
    colnames(M) = paste0(colnames(M),gcTag)

    return(M)
}

## ####################################################
#' This converts an average gene expression matrix to a
#' data frame.
#'
#' @param M - An average gene expression matrix.
#' @return A data frame with columns cellCluster, geneCluster
#' and average expression
#' @export
getAverageExpressionDF = function(M)
{
    N = nrow(M) * ncol(M)
    df = data.frame(cellCluster=character(N),
                    geneCluster=character(N),
                    expression=numeric(N))

    finger = 1
    for(i in 1:nrow(M))
    {
        for(j in 1:ncol(M))
        {
            df$cellCluster[finger] = rownames(M)[i]
            df$geneCluster[finger] = colnames(M)[j]
            df$expression[finger] = M[i,j]

            finger = finger + 1
        }
    }

    return(df)
}

## ####################################################
#' This function strips out non-gene information from
#' the beginning of GO sets, etc.
#'
#' @param geneSet - a list of gene sets
#' @return a named list of gene sets
#' @export
stripGeneSet = function(geneSet)
{
    names = character(length(geneSet))
    for(i in 1:length(geneSet))
    {
        names[i] = geneSet[[i]][1]
        geneSet[[i]] = geneSet[[i]][3:length(geneSet[[i]])]
    }

    return(geneSet)
}

## ####################################################
#' This compares the gene clusters to other gene sets
#' e.g., GO, Hallmark, and determines the p-value for
#' their overlaps when compared to a set of background
#' genes.
#'
#' @param geneSets - a named list of gene sets
#' @param clusterDF - a data frame giving the cluster
#' membership of each gene with columns gene and geneCluster
#' @param backgroundGenes - a character vector of genes
#' 
geneSetsVsGeneClustersPValueMatrix = function(geneSets,
                                              clusterDF,
                                              backgroundGenes)
{
    background = length(backgroundGenes)
    clusters = as.numeric(unique(clusterDF$geneCluster))
    clusters = clusters[order(clusters)]
    NClusters = length(clusters)
    NGeneSets = length(geneSets)

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

            ## Maybe there's a better way to do this:
            M[i,j] = geneListPValue(A,B,C,background)
        }
    }
    
    return(M)
}

## ####################################################
#' This orders the gene set p-values (or -log10 p-values) and
#' applies a cutoff (if given) to show only the significant
#' gene sets for each gene cluster
#'
#' @param M - A matrix of gene set p-values (or their logs)
#' to be ordered by their significance
#' @param ascending - Direction in which to order the columns.
#' Defaults to TRUE, so that p-values will be ordered according
#' to decreasing significance, should be set to FALSE if ordering
#' -log p-value
#' @param cutoff - Used to suppress non-significant gene sets in
#' each column.  The value of ascending determines the direction
#' of the cutoff
#' @return This returns a data frame where each gene cluster has
#' two columns, one giving the gene set name and the other giving
#' its significance
orderGeneSetPValues = function(M,ascending=TRUE,cutoff=NULL)
{
    for(i in 1:ncol(M))
    {
        a = data.frame(geneSet=rownames(M),
                       value=M[,i],
                       stringsAsFactors=FALSE)

        ## Order by significance:
        if(ascending)
            a = a[order(a$value),]
        else
            a = a[order(-a$value),]

        ## Supress insignificance:
        if(ascending)
            idx = a$value <= cutoff
        else
            idx = a$value >= cutoff

        a$value = as.character(a$value)
        a$geneSet[!idx] = ''
        a$value[!idx]

        ## Accumulate these:
        if(i == 1)
            df = a
        else
            df = cbind(df,a)
    }
        
    return(df)
}

## ####################################################
sankeyFromMatrix = function(M,disambiguation=c('R_','C_'),
                            fontSize=20,minus='red',plus='blue')
{
    ## Maybe the matrix doesn't have row and column names:
    if(is.null(rownames(M)))
        rownames(M) = paste0(disambiguation[1],1:nrow(M))
    if(is.null(colnames(M)))
        colnames(M) = paste0(disambiguation[1],1:ncol(M))
 
    ## Create the links DF:
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

    ## Color the links DF:
    idx = links$value > 0
    links$group = ''
    links$group[idx] = 'plus'
    links$group[!idx] = 'minus'

    links$value = abs(links$value)

    ## Create the nodes DF:
    nodes = unique(c(links$source,links$target))
    nodes = data.frame(name=nodes,
                       stringsAsFactors=FALSE)

    links$IDsource = match(links$source,nodes$name) - 1
    links$IDtarget = match(links$target,nodes$name) - 1

    linkColor = 'd3.scaleOrdinal() .domain(["minus","plus"]) .range(["X", "Y"])'
    linkColor = str_replace(linkColor,'X',minus)
    linkColor = str_replace(linkColor,'Y',plus)    

    p = sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget", 
                      Value = "value", NodeID = "name", 
                      colourScale=linkColor, LinkGroup="group")

    return(p)
}


## ####################################################
#' This function extracts a nearest neighbor network
#' from a Seurat object
#'
#' @param f - a Seurat object
#' @return - This returns dataframe of neighbors:
#' nodeA - node names for node A 
#' nodeB - node names for node B
#' weight - edge weight
getNearestNeighborListsSeurat = function(f){
    
    #convert to dgTMatrix and extract relevant information
    graph = as(f@graphs$RNA_snn, "dgTMatrix") 
    neighborListDf = data.frame("nodeA" = graph@Dimnames[[1]][graph@i+1],
                                "nodeB" =  graph@Dimnames[[2]][graph@j+1], 
                                "weight" = graph@x)
    
    #remove self-loops
    neighborListDf = 
        neighborListDf[neighborListDf$nodeA != neighborListDf$nodeB,]
    return(neighborListDf)
}


## ####################################################
#' This function generates random indices for node B
#'
#' @param neighborListDf - a dataframe containing the neighborlist
#' @param n - the number of times to randomise indices
#' @param useWeights - whether to preserve edgeweights.
#' @return - a matrix with randomised indices for node B
randomiseNodeIndices = function(neighborListDf, n = 100, useWeights = F){
    
    #determine number of edges and create empty matrix for randomised indices
    nEdges = nrow(neighborListDf)
    indices = 1:nEdges
    randomIndices = matrix(, nrow = nEdges, ncol = n)
    
    #check if weights are to be used
    if (useWeights){
        
        #determine unique weights
        weights = unique(neighborListDf$weight)
        for (i in 1:n){
            
            #randomise indices within each weight category
            for (weight in weights){
                selected = neighborListDf$weight == weight
                randomIndices[selected,i] =
                    sample(indices[selected])
            }
        }
    }
    
    #otherwise ignore weights and randomise indices
    else {
        for (i in 1:n){
            randomIndices[,i] = sample(indices)
        }
    }
    return(randomIndices)
}


## ####################################################
#' This function reads in gene sets in .gmt format
#'
#' @param gmtFile - a .gmt file containing gene sets 
#' @return - A named list of gene sets 
readGmt = function(gmtFile){
 lines = readLines(gmtFile)
 geneSets = list()
 for (line in lines){
     info = strsplit(line, "\t")[[1]]
     geneSets[paste(info[1],info[2])] = list(info[3:length(info)])
 } 
 return(geneSets)
}


## ####################################################
#' This function annotates genes with terms
#'
#' @param gmtFile - a .gmt file containing gene sets to annotate genes with
#' @return - A list where names are genes and values are lists of terms
annotateGenes = function(geneSets){
  #this is a bit slow for large collections of gene sets (a few minutes)
  #think about whether it's necessary to speed up
  genesAnno = list()
  for (geneSet in names(geneSets)){
    for (gene in geneSets[[geneSet]]){
        if (gene %in% names(genesAnno)){
            genesAnno[[gene]] = c(genesAnno[[gene]], geneSet)
        }
       else {
           genesAnno[[gene]] = c(geneSet)
       }
      }
    }
    return(genesAnno)
}


## ####################################################
#' This function predicts the functions of genes based on the functions of their
#' neighbours.
#'
#' @param fPrime - a transposed Seurat object (generated with 
#' transposeSeuratObject())
#' @param genesAnno - genes annotated with gene sets
#' @return - A list where names are genes and values are lists of terms.
#' The values of the lists of terms are calculated according to the weights
#' of the edges connecting the neighbors. 
neighborTerms = function(fPrime,genesAnno, normalise = T){
    
    #determine genes
    genes = rownames(fPrime@graphs$RNA_snn)
    genesPredictedTerms = list()
    i = 1
    
    #iterate through genes
    for (gene in genes){
        if (i %% 100 == 0){
            print(paste(i, "genes processed"))
        }
        i = i + 1
        genesPredictedTerms[[gene]] = list()
        
        #determine neighbors
        neighbors = genes[fPrime@graphs$RNA_snn[gene,] > 0]
        neighbors = neighbors[neighbors != gene]
        
        #calculate total weight of edges to neighbors
        total = sum(fPrime@graphs$RNA_snn[gene,])
        
        #iterate through neighbors
        for (neighbor in neighbors){
            
            #determine weight of connecting edge and normalise if T
            weight = fPrime@graphs$RNA_snn[gene,neighbor]
            if (normalise){
                weight = weight/total
            }
            if (!(neighbor %in% names(genesAnno))){
                next
            }
            
            #extract terms for neighbor
            terms = genesAnno[[neighbor]]
            
            #add these to predicted terms with appropriate weights
            for (term in terms){
                if (term %in% names(genesPredictedTerms[[gene]])){
                    genesPredictedTerms[[gene]][[term]] = 
                        genesPredictedTerms[[gene]][[term]] + weight
                } else {
                    genesPredictedTerms[[gene]][[term]] = weight 
                }
            }
        }
        termNames = names(genesPredictedTerms[[gene]])
        genesPredictedTerms[[gene]] = as.numeric(genesPredictedTerms[[gene]])
        names(genesPredictedTerms[[gene]]) = termNames
        genesPredictedTerms[[gene]] =  
            genesPredictedTerms[[gene]][order(-genesPredictedTerms[[gene]])]
    }
    return(genesPredictedTerms)
}

